#Script takes as its sole parameter a model ID & produces a fit model
#and some raw benchmarking data (to be cooked into a report)

#Note the use of wait(): dask's async is mostly obviated to faciliatate profiling
#note that this doesn't slow computation: wait is only used when the next step would
#require the previous step to finish anyway. It's not lazy evaluation if you need the data immediately. 

#general imports
import time
import sys
import pandas as pd
import numpy as np
from formulaic import Formula
from datetime import datetime
import importlib
import pickle
import os
import functools
import json
import psutil

#import time



#import datetime
#import uuid


#dask imports
from dask_jobqueue import SLURMCluster
from dask.distributed import Client, wait
from dask.distributed import get_worker

from dask import delayed


### Utility functions

def abort_on_failure(future,client):
    """
    Call this function with a completed future that is strictly necessary for the task at hand.
    If it didn't work, we'll crash. 
    """
    if future.status=="error":
        fprint("[!] Computation failed. Aborting.")
        fprint("[!] Check slave node logs for details.")
        client.shutdown()
        assert 1==2


def log_usage(log_dir="dask_task_logs"):
    """Decorator factory that profiles functions submitted to dask"""
    os.makedirs(log_dir, exist_ok=True)

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # Try to get the worker name/address
            try:
                worker = get_worker()
                worker_name = worker.name or worker.address
            except ValueError:
                worker_name = "unknown"

            start_wall = time.time()
            start_cpu = psutil.Process().cpu_times()
            process = psutil.Process()
            mem_start = process.memory_info().rss

            result = func(*args, **kwargs)

            end_wall = time.time()
            end_cpu = process.cpu_times()
            mem_end = process.memory_info().rss
            mem_peak = max(mem_start, mem_end)

            stats = {
                "function": func.__name__,
                "worker": worker_name,
                "wall_time_sec": end_wall - start_wall,
                "cpu_user_time": end_cpu.user - start_cpu.user,
                "cpu_system_time": end_cpu.system - start_cpu.system,
                "mem_peak_rss": mem_peak,
            }

            safe_worker_name = str(worker_name).replace(":", "_").replace("/", "_")
            fname = f"{func.__name__}_{safe_worker_name}_{uuid.uuid4().hex[:8]}.json"
            fpath = os.path.join(log_dir, fname)
            with open(fpath, "w") as f:
                json.dump(stats, f, indent=2)

            return result
        return wrapper
    return decorator

def fprint(string):
    """Wraps print with flush=True"""
    print(string,flush=True)

def stamp():
    """returns a formatted string with time since test start"""
    return f'@ {time.time()-start_time:.2f}'

### Functions to be passed as jobs to dask

def statsmodels_fit(p,method="bfgs"):
    import statsmodels.discrete.count_model as smdc
    X,y,Z=p
    zinb_model = smdc.ZeroInflatedNegativeBinomialP(y, X, exog_infl=Z)

    n_count_params = zinb_model.exog.shape[1]      # Count model parameters
    n_infl_params = zinb_model.exog_infl.shape[1]    # Inflation model parameters
    n_total = n_count_params + n_infl_params + 1 # adding 1 for alpha
    start_params = np.full(n_total, 0.1)

    zinb_result = zinb_model.fit(start_params=start_params,maxiter=1000,method=method)

    return zinb_result



def load_csv(path):
    return pd.read_csv(path)

def load_tsv(path):
    return pd.read_csv(path,sep="\t")

@log_usage()
def create_matricies(main_form,zin_form,data):
    y, X=Formula(main_form).get_model_matrix(data,output='pandas')
    Z=Formula(zin_form).get_model_matrix(data,output='pandas')
    return(X, y, Z)



### Main testing routine
start_time=-1
def main():
    global start_time
    start_time=time.time()

    #load model specification
    if len(sys.argv) !=2:
        fprint("[x] Wrong number of arguments. Aborting.")
        return -1
    
    #useful paths
    data_root="/gpfs/gibbs/pi/reilly/tabula_data"

    model_code=sys.argv[1]
    modelspecs=pd.read_csv("modelspecs.tsv",sep="\t",index_col=0)




    fprint(f'[+] Creating cluster {stamp()}')
    
    now=datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    #Make sure memory per slurm job is large enough to hold the data...
    cluster=SLURMCluster(
            cores=2,#cores per slurm job
            memory="16G",#memory per slurm job
            processes=1,#dask workers per slurm job
            job_extra_directives=["-p ycga", 
                f"--job-name=simclust_worker",
                "--time=08:00:00",
                f"--output=slave_{model_code}_{now}_%j.out"]
        )

    client = Client(cluster)

    fprint(f"[+] Cluster started {stamp()}. Monitor on {cluster.dashboard_link}")

    
    #cluster.adapt(minimum_jobs=1, maximum_jobs=10)

    fprint(f'[+] Loading data {stamp()}')
    
    cluster.scale(jobs=1)

    dat_future=None

    if "Fake" in modelspecs.loc[model_code, 'dataset']:
        dat_future=client.submit(load_csv,f"{data_root}/simulated/fake_cres.csv")
    elif "Shendure" in modelspecs.loc[model_code, 'dataset']:
        dat_future=client.submit(load_tsv,f"{data_root}/shendure/shendure_counts_grouped.txt")
    else:
        fprint(f'[x] Unknown dataset. Aborting. {stamp()}')
        return -1


    wait(dat_future)
    abort_on_failure(dat_future,client)

    fprint(f'[+] Creating matricies {stamp()}')

    mats_future = client.submit(create_matricies,
        data=dat_future,
        zin_form=modelspecs.loc[model_code, "z_equ"],
        main_form=modelspecs.loc[model_code, "main_equ"]
    )
    
    wait(mats_future)
    abort_on_failure(mats_future,client)


    fprint(f'[+] Fitting {stamp()}')

    model_future=client.submit(statsmodels_fit,mats_future)

    

    wait(model_future)
    abort_on_failure(model_future,client)
    

    fprint(f'[+] Dumping model {stamp()}')

    #'now' isn't now anymore, but this is easier for pairity/lookups...
    with open(f"{data_root}/speed_test/models/{model_code}_{now}.pkl", "wb") as f:
        pickle.dump(model_future.result(), f)
    

    fprint(f'[+] Shutting down {stamp()}')

    client.shutdown()

if __name__=="__main__":
    main()