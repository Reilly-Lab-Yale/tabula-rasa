#Script takes as its a model ID & produces a fit model
#and some raw benchmarking data


#stock imports
import time
import sys
import os
from datetime import datetime
import pickle
import importlib
import random
import string

# installed imports
import pandas as pd
import numpy as np
import tensorflow as tf
from formulaic import Formula
import dill

#dask imports
from dask_jobqueue import SLURMCluster
from dask.distributed import Client, wait
from dask.distributed import get_worker
from dask.distributed import performance_report
from dask import delayed


### constants

data_root="/gpfs/gibbs/pi/reilly/tabula_data"

#not specified in model table since will req. optimization here

#hard-coded max times for workers used to 
#fit different models.
WORKER_TIMES={
    "c9100090":1,
    "c9100010":1,
    "c9200090":1,
    "c9200010":1,
    "c9010090":1,
    "c9011090":1,
    "c0100000":12,
    "c0100010":12,
    "c0200000":12,
    "c0200010":12,
    "c0011000":12
}

STATSMODELS_MAXITER=1000

MAX_PARALLEL=10

### Utility functions

def crash_if_no_gpu():
    gpus = tf.config.list_physical_devices('GPU')
    if not gpus:
        raise RuntimeError("No GPU devices found, despite the fact that you asked for one.")

def abort_on_failure(future,client):
    """
    Call this function with a completed future that is strictly necessary for the task at hand.
    If it didn't work, we'll crash. 
    """
    status=None
    if type(future)==type([3]):
        #if it's a list
        if any(a_future.status=="error" for a_future in future):
            status="error"
    else:
        #it's not a list, presumably just one future
        status=future.status
    if future.status=="error":
        sprint("[!] Computation failed. Aborting.")
        sprint("[!] Check slave node logs for details.")
        client.shutdown()
        assert 1==2

#pesuto-logging functions for quick debugging
def rand_tag(length=6):
    """
    Some of my print statements seem to be getting duplicated
    so I'm going to stamp them with a random string...
    """
    return ''.join(random.choices(string.ascii_letters + string.digits, k=length))

def stamp():
    """returns a formatted string with time since test start"""
    return f'@ {time.time()-start_time:.4f}'

def fprint(string):
    """Wraps print with flush=True"""
    print(string,flush=True)

def sprint(string):
    """
    stamp-print
    """
    fprint(f"{string} {rand_tag()} {stamp()}")



### Functions to be passed as jobs to dask

def statsmodels_fit(p,method,name):
    sprint(f"[W] [+] Beginning fitting {name}")
    
    import statsmodels.discrete.count_model as smdc
    X,y,Z=p
    zinb_model = smdc.ZeroInflatedNegativeBinomialP(y, X, exog_infl=Z)

    n_count_params = zinb_model.exog.shape[1]      # Count model parameters
    n_infl_params = zinb_model.exog_infl.shape[1]    # Inflation model parameters
    n_total = n_count_params + n_infl_params + 1 # adding 1 for alpha
    start_params = np.full(n_total, 0.1)

    zinb_result = zinb_model.fit(start_params=start_params,maxiter=STATSMODELS_MAXITER,method=method)

    sprint(f"[W] [+] Done fitting {name}. Serializing result for transfer")

    return dill.dumps(zinb_result)

def tensorzinb_fit(p,method,name):
    sprint(f"[W] [+] Beginning tensor fitting {name}")
    
    from tensorzinb.tensorzinb import TensorZINB
    X,y,Z=p

    zinbo=TensorZINB(y["umis_mpra_bc"].to_numpy().reshape((-1,1)),X,exog_infl=Z.to_numpy())#,same_dispersion=True
    zinb_result=zinbo.fit(init_method="nb")

    sprint(f"[W] [+] Done fitting {name}. Serializing result for transfer")
    
    fprint(zinb_result)

    return dill.dumps(zinb_result)

def load_csv(path):
    sprint(f"[W] [+] Loading {path}")
    return pd.read_csv(path)

def load_tsv(path):
    sprint(f"[W] [+] Loading {path}")
    return pd.read_csv(path,sep="\t")

def create_matricies(main_form,zin_form,data):
    sprint("[W] [+] Creating matrix")
    y, X=Formula(main_form).get_model_matrix(data,output='pandas')
    Z=Formula(zin_form).get_model_matrix(data,output='pandas')
    return(X, y, Z)

def dump_unified_model(model_future,filename):
    sprint("[W] [+] Dumping unified model")
    with open(filename, "wb") as f:
        dill.dump(dill.loads(model_future), f)
        f.flush()
        os.fsync(f.fileno())

def dump_broken_model(model_future,types,filename):
    sprint(f"[W] [+] Materalizing & de-seralizing model")
    
    materalized_models = {
        t:dill.loads(model_future[t])
        for t in types
    }

    print(f"[W] [+] Dumping to disc! {stamp()}")
    with open(filename, "wb") as f:
        dill.dump(materalized_models, f)
        f.flush()
        os.fsync(f.fileno())

### Main testing routine
start_time=-1
def main():
    
    global start_time
    start_time=time.time()

    #load model specification
    if len(sys.argv) !=2:
        fprint("[!] Wrong number of arguments. Aborting.")
        return -1
    
    #useful paths
    model_code=sys.argv[1]
    modelspecs=pd.read_csv(f"{data_root}/speed_test/modelspecs.tsv",sep="\t",index_col=0)

    sprint(f'[+] Creating cluster')
    
    now=datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    
    cluster=None

    #Key decision in cluster creation is ± GPU
    #Make sure memory per slurm job is large enough to hold the data...

    hardware=modelspecs.loc[model_code,"hardware"].split("_")[1]

    if hardware=="CPU":
        cluster=SLURMCluster(
                cores=2,#cores per slurm job
                memory="16G",#memory per slurm job
                processes=1,#dask workers per slurm job
                job_extra_directives=["-p ycga", 
                    f"--job-name=simclust_worker",
                    f"--time={WORKER_TIMES[model_code]}:00:00",
                    f"--output=slave_{model_code}_{now}_%j.out"]
            )
    elif hardware=="A100":
        cluster=SLURMCluster(
                cores=2,#cores per slurm job
                memory="16G",#memory per slurm job
                processes=1,#dask workers per slurm job
                job_extra_directives=["-p gpu", 
                    "--gpus=a100:1",
                    "--job-name=simclust_worker",
                    f"--time={WORKER_TIMES[model_code]}:00:00",
                    f"--output=slave_{model_code}_{now}_%j.out"]
            )
    

    client = Client(cluster)

    sprint(f"[+] Cluster started. Monitor on {cluster.dashboard_link}")

    method=modelspecs.loc[model_code,"sm_optimizer"].split("_")[1]

    #start logging!
    with performance_report(filename=f"{model_code}_{now}_report.html"):
        
        #cluster.adapt(minimum_jobs=1, maximum_jobs=10)

        #initial scaling for init tasks : loading data
        cluster.scale(jobs=1)

        dat_future=None

        if "Fake" in modelspecs.loc[model_code, 'dataset']:
            dat_future=client.submit(load_csv,f"{data_root}/simulated/fake_cres.csv")
        elif "Shendure" in modelspecs.loc[model_code, 'dataset']:
            dat_future=client.submit(load_tsv,f"{data_root}/shendure/shendure_counts_grouped.txt")
        else:
            fprint(f'[!] Unknown dataset. Aborting. {stamp()}')
            return -1


        model_future=None

        #for a unified model, we just fit on one node. For a broken model, we will run in parallel.
        if pd.isna(modelspecs.loc[model_code, 'broken_by']):
            sprint(f'[i] Unified model : executing in one job.')

            mats_future = client.submit(create_matricies,
                data=dat_future,
                zin_form=modelspecs.loc[model_code, "z_equ"],
                main_form=modelspecs.loc[model_code, "main_equ"]
            )
            
            model_future=None
            #two options: statsmodels & tensorzinb. pick one & proceed...
            if modelspecs.loc[model_code, "lib"] == "statsmodels (0)":
                model_future=client.submit(statsmodels_fit,mats_future,method,"UNIFIED")
            elif modelspecs.loc[model_code, "lib"] == "tensorzinb (1)":
                model_future=client.submit(tensorzinb_fit,mats_future,method,"UNIFIED")
            else:
                sprint("[!] Unknown modeling library. Aborting.")
                return -1
            


            #'now' isn't now anymore, but this is easier for pairity/lookups...
            
            dump_ret=client.submit(dump_unified_model,model_future,f"{data_root}/speed_test/models/{model_code}_{now}.pkl")
            wait(dump_ret)
            
        else:#end unified model, begin broken model
            sprint(f'[i] Broken model, parallelizing')
            #broken model. 
            #First, get unique cell-types
            types_future = client.submit(lambda df: df[modelspecs.loc[model_code, 'broken_by']].unique(), dat_future)

            #wait & grab result (small enough for master node)
            types = types_future.result()
            abort_on_failure(types_future,client)
            
            #now scale up the cluster...
            num_workers=min(MAX_PARALLEL,len(types))
            sprint(f"[+] scaling up the cluster to {num_workers} workers")
            cluster.scale(jobs=num_workers)

            sprint(f'[+] Proceeding for one model for each of {types}')

            #timestamp at future creation & execution separately so that
            #we can see if there is substantial slowdown in the for loop
            #or in subsetting
            #which I am suspicious of
            

            mats_futures = {
                t: client.submit(
                    create_matricies,
                    data=client.submit(lambda df, t=t: df[df[modelspecs.loc[model_code, 'broken_by']] == t], dat_future, t),
                    zin_form=modelspecs.loc[model_code, "z_equ"],
                    main_form=modelspecs.loc[model_code, "main_equ"]
                )
                for t in types
            }

            
            model_future = {
                t: client.submit(
                        statsmodels_fit,
                        mats_futures[t],
                        method,
                        t
                    )
                for t in types
            }
            
            
            dump_ret=client.submit(dump_broken_model,model_future,types,f"{data_root}/speed_test/models/{model_code}_{now}.pkl")

            wait(dump_ret)
        

        #end broken & unif modeling

    #end logging performance report (end of with block).


    sprint(f'[+] Done with all tasks. Shutting down')

    client.shutdown()

if __name__=="__main__":
    main()