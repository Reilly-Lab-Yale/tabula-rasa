#Script takes as its sole parameter a model ID & produces a fit model
#And some raw benchmarking data (to be cooked into a report)

#Note the use of wait(): dask's async is mostly obviated to faciliatate profiling
#(note that this doesn't slow computation: wait is only used when the next step would)
#require the previous step to finish anyway. It's not lazy evaluation if you need the data immediately. 

#general imports
import time
import sys
import pandas as pd
import numpy as np
from formulaic import Formula
from datetime import datetime

#dask imports
from dask_jobqueue import SLURMCluster
from dask.distributed import Client, wait
from dask import delayed

### Functions to be passed as jobs to dask

def statsmodels_fit(X,y,Z,method="bfgs"):
    zinb_model = smdc.ZeroInflatedNegativeBinomialP(y, X, exog_infl=Z)

    n_count_params = zinb_model.exog.shape[1]      # Count model parameters
    n_infl_params = zinb_model.exog_infl.shape[1]    # Inflation model parameters
    n_total = n_count_params + n_infl_params + 1 # adding 1 for alpha
    start_params = np.full(n_total, 0.1)

    zinb_result = zinb_model.fit(start_params=start_params,maxiter=1000,method=method)

    return zinb_result




def load_data(path):
    return pd.read_csv(path,sep="\t")

def create_matricies(main_form,zin_form,data):
    y, X=Formula(main_form).get_model_matrix(data,output='pandas')
    Z=Formula(zin_form).get_model_matrix(data,output='pandas')
    return(X, y, Z)

### Utility functions

def fprint(string):
    """Wraps print with flush=True"""
    print(string,flush=True)

def stamp():
    """returns a formatted string with time since test start"""
    return f'@ {time.time()-start_time:.2f}'

start_time=-1

### Main testing routine

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

    fprint(modelspecs.loc[model_code, 'lib'])

    fprint(f'[+] Creating cluster {stamp()}')

    if "statsmodels" in modelspecs.loc[model_code, 'lib']:
        fprint("[+] Importing statsmodels.")
        import statsmodels.discrete.count_model as smdc
    elif "tensorzinb" in modelspecs.loc[model_code, 'lib']:
        fprint("[+] Importing tensorzinb.")
        from tensorzinb.tensorzinb import TensorZINB
    else:
        fprint("[x] Unknown library spec. Aborting.")
        return -1

    fprint(f'[+] Creating cluster {stamp()}')
    
    now=datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    #Make sure memory per slurm job is large enough to hold the data...
    cluster=SLURMCluster(
            cores=2,#cores per slurm job
            memory="16G",#memory per slurm job
            processes=1,#dask workers per slurm job
            job_extra_directives=["-p ycga", 
                f"--job-name=simclust",
                "--time=08:00:00",
                f"--output={model_code}_{now}_%j.out"]
        )

    client = Client(cluster)

    fprint(f"[!] Cluster started {stamp()}. Monitor on {cluster.dashboard_link}")

    
    #cluster.adapt(minimum_jobs=1, maximum_jobs=10)

    fprint(f'[+] Loading data {stamp()}')
    
    cluster.scale(jobs=1)

    dat_future=None

    if "Fake" in modelspecs.loc[model_code, 'dataset']:
        dat_future=client.submit(load_data,f"{data_root}/simulated/fake_cres.csv")
    elif "Shendure" in modelspecs.loc[model_code, 'dataset']:
        dat_future=client.submit(load_data,f"{data_root}/shendure/shendure_counts_grouped.txt")
    else:
        fprint(f'[x] Unknown dataset. Aborting. {stamp()}')
        return -1


    wait(dat_future)

    fprint(f'[+] Creating matricies {stamp()}')

    mats_future = client.submit(create_matricies,
        data=dat_future,
        zin_form=modelspecs.loc[model_code, "z_equ"],
        main_form=modelspecs.loc[model_code, "main_equ"]
    )
    
    wait(mats_future)

    #fprint(f'[+] Fitting {stamp()}')

    #fprint(f'[+] Gathering statistics {stamp()}')

    fprint(f'[+] Shutting down {stamp()}')

    client.shutdown()

if __name__=="__main__":
    main()