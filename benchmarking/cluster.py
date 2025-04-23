#Script takes as its a model ID & produces a fit model
#and some raw benchmarking data

#Note the use of wait(): dask's async is mostly obviated
#note that this doesn't slow computation: wait is only used when the next step would
#require the previous step to finish anyway. It's not lazy evaluation if you need the data immediately. 
#and it's all ufunc so there's no computation graph optimization

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
        "c900090":1,
        "c900010":1,
        "c000000":12,
        "c000010":12,
    }

STATSMODELS_MAXITER=5000

OPTIMIZERS={
    'c000010':'cg'
}

### Utility functions

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
        fprint("[!] Computation failed. Aborting.")
        fprint("[!] Check slave node logs for details.")
        client.shutdown()
        assert 1==2


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
        fprint("[!] Wrong number of arguments. Aborting.")
        return -1
    
    #useful paths
    

    model_code=sys.argv[1]
    modelspecs=pd.read_csv(f"{data_root}/speed_test/modelspecs.tsv",sep="\t",index_col=0)




    fprint(f'[+] Creating cluster {stamp()}')
    
    now=datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


    

    #Make sure memory per slurm job is large enough to hold the data...
    cluster=SLURMCluster(
            cores=2,#cores per slurm job
            memory="16G",#memory per slurm job
            processes=1,#dask workers per slurm job
            job_extra_directives=["-p ycga", 
                f"--job-name=simclust_worker",
                f"--time={WORKER_TIMES[model_code]}:00:00",
                f"--output=slave_{model_code}_{now}_%j.out"]
        )
    MAX_PARALLEL=10

    client = Client(cluster)

    fprint(f"[+] Cluster started {stamp()}. Monitor on {cluster.dashboard_link}")

    with performance_report(filename=f"{model_code}_{now}_report.html"):
        
        #cluster.adapt(minimum_jobs=1, maximum_jobs=10)

        fprint(f'[+] Loading data {stamp()}')
        
        cluster.scale(jobs=1)

        dat_future=None

        if "Fake" in modelspecs.loc[model_code, 'dataset']:
            dat_future=client.submit(load_csv,f"{data_root}/simulated/fake_cres.csv")
        elif "Shendure" in modelspecs.loc[model_code, 'dataset']:
            dat_future=client.submit(load_tsv,f"{data_root}/shendure/shendure_counts_grouped.txt")
        else:
            fprint(f'[!] Unknown dataset. Aborting. {stamp()}')
            return -1


        wait(dat_future)
        abort_on_failure(dat_future,client)

        model_future=None

        #for a unified model, we just fit on one node. For a broken model, we will run in parallel.
        if pd.isna(modelspecs.loc[model_code, 'broken_by']):
            fprint(f'[i] Unified model : executing in one job.')
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
            #move this to a function
            with open(f"{data_root}/speed_test/models/{model_code}_{now}.pkl", "wb") as f:
                pickle.dump(model_future.result(), f)
            
        else:
            fprint(f'[i] Broken model, parallelizing {stamp()}')
            #broken model. 
            #First, get unique cell-types
            types_future = client.submit(lambda df: df[modelspecs.loc[model_code, 'broken_by']].unique(), dat_future)

            #wait & grab result (small enough for master node)
            types = types_future.result()
            abort_on_failure(types_future,client)
            
            #now scale up the cluster...
            num_workers=min(MAX_PARALLEL,len(types))
            fprint(f"[+] scaling up the cluster to {num_workers} workers")
            cluster.scale(jobs=num_workers)

            fprint(f'[+] Proceeding for one model for each of {types} {stamp()}')

            #timestamp at future creation & execution separately so that
            #we can see if there is substantial slowdown in the for loop
            #or in subsetting
            #which I am suspicious of
            
            fprint(f'[+] Creating matricies : creating futures {stamp()}')

            mats_futures = {
                t: client.submit(
                    create_matricies,
                    data=client.submit(lambda df, t=t: df[df[modelspecs.loc[model_code, 'broken_by']] == t], dat_future, t),
                    zin_form=modelspecs.loc[model_code, "z_equ"],
                    main_form=modelspecs.loc[model_code, "main_equ"]
                )
                for t in types
            }

            

            fprint(f'[+] Creating matricies : execution {stamp()}')
            wait(list(mats_futures.values()))

            fprint(f'[+] Fitting : creating futures {stamp()}')
            model_future = {
                t: client.submit(
                        statsmodels_fit,
                        mats_futures[t]
                    )
                for t in types
            }
            
            fprint(f'[+] Fitting : execution {stamp()}')
            wait(list(model_future.values()))

            fprint(f'[+] Dumping model {stamp()}')
            
            materalized_models = {
                t:model_future[t].result()
                for t in types
            }
            
            with open(f"{data_root}/speed_test/models/{model_code}_{now}.pkl", "wb") as f:
                pickle.dump(materalized_models, f)

        #end broken & unif modeling
        
            

    #end performance report (end of with block).


    fprint(f'[+] Shutting down {stamp()}')

    client.shutdown()

if __name__=="__main__":
    main()