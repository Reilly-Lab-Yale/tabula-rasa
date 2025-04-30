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
from formulaic import Formula

#dask imports
from dask_jobqueue import SLURMCluster
from dask.distributed import Client, wait
from dask.distributed import get_worker
from dask.distributed import performance_report
#from dask import delayed


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
    "c9011010":1,
    "c9200020":1,
    "c020020":4,
    "c0100000":4,
    "c0100010":4,
    "c0200000":4,
    "c0200010":4,
    "c0011000":4,
    "c0012010":4,
    "c0010010":4,
    "c0100020":4,
    "c0011020":4,
    "c0010020":4
}

STATSMODELS_MAXITER=1000

MAX_PARALLEL=10
MAX_A100=3

MIN_PTS=3

### Utility functions

def crash_if_no_gpu():
    import tensorflow as tf
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




#pseudo-logging functions for quick debugging
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

def strip_training_data(model_result):
    """Fix me"""
    model = model_result.model
    model.endog = None
    model.exog = None
    model.exog_infl = None
    return model_result

def toosmall(y):
    if sum(y["umis_mpra_bc"].to_numpy()>0)<MIN_PTS:
        return True
    else:
        return False


### Functions to be passed as jobs to dask

def statsmodels_fit(p,method,name):
    sprint(f"[W] [+] Beginning fitting {name}")
    X,y,Z=p
    if toosmall(y):
        sprint(f"[W] [+] Not enough data for {name}")
        return "too_small"
    
    import statsmodels.discrete.count_model as smdc
    
    zinb_model = smdc.ZeroInflatedNegativeBinomialP(y, X, exog_infl=Z)

    n_count_params = zinb_model.exog.shape[1]      # Count model parameters
    n_infl_params = zinb_model.exog_infl.shape[1]    # Inflation model parameters
    n_total = n_count_params + n_infl_params + 1 # adding 1 for alpha
    start_params = np.full(n_total, 0.1)

    zinb_result = zinb_model.fit(start_params=start_params,maxiter=STATSMODELS_MAXITER,method=method)

    sprint(f"[W] [+] Done fitting {name}.")

    return strip_training_data(zinb_result)

def tensorzinb_fit(p,method,name):
    sprint(f"[W] [+] Beginning tensor fitting {name}")
    X,y,Z=p
    if toosmall(y):
        sprint(f"[W] [+] Not enough data for {name}")
        return "too_small"

    from tensorzinb.tensorzinb import TensorZINB

    zinbo=TensorZINB(y["umis_mpra_bc"].to_numpy().reshape((-1,1)),X,exog_infl=Z.to_numpy())#,same_dispersion=True
    zinb_result=None
    try:
        zinb_result=zinbo.fit(init_method="nb")
    except InvalidIndexError:
        sprint(f"[W] [!] Index error in {name}") 
        return "index_error"

    sprint(f"[W] [+] Done fitting {name}.ls")
    

    return zinb_result


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

## To be used locally

def dump_unified_model(model_future,filename):
    sprint("[+] Called to dump unified model")

    with open(filename, "wb") as f:
        pickle.dump(model_future.result(), f)
        f.flush()
        os.fsync(f.fileno())

def dump_broken_model(model_future,types,filename,sent=False):
    sprint("[+] Called to dump separated model")
    
    materalized_models = {
        t:model_future[t].result()
        for t in types
    }

    if sent:
        materalized_models["multi_sent"]=True

    print(f"[+] Dumping to disc! {stamp()}")
    with open(filename, "wb") as f:
        pickle.dump(materalized_models, f)
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

    hardware_code,hardware=modelspecs.loc[model_code,"hardware"].split("_")

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
    

    client = Client(cluster,
        timeout=f"{5*60}s",   # Client <-> scheduler timeout 
        heartbeat_interval="20s"  # Worker heartbeat interval
    )

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
            
            
            dump_unified_model(model_future,f"{data_root}/speed_test/models/{model_code}_{now}.pkl")
        
            
        else:#end unified model, begin broken model
            sprint(f'[i] Broken model, parallelizing')
            #broken model. 
            #First, get unique cell-types
            types_future = client.submit(lambda df: df[modelspecs.loc[model_code, 'broken_by']].unique(), dat_future)

            #wait & grab result (small enough for master node)
            types = types_future.result()
            abort_on_failure(types_future,client)
            
            #now scale up the cluster...
            #how many nodes we request will depend on particular hardware
            #(we request less of scarcer hardware....)
            
            num_workers=-1

            hardware_choice=modelspecs.loc[model_code,"hardware"]

            #given initial profiling, I think video memory is the limiting factor
            #so I'm not sure if pooling multiple workers per task will really make a difference..?

            if hardware_choice=="0_CPU":
                #CPUs are cheap. Just make as many workers as there are models to fit,
                #up to MAX_PARALLEL
                num_workers=min(MAX_PARALLEL,len(types))
            elif hardware_choice=="1_A100":
                #1_A100 is the "get a job fast for debugging" setting
                num_workers=1
            elif hardware_choice=="2_A100":
                #this is the "let's be serious and allocate a bunch of A100s" mode
                num_workers=min(MAX_A100,len(types))
            elif hardware_choice=="3_A100":
                #this is "try pooling" mode. 
                sprint("[!] Unimplemented hardware choice, aborting.")
                return -1
            else:
                sprint("[!] Unrecognized hardware choice, aborting.")
                return -1

            
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

            #split control flow here on the basis of library
            sent=False
            model_future=None
            if modelspecs.loc[model_code,"lib"]=="statsmodels (0)":
                model_future = {
                    t: client.submit(
                            statsmodels_fit,
                            mats_futures[t],
                            method,
                            t
                        )
                    for t in types
                }
            elif modelspecs.loc[model_code,"lib"]=="tensorzinb (1)":
                sent=True
                #broken model + tensorzinb = requires sentinel value
                model_future = {
                    t: client.submit(
                            tensorzinb_fit,
                            mats_futures[t],
                            method,
                            t
                        )
                    for t in types
                }#p method name
            else:
                sprint("[!] Unrecognized library choice, aborting.")
                return -1

            
            dump_broken_model(model_future,types,f"{data_root}/speed_test/models/{model_code}_{now}.pkl",sent=sent)

        

        #end broken & unif modeling

    #end logging performance report (end of with block).


    sprint(f'[+] Done with all tasks. Shutting down')

    client.shutdown()

if __name__=="__main__":
    main()