### imports ###
import pandas as pd
import numpy as np
import time
import pickle
from formulaic import Formula
import seaborn as sns
import matplotlib.pyplot as plt
import os
from tensorzinb.tensorzinb import TensorZINB
import scMPRAforge as scm


### create dask cluster ###

from dask_jobqueue import SLURMCluster
from dask.distributed import Client

local=False
if not local:
    cluster=SLURMCluster(
        cores=4,#cores per slurm job
        memory="128G",#memory per slurm job
        processes=4,#dask workers per slurm job
        job_extra_directives=["-p ycga", 
            f"--job-name=simclust_worker",
            f"--time=36:00:00",
            f"--output=worker_%j.out"]
    )

    cluster.scale(jobs=4)

    client = Client(cluster,
            timeout=f"{5*60}s",   # Client <-> scheduler timeout 
            heartbeat_interval="20s"  # Worker heartbeat interval
        )
else:
    from dask.distributed import Client, LocalCluster
    cluster=LocalCluster()
    client = Client(cluster)

print(client.dashboard_link,flush=True)

#data_root="/nfs/roberts/project/pi_skr2/shared/tabula_data"
data_root="/vast/palmer/pi/reilly/tabula_data"
path=f"{data_root}/shendure"
name="shendure_ortho_consider_missing_20260310"

print("! fitting",flush=True)
### fit ###

if os.path.isdir(path+"/"+name):
    print("[+] Model found. Loading...")
    primordial=scm.ortho.load(client,path,name)
    shendure=primordial.training_data
else:
    print("[+] Model not found. Creating...")

    #load data
    shendure=scm.scMPRA_data.from_tsv(f"{path}/shendure_processed.tsv")
    shendure.set_negative_controls(["minP","noP"])
    shendure.set_reference_cell("Pluripotent")
    shendure.ortho_filter()
    
    shendure.set_consider_missing(True)

    primordial=scm.ortho()
    primordial.criss_cross(client=client,
                       dat=shendure)
    primordial.extract_params(client)
    primordial.save(path,name)

print("! Done, shutting down",flush=True)

import time
time.sleep(10)

client.close()
cluster.close()

time.sleep(10)
