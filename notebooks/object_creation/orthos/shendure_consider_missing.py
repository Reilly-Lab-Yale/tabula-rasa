### imports ###
import pandas as pd
import numpy as np
import time
import pickle
from pathlib import Path
from formulaic import Formula
import seaborn as sns
import matplotlib.pyplot as plt
import os
import signal
from tensorzinb.tensorzinb import TensorZINB
import scMPRAforge as scm

# Allow graceful shutdown via SIGTERM (e.g. scancel) so finally/context managers run
signal.signal(signal.SIGTERM, lambda s, f: (_ for _ in ()).throw(KeyboardInterrupt()))

### create dask cluster ###

from dask_jobqueue import SLURMCluster
from dask.distributed import Client, performance_report

local=False
if not local:
    worker_cores = 1
    worker_memory = "128G"
    worker_processes = 1
    cluster_jobs = 16

    cluster=SLURMCluster(
        cores=worker_cores,#cores per slurm job
        memory=worker_memory,#memory per slurm job
        processes=worker_processes,#dask workers per slurm job
        job_extra_directives=["-p ycga", 
            f"--job-name=simclust_worker",
            f"--time=36:00:00",
            f"--output=worker_%j.out"]
    )

    cluster.scale(jobs=cluster_jobs)

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
report_path = Path(__file__).with_name(f"{name}_dask_performance_report.html")

print("! fitting",flush=True)
print(f"[+] Writing Dask performance report to {report_path}", flush=True)
### fit ###

try:
    with performance_report(filename=str(report_path)):
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
finally:
    print("! Done, shutting down",flush=True)
    time.sleep(10)
    client.close()
    cluster.close()
    time.sleep(10)
