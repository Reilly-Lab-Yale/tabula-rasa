import scMPRAforge as scm
from dask_jobqueue import SLURMCluster
from dask.distributed import Client, LocalCluster
from pathlib import Path

cluster=LocalCluster(worker_extra_args=["--resources", "FIT=1"])
client=Client(cluster)

#We want to randomly assign one of the workers to 
#take the high memory tasks, so that no more than one such task is
#executing at a time.
workers = list(client.scheduler_info()["workers"])

# choose first worker
target = workers[0]

client.run(
    lambda dask_worker: dask_worker.resources.update({"CELL_DESIGN": 1}),
    workers=[target],
)

#    cluster=SLURMCluster(
#        cores=3,#cores per slurm job
#        memory="512G",#memory per slurm job
#        processes=3,#dask workers per slurm job,
#        worker_extra_args=["--resources", "FIT=1"],
#        job_extra_directives=["-p week", 
#            f"--job-name=simclust_worker",
#            f"--time=72:00:00",
#            f"--output=worker_%j.out"]
#    )
#    cluster.scale(jobs=2)
client = Client(cluster,
        timeout=f"{5*60}s",   # Client <-> scheduler timeout 
        heartbeat_interval="20s"  # Worker heartbeat interval
    )

data_root=Path("/nfs/roberts/project/pi_skr2/shared/tabula_data")

print(client.dashboard_link,flush=True)

sim=scm.de_novo_simulation(location=data_root,
                            name="twothird_pow_sim_2026-02-05",
                            client=client)


print("LOADED",flush=True)

sim.fit_orthos(serial_orthos=False,direction="by_cell_type")

print("SUBMITTED",flush=True)

sim.save()

print("DONE",flush=True)

import time
time.sleep(60)

client.close()
cluster.close()

time.sleep(60)