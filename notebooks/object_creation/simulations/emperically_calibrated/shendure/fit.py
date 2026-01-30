import scMPRAforge as scm
from dask_jobqueue import SLURMCluster
from dask.distributed import Client, LocalCluster
from pathlib import Path

data_root=Path("/nfs/roberts/project/pi_skr2/shared/tabula_data")

local=False

if local:
    cluster=LocalCluster(memory_limit='48G')
    client=Client(cluster)
else:
    cluster=SLURMCluster(
        cores=2,#cores per slurm job
        memory="32G",#memory per slurm job
        processes=1,#dask workers per slurm job,
        job_extra_directives=["-p day", 
            f"--job-name=simclust_worker",
            f"--time=8:00:00",
            f"--output=worker_%j.out"]
    )
    #cluster.scale(jobs=20)
    cluster.scale(jobs=4)
    client = Client(cluster,
            timeout=f"{5*60}s",   # Client <-> scheduler timeout 
            heartbeat_interval="20s"  # Worker heartbeat interval
        )

print(client.dashboard_link,flush=True)

sim=scm.de_novo_simulation(location=data_root,
                            name="twothird_pow_sim_2026-01-30",
                            client=client)

sim.fit_orthos(serial_orthos=True)
sim.save()
print("DONE FIT",flush=True)

sim.precompute_wald(cov_method="sandwich")
sim.save()
print("DONE SANDWICH",flush=True)

sim.precompute_wald(cov_method="opg")
sim.save()
print("DONE OPG",flush=True)


client.close()
cluster.close()

import time
time.sleep(30)