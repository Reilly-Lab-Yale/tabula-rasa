### patch_design_matrices.py
#
# One-time script to regenerate and save design matrices for a partially-saved
# ortho (e.g. after a driver OOM that killed the process during save()).
#
# The 4 model/parameter pkl files must already exist in the save directory.
# Design matrices are recomputed from the original data and saved as sharded
# directories (one pkl per level) without pulling anything to the driver.
###

import time
import signal
from pathlib import Path
import scMPRAforge as scm

# Allow graceful shutdown via SIGTERM (e.g. scancel)
signal.signal(signal.SIGTERM, lambda s, f: (_ for _ in ()).throw(KeyboardInterrupt()))

### create dask cluster ###

from dask_jobqueue import SLURMCluster
from dask.distributed import Client, performance_report

local=False
if not local:
    worker_cores = 1
    worker_memory = "96G"
    worker_processes = 1
    cluster_jobs = 16

    cluster=SLURMCluster(
        cores=worker_cores,
        memory=worker_memory,
        processes=worker_processes,
        job_extra_directives=["-p ycga",
            f"--job-name=simclust_worker",
            f"--time=2-23:40:00",
            f"--output=worker_%j.out"]
    )

    cluster.scale(jobs=cluster_jobs)

    client = Client(cluster,
            timeout=f"{5*60}s",
            heartbeat_interval="20s"
        )
else:
    from dask.distributed import Client, LocalCluster
    cluster=LocalCluster()
    client = Client(cluster)

print(client.dashboard_link, flush=True)

data_root="/vast/palmer/pi/reilly/tabula_data"
path=f"{data_root}/shendure"
name="shendure_ortho_consider_missing_20260320"

print("! Loading data", flush=True)
try:
    shendure=scm.scMPRA_data.from_tsv(f"{path}/shendure_processed.tsv")
    shendure.set_negative_controls(["minP","noP"])
    shendure.set_reference_cell("Pluripotent")
    shendure.ortho_filter()
    shendure.set_consider_missing(True)

    print("! Loading partial ortho", flush=True)
    primordial=scm.ortho.load(client, path, name)
    print(f"  by_cre_design:       {primordial.by_cre_design}", flush=True)
    print(f"  by_cell_type_design: {primordial.by_cell_type_design}", flush=True)

    print("! Recomputing design matrices", flush=True)
    primordial.recompute_design_matrices(client, shendure)

    print("! Saving design matrices", flush=True)
    full_path = Path(path) / name
    from scMPRAforge.core import _save_design_dict
    _save_design_dict(client, primordial.by_cre_design,       full_path / "by_cre_design")
    _save_design_dict(client, primordial.by_cell_type_design, full_path / "by_cell_type_design")

    print("! Done", flush=True)
finally:
    print("! Shutting down", flush=True)
    time.sleep(10)
    client.close()
    cluster.close()
    time.sleep(10)
