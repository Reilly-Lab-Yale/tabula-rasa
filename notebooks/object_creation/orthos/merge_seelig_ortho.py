#!/usr/bin/env python
"""
Merges seelig_ortho_20260320_by_cre + seelig_ortho_20260320_by_cell_type
into the final seelig_ortho_20260320.

Run after both half-jobs complete:
    ipython merge_seelig_ortho.py
"""
import scMPRAforge as scm
from dask_jobqueue import SLURMCluster
from dask.distributed import Client
from pathlib import Path

data_root = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data")
path = data_root / "seelig"
name_cre = "seelig_ortho_20260320_by_cre"
name_ct  = "seelig_ortho_20260320_by_cell_type"
name_out = "seelig_ortho_20260320"

cluster = SLURMCluster(
    cores=1, memory="64G", processes=1,
    job_extra_directives=["-p priority", "-A prio_skr2", "--time=4:00:00"]
)
cluster.scale(jobs=4)
client = Client(cluster)
client.wait_for_workers(n_workers=4, timeout=300)

cre = scm.ortho.load(client, path, name_cre)
ct  = scm.ortho.load(client, path, name_ct)

# Combine: take by_cre side from cre, by_cell_type side from ct
cre.by_cell_type          = ct.by_cell_type
cre.by_cell_type_design   = ct.by_cell_type_design
cre.by_cell_type_parameters = ct.by_cell_type_parameters

cre.save(path, name_out, client=client)
print(f"[+] Saved merged ortho to {path / name_out}")

client.close()
cluster.close()
