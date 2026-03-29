#!/usr/bin/env python
"""
Merge cohen NB+cm by_cre and by_cell_type halves into a single ortho.
Run after both fit_by_cre.py and fit_by_cell_type.py have completed.
"""
import scMPRAforge as scm
from dask.distributed import Client, LocalCluster
from pathlib import Path

data_root = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data")
path = data_root / "cohen"
name_cre = "cohen_cm_nb_20260329_by_cre"
name_ct  = "cohen_cm_nb_20260329_by_cell_type"
name_out = "cohen_cm_nb_20260329"

cluster = LocalCluster(n_workers=2, threads_per_worker=1, memory_limit="30GB")
client = Client(cluster)

print("[+] Loading by_cre half...", flush=True)
cre = scm.ortho.load(client, path, name_cre)

print("[+] Loading by_cell_type half...", flush=True)
ct = scm.ortho.load(client, path, name_ct)

cre.by_cell_type            = ct.by_cell_type
cre.by_cell_type_design     = ct.by_cell_type_design
cre.by_cell_type_parameters = ct.by_cell_type_parameters

print("[+] Saving merged ortho...", flush=True)
cre.save(path, name_out, client=client)
print("[+] Done.", flush=True)

client.close()
cluster.close()
