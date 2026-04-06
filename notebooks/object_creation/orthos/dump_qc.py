#!/usr/bin/env python
"""
Dump raw QC dicts (by_cell_qc, by_cre_qc) to pickle for each takeshi ortho.
Output: qc/{ortho_name}/by_cell_qc.pkl, qc/{ortho_name}/by_cre_qc.pkl
"""
import sys, pickle
from pathlib import Path
sys.path.insert(0, '/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa-takeshi')

from dask.distributed import Client, LocalCluster
import scMPRAforge.core as scm

DATA_ROOT = "/nfs/roberts/project/pi_skr2/shared/tabula_data_new/takeshi"
QC_DIR = Path(__file__).resolve().parent / "qc"

ORTHOS = [
    "takeshi_obs_nb_phantom",
    "takeshi_obs_zinb_phantom",
    "takeshi_cm_nb_phantom",
    "takeshi_cm_zinb_phantom",
]

cluster = LocalCluster(n_workers=1, threads_per_worker=2, processes=False, memory_limit="60GB")
client = Client(cluster)
print(f"Dashboard: {client.dashboard_link}", flush=True)

for name in ORTHOS:
    print(f"\n=== {name} ===", flush=True)
    ortho_obj = scm.ortho.load(client, DATA_ROOT, name)
    ortho_obj.compute_model_qc()

    out_dir = QC_DIR / name
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(out_dir / "by_cell_qc.pkl", "wb") as f:
        pickle.dump(ortho_obj.by_cell_qc, f)
    with open(out_dir / "by_cre_qc.pkl", "wb") as f:
        pickle.dump(ortho_obj.by_cre_qc, f)

    print(f"  by_cell_qc: {len(ortho_obj.by_cell_qc)} entries", flush=True)
    print(f"  by_cre_qc:  {len(ortho_obj.by_cre_qc)} entries", flush=True)
    print(f"  Saved to {out_dir}", flush=True)

client.close()
cluster.close()
print("\nDone.", flush=True)
