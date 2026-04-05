"""
Smoke test: precompute_wald on a real phantom ortho.
Run via sbatch (needs memory for Dask graph serialization).

Usage:
    sbatch --mem=64G -c 4 -t 1:00:00 -p priority -A prio_skr2 \
        --wrap="module load miniconda; conda activate tz; python3 scMPRAforge/tests/smoke_test_wald.py"
"""
import sys
sys.path.insert(0, '/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa')

from dask.distributed import Client, LocalCluster
import scMPRAforge.core as scm
import numpy as np

DATA_DIR = "/nfs/roberts/project/pi_skr2/shared/tabula_data_new/cohen"
ORTHO = "cohen_obs_nb_phantom_20260401"

cluster = LocalCluster(n_workers=1, threads_per_worker=2, processes=False, memory_limit="50GB")
client = Client(cluster)
print(f"Dashboard: {client.dashboard_link}", flush=True)

print(f"Loading {ORTHO}...", flush=True)
ortho = scm.ortho.load(client, DATA_DIR, ORTHO)
print(f"  by_cre: {len(ortho.by_cre.model)} models", flush=True)
print(f"  by_cell_type: {len(ortho.by_cell_type.model)} models", flush=True)

print("Running precompute_wald (sandwich)...", flush=True)
ortho.precompute_wald(client, cov_method="sandwich")
print("precompute_wald complete.", flush=True)

wp = ortho.wald_precomp.flattened_copy()

n_ok_ct = 0
n_fail_ct = 0
for ct, entry in wp.by_cell_type.items():
    ok = np.all(np.isfinite(entry.se_x_mu)) and np.all(entry.se_x_mu > 0)
    if ok:
        n_ok_ct += 1
    else:
        n_fail_ct += 1
        print(f"  ISSUE by_ct {ct}: debug={entry.debug_msg}", flush=True)

n_ok_cr = 0
n_fail_cr = 0
for cr, entry in wp.by_cre.items():
    ok = np.all(np.isfinite(entry.se_x_mu)) and np.all(entry.se_x_mu > 0)
    if ok:
        n_ok_cr += 1
    else:
        n_fail_cr += 1

print(f"\nResults:", flush=True)
print(f"  by_cell_type: {n_ok_ct} ok, {n_fail_ct} failed (of {n_ok_ct+n_fail_ct})", flush=True)
print(f"  by_cre: {n_ok_cr} ok, {n_fail_cr} failed (of {n_ok_cr+n_fail_cr})", flush=True)

# Print a few examples
for ct, entry in list(wp.by_cell_type.items())[:2]:
    print(f"\n  by_ct '{ct}':", flush=True)
    print(f"    se_x_mu = {entry.se_x_mu}", flush=True)
    print(f"    debug = {entry.debug_msg[:80]}", flush=True)

client.close()
cluster.close()
print("\nSmoke test done.", flush=True)
