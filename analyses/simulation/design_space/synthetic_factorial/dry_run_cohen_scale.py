#!/usr/bin/env python3
"""Dry-run: simulate one cohen-scale sample to measure wall time + memory
before committing to a bounds-expanded top-up sweep.

Cohen empirical: n_cells=18633, n_cres=116, bcs_per_cre=17244,
                 moi=149, alpha=1.39, minP=0.94.

This is way outside the original LHS box (bcs_per_cre 34x past LHS max,
moi 5x past). We expect this single sample to take much longer than a
typical original-sweep sample. Record the timing and use it to size the
top-up sweep.
"""

import sys, time, gc
from pathlib import Path
import numpy as np
import pandas as pd

_repo_root = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(_repo_root))

import scMPRAforge as scm
from dask.distributed import Client
from dask_jobqueue import SLURMCluster

OUT_ROOT = Path("/nfs/roberts/scratch/pi_skr2/mcn26/synthetic_factorial_sims/dry_run_cohen_scale")
OUT_ROOT.mkdir(parents=True, exist_ok=True)

# Cohen-scale single point.
COHEN = dict(
    n_cells=18633,
    n_cres=116,
    bcs_per_cre=17244.5,
    moi=149.15,
    lib_alpha_nb=1.3873,
    minP=0.9363,
    activity_max_mult=4.0,
)

CELL_TYPE = "reference"
N_SIMS = 5  # match full-run setting
N_REPS = 1  # just one rep for the timing estimate; multiplies cleanly


def make_bounds(p):
    b = scm.SHENDURE_BOUNDS.copy()
    if CELL_TYPE not in b.cells_per_cell_type.index:
        ct = b.cells_per_cell_type.index[0]
        b.cells_per_cell_type = b.cells_per_cell_type.loc[[ct]].copy()
        b.cells_per_cell_type.index = [CELL_TYPE]
    else:
        b.cells_per_cell_type = b.cells_per_cell_type.loc[[CELL_TYPE]].copy()
    b.cells_per_cell_type.loc[CELL_TYPE] = int(p["n_cells"])
    b.library_model.mu_nb = float(p["bcs_per_cre"])
    b.library_model.alpha_nb = float(p["lib_alpha_nb"])
    b.library_model.update_alt_nb_param()
    b.set_effective_moi(float(p["moi"]))
    b.reference_activity = float(p["minP"])
    b.min_mpra_umi = 0.05 * float(p["minP"])
    return b


def main():
    cluster = SLURMCluster(
        cores=1,
        memory="48G",
        processes=1,
        job_script_prologue=[f"export PYTHONPATH={_repo_root}:$PYTHONPATH"],
        job_extra_directives=[
            "-p priority",
            "--account=prio_skr2",
            "--job-name=dryrun_worker",
            "--time=8:00:00",
            "--output=worker_%j.out",
        ],
    )
    cluster.adapt(minimum=1, maximum=50, interval="30s")
    client = Client(cluster, timeout="600s", heartbeat_interval="20s")
    print(f"dashboard: {client.dashboard_link}", flush=True)

    bound = make_bounds(COHEN)
    minP = COHEN["minP"]
    max_act = COHEN["activity_max_mult"] * minP
    min_act = 0.05 * minP

    overall = time.time()
    rep_dir = OUT_ROOT
    rep_dir.mkdir(exist_ok=True)
    print(f"start: {time.strftime('%H:%M:%S')}", flush=True)

    t0 = time.time()
    print("[T+0] one_library_replicate ...", flush=True)
    _, sim = scm.one_library_replicate(
        root=rep_dir, n_sims=N_SIMS, client=client,
        flatten_overtransfection=True, bound=bound,
        n_cres=COHEN["n_cres"], min=min_act, max=max_act,
        minP=minP, cell_type=CELL_TYPE,
    )
    t_olr = time.time() - t0
    print(f"  one_library_replicate: {t_olr:.1f}s", flush=True)

    t1 = time.time()
    print("[T+%.0f] sim.save ..." % (time.time() - overall), flush=True)
    sim.save()
    t_save1 = time.time() - t1
    print(f"  save: {t_save1:.1f}s", flush=True)

    t2 = time.time()
    print("[T+%.0f] make_all_by_celltype_hypotheses ..." % (time.time() - overall), flush=True)
    example = scm.scMPRA_data.from_parquet(sim.scmpradatp / "0.scmpra")
    hs = scm.make_all_by_celltype_hypotheses(counts=example, reference_cre="reference")
    t_hs = time.time() - t2
    print(f"  hypothesis set ({len(hs.df)} hypotheses): {t_hs:.1f}s", flush=True)

    t3 = time.time()
    print("[T+%.0f] add_hypothesis_set + mwu ..." % (time.time() - overall), flush=True)
    sim.add_hypothesis_set("hs_act", hs)
    sim.mwu("hs_act")
    t_mwu_submit = time.time() - t3
    print(f"  mwu submitted: {t_mwu_submit:.1f}s", flush=True)

    t4 = time.time()
    print("[T+%.0f] sim.save (final, blocks on mwu) ..." % (time.time() - overall), flush=True)
    sim.save()
    t_mwu_total = time.time() - t4
    print(f"  mwu total wall (incl block): {t_mwu_total:.1f}s", flush=True)

    total = time.time() - overall
    print(f"\n=== TOTAL one cohen-scale rep: {total:.1f}s ({total/60:.1f} min) ===", flush=True)
    print(f"  one_library_replicate:  {t_olr:.1f}s ({t_olr/total*100:.0f}%)")
    print(f"  pre-mwu save:           {t_save1:.1f}s ({t_save1/total*100:.0f}%)")
    print(f"  hypothesis_set build:   {t_hs:.1f}s ({t_hs/total*100:.0f}%)")
    print(f"  mwu (submit + block):   {t_mwu_submit + t_mwu_total:.1f}s ({(t_mwu_submit+t_mwu_total)/total*100:.0f}%)")

    # Disk usage
    import subprocess
    du = subprocess.check_output(["du", "-sh", str(OUT_ROOT)]).decode().split()[0]
    print(f"\ndisk used: {du}", flush=True)

    client.close()
    cluster.close()


if __name__ == "__main__":
    main()
