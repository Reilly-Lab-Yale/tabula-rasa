#!/usr/bin/env python3
"""
MWU retest: walk all design-space sims and run MWU on each.

This script loads each sim that already has ttest results, runs
sim.mwu("hs_act") to produce MWU results alongside the existing ttest
results, then saves. Sims missing scmpra count data are skipped.

Usage:
    python mwu_retest.py           # run MWU retest on all axes
    python mwu_retest.py bcs_per_cre  # single axis
"""

import sys
import os
import gc
from pathlib import Path

_repo_root = Path(__file__).resolve().parents[3]
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

import scMPRAforge as scm
from dask.distributed import Client, LocalCluster
from dask_jobqueue import SLURMCluster

SIM_ROOT = Path("/nfs/roberts/scratch/pi_skr2/mcn26/design_space_sims/2026-04-11_design_space")
N_SIMS = 5
AXES = ["bcs_per_cre", "n_cells", "moi", "n_cres"]


def retest_axis(axis: str, client):
    axis_dir = SIM_ROOT / axis
    if not axis_dir.exists():
        print(f"  {axis}: no sim dir, skipping", flush=True)
        return

    done = 0
    skipped = 0
    already = 0
    for pt_dir in sorted(axis_dir.iterdir()):
        if not pt_dir.is_dir():
            continue
        print(f"\n  {pt_dir.name}", flush=True)
        for d in sorted(pt_dir.iterdir()):
            if not d.is_dir():
                continue
            # Skip sims that already have MWU results
            mwu_dir = d / "tests" / "hs_act" / "mwu"
            if mwu_dir.exists() and any(mwu_dir.iterdir()):
                already += 1
                continue
            # Must have ttest results (confirms hypothesis set exists)
            tt_dir = d / "tests" / "hs_act" / "ttest"
            if not (tt_dir.exists() and any(tt_dir.iterdir())):
                skipped += 1
                continue
            # Must have all essential files
            required = [d / "ground_truth.tsv.gz", d / "state.parquet"]
            required += [tt_dir / f"{i}_results.tsv" for i in range(N_SIMS)]
            for i in range(N_SIMS):
                required.append(d / f"simulated_scmpra/{i}.scmpra/data.parquet")
            if not all(p.exists() for p in required):
                skipped += 1
                continue
            sim = None
            try:
                sim = scm.de_novo_simulation(
                    location=pt_dir, name=d.name, client=client
                )
                sim.mwu("hs_act")
                sim.save()
                done += 1
                if done % 10 == 0:
                    print(f"    {done} retested...", flush=True)
            except Exception as e:
                print(f"    skip {d.name}: {type(e).__name__}: {e}", flush=True)
                skipped += 1
            finally:
                del sim
                gc.collect()

    print(f"\n  {axis}: {done} retested, {already} already had MWU, {skipped} skipped", flush=True)


def _make_slurm_client():
    cluster = SLURMCluster(
        cores=1,
        memory="24G",
        processes=1,
        env_extra=[f"export PYTHONPATH={_repo_root}:$PYTHONPATH"],
        job_extra_directives=[
            "-p priority",
            "--account=prio_skr2",
            "--job-name=mwu_worker",
            "--time=8:00:00",
            "--output=worker_%j.out",
        ],
    )
    cluster.scale(jobs=20)
    client = Client(cluster, timeout="300s", heartbeat_interval="20s")
    print(f"Dask dashboard: {client.dashboard_link}", flush=True)
    return cluster, client


if __name__ == "__main__":
    axes = [sys.argv[1]] if len(sys.argv) > 1 else AXES
    for a in axes:
        if a not in AXES:
            print(f"Unknown axis: {a!r}. Choose from {AXES}")
            sys.exit(1)

    cluster, client = _make_slurm_client()
    try:
        for axis in axes:
            print(f"\n{'='*60}\nMWU retest: {axis}\n{'='*60}", flush=True)
            retest_axis(axis, client)
    finally:
        client.close()
        cluster.close()

    print("\nDone.", flush=True)
