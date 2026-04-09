#!/usr/bin/env python3
"""
Cohen t-test power -- all cell types

Assesses Welch's t-test power (TPR) as a function of CRE activity relative to
`minP`, for every cell type in the Cohen dataset.

Cohen is an **episomal** MPRA: all CREs are transfected into all cell types
(no integrating bottleneck). Therefore we run **one set of simulations** with
all cell types included, rather than separate per-cell-type runs.
Each simulation draws the same synthetic CREs into all cell types with
identical activities (`n_cres` confirmed uniform across cell types in the
empirical data).

Steps:
1. Draw `n_cres` synthetic CREs from uniform(min_activity, max_activity),
   one fixed at `minP` as reference
2. Replicate the same CREs across all cell types in the ground truth
3. Simulate `n_sims` replicates x `n_library_reps` independent libraries
   using empirical Cohen bounds
4. Run Welch's t-test for all CREs vs. reference -- both with reporter
   (standard) and without reporter (zero-deflated)
5. Aggregate and plot per-cell-type power curves, overlaying both conditions

Usage:
    python cohen_power_ttest_all_cell_types.py [phase]

Phases:
    simulate  -- generate simulations and run t-test (both conditions)
    plot      -- aggregate results and produce SVG
    all       -- run both phases sequentially (default)
"""

import sys
import os
import math
import uuid
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# Ensure the repo root is on PYTHONPATH so scMPRAforge is importable
# when running from the script's own directory.
_repo_root = Path(__file__).resolve().parents[4]
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

import scMPRAforge as scm
from dask.distributed import Client, LocalCluster
from dask_jobqueue import SLURMCluster

# ---------------------------------------------------------------------------
# Paths and parameters
# ---------------------------------------------------------------------------

DATA_ROOT = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new")
OUTPUT_DIR = Path(__file__).resolve().parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

SIM_DATE = "2026-04-08"
SIM_DIR = DATA_ROOT / "simulated" / f"{SIM_DATE}_cohen_pow"

# Cohen is episomal: all CREs present in all cell types. Verified uniform
# n_cres=115 across cell types (Bipolar, Interneuron, Mueller Glia, Rod)
# from retina_single_counting_u6.scmpra. Hardcoded to avoid loading the
# 7.3 GB training data parquet at import time.
N_CRES = 115

CELL_TYPES = list(scm.COHEN_BOUNDS.cells_per_cell_type.index)

N_LIBRARY_REPS = 100
N_SIMS = 5
MINP = scm.COHEN_BOUNDS.reference_activity

# Activity range chosen to capture full power curve on both sides of FC=1.
# Left side (FC<1) reaches ~0.8 power around FC~0.75; right side mirrors at
# FC~1.25. Range 0.5*minP--1.6*minP gives comfortable margin to see full
# S-curve.
# min_activity = scm.COHEN_BOUNDS.min_mpra_umi  # old: too wide, squashes
#                                                  transition zone
MIN_ACTIVITY = MINP * 0.5
MAX_ACTIVITY = MINP * 1.6

BATCH_SIZE = 10


# ---------------------------------------------------------------------------
# Phase: simulate -- generate data and run MWU (both conditions)
# ---------------------------------------------------------------------------

def phase_simulate(client):
    print(f"minP={MINP:.6f}, min={MIN_ACTIVITY:.6f}, max={MAX_ACTIVITY:.6f}")
    print(f"n_cres={N_CRES}, cell_types={CELL_TYPES}")
    print(f"library reps: {N_LIBRARY_REPS}, sims per rep: {N_SIMS}")
    print()

    for ct in CELL_TYPES:
        cells = scm.COHEN_BOUNDS.cells_per_cell_type.get(ct, "N/A")
        print(f"  {ct}: cells={cells}")

    # --- Helper: process a batch of sims (gamut -> save -> mwu -> save) ---
    # Cohen is episomal: same CREs transfected into all cell types.
    #
    # Each batch is fully completed (gamut -> save -> mwu -> save) before the
    # next starts. Sim objects are discarded after each batch; only
    # (location, name) tuples are retained. This keeps main-process memory
    # flat regardless of how many batches have run.

    def _process_batch(batch, hs_all_ct):
        """Block until batch is computed, run both MWU conditions, save."""
        for s in batch:
            s.save()

        if hs_all_ct is None:
            example_data = scm.scMPRA_data.from_parquet(
                batch[0].scmpradatp / "0.scmpra"
            )
            hs_all_ct = scm.make_all_by_celltype_hypotheses(
                counts=example_data, reference_cre="reference"
            )

        for s in batch:
            s.add_hypothesis_set("hs_all_ct", hs_all_ct)

            # +reporter (standard)
            s.ttest("hs_all_ct")

            # -reporter (zero-deflated): drops all zero-count observations
            # before testing. This simulates what a real no-reporter experiment
            # looks like on disc: only non-zero MPRA signal is observed.
            s.ttest("hs_all_ct", has_reporter=False)

        for s in batch:
            s.save()

        return hs_all_ct, [(s.location, s.name) for s in batch]

    # --- Resume: collect already-completed sims from disk ---
    # A sim is "complete" if it has both mwu/ and mwu_deflated/ results.
    sim_paths = []
    if SIM_DIR.exists():
        for d in sorted(SIM_DIR.iterdir()):
            if not d.is_dir():
                continue
            tt_dir = d / "tests" / "hs_all_ct" / "ttest"
            defl_dir = d / "tests" / "hs_all_ct" / "ttest_deflated"
            if (tt_dir.exists() and any(tt_dir.iterdir())
                    and defl_dir.exists() and any(defl_dir.iterdir())):
                sim_paths.append((SIM_DIR, d.name))
    print(f"Resuming: found {len(sim_paths)} complete sims on disk.")

    n_remaining = N_LIBRARY_REPS - len(sim_paths)
    print(f"Will run {n_remaining} fresh sims to reach {N_LIBRARY_REPS} total.")

    hs_all_ct = None
    batch = []

    for i in range(n_remaining):
        rng = np.random.default_rng()

        cre_gt = rng.uniform(MIN_ACTIVITY, MAX_ACTIVITY, size=N_CRES - 1)
        cre_gt = np.append(cre_gt, MINP)
        names = [f"synthcre_{j}" for j in range(N_CRES - 1)] + ["reference"]

        gt_df = pd.concat([
            pd.DataFrame({"cre_id": names, "mu": cre_gt, "cell_type": ct})
            for ct in CELL_TYPES
        ], ignore_index=True)

        libraries = [
            scm.simulate_library(
                CREs=pd.Series(names),
                library_model=scm.COHEN_BOUNDS.library_model,
            )
            for _ in range(N_SIMS)
        ]

        sim = scm.de_novo_simulation(
            location=SIM_DIR,
            name=f"sim_{uuid.uuid4().hex[:8]}",
            client=client,
            libraries=libraries,
            library_mapping="corresponding",
            flatten_overtransfection=True,
            n_sims=N_SIMS,
            experiment_bounds=scm.COHEN_BOUNDS,
            ground_truth=gt_df,
        )
        sim.gamut()
        batch.append(sim)

        if len(batch) == BATCH_SIZE:
            hs_all_ct, paths = _process_batch(batch, hs_all_ct)
            sim_paths.extend(paths)
            batch = []
            print(
                f"Batch done -- {len(sim_paths)}/{N_LIBRARY_REPS} complete.",
                flush=True,
            )

    if batch:
        hs_all_ct, paths = _process_batch(batch, hs_all_ct)
        sim_paths.extend(paths)
        print(
            f"Batch done -- {len(sim_paths)}/{N_LIBRARY_REPS} complete.",
            flush=True,
        )

    print(f"All {len(sim_paths)} sims done and tested.", flush=True)
    return sim_paths


# ---------------------------------------------------------------------------
# Phase: plot -- aggregate and produce SVG
# ---------------------------------------------------------------------------

def phase_plot(client):
    # Reload sim handles from disk.
    sim_paths = []
    if SIM_DIR.exists():
        for d in sorted(SIM_DIR.iterdir()):
            if not d.is_dir():
                continue
            tt_dir = d / "tests" / "hs_all_ct" / "ttest"
            if tt_dir.exists() and any(tt_dir.iterdir()):
                sim_paths.append((SIM_DIR, d.name))
    print(f"Loading {len(sim_paths)} sims from disk.", flush=True)

    sims = [
        scm.de_novo_simulation(location=loc, name=name, client=client)
        for loc, name in sim_paths
    ]

    # Aggregate results, preserving cell type (sum_pow drops it, so we roll
    # our own).
    def sum_pow_by_ct(sims, hypothesis_set_name, test_type):
        """Like sum_pow but returns dict {cell_type -> DataFrame}."""
        rows = []
        for sim in sims:
            for i in range(sim.get_state_field("n_sims")):
                m = sim._merge_in_ground_truth(
                    hypothesis_set_name, test_type, i
                )
                m["fc"] = m["comparison_truth"] / m["reference_truth"]
                rows.append(
                    m[["comparison_cell_type", "reject_null", "fc"]]
                )
        combined = pd.concat(rows, ignore_index=True)
        return {
            ct: combined[combined["comparison_cell_type"] == ct][
                ["reject_null", "fc"]
            ].reset_index(drop=True)
            for ct in CELL_TYPES
        }

    all_mergy_reporter = sum_pow_by_ct(sims, "hs_all_ct", "ttest")
    all_mergy_deflated = sum_pow_by_ct(sims, "hs_all_ct", "ttest_deflated")

    # Save aggregated dataframes
    for label, mergy_dict in [
        ("reporter", all_mergy_reporter),
        ("deflated", all_mergy_deflated),
    ]:
        rows = []
        for ct, df in mergy_dict.items():
            df = df.copy()
            df["cell_type"] = ct
            rows.append(df)
        combined = pd.concat(rows, ignore_index=True)
        out_path = OUTPUT_DIR / f"cohen_power_df_{label}.parquet"
        combined.to_parquet(out_path)
        print(f"Saved: {out_path}", flush=True)

    # --- Plot: overlay +reporter and deflated on each subplot ---
    def _pow_curve_data(mergy, n_bins=100):
        df = mergy.copy()
        df["fc"] = pd.cut(df["fc"], bins=n_bins)
        binned = (
            df.groupby("fc", observed=True)["reject_null"]
            .mean()
            .reset_index(name="reject_frac")
        )
        binned["bin_center"] = binned["fc"].apply(lambda x: x.mid)
        return binned

    ncols = 2
    nrows = math.ceil(len(CELL_TYPES) / ncols)
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(5 * ncols, 4 * nrows), sharey=True
    )
    axes = axes.flatten()

    for i, ct in enumerate(CELL_TYPES):
        ct_display = "Rod" if ct == "reference" else ct
        ax = axes[i]

        # +reporter curve
        binned_rep = _pow_curve_data(all_mergy_reporter[ct])
        ax.plot(
            binned_rep["bin_center"], binned_rep["reject_frac"],
            color="steelblue", marker="o", markersize=2, linewidth=1,
            label="+reporter",
        )

        # deflated curve
        binned_def = _pow_curve_data(all_mergy_deflated[ct])
        ax.plot(
            binned_def["bin_center"], binned_def["reject_frac"],
            color="coral", marker="s", markersize=2, linewidth=1,
            label="-reporter (deflated)",
        )

        ax.axhline(0.8, color="black", linestyle="--", lw=0.8)
        ax.axvline(1.0, color="grey",  linestyle=":",  lw=0.8)
        ax.set_title(ct_display, fontsize=9)
        ax.set_xlabel("FC (activity / minP)", fontsize=8)
        ax.set_ylabel("Power (TPR)", fontsize=8)
        ax.set_ylim(0, 1)
        cells = scm.COHEN_BOUNDS.cells_per_cell_type.get(ct, "?")
        ax.text(
            0.97, 0.05,
            f"n_cres={N_CRES}, cells={cells}",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=6, color="grey",
        )
        if i == 0:
            ax.legend(fontsize=7, loc="upper left")

    for j in range(len(CELL_TYPES), len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(
        "Welch's t-test Power by Cell Type -- Cohen (episomal)\n"
        "+reporter (obs condition) vs -reporter (zero-deflated)",
        fontsize=12,
    )
    plt.tight_layout()
    svg_path = OUTPUT_DIR / "power_ttest_reporter_comparison.svg"
    fig.savefig(svg_path, format="svg", bbox_inches="tight")
    print(f"Saved: {svg_path}", flush=True)


# ---------------------------------------------------------------------------
# Main dispatch
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    phase = sys.argv[1] if len(sys.argv) > 1 else "all"

    def _make_slurm_client():
        """SLURMCluster for heavy simulate phase."""
        cluster = SLURMCluster(
            cores=1,
            memory="64G",
            processes=1,
            env_extra=[
                f"export PYTHONPATH={_repo_root}:$PYTHONPATH",
            ],
            job_extra_directives=[
                "-p priority",
                "--account=prio_skr2",
                "--job-name=cohen_pow_worker",
                "--time=6:00:00",
                "--output=worker_%j.out",
            ],
        )
        cluster.scale(jobs=10)
        client = Client(
            cluster,
            timeout=f"{5*60}s",
            heartbeat_interval="20s",
        )
        print(f"Dask dashboard: {client.dashboard_link}", flush=True)
        return cluster, client

    def _make_local_client():
        """Lightweight LocalCluster for plot/aggregation phase."""
        import psutil
        _slurm_mem_mb = os.environ.get("SLURM_MEM_PER_NODE")
        if _slurm_mem_mb:
            _mem_limit = f"{int(int(_slurm_mem_mb) * 0.9)}MB"
        else:
            _mem_limit = f"{int(psutil.virtual_memory().total * 0.9 / 1e9)}GB"
        cluster = LocalCluster(
            n_workers=1,
            threads_per_worker=2,
            processes=False,
            memory_limit=_mem_limit,
        )
        client = Client(cluster)
        print(f"Dask dashboard: {client.dashboard_link}", flush=True)
        return cluster, client

    if phase == "all":
        # Simulate with SLURM workers, then plot locally
        print(f"\n{'='*60}", flush=True)
        print("Phase: simulate", flush=True)
        print(f"{'='*60}", flush=True)
        cluster, client = _make_slurm_client()
        phase_simulate(client)
        client.close()
        cluster.close()

        print(f"\n{'='*60}", flush=True)
        print("Phase: plot", flush=True)
        print(f"{'='*60}", flush=True)
        cluster, client = _make_local_client()
        phase_plot(client)
        client.close()
        cluster.close()

    elif phase == "simulate":
        cluster, client = _make_slurm_client()
        phase_simulate(client)
        client.close()
        cluster.close()

    elif phase == "plot":
        cluster, client = _make_local_client()
        phase_plot(client)
        client.close()
        cluster.close()

    else:
        print(
            f"Unknown phase: {phase!r}. "
            f"Choose from: simulate, plot, or all"
        )
        sys.exit(1)

    print("\nDone.", flush=True)
