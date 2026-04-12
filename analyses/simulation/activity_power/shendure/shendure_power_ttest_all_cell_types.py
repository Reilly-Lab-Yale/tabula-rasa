#!/usr/bin/env python3
"""
Shendure t-test power -- all cell types

Assesses Welch's t-test power (TPR) as a function of CRE_oi strength relative
to `minP`, for every cell type in the Shendure dataset.

For each cell type:
1. Draw `n_cres` synthetic CREs from uniform(min, max), with one fixed at
   `minP` as the reference
2. Simulate `cells_per_cell_type` cells (matching the real count for that cell
   type) across `n_sims` replicates, repeated 100 times
3. Run Welch's t-test for all CREs vs. reference -- both with reporter
   (standard) and without reporter (zero-deflated, simulating a no-reporter
   experiment)
4. Aggregate and plot per-cell-type power curves, overlaying both conditions

Each cell type is simulated independently with its own empirical
`(n_cres, cell_count)` pair.
`one_library_replicate` now accepts a `cell_type` parameter so real names flow
through cleanly.

Usage:
    python shendure_power_ttest_all_cell_types.py [phase]

Phases:
    simulate  -- generate simulations and run t-test (both conditions)
    plot      -- aggregate results and produce SVG
    all       -- run both phases sequentially (default)
"""

import sys
import os
import pickle
import math
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

SIM_DATE = "2026-04-10"
SIM_DIR = DATA_ROOT / "simulated" / f"{SIM_DATE}_shendure_pow"

# Load per-cell-type n_cres from the canonical phantom ortho parameters
# (avoids loading the full ortho object, which has lazy file references to the
# training data).
ORTHO_DIR = DATA_ROOT / "shendure" / "shendure_obs_nb_phantom"
with open(ORTHO_DIR / "by_cell_type_parameters.pkl", "rb") as f:
    _params = pickle.load(f)

CELL_TYPES = sorted(_params.nb.keys())
N_CRES_PER_CT = {ct: len(_params.nb[ct]) for ct in CELL_TYPES}
del _params

# Original notebook used FC range ~[0.18, 2.59] (with old minP=0.019,
# MAX_ACTIVITY=0.05). FC=2.59 turned out to be too narrow -- the deflated
# (no-reporter) condition fails to reach 80% power for most cell types within
# that range. Extended to FC=4 to capture the threshold for the long-tail cell
# types. N_LIBRARY_REPS bumped 100 -> 156 (1.56x) to maintain sim density per
# unit FC across the wider range.
MINP = scm.SHENDURE_BOUNDS.reference_activity
MAX_ACTIVITY = 4.0 * MINP
MIN_ACTIVITY = scm.SHENDURE_BOUNDS.min_mpra_umi
N_LIBRARY_REPS = 156
N_SIMS = 5


# ---------------------------------------------------------------------------
# Phase: simulate -- generate data and run t-test (both +reporter and deflated)
# ---------------------------------------------------------------------------

def phase_simulate(client):
    print(f"minP={MINP:.6f}, min={MIN_ACTIVITY:.6f}, max={MAX_ACTIVITY}")
    print(f"Cell types: {len(CELL_TYPES)}, library reps: {N_LIBRARY_REPS}, "
          f"sims per rep: {N_SIMS}")
    print()

    for ct in CELL_TYPES:
        cells = scm.SHENDURE_BOUNDS.cells_per_cell_type.get(ct, "N/A")
        print(f"  {ct}: n_cres={N_CRES_PER_CT[ct]}, cells={cells}")

    all_sims = {}

    for cell_type in CELL_TYPES:
        print(f"\n=== {cell_type} ===", flush=True)

        bound_ct = scm.SHENDURE_BOUNDS.copy()
        bound_ct.cells_per_cell_type = (
            scm.SHENDURE_BOUNDS.cells_per_cell_type.loc[[cell_type]]
        )
        ct_dir = SIM_DIR / cell_type

        # Resume: reload already-complete sims from disk.
        # A sim is "complete" if it has both ttest/ and ttest_deflated/ results.
        sims_ct = []
        if ct_dir.exists():
            for d in sorted(ct_dir.iterdir()):
                tt_dir = d / "tests" / "hs_all_ct" / "ttest"
                defl_dir = d / "tests" / "hs_all_ct" / "ttest_deflated"
                if (tt_dir.exists() and any(tt_dir.iterdir())
                        and defl_dir.exists() and any(defl_dir.iterdir())):
                    sims_ct.append(
                        scm.de_novo_simulation(
                            location=ct_dir, name=d.name, client=client
                        )
                    )
        print(f"  Resumed {len(sims_ct)} complete sims from disk.", flush=True)

        # Run any remaining reps
        n_remaining = N_LIBRARY_REPS - len(sims_ct)
        print(f"  Running {n_remaining} fresh sims.", flush=True)
        for i in range(n_remaining):
            _, sim = scm.one_library_replicate(
                root=ct_dir,
                n_sims=N_SIMS,
                client=client,
                flatten_overtransfection=True,
                bound=bound_ct,
                n_cres=N_CRES_PER_CT[cell_type],
                min=MIN_ACTIVITY,
                max=MAX_ACTIVITY,
                minP=MINP,
                cell_type=cell_type,
            )
            sims_ct.append(sim)

        all_sims[cell_type] = sims_ct
        print(f"  {len(sims_ct)} replicates ready.", flush=True)

    # Save all sims
    for sims in all_sims.values():
        for sim in sims:
            sim.save()
    print("\nAll sims saved.", flush=True)

    # --- t-test: +reporter (standard) ---
    # Build a hypothesis set per cell type and run t-test (skip if already done).
    for cell_type, sims in all_sims.items():
        print(f"t-test (+reporter): {cell_type}", flush=True)
        hs = None
        for sim in sims:
            tt_dir = sim.location / sim.name / "tests" / "hs_all_ct" / "ttest"
            if tt_dir.exists() and any(tt_dir.iterdir()):
                continue  # already done
            if hs is None:
                example = scm.scMPRA_data.from_parquet(
                    sim.scmpradatp / "0.scmpra"
                )
                hs = scm.make_all_by_celltype_hypotheses(
                    counts=example, reference_cre="reference"
                )
            sim.add_hypothesis_set("hs_all_ct", hs)
            sim.ttest("hs_all_ct")
            sim.save()

    for sims in all_sims.values():
        for sim in sims:
            sim.save()
    print("t-test (+reporter) done and saved.", flush=True)

    # --- t-test: deflated (no reporter) ---
    # Same simulated data, but all zero-count observations are dropped before
    # testing. This simulates what a real no-reporter experiment looks like on
    # disc: only non-zero MPRA signal is observed. Without this, simulated data
    # gives the test an unfair advantage because the simulation generates
    # explicit zeros for transfected-but-silent events, which a no-reporter
    # experiment would never observe.
    for cell_type, sims in all_sims.items():
        print(f"t-test (deflated): {cell_type}", flush=True)
        hs = None
        for sim in sims:
            defl_dir = (
                sim.location / sim.name / "tests" / "hs_all_ct" / "ttest_deflated"
            )
            if defl_dir.exists() and any(defl_dir.iterdir()):
                continue  # already done
            if hs is None:
                example = scm.scMPRA_data.from_parquet(
                    sim.scmpradatp / "0.scmpra"
                )
                hs = scm.make_all_by_celltype_hypotheses(
                    counts=example, reference_cre="reference"
                )
            # Ensure hypothesis set exists (may have been added in +reporter
            # pass, but if we resumed from disk it won't be).
            hypo_file = (
                sim.location / sim.name / "tests" / "hs_all_ct" / "hypotheses.tsv"
            )
            if not hypo_file.exists():
                sim.add_hypothesis_set("hs_all_ct", hs)
            sim.ttest("hs_all_ct", has_reporter=False)
            sim.save()

    for sims in all_sims.values():
        for sim in sims:
            sim.save()
    print("t-test (deflated) done and saved.", flush=True)


# ---------------------------------------------------------------------------
# Phase: plot -- aggregate and produce SVG
# ---------------------------------------------------------------------------

def phase_plot(client):
    # Reload all sims from disk
    all_sims = {}
    for cell_type in CELL_TYPES:
        ct_dir = SIM_DIR / cell_type
        sims_ct = []
        if ct_dir.exists():
            for d in sorted(ct_dir.iterdir()):
                if not d.is_dir():
                    continue
                tt_dir = d / "tests" / "hs_all_ct" / "ttest"
                if tt_dir.exists() and any(tt_dir.iterdir()):
                    sims_ct.append(
                        scm.de_novo_simulation(
                            location=ct_dir, name=d.name, client=client
                        )
                    )
        all_sims[cell_type] = sims_ct
        print(f"  {cell_type}: loaded {len(sims_ct)} sims", flush=True)

    # Aggregate results for both conditions
    all_mergy_reporter = {
        ct: scm.sum_pow(sims, hypothesis_set_name="hs_all_ct", test_type="ttest")
        for ct, sims in all_sims.items()
    }
    all_mergy_deflated = {
        ct: scm.sum_pow(
            sims, hypothesis_set_name="hs_all_ct", test_type="ttest_deflated"
        )
        for ct, sims in all_sims.items()
    }

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
        out_path = OUTPUT_DIR / f"power_df_{label}.parquet"
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

    ncols = 3
    nrows = math.ceil(len(CELL_TYPES) / ncols)
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(5 * ncols, 4 * nrows), sharey=True
    )
    axes = axes.flatten()

    for i, cell_type in enumerate(CELL_TYPES):
        cell_type_display = (
            "Pluripotent" if cell_type == "reference" else cell_type
        )
        ax = axes[i]

        # +reporter curve
        binned_rep = _pow_curve_data(all_mergy_reporter[cell_type])
        ax.plot(
            binned_rep["bin_center"], binned_rep["reject_frac"],
            color="steelblue", marker="o", markersize=2, linewidth=1,
            label="+reporter",
        )

        # deflated curve
        binned_def = _pow_curve_data(all_mergy_deflated[cell_type])
        ax.plot(
            binned_def["bin_center"], binned_def["reject_frac"],
            color="coral", marker="s", markersize=2, linewidth=1,
            label="-reporter (deflated)",
        )

        ax.axhline(0.8, color="black", linestyle="--", lw=0.8)
        ax.axvline(1.0, color="grey", linestyle=":", lw=0.8)
        ax.set_title(cell_type_display, fontsize=9)
        ax.set_xlabel("FC (activity / minP)", fontsize=8)
        ax.set_ylabel("Power (TPR)", fontsize=8)
        ax.set_ylim(0, 1)
        cells = scm.SHENDURE_BOUNDS.cells_per_cell_type.get(cell_type, "?")
        ax.text(
            0.97, 0.05,
            f"n_cres={N_CRES_PER_CT[cell_type]}, cells={cells}",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=6, color="grey",
        )
        if i == 0:
            ax.legend(fontsize=7, loc="upper left")

    for j in range(len(CELL_TYPES), len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(
        "Welch's t-test Power by Cell Type -- Shendure\n"
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
        # Worker memory: peak RSS observed in the FC=2.59 run was ~17 GB on
        # the high-cell-count cell types (Pluripotent, NeuroectodermBrain).
        # 24 GB gives ~1.4x headroom. Per-rep memory does not scale with
        # N_LIBRARY_REPS (reps are sequential), so this is sufficient for the
        # FC=4 run too. Worker count doubled 10 -> 20 to halve wall-clock on
        # the t-test phase, which fans out wider than the inner-rep graph.
        cluster = SLURMCluster(
            cores=1,
            memory="24G",
            processes=1,
            env_extra=[
                f"export PYTHONPATH={_repo_root}:$PYTHONPATH",
            ],
            job_extra_directives=[
                "-p priority",
                "--account=prio_skr2",
                "--job-name=shendure_pow_worker",
                "--time=8:00:00",
                "--output=worker_%j.out",
            ],
        )
        cluster.scale(jobs=20)
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
