#!/usr/bin/env python3
"""
Seelig t-test power -- all cell types.

Single-condition power analysis for the Seelig dataset, which has no
transfection reporter. The only realistic test condition is the naive
drop-zeros (has_reporter=False) Welch t-test on observed nonzero counts:
there is no reporter signal to anchor zero imputation, and MOIB-at-test-time
was empirically rejected (see moib_test_pilot.py for the negative result --
inflated Type I error). Fit-time MOIB Wald (via the seelig_cm_moib_nb_phantom
ortho) is the principled alternative but is not what this script measures.

For each cell type:
1. Draw `n_cres` synthetic CREs from uniform(min, max), with one fixed at
   `minP` as the reference. n_cres comes from the canonical
   seelig_cm_nb_phantom ortho.
2. Simulate `cells_per_cell_type` cells across N_SIMS replicates, repeated
   N_LIBRARY_REPS times.
3. Run Welch's t-test in the deflated condition only (has_reporter=False --
   drop zeros before testing). The simulator generates the full
   transfected-but-silent zero set; this condition reflects what a no-reporter
   experiment would actually observe.
4. Aggregate and plot per-cell-type power curves.

Usage:
    python seelig_power_ttest_all_cell_types.py [phase]

Phases:
    simulate  -- generate simulations and run t-test
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
SIM_DIR = DATA_ROOT / "simulated" / f"{SIM_DATE}_seelig_pow"

# Use the canonical phantom ortho for n_cres per cell type. Mu values are not
# read here -- the simulation draws synthetic CRE strengths from
# uniform(MIN_ACTIVITY, MAX_ACTIVITY) bounded by SEELIG_BOUNDS.
ORTHO_DIR = DATA_ROOT / "seelig" / "seelig_cm_nb_phantom"
with open(ORTHO_DIR / "by_cell_type_parameters.pkl", "rb") as f:
    _params = pickle.load(f)

CELL_TYPES = sorted(_params.nb.keys())
N_CRES_PER_CT = {ct: len(_params.nb[ct]) for ct in CELL_TYPES}
del _params

MINP = scm.SEELIG_BOUNDS.reference_activity
MAX_ACTIVITY = 4.0 * MINP
MIN_ACTIVITY = scm.SEELIG_BOUNDS.min_mpra_umi
N_LIBRARY_REPS = 156
N_SIMS = 5


# ---------------------------------------------------------------------------
# Phase: simulate
# ---------------------------------------------------------------------------

def phase_simulate(client):
    print(f"minP={MINP:.6f}, min={MIN_ACTIVITY:.6e}, max={MAX_ACTIVITY:.6f}")
    print(f"Cell types: {len(CELL_TYPES)}, library reps: {N_LIBRARY_REPS}, "
          f"sims per rep: {N_SIMS}")
    print()

    for ct in CELL_TYPES:
        cells = scm.SEELIG_BOUNDS.cells_per_cell_type.get(ct, "N/A")
        print(f"  {ct}: n_cres={N_CRES_PER_CT[ct]}, cells={cells}")

    all_sims = {}

    for cell_type in CELL_TYPES:
        print(f"\n=== {cell_type} ===", flush=True)

        bound_ct = scm.SEELIG_BOUNDS.copy()
        bound_ct.cells_per_cell_type = (
            scm.SEELIG_BOUNDS.cells_per_cell_type.loc[[cell_type]]
        )
        ct_dir = SIM_DIR / cell_type

        # Resume: a sim is "complete" if its ttest_deflated/ dir exists and is
        # non-empty.
        sims_ct = []
        if ct_dir.exists():
            for d in sorted(ct_dir.iterdir()):
                defl_dir = d / "tests" / "hs_all_ct" / "ttest_deflated"
                if defl_dir.exists() and any(defl_dir.iterdir()):
                    sims_ct.append(
                        scm.de_novo_simulation(
                            location=ct_dir, name=d.name, client=client
                        )
                    )
        print(f"  Resumed {len(sims_ct)} complete sims from disk.", flush=True)

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

    for sims in all_sims.values():
        for sim in sims:
            sim.save()
    print("\nAll sims saved.", flush=True)

    # --- t-test: deflated (no reporter) ---
    # Drop zeros before testing. This matches what a no-reporter experiment
    # like Seelig actually observes -- only nonzero MPRA counts make it onto
    # disk.
    for cell_type, sims in all_sims.items():
        print(f"t-test (deflated): {cell_type}", flush=True)
        hs = None
        for sim in sims:
            defl_dir = (
                sim.location / sim.name / "tests" / "hs_all_ct" / "ttest_deflated"
            )
            if defl_dir.exists() and any(defl_dir.iterdir()):
                continue
            if hs is None:
                example = scm.scMPRA_data.from_parquet(
                    sim.scmpradatp / "0.scmpra"
                )
                hs = scm.make_all_by_celltype_hypotheses(
                    counts=example, reference_cre="reference"
                )
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
# Phase: plot
# ---------------------------------------------------------------------------

def phase_plot(client):
    all_sims = {}
    for cell_type in CELL_TYPES:
        ct_dir = SIM_DIR / cell_type
        sims_ct = []
        if ct_dir.exists():
            for d in sorted(ct_dir.iterdir()):
                if not d.is_dir():
                    continue
                defl_dir = d / "tests" / "hs_all_ct" / "ttest_deflated"
                if defl_dir.exists() and any(defl_dir.iterdir()):
                    sims_ct.append(
                        scm.de_novo_simulation(
                            location=ct_dir, name=d.name, client=client
                        )
                    )
        all_sims[cell_type] = sims_ct
        print(f"  {cell_type}: loaded {len(sims_ct)} sims", flush=True)

    all_mergy_deflated = {
        ct: scm.sum_pow(
            sims, hypothesis_set_name="hs_all_ct", test_type="ttest_deflated"
        )
        for ct, sims in all_sims.items()
    }

    rows = []
    for ct, df in all_mergy_deflated.items():
        df = df.copy()
        df["cell_type"] = ct
        rows.append(df)
    combined = pd.concat(rows, ignore_index=True)
    out_path = OUTPUT_DIR / "seelig_power_df_deflated.parquet"
    combined.to_parquet(out_path)
    print(f"Saved: {out_path}", flush=True)

    # --- Plot ---
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

    ncols = min(2, len(CELL_TYPES))
    nrows = math.ceil(len(CELL_TYPES) / ncols)
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(5 * ncols, 4 * nrows), sharey=True
    )
    if len(CELL_TYPES) == 1:
        axes = np.array([axes])
    axes = np.array(axes).flatten()

    for i, cell_type in enumerate(CELL_TYPES):
        # In the Seelig bounds, "reference" is HepG2 (chosen at preprocessing
        # time as the higher-cell-count of the two CTs: 5858 HepG2 vs 4478
        # K562). Both are cell lines, neither is pluripotent.
        ct_display = "HepG2" if cell_type == "reference" else cell_type
        ax = axes[i]

        binned = _pow_curve_data(all_mergy_deflated[cell_type])
        ax.plot(
            binned["bin_center"], binned["reject_frac"],
            color="coral", marker="s", markersize=2, linewidth=1,
            label="-reporter (deflated)",
        )

        ax.axhline(0.8, color="black", linestyle="--", lw=0.8)
        ax.axvline(1.0, color="grey", linestyle=":", lw=0.8)
        ax.set_title(ct_display, fontsize=10)
        ax.set_xlabel("FC (activity / minP)", fontsize=8)
        ax.set_ylabel("Power (TPR)", fontsize=8)
        ax.set_ylim(0, 1)
        cells = scm.SEELIG_BOUNDS.cells_per_cell_type.get(cell_type, "?")
        ax.text(
            0.97, 0.05,
            f"n_cres={N_CRES_PER_CT[cell_type]}, cells={cells}",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=7, color="grey",
        )
        if i == 0:
            ax.legend(fontsize=7, loc="upper left")

    for j in range(len(CELL_TYPES), len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(
        "Welch's t-test Power by Cell Type -- Seelig\n"
        "no transfection reporter (deflated test on nonzero counts)",
        fontsize=12,
    )
    plt.tight_layout()
    svg_path = OUTPUT_DIR / "seelig_power_ttest_deflated.svg"
    fig.savefig(svg_path, format="svg", bbox_inches="tight")
    print(f"Saved: {svg_path}", flush=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    phase = sys.argv[1] if len(sys.argv) > 1 else "all"

    def _make_slurm_client():
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
                "--job-name=seelig_pow_worker",
                "--time=8:00:00",
                "--exclude=a1132u18n02",
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
