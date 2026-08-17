#!/usr/bin/env python3
"""
Seelig activity power -- all cell types, both deflated test arms.

Single-condition power analysis for the Seelig dataset, which has no
transfection reporter. The only realistic test condition is the naive
drop-zeros (has_reporter=False) test on observed nonzero counts: there is no
reporter signal to anchor zero imputation, and MOIB-at-test-time was
empirically rejected (see moib_test_pilot.py for the negative result --
inflated Type I error). Fit-time MOIB Wald (via the seelig_cm_moib_nb_phantom
ortho) is the principled alternative but is not what this script measures.

Because there is no reporter, only the deflated arms apply; the +reporter arms
run for Shendure have no counterpart here. Both Welch's t-test and MWU are run
on identical simulated data, so the two arms are exactly paired: the t-test arm
is the reproduction check, MWU is the test Results 2.2 adopts.

For each cell type:
1. Draw `n_cres` synthetic CREs from uniform(min, max), with one fixed at
   `minP` as the reference. n_cres comes from the canonical
   seelig_cm_nb_phantom ortho.
2. Simulate `cells_per_cell_type` cells across N_SIMS replicates, repeated
   N_LIBRARY_REPS times.
3. Run every arm in ARMS against the same sims.
4. Aggregate and plot per-cell-type power curves, overlaying both arms.

Usage:
    python seelig_power_ttest_all_cell_types.py [phase]

Phases:
    simulate  -- generate simulations and run both test arms
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

SIM_DATE = "2026-08-17"
SIM_DIR = DATA_ROOT / "simulated" / f"{SIM_DATE}_seelig_pow"

# {t-test, MWU} x {deflated}, as (test directory, method, reporter). Yin ran no
# transfection reporter, so the +reporter arms that Shendure carries have no
# counterpart and are omitted. Aggregates are written per arm with SIM_DATE in
# the filename.
ARMS = [
    ("ttest_deflated", "ttest", False),
    ("mwu_deflated",   "mwu",   False),
]

# Use the canonical phantom ortho for n_cres per cell type. Mu values are not
# read here -- the simulation draws synthetic CRE strengths from
# uniform(MIN_ACTIVITY, MAX_ACTIVITY) bounded by SEELIG_BOUNDS.
ORTHO_DIR = DATA_ROOT / "seelig" / "seelig_cm_nb_phantom"
with open(ORTHO_DIR / "by_cell_type_parameters.pkl", "rb") as f:
    _params = pickle.load(f)

CELL_TYPES = sorted(_params.nb.keys())
N_CRES_PER_CT = {ct: len(_params.nb[ct]) for ct in CELL_TYPES}
del _params

assert CELL_TYPES, "ortho yielded no cell types"
assert all(N_CRES_PER_CT[ct] > 1 for ct in CELL_TYPES), (
    f"a cell type has too few CREs to hold out a reference: {N_CRES_PER_CT}"
)

MINP = scm.SEELIG_BOUNDS.reference_activity
MAX_ACTIVITY = 4.0 * MINP
MIN_ACTIVITY = scm.SEELIG_BOUNDS.min_mpra_umi
N_LIBRARY_REPS = 156
N_SIMS = 5

assert MIN_ACTIVITY < MINP < MAX_ACTIVITY, (
    f"activity bounds must bracket the reference: "
    f"min={MIN_ACTIVITY:.6e}, minP={MINP:.6e}, max={MAX_ACTIVITY:.6e}"
)


def _has_all_data(sim_dir):
    """Whether a sim directory holds a complete set of simulated datasets.

    Resume keys on the data rather than on the test results: the arm loop
    already skips arms it finds on disk, so a run interrupted partway through
    testing can pick up its existing sims and compute only the missing arms.
    Keying on test results instead would discard hours of simulation whenever
    a run died between arms.
    """
    if not (sim_dir / "state.parquet").exists():
        return False
    return all(
        (sim_dir / "simulated_scmpra" / f"{i}.scmpra" / "data.parquet").exists()
        for i in range(N_SIMS)
    )


def _has_all_arms(sim_dir):
    arm_dirs = [sim_dir / "tests" / "hs_all_ct" / a for a, _, _ in ARMS]
    return all(p.exists() and any(p.iterdir()) for p in arm_dirs)


# ---------------------------------------------------------------------------
# Phase: simulate -- generate data and run every arm in ARMS
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
        assert len(bound_ct.cells_per_cell_type) == 1, (
            f"per-cell-type bound must hold exactly one cell type, got "
            f"{len(bound_ct.cells_per_cell_type)} for {cell_type}"
        )
        ct_dir = SIM_DIR / cell_type

        # Resume: reload every sim whose simulated data is complete.
        sims_ct = []
        n_untested = 0
        if ct_dir.exists():
            for d in sorted(ct_dir.iterdir()):
                if not d.is_dir() or not _has_all_data(d):
                    continue
                sims_ct.append(
                    scm.de_novo_simulation(
                        location=ct_dir, name=d.name, client=client
                    )
                )
                if not _has_all_arms(d):
                    n_untested += 1
        print(
            f"  Resumed {len(sims_ct)} sims from disk "
            f"({n_untested} still missing at least one arm).",
            flush=True,
        )

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
        # An interrupted run can leave more sims on disk than requested; extra
        # replicates only sharpen the estimate, too few means data went missing.
        assert len(sims_ct) >= N_LIBRARY_REPS, (
            f"{cell_type}: expected at least {N_LIBRARY_REPS} replicates, have "
            f"{len(sims_ct)}"
        )
        print(f"  {len(sims_ct)} replicates ready.", flush=True)

    for sims in all_sims.values():
        for sim in sims:
            sim.save()
    print("\nAll sims saved.", flush=True)

    # --- run every arm on the same simulated data ---
    # The deflated arms drop all zero-count observations before testing, which
    # is what a no-reporter experiment actually sees on disc. Without that the
    # test gets an unfair advantage: the simulation writes explicit zeros for
    # transfected-but-silent events that such an experiment could never
    # observe. Both arms read the same sims, so they are exactly paired.
    for arm, method, has_reporter in ARMS:
        for cell_type, sims in all_sims.items():
            print(f"{arm}: {cell_type}", flush=True)
            hs = None
            for sim in sims:
                arm_dir = sim.location / sim.name / "tests" / "hs_all_ct" / arm
                if arm_dir.exists() and any(arm_dir.iterdir()):
                    continue  # already done
                if hs is None:
                    example = scm.scMPRA_data.from_parquet(
                        sim.scmpradatp / "0.scmpra"
                    )
                    hs = scm.make_all_by_celltype_hypotheses(
                        counts=example, reference_cre="reference"
                    )
                # May already exist from an earlier arm, but not if we resumed
                # from disk.
                hypo_file = (
                    sim.location / sim.name / "tests" / "hs_all_ct"
                    / "hypotheses.tsv"
                )
                if not hypo_file.exists():
                    sim.add_hypothesis_set("hs_all_ct", hs)
                getattr(sim, method)("hs_all_ct", has_reporter=has_reporter)
                sim.save()

        for sims in all_sims.values():
            for sim in sims:
                sim.save()
        print(f"{arm} done and saved.", flush=True)


# ---------------------------------------------------------------------------
# Phase: plot -- aggregate and produce SVG
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
                if _has_all_arms(d):
                    sims_ct.append(
                        scm.de_novo_simulation(
                            location=ct_dir, name=d.name, client=client
                        )
                    )
        all_sims[cell_type] = sims_ct
        print(f"  {cell_type}: loaded {len(sims_ct)} sims", flush=True)

    assert all(all_sims[ct] for ct in CELL_TYPES), (
        f"a cell type has no complete sims: "
        f"{ {ct: len(s) for ct, s in all_sims.items()} }"
    )

    per_arm = {}
    for arm, _, _ in ARMS:
        per_arm[arm] = {
            ct: scm.sum_pow(
                sims, hypothesis_set_name="hs_all_ct", test_type=arm
            )
            for ct, sims in all_sims.items()
        }
        rows = []
        n_expected = 0
        for ct, df in per_arm[arm].items():
            n_expected += len(df)
            df = df.copy()
            df["cell_type"] = ct
            rows.append(df)
        combined = pd.concat(rows, ignore_index=True)
        assert len(combined) == n_expected, (
            f"{arm}: concat changed row count, {n_expected} -> "
            f"{len(combined)}"
        )
        assert combined["reject_null"].notna().all(), (
            f"{arm}: {combined['reject_null'].isna().sum()} null test outcomes"
        )
        assert (combined["fc"] > 0).all(), (
            f"{arm}: {(combined['fc'] <= 0).sum()} non-positive fold changes"
        )
        assert set(combined["cell_type"]) == set(CELL_TYPES), (
            f"{arm}: cell types in aggregate {sorted(set(combined['cell_type']))} "
            f"do not match {CELL_TYPES}"
        )
        out_path = OUTPUT_DIR / f"power_df_{arm}_{SIM_DATE}.parquet"
        combined.to_parquet(out_path)
        print(f"Saved: {out_path} ({len(combined):,} rows)", flush=True)

    # Power is monotone in effect size, so the strongest CREs must be easier to
    # detect than the weakest. A violation means the arms or the ground truth
    # got crossed somewhere upstream.
    for arm, _, _ in ARMS:
        for ct, df in per_arm[arm].items():
            lo = df.loc[df["fc"] < 1.05, "reject_null"].mean()
            hi = df.loc[df["fc"] > 3.5, "reject_null"].mean()
            assert lo < hi, (
                f"{arm}/{ct}: power does not increase with effect size, "
                f"reject rate near FC=1 is {lo:.3f} but {hi:.3f} at FC>3.5"
            )

    # --- Plot: overlay both arms on each subplot ---
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
    axes = np.array(axes).flatten()

    arm_style = {
        "ttest_deflated": ("coral", "s", "t-test (deflated)"),
        "mwu_deflated": ("steelblue", "o", "MWU (deflated)"),
    }

    for i, cell_type in enumerate(CELL_TYPES):
        # In the Seelig bounds, "reference" is HepG2 (chosen at preprocessing
        # time as the higher-cell-count of the two CTs: 5858 HepG2 vs 4478
        # K562). Both are cell lines, neither is pluripotent.
        ct_display = "HepG2" if cell_type == "reference" else cell_type
        ax = axes[i]

        for arm, _, _ in ARMS:
            color, marker, label = arm_style[arm]
            binned = _pow_curve_data(per_arm[arm][cell_type])
            ax.plot(
                binned["bin_center"], binned["reject_frac"],
                color=color, marker=marker, markersize=2, linewidth=1,
                label=label,
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
        "Power by Cell Type -- Seelig\n"
        "no transfection reporter (deflated tests on nonzero counts)",
        fontsize=12,
    )
    plt.tight_layout()
    svg_path = OUTPUT_DIR / f"power_arm_comparison_{SIM_DATE}.svg"
    fig.savefig(svg_path, format="svg", bbox_inches="tight")
    print(f"Saved: {svg_path}", flush=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    phase = sys.argv[1] if len(sys.argv) > 1 else "all"

    def _make_slurm_client():
        # Worker memory: the Shendure power run peaked at ~9.2 GB per worker.
        # A Seelig sim object carries n_cres x n_cells about 4.6x larger, so
        # the per-worker footprint is sized up in proportion with headroom.
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
                "--job-name=seelig_pow_worker",
                "--time=12:00:00",
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
