#!/usr/bin/env python3
"""
Design-space 1D sweeps: vary one experimental-design axis at a time, hold
others at the anchor (shendure Pluripotent), measure Welch t-test power.

Axes:
    bcs_per_cre  -- library_model.mu_nb (NB mean for per-CRE barcode count)
    n_cells      -- cells_per_cell_type[ct]
    moi          -- transfection_model.mu_nb (via bounds.set_effective_moi)
    n_cres       -- one_library_replicate(n_cres=...) kwarg

Rationale: the naive cross-dataset regression
(analyses/simulation/cross_dataset_power_regression) identified bcs_per_cre
as the largest-partial-R^2 predictor of FC@50%-power, but that regression is
only 16 rows and confounded across datasets. This sweep tests the claim in
simulation by holding everything else fixed.

Usage:
    python design_space_sweep.py simulate bcs_per_cre
    python design_space_sweep.py plot bcs_per_cre

Phases:
    simulate  -- generate simulations and run t-test for the chosen axis
    plot      -- aggregate power_df, make SVG
    all       -- both sequentially (default)
"""

import sys
import os
import json
import math
import copy
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Put repo root on PYTHONPATH so scMPRAforge imports cleanly.
_repo_root = Path(__file__).resolve().parents[3]
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

import scMPRAforge as scm
from dask.distributed import Client, LocalCluster
from dask_jobqueue import SLURMCluster

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

DATA_ROOT = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new")
OUT = Path(__file__).resolve().parent / "output"
OUT.mkdir(exist_ok=True)
SWEEPS_DIR = OUT / "sweeps"
SWEEPS_DIR.mkdir(exist_ok=True)

# Sims live on scratch: project partition hit 100% full during the initial
# bcs_per_cre runs, causing slurm to cancel jobs before they could even write
# slurm-*.out (0:53 exit signature with 0 elapsed time). Scratch has 5 TB free.
SIM_DATE = "2026-04-11"
SIM_ROOT = Path("/nfs/roberts/scratch/pi_skr2/mcn26/design_space_sims") / f"{SIM_DATE}_design_space"

# ---------------------------------------------------------------------------
# Anchor: shendure, Pluripotent (stored as "reference" in bounds + orthos)
# ---------------------------------------------------------------------------

ANCHOR_CT = "reference"
ANCHOR_BOUNDS = scm.SHENDURE_BOUNDS

ANCHOR_N_CELLS = int(ANCHOR_BOUNDS.cells_per_cell_type.loc[ANCHOR_CT])
ANCHOR_BCS_PER_CRE = float(ANCHOR_BOUNDS.library_model.mu_nb)
ANCHOR_MOI = float(ANCHOR_BOUNDS.get_effective_moi())

# n_cres anchor: match what the empirical shendure obs fit used for this CT.
_PARAMS_PATH = DATA_ROOT / "shendure/shendure_obs_nb_phantom/by_cell_type_parameters.pkl"
with open(_PARAMS_PATH, "rb") as _f:
    _params = pickle.load(_f)
ANCHOR_N_CRES = len(_params.nb[ANCHOR_CT])
del _params

# Activity range: copy shendure power-sim defaults so results are comparable.
MINP = float(ANCHOR_BOUNDS.reference_activity)
MIN_ACT = float(ANCHOR_BOUNDS.min_mpra_umi)
MAX_ACT = 4.0 * MINP

# Compute budget: 50 reps x 5 sims = 250 realisations per grid point.
# "boil the seas" -- worker count bumped accordingly in _make_slurm_client.
N_SIMS = 5
N_LIBRARY_REPS = 50

# ---------------------------------------------------------------------------
# Axis grids (log-spaced where it makes sense)
# ---------------------------------------------------------------------------

AXES = {
    # bcs_per_cre: NB mean of the per-CRE barcode count draw. Shendure anchor
    # sits at ~136 (piggyBac + dense library). Bracket two decades around it.
    # mu_nb=3 will produce ~30% zero-barcode CREs, which cleanly represents a
    # low-complexity library (those CREs just drop out of downstream tests).
    "bcs_per_cre": np.array([3.0, 10.0, 30.0, 100.0, 300.0, 1000.0]),
    # n_cells: anchor ~8200. Span ~16x below to ~16x above.
    "n_cells":     np.array([500, 2000, 8000, 32000, 128000], dtype=int),
    # moi: plasmid-dose knob. Anchor ~18 (piggyBac). Bracket widely -- MOI<1
    # is sparse-transfection regime, MOI>>18 is collision-dominated.
    "moi":         np.array([0.5, 2.0, 6.0, 18.0, 54.0]),
    # n_cres: library breadth. Anchor 207. Spans from single-panel to genome-
    # scale at fixed library bcs_per_cre (so the total barcode count grows
    # with n_cres).
    "n_cres":      np.array([50, 200, 1000, 5000], dtype=int),
}


def predicted_collision_rate(n_cres: float, bcs_per_cre: float, moi: float) -> float:
    """Forward coupon-collector estimate of the per-event collision rate.

    Mirrors the post-hoc estimator in Bounds.from_ortho (core.py ~1212), which
    inverts observed unique-barcode counts to get total transfection events.
    Here we go forward: given a library of N unique barcodes and a mean of
    `moi` independent transfection draws per cell, the expected number of
    unique barcodes per cell is N*(1 - (1-1/N)**moi). Anything in excess of
    that count was a collision.

    Returns a fraction in [0, 1). Note that "high collision" is conventionally
    flagged at >2% (WARN_MULTI_TRANSFECTION_PERCENT in core.py).
    """
    N = float(n_cres) * float(bcs_per_cre)
    if N <= 0 or moi <= 0:
        return 0.0
    expected_unique = N * (1.0 - (1.0 - 1.0 / N) ** moi)
    collisions = max(0.0, moi - expected_unique)
    return collisions / moi


def collision_rate_for_point(axis: str, value) -> float:
    """Compute the predicted collision rate at a sweep point with everything
    else held at the anchor."""
    n_cres = int(value) if axis == "n_cres" else ANCHOR_N_CRES
    bpc = float(value) if axis == "bcs_per_cre" else ANCHOR_BCS_PER_CRE
    moi = float(value) if axis == "moi" else ANCHOR_MOI
    return predicted_collision_rate(n_cres, bpc, moi)


# Conventional "this experiment is collision-dominated" threshold; matches
# WARN_MULTI_TRANSFECTION_PERCENT in core.py.
COLLISION_WARN_PCT = 2.0


def make_bound_for_point(axis: str, value):
    """Return a mutated deep copy of ANCHOR_BOUNDS for one sweep point.

    Only the single axis under test is perturbed; everything else stays at
    the shendure Pluripotent anchor. n_cres is handled by the caller (it's a
    one_library_replicate kwarg, not a bounds field).
    """
    b = ANCHOR_BOUNDS.copy()
    # Collapse cells_per_cell_type to just the anchor CT so the sim runs a
    # single cell type (matches shendure power-sim pattern).
    b.cells_per_cell_type = ANCHOR_BOUNDS.cells_per_cell_type.loc[[ANCHOR_CT]].copy()

    if axis == "n_cells":
        b.cells_per_cell_type.loc[ANCHOR_CT] = int(value)
    elif axis == "bcs_per_cre":
        b.library_model.mu_nb = float(value)
        b.library_model.update_alt_nb_param()
    elif axis == "moi":
        b.set_effective_moi(float(value))
    elif axis == "n_cres":
        pass  # applied as one_library_replicate kwarg
    else:
        raise ValueError(f"Unknown axis: {axis!r}")
    return b


def _point_dirname(axis: str, value) -> str:
    # Keep filesystem-safe names even for floats.
    if isinstance(value, (int, np.integer)):
        return f"{axis}__{int(value)}"
    return f"{axis}__{float(value):g}".replace(".", "p")


# ---------------------------------------------------------------------------
# Phase: simulate
# ---------------------------------------------------------------------------

def phase_simulate(axis: str, client):
    grid = AXES[axis]
    sweep_dir = SIM_ROOT / axis
    sweep_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n=== axis: {axis}, grid: {list(grid)} ===", flush=True)
    print(f"  anchor n_cells={ANCHOR_N_CELLS}, bcs/cre={ANCHOR_BCS_PER_CRE:.3f}, "
          f"moi={ANCHOR_MOI:.3f}, n_cres={ANCHOR_N_CRES}", flush=True)

    for value in grid:
        pt_name = _point_dirname(axis, value)
        pt_dir = sweep_dir / pt_name
        pt_dir.mkdir(exist_ok=True)
        n_cres = int(value) if axis == "n_cres" else ANCHOR_N_CRES
        coll_pct = collision_rate_for_point(axis, value) * 100.0
        print(f"\n-- {pt_name}  (n_cres={n_cres}, predicted collision={coll_pct:.2f}%) --", flush=True)

        try:
            import gc
            bound = make_bound_for_point(axis, value)

            # A sim is "done" iff it has ttest results on disk. Anything else
            # (including dirs from interrupted prior runs that wrote
            # state.parquet but not the transcription parquets) is treated as
            # absent. We never load done sims eagerly -- holding 50
            # de_novo_simulation objects in driver memory OOMs the coordinator
            # for high-volume grid points. The aggregator re-reads from disk.
            n_done = 0
            for d in sorted(pt_dir.iterdir()):
                if not d.is_dir():
                    continue
                tt_dir = d / "tests" / "hs_act" / "ttest"
                if tt_dir.exists() and any(tt_dir.iterdir()):
                    n_done += 1
            print(f"   on-disk: {n_done} done", flush=True)

            n_remaining = N_LIBRARY_REPS - n_done
            hs = None  # built lazily from the first sim's scmpra
            for i in range(n_remaining):
                _, sim = scm.one_library_replicate(
                    root=pt_dir,
                    n_sims=N_SIMS,
                    client=client,
                    flatten_overtransfection=True,
                    bound=bound,
                    n_cres=n_cres,
                    min=MIN_ACT,
                    max=MAX_ACT,
                    minP=MINP,
                    cell_type=ANCHOR_CT,
                )
                # save() blocks on the transcription dask futures (gamut just
                # submits them) and writes the .scmpra parquet trees that
                # from_parquet expects. MUST come before any read.
                sim.save()
                if hs is None:
                    example = scm.scMPRA_data.from_parquet(sim.scmpradatp / "0.scmpra")
                    hs = scm.make_all_by_celltype_hypotheses(
                        counts=example, reference_cre="reference"
                    )
                sim.add_hypothesis_set("hs_act", hs)
                sim.ttest("hs_act")
                sim.save()
                print(f"   rep {i+1}/{n_remaining} done", flush=True)
                # Drop sim refs so the next iteration starts clean.
                del sim
                gc.collect()
            print(f"   axis point {pt_name} complete", flush=True)
        except Exception as e:
            # Ragged-edge tolerance: log and move on so one bad point can't
            # kill an entire axis. The aggregator will simply skip empty
            # point directories.
            print(f"   FAILED at {pt_name}: {type(e).__name__}: {e}", flush=True)
            import traceback; traceback.print_exc()


# ---------------------------------------------------------------------------
# Phase: plot
# ---------------------------------------------------------------------------

def _aggregate(axis: str, client, test_type: str = "ttest") -> pd.DataFrame:
    """Walk scratch for this axis, load each surviving sim defensively, and
    build a (reject_null, fc, value) DataFrame for plotting.

    The 2026-04-11 rsync move from project -> scratch lost files from many
    sim dirs in non-uniform ways: some are missing ground_truth.tsv.gz,
    some are missing futures.tsv.gz, some have partial ttest/*_results.tsv.
    This function treats every file read as fallible and skips any sim that
    can't be fully consumed, reporting a count at the end.
    """
    grid = AXES[axis]
    sweep_dir = SIM_ROOT / axis
    rows = []
    skipped_load = 0
    skipped_merge = 0
    for value in grid:
        pt_dir = sweep_dir / _point_dirname(axis, value)
        if not pt_dir.exists():
            continue
        for d in sorted(pt_dir.iterdir()):
            if not d.is_dir():
                continue
            res_dir = d / "tests" / "hs_act" / test_type
            if not (res_dir.exists() and any(res_dir.iterdir())):
                continue
            # Pre-check files that de_novo_simulation.__init__ loads eagerly
            # and that _merge_in_ground_truth needs per-index.
            required = [
                d / "ground_truth.tsv.gz",
                d / "futures.tsv.gz",
                d / "state.parquet",
            ] + [res_dir / f"{i}_results.tsv" for i in range(N_SIMS)]
            if not all(p.exists() for p in required):
                skipped_load += 1
                continue
            try:
                sim = scm.de_novo_simulation(location=pt_dir, name=d.name, client=client)
            except Exception as e:
                print(f"  skip load {d.name}: {type(e).__name__}: {e}", flush=True)
                skipped_load += 1
                continue
            # Per-sim merge (replicates scm.sum_pow inline so a failure on one
            # sim does not poison the whole grid point).
            try:
                reps = sim.get_state_field("n_sims")
                for i in range(reps):
                    mergy = sim._merge_in_ground_truth(
                        hypothesis_set_name="hs_act", test_type=test_type, index=i
                    )
                    mergy["comparison_truth"] = mergy["comparison_truth"].astype(float)
                    mergy["reference_truth"] = mergy["reference_truth"].astype(float)
                    mergy["fc"] = mergy["comparison_truth"] / mergy["reference_truth"]
                    mergy = mergy[["reject_null", "fc"]].copy()
                    mergy["axis"] = axis
                    mergy["value"] = float(value)
                    rows.append(mergy)
            except Exception as e:
                print(f"  skip merge {d.name}: {type(e).__name__}: {e}", flush=True)
                skipped_merge += 1
                continue
    if skipped_load or skipped_merge:
        print(
            f"  aggregate ({test_type}): skipped {skipped_load} un-loadable + "
            f"{skipped_merge} un-mergeable sim dirs",
            flush=True,
        )
    if not rows:
        raise RuntimeError(f"No usable sims found for axis {axis!r} with test_type={test_type!r}")
    return pd.concat(rows, ignore_index=True)


def phase_plot(axis: str, client, test_type: str = "ttest"):
    suffix = f"_{test_type}" if test_type != "ttest" else ""
    combined = _aggregate(axis, client, test_type=test_type)
    sweep_out = SWEEPS_DIR / axis
    sweep_out.mkdir(exist_ok=True)
    combined.to_parquet(sweep_out / f"power_df{suffix}.parquet")
    print(f"Saved: {sweep_out / f'power_df{suffix}.parquet'}", flush=True)

    # Power curves per grid point
    fig, ax = plt.subplots(figsize=(6, 4))
    cmap = plt.get_cmap("viridis")
    values = sorted(combined["value"].unique())
    for i, v in enumerate(values):
        sub = combined[combined["value"] == v].copy()
        sub = sub[sub["fc"].notna()]
        sub["bin"] = pd.cut(sub["fc"], bins=40)
        binned = (
            sub.groupby("bin", observed=True)
            .agg(rej=("reject_null", "mean"), mid=("fc", "mean"))
            .reset_index()
            .dropna()
        )
        coll = collision_rate_for_point(axis, v) * 100.0
        # Mark collision-dominated points with a dashed linestyle.
        ls = "--" if coll > COLLISION_WARN_PCT else "-"
        ax.plot(binned["mid"], binned["rej"],
                marker="o", markersize=2, lw=1, linestyle=ls,
                color=cmap(i / max(1, len(values) - 1)),
                label=f"{axis}={v:g}  (coll={coll:.2f}%)")
    ax.axhline(0.8, color="k", ls="--", lw=0.6)
    ax.axvline(1.0, color="grey", ls=":", lw=0.6)
    ax.set_xlabel("FC (activity / minP)")
    ax.set_ylabel("power (reject rate)")
    ax.set_ylim(0, 1)
    ax.set_title(f"Design-space sweep: {axis} ({test_type})\nshendure anchor, {ANCHOR_CT} "
                 f"(dashed = predicted collision >{COLLISION_WARN_PCT}%)")
    ax.legend(fontsize=7, loc="best")
    plt.tight_layout()
    fig.savefig(sweep_out / f"power_curves{suffix}.svg", format="svg", bbox_inches="tight")
    print(f"Saved: {sweep_out / f'power_curves{suffix}.svg'}", flush=True)

    # Summary: power@FC=2 per grid point
    def _power_at_fc(df, target=2.0, half=0.15):
        sub = df[(df["fc"] >= target - half) & (df["fc"] <= target + half)]
        if len(sub) == 0:
            return np.nan
        return float(sub["reject_null"].mean())

    summary = (
        combined.groupby("value")
        .apply(_power_at_fc)
        .reset_index(name="power_at_fc2")
    )
    summary["collision_pct"] = [
        collision_rate_for_point(axis, v) * 100.0 for v in summary["value"]
    ]
    summary.to_parquet(sweep_out / f"power_at_fc2{suffix}.parquet")
    print(f"\nPower@FC=2 by {axis}:")
    print(summary.to_string(index=False))

    # Two-panel summary: power@FC=2 (top) and predicted collision% (bottom),
    # sharing the x axis. Collision-dominated zone shaded.
    fig, (ax_pow, ax_col) = plt.subplots(
        2, 1, figsize=(5, 5), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )

    safe = summary["collision_pct"] <= COLLISION_WARN_PCT
    ax_pow.plot(summary["value"], summary["power_at_fc2"],
                color="steelblue", lw=1)
    ax_pow.scatter(summary["value"][safe], summary["power_at_fc2"][safe],
                   color="steelblue", s=40, label="ok")
    ax_pow.scatter(summary["value"][~safe], summary["power_at_fc2"][~safe],
                   facecolor="white", edgecolor="steelblue", s=40,
                   label=f"collision >{COLLISION_WARN_PCT}%")
    ax_pow.axhline(0.8, color="k", ls="--", lw=0.6)
    ax_pow.set_ylabel("power @ FC=2")
    ax_pow.set_ylim(0, 1)
    ax_pow.set_title(f"power@FC=2 vs {axis}  (shendure anchor)")
    ax_pow.legend(fontsize=7, loc="best")

    ax_col.plot(summary["value"], summary["collision_pct"],
                color="firebrick", marker="o", lw=1)
    ax_col.axhline(COLLISION_WARN_PCT, color="k", ls="--", lw=0.6)
    ax_col.set_ylabel("predicted\ncollision %")
    ax_col.set_xlabel(axis)
    if axis in ("bcs_per_cre", "n_cells", "moi"):
        ax_col.set_xscale("log")
    plt.tight_layout()
    fig.savefig(sweep_out / f"power_at_fc2{suffix}.svg", format="svg", bbox_inches="tight")
    print(f"Saved: {sweep_out / f'power_at_fc2{suffix}.svg'}", flush=True)


# ---------------------------------------------------------------------------
# Anchor snapshot (write once, on any invocation)
# ---------------------------------------------------------------------------

def write_anchor():
    anchor = dict(
        bounds="SHENDURE_BOUNDS",
        cell_type=ANCHOR_CT,
        n_cells=ANCHOR_N_CELLS,
        bcs_per_cre=ANCHOR_BCS_PER_CRE,
        moi=ANCHOR_MOI,
        n_cres=ANCHOR_N_CRES,
        minP=MINP,
        min_activity=MIN_ACT,
        max_activity=MAX_ACT,
        n_sims=N_SIMS,
        n_library_reps=N_LIBRARY_REPS,
        axes={k: list(map(float, v)) for k, v in AXES.items()},
    )
    with open(OUT / "anchor.json", "w") as f:
        json.dump(anchor, f, indent=2)
    print(f"Wrote anchor: {OUT / 'anchor.json'}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _make_slurm_client():
    cluster = SLURMCluster(
        cores=1,
        memory="24G",
        processes=1,
        env_extra=[f"export PYTHONPATH={_repo_root}:$PYTHONPATH"],
        job_extra_directives=[
            "-p priority",
            "--account=prio_skr2",
            "--job-name=dsweep_worker",
            "--time=8:00:00",
            "--output=worker_%j.out",
        ],
    )
    cluster.scale(jobs=60)
    client = Client(cluster, timeout="300s", heartbeat_interval="20s")
    print(f"Dask dashboard: {client.dashboard_link}", flush=True)
    return cluster, client


def _make_local_client():
    import psutil
    mem_mb = os.environ.get("SLURM_MEM_PER_NODE")
    if mem_mb:
        mem_limit = f"{int(int(mem_mb) * 0.9)}MB"
    else:
        mem_limit = f"{int(psutil.virtual_memory().total * 0.9 / 1e9)}GB"
    cluster = LocalCluster(
        n_workers=1, threads_per_worker=2, processes=False, memory_limit=mem_limit
    )
    client = Client(cluster)
    return cluster, client


if __name__ == "__main__":
    phase = sys.argv[1] if len(sys.argv) > 1 else "all"
    axis = sys.argv[2] if len(sys.argv) > 2 else "bcs_per_cre"
    # Optional 3rd arg: test type for plot phase (ttest or mwu)
    test_type = sys.argv[3] if len(sys.argv) > 3 else "ttest"
    if axis not in AXES:
        print(f"Unknown axis {axis!r}. Choose from: {list(AXES)}")
        sys.exit(1)

    write_anchor()

    if phase in ("simulate", "all"):
        print(f"\n{'='*60}\nPhase: simulate ({axis})\n{'='*60}", flush=True)
        cluster, client = _make_slurm_client()
        try:
            phase_simulate(axis, client)
        finally:
            client.close()
            cluster.close()

    if phase in ("plot", "all"):
        print(f"\n{'='*60}\nPhase: plot ({axis}, {test_type})\n{'='*60}", flush=True)
        cluster, client = _make_local_client()
        try:
            phase_plot(axis, client, test_type=test_type)
        finally:
            client.close()
            cluster.close()

    print("\nDone.", flush=True)
