#!/usr/bin/env python3
"""
Cohen null calibration -- t-test

Measures per-anchor false positive rate of Welch's t-test under the null
(two independent CREs drawn at the same baseline mu) across the Cohen
(episomal, retinal) dataset's activity range. No fold change offset, no
multiple testing correction -- raw p<0.05 rejection rate per anchor baseline.

Usage:
    python cohen_pairwise_null_calibration_ttest.py [phase]

Phases:
    simulate  -- generate simulations and run t-test
    plot      -- aggregate results and produce FPR calibration plots
    all       -- run both phases sequentially (default)
"""

import sys
import os
import json
import pickle
import uuid
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

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
SIM_DIR = DATA_ROOT / "simulated" / f"{SIM_DATE}_cohen_null_ttest"

ORTHO_DIR = DATA_ROOT / "cohen" / "cohen_obsingle_nb_phantom"
with open(ORTHO_DIR / "by_cell_type_parameters.pkl", "rb") as f:
    _params = pickle.load(f)

CELL_TYPES = sorted(_params.nb.keys())
N_CRES = max(len(_params.nb[ct]) for ct in CELL_TYPES)
_all_mus = np.concatenate([_params.nb[ct]["mu"].values for ct in CELL_TYPES])
MU_P5 = float(np.percentile(_all_mus, 5))
MU_P95 = float(np.percentile(_all_mus, 95))
del _params, _all_mus

MINP = scm.COHEN_BOUNDS.reference_activity

N_TOTAL_ANCHORS = 24
ANCHOR_BASELINES = np.geomspace(MU_P5, MU_P95, N_TOTAL_ANCHORS)

CRES_PER_ANCHOR = 2  # one null pair per anchor
N_LIBRARY_REPS = 20
N_SIMS = 5
ALPHA = 0.05


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def build_ground_truth(cell_types, n_cres, anchor_baselines, minP):
    """Per-anchor null pairs (both CREs at the same mu) + fillers + reference."""
    names = []
    mus = []

    for a_idx, baseline in enumerate(anchor_baselines):
        names.append(f"anchor_{a_idx}_a")
        mus.append(baseline)
        names.append(f"anchor_{a_idx}_b")
        mus.append(baseline)

    n_used = len(names) + 1  # +1 for reference
    n_filler = n_cres - n_used
    for j in range(n_filler):
        names.append(f"filler_{j}")
        mus.append(minP)

    names.append("reference")
    mus.append(minP)

    assert len(names) == n_cres, f"Expected {n_cres} CREs, got {len(names)}"

    return pd.concat([
        pd.DataFrame({"cre_id": names, "mu": mus, "cell_type": ct})
        for ct in cell_types
    ], ignore_index=True)


def build_null_hypotheses(cell_types, n_anchors):
    rows = []
    for ct in cell_types:
        for a_idx in range(n_anchors):
            rows.append({
                "comparison_CRE": f"anchor_{a_idx}_a",
                "comparison_cell_type": ct,
                "reference_CRE": f"anchor_{a_idx}_b",
                "reference_cell_type": ct,
                "meta": f"null_anchor_{a_idx}",
            })
    return scm.HypothesisSet.from_dataframe(pd.DataFrame(rows))


def _sim_is_complete(sim_dir):
    tt_dir = sim_dir / "tests" / "hs_null" / "ttest"
    return tt_dir.exists() and any(tt_dir.iterdir())


# ---------------------------------------------------------------------------
# Phase: simulate
# ---------------------------------------------------------------------------

def phase_simulate(client):
    print(f"minP = {MINP:.6f}")
    print(f"n_cres = {N_CRES}")
    print(f"Cell types: {CELL_TYPES}")
    print(f"mu p5 = {MU_P5:.4f}, p95 = {MU_P95:.4f}")
    print(f"Anchor baselines ({N_TOTAL_ANCHORS}): "
          f"{[f'{x:.4f}' for x in ANCHOR_BASELINES]}")
    print(f"Library reps: {N_LIBRARY_REPS}, sims per rep: {N_SIMS}")
    print(f"Total sims: {N_LIBRARY_REPS}", flush=True)

    hs = build_null_hypotheses(CELL_TYPES, N_TOTAL_ANCHORS)
    cohort_dir = SIM_DIR / "cohort_0"

    existing_complete = []
    if cohort_dir.exists():
        for d in sorted(cohort_dir.iterdir()):
            if d.is_dir() and _sim_is_complete(d):
                existing_complete.append(d.name)
    print(f"Resumed {len(existing_complete)} complete sims.", flush=True)

    n_remaining = N_LIBRARY_REPS - len(existing_complete)
    print(f"Running {n_remaining} fresh sims.", flush=True)

    for lib_rep in range(n_remaining):
        gt_df = build_ground_truth(
            CELL_TYPES, N_CRES, ANCHOR_BASELINES, MINP
        )

        libraries = [
            scm.simulate_library(
                CREs=pd.Series(gt_df["cre_id"].unique()),
                library_model=scm.COHEN_BOUNDS.library_model,
            )
            for _ in range(N_SIMS)
        ]

        sim_name = f"sim_{uuid.uuid4().hex[:8]}"
        sim = scm.de_novo_simulation(
            location=cohort_dir,
            name=sim_name,
            client=client,
            libraries=libraries,
            library_mapping="corresponding",
            flatten_overtransfection=True,
            n_sims=N_SIMS,
            experiment_bounds=scm.COHEN_BOUNDS,
            ground_truth=gt_df,
        )
        sim.gamut()
        sim.add_hypothesis_set("hs_null", hs)
        sim.ttest("hs_null")
        sim.save()
        print(f"  rep {lib_rep+1}/{n_remaining}: {sim_name}", flush=True)


# ---------------------------------------------------------------------------
# Phase: plot
# ---------------------------------------------------------------------------

def phase_plot(client):
    sim_registry = []
    cohort_dir = SIM_DIR / "cohort_0"
    if cohort_dir.exists():
        for d in sorted(cohort_dir.iterdir()):
            if d.is_dir() and _sim_is_complete(d):
                sim_registry.append((cohort_dir, d.name))
    print(f"Loaded {len(sim_registry)} sims from disk.", flush=True)

    with open(OUTPUT_DIR / "cohen_null_ttest_sim_registry.json", "w") as f:
        json.dump([(str(p), n) for p, n in sim_registry], f, indent=2)

    null_rows = []
    for cd, sim_name in sim_registry:
        sim = scm.de_novo_simulation(location=cd, name=sim_name, client=client)
        gt = sim.ground_truth

        for i in range(sim.get_state_field("n_sims")):
            m = sim._merge_in_ground_truth(
                hypothesis_set_name="hs_null", test_type="ttest", index=i
            )
            for _, row in m.iterrows():
                cre_name = row["comparison_CRE"]
                ct = row["comparison_cell_type"]
                if not cre_name.startswith("anchor_"):
                    continue
                a_idx = int(cre_name.split("_")[1])
                anchor_mu = gt[
                    (gt["cre_id"] == cre_name)
                    & (gt["cell_type"] == ct)
                ]["mu"].iloc[0]
                null_rows.append({
                    "cell_type": ct,
                    "anchor_idx": a_idx,
                    "baseline_mu": anchor_mu,
                    "p_value": row["p_value"],
                    "reject_null": row["p_value"] < ALPHA,
                })

    null_df = pd.DataFrame(null_rows)
    null_df.to_parquet(OUTPUT_DIR / "cohen_null_ttest_df.parquet", index=False)
    print(f"Saved null_df ({len(null_df)}).", flush=True)

    fpr_grid = (
        null_df
        .groupby(["cell_type", "anchor_idx", "baseline_mu"])
        .agg(fpr=("reject_null", "mean"),
             n_tests=("reject_null", "count"))
        .reset_index()
    )
    fpr_grid.to_parquet(
        OUTPUT_DIR / "cohen_null_ttest_fpr_grid.parquet", index=False
    )
    print(f"Grid cells: {len(fpr_grid)}", flush=True)

    for ct in CELL_TYPES:
        ct_data = null_df[null_df["cell_type"] == ct]
        if len(ct_data) > 0:
            print(f"  {ct}: overall FPR@{ALPHA} = "
                  f"{ct_data['reject_null'].mean():.4f} "
                  f"(n={len(ct_data)})", flush=True)

    # --- Heatmap: rows=cell_types, cols=anchors, color = FPR (centered on alpha) ---
    pivot = fpr_grid.pivot_table(
        index="cell_type", columns="anchor_idx", values="fpr"
    ).reindex(CELL_TYPES)

    fig, ax = plt.subplots(figsize=(14, max(2.5, 0.6 * len(CELL_TYPES) + 1)))
    vmax = max(0.10, float(np.nanmax(pivot.values)) if pivot.size else ALPHA)
    sns.heatmap(
        pivot, ax=ax,
        cmap="RdYlGn_r", center=ALPHA, vmin=0, vmax=vmax,
        linewidths=0.3, linecolor="white",
        annot=True, fmt=".3f", annot_kws={"size": 7},
        cbar_kws={"label": f"FPR (alpha={ALPHA})"},
        xticklabels=[f"{ANCHOR_BASELINES[i]:.3f}" for i in pivot.columns],
        yticklabels=[
            ("Rod" if ct == "reference" else ct) for ct in pivot.index
        ],
    )
    ax.set_xlabel("Baseline activity (mu)", fontsize=10)
    ax.set_ylabel("Cell type", fontsize=10)
    ax.set_title(
        f"Cohen null FPR calibration -- Welch t-test (alpha={ALPHA})\n"
        f"green = well-calibrated, red = inflated FPR "
        f"({N_LIBRARY_REPS}x{N_SIMS} null pairs per cell)",
        fontsize=11,
    )
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    plt.tight_layout()
    svg_path = OUTPUT_DIR / "cohen_null_ttest_fpr_heatmap.svg"
    fig.savefig(svg_path, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {svg_path}", flush=True)

    # --- Line plot: FPR vs baseline mu, one line per cell type ---
    fig, ax = plt.subplots(figsize=(9, 5))
    palette = sns.color_palette("tab10", n_colors=len(CELL_TYPES))
    for i, ct in enumerate(CELL_TYPES):
        ct_data = fpr_grid[fpr_grid["cell_type"] == ct].sort_values(
            "baseline_mu"
        )
        label = "Rod" if ct == "reference" else ct
        ax.plot(ct_data["baseline_mu"], ct_data["fpr"], "o-",
                color=palette[i], label=label, markersize=4, lw=1.5)
    ax.axhline(y=ALPHA, color="black", linestyle="--", lw=1, alpha=0.7,
               label=f"alpha={ALPHA}")
    ax.set_xscale("log")
    ax.set_xlabel("Baseline activity (mu)", fontsize=11)
    ax.set_ylabel(f"FPR at p<{ALPHA}", fontsize=11)
    ax.set_title(
        f"Cohen null FPR calibration -- Welch t-test", fontsize=12,
    )
    ax.legend(fontsize=8, loc="upper right", ncol=2)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    svg_path2 = OUTPUT_DIR / "cohen_null_ttest_fpr_lines.svg"
    fig.savefig(svg_path2, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {svg_path2}", flush=True)


# ---------------------------------------------------------------------------
# Main dispatch
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
                "--job-name=cohen_null_worker",
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
        print(f"Unknown phase: {phase!r}. Choose from: simulate, plot, all")
        sys.exit(1)

    print("\nDone.", flush=True)
