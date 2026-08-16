#!/usr/bin/env python3
"""
Shendure pairwise CRE power analysis -- MWU

Measures MWU test power for pairwise CRE comparisons as a function of
(baseline activity, log2 fold change) in the Shendure (integrating, iPSC)
dataset.

Shendure is integrating: each cell type has its own fitted CRE set. Each cell
type is simulated independently with its own n_cres and cells_per_cell_type.
Cell types with fewer CREs are split into multiple cohorts that tile the
anchor space; results are aggregated across cohorts.

Usage:
    python shendure_pairwise_power_mwu.py [phase]

Phases:
    simulate  -- generate simulations and run MWU test
    plot      -- aggregate results and produce SVG plots
    all       -- run both phases sequentially (default)
"""

import sys
import os
import json
import pickle
import uuid
import hashlib
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
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import scMPRAforge as scm
import kircher_reference
from dask.distributed import Client, LocalCluster
from dask_jobqueue import SLURMCluster

# ---------------------------------------------------------------------------
# Paths and parameters
# ---------------------------------------------------------------------------

DATA_ROOT = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new")
OUTPUT_DIR = Path(__file__).resolve().parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

SIM_DATE = "2026-08-13"
SIM_DIR = DATA_ROOT / "simulated" / f"{SIM_DATE}_shendure_pw"

ORTHO_DIR = DATA_ROOT / "shendure" / "shendure_obs_nb"
with open(ORTHO_DIR / "by_cell_type_parameters.pkl", "rb") as f:
    _params = pickle.load(f)

CELL_TYPES = sorted(_params.nb.keys())
N_CRES_PER_CT = {ct: len(_params.nb[ct]) for ct in CELL_TYPES}
_all_mus = np.concatenate([_params.nb[ct]["mu"].values for ct in CELL_TYPES])
MU_P5 = float(np.percentile(_all_mus, 5))
MU_P95 = float(np.percentile(_all_mus, 95))
del _params, _all_mus

MINP = scm.SHENDURE_BOUNDS.reference_activity

SIM_SEED = 20260813
N_TOTAL_ANCHORS = 24
ANCHOR_BASELINES = np.geomspace(MU_P5, MU_P95, N_TOTAL_ANCHORS)
FC_LOG2_OFFSETS = np.array([0.02, 0.03, 0.05, 0.07, 0.10, 0.15, 0.20, 0.30,
                            0.40, 0.45])

CRES_PER_ANCHOR = 1 + len(FC_LOG2_OFFSETS)
N_LIBRARY_REPS = 20
N_SIMS = 5

N_FILLER_NULL = 5


def _cohorts_for_ct(n_cres):
    """Round-robin tile the 24 anchors into cohorts that fit in n_cres."""
    max_anchors = (n_cres - 1) // CRES_PER_ANCHOR  # -1 for reference
    n_cohorts = int(np.ceil(N_TOTAL_ANCHORS / max_anchors))
    cohorts = [[] for _ in range(n_cohorts)]
    for i in range(N_TOTAL_ANCHORS):
        cohorts[i % n_cohorts].append(i)
    return cohorts


COHORTS_BY_CT = {ct: _cohorts_for_ct(N_CRES_PER_CT[ct]) for ct in CELL_TYPES}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def build_ground_truth(cell_type, anchor_indices, n_cres,
                       anchor_baselines, fc_log2_offsets, minP):
    """Build ground truth DataFrame for one (cell_type, cohort)."""
    names = []
    mus = []

    for a_idx in anchor_indices:
        baseline = anchor_baselines[a_idx]
        names.append(f"anchor_{a_idx}")
        mus.append(baseline)
        for fc_l2 in fc_log2_offsets:
            names.append(f"anchor_{a_idx}_fc_{fc_l2:.2f}")
            mus.append(baseline * 2**fc_l2)

    n_used = len(names) + 1
    n_filler = n_cres - n_used
    for j in range(n_filler):
        names.append(f"filler_{j}")
        mus.append(minP)

    names.append("reference")
    mus.append(minP)

    assert len(names) == n_cres, f"Expected {n_cres} CREs, got {len(names)}"

    return pd.DataFrame({
        "cre_id": names, "mu": mus, "cell_type": cell_type,
    })


def build_pairwise_hypotheses(anchor_indices, cell_type, fc_log2_offsets,
                              n_filler_null=N_FILLER_NULL):
    rows = []

    for a_idx in anchor_indices:
        anchor_name = f"anchor_{a_idx}"
        for fc_l2 in fc_log2_offsets:
            variant_name = f"anchor_{a_idx}_fc_{fc_l2:.2f}"
            rows.append({
                "comparison_CRE": variant_name,
                "comparison_cell_type": cell_type,
                "reference_CRE": anchor_name,
                "reference_cell_type": cell_type,
                "meta": f"power_fc_{fc_l2:.2f}",
            })

    for j in range(0, n_filler_null * 2, 2):
        rows.append({
            "comparison_CRE": f"filler_{j}",
            "comparison_cell_type": cell_type,
            "reference_CRE": f"filler_{j+1}",
            "reference_cell_type": cell_type,
            "meta": "null_check",
        })

    return scm.HypothesisSet.from_dataframe(pd.DataFrame(rows))


def _sim_is_complete(sim_dir):
    mwu_dir = sim_dir / "tests" / "hs_pairwise" / "mwu"
    return mwu_dir.exists() and any(mwu_dir.iterdir())


# ---------------------------------------------------------------------------
# Phase: simulate
# ---------------------------------------------------------------------------

def phase_simulate(client):
    print(f"minP = {MINP:.6f}")
    print(f"mu p5 = {MU_P5:.4f}, p95 = {MU_P95:.4f}")
    print(f"Anchor baselines ({N_TOTAL_ANCHORS}): "
          f"{[f'{x:.4f}' for x in ANCHOR_BASELINES]}")
    print(f"FC offsets (log2): {list(FC_LOG2_OFFSETS)}")
    print(f"Library reps: {N_LIBRARY_REPS}, sims per rep: {N_SIMS}")
    print()
    for ct in CELL_TYPES:
        cohorts = COHORTS_BY_CT[ct]
        cells = scm.SHENDURE_BOUNDS.cells_per_cell_type.get(ct, "?")
        print(f"  {ct}: n_cres={N_CRES_PER_CT[ct]}, cells={cells}, "
              f"cohorts={len(cohorts)}, anchors/cohort="
              f"{[len(c) for c in cohorts]}")
    total_sims = sum(
        len(COHORTS_BY_CT[ct]) for ct in CELL_TYPES
    ) * N_LIBRARY_REPS
    print(f"\nTotal sims: {total_sims}", flush=True)

    for ct in CELL_TYPES:
        cohorts = COHORTS_BY_CT[ct]
        n_cres = N_CRES_PER_CT[ct]

        bound_ct = scm.SHENDURE_BOUNDS.copy()
        bound_ct.cells_per_cell_type = (
            scm.SHENDURE_BOUNDS.cells_per_cell_type.loc[[ct]]
        )

        for cohort_idx, anchor_indices in enumerate(cohorts):
            print(f"\n=== {ct} cohort {cohort_idx}/{len(cohorts)-1} "
                  f"({len(anchor_indices)} anchors) ===", flush=True)

            hs = build_pairwise_hypotheses(
                anchor_indices, ct, FC_LOG2_OFFSETS
            )

            cohort_dir = SIM_DIR / ct / f"cohort_{cohort_idx}"

            existing_complete = []
            if cohort_dir.exists():
                for d in sorted(cohort_dir.iterdir()):
                    if d.is_dir() and _sim_is_complete(d):
                        existing_complete.append(d.name)
            print(f"  Resumed {len(existing_complete)} complete sims.",
                  flush=True)

            n_remaining = N_LIBRARY_REPS - len(existing_complete)
            print(f"  Running {n_remaining} fresh sims.", flush=True)

            for lib_rep in range(n_remaining):
                rep_idx = len(existing_complete) + lib_rep
                tag = f"{SIM_SEED}:{ct}:{cohort_idx}:{rep_idx}"
                np.random.seed(int(hashlib.sha256(tag.encode()).hexdigest()[:8], 16))
                gt_df = build_ground_truth(
                    ct, anchor_indices, n_cres,
                    ANCHOR_BASELINES, FC_LOG2_OFFSETS, MINP
                )

                libraries = [
                    scm.simulate_library(
                        CREs=pd.Series(gt_df["cre_id"].values),
                        library_model=scm.SHENDURE_BOUNDS.library_model,
                    )
                    for _ in range(N_SIMS)
                ]

                sim_name = "sim_" + hashlib.sha256(tag.encode()).hexdigest()[:8]
                sim = scm.de_novo_simulation(
                    location=cohort_dir,
                    name=sim_name,
                    client=client,
                    libraries=libraries,
                    library_mapping="corresponding",
                    flatten_overtransfection=True,
                    n_sims=N_SIMS,
                    experiment_bounds=bound_ct,
                    ground_truth=gt_df,
                )
                sim.gamut()
                sim.add_hypothesis_set("hs_pairwise", hs)
                sim.mwu("hs_pairwise")
                sim.save()
                print(f"    rep {lib_rep+1}/{n_remaining}: {sim_name}",
                      flush=True)

            print(f"  {ct} cohort {cohort_idx} done.", flush=True)


# ---------------------------------------------------------------------------
# Phase: plot
# ---------------------------------------------------------------------------

def phase_plot(client):
    sim_registry = []  # (cell_type, cohort_idx, cohort_dir, sim_name)
    for ct in CELL_TYPES:
        cohorts = COHORTS_BY_CT[ct]
        for cohort_idx in range(len(cohorts)):
            cohort_dir = SIM_DIR / ct / f"cohort_{cohort_idx}"
            if not cohort_dir.exists():
                continue
            for d in sorted(cohort_dir.iterdir()):
                if d.is_dir() and _sim_is_complete(d):
                    sim_registry.append((ct, cohort_idx, cohort_dir, d.name))
    print(f"Loaded {len(sim_registry)} sims from disk.", flush=True)

    with open(OUTPUT_DIR / "shendure_pairwise_power_sim_registry.json",
              "w") as f:
        json.dump(
            [(ct, ci, str(p), n) for ct, ci, p, n in sim_registry],
            f, indent=2,
        )

    power_rows = []
    null_rows = []

    for ct, cohort_idx, cohort_dir, sim_name in sim_registry:
        sim = scm.de_novo_simulation(
            location=cohort_dir, name=sim_name, client=client
        )
        gt = sim.ground_truth

        for i in range(sim.get_state_field("n_sims")):
            m = sim._merge_in_ground_truth(
                hypothesis_set_name="hs_pairwise", test_type="mwu", index=i
            )

            for _, row in m.iterrows():
                cre_name = row["comparison_CRE"]

                if "_fc_" in cre_name:
                    parts = cre_name.split("_fc_")
                    anchor_name = parts[0]
                    log2_fc = float(parts[1])
                    anchor_mu = gt[
                        gt["cre_id"] == anchor_name
                    ]["mu"].iloc[0]

                    power_rows.append({
                        "cell_type": ct,
                        "baseline_mu": anchor_mu,
                        "log2_fc": log2_fc,
                        "p_value": row["p_value"],
                        "reject_null": row["p_value"] < 0.05,
                    })
                elif cre_name.startswith("filler_"):
                    null_rows.append({
                        "cell_type": ct,
                        "p_value": row["p_value"],
                    })

    power_df = pd.DataFrame(power_rows)
    null_df = pd.DataFrame(null_rows)

    power_df.to_parquet(
        OUTPUT_DIR / "shendure_pairwise_power_df.parquet", index=False
    )
    null_df.to_parquet(
        OUTPUT_DIR / "shendure_pairwise_null_df.parquet", index=False
    )
    print(f"Saved power_df ({len(power_df)}) and null_df ({len(null_df)}).",
          flush=True)

    for ct in CELL_TYPES:
        ct_null = null_df[null_df["cell_type"] == ct]
        if len(ct_null) > 0:
            fpr = ct_null["p_value"].lt(0.05).mean()
            print(f"  {ct}: null FPR@0.05 = {fpr:.4f} (n={len(ct_null)})",
                  flush=True)

    power_grid = (
        power_df
        .groupby(["cell_type", "baseline_mu", "log2_fc"])
        .agg(power=("reject_null", "mean"),
             n_tests=("reject_null", "count"))
        .reset_index()
    )
    print(f"Grid cells: {len(power_grid)}", flush=True)

    # --- Heatmap grid ---
    ncols = 2
    nrows = int(np.ceil(len(CELL_TYPES) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 5 * nrows))
    axes = np.atleast_1d(axes).flatten()

    for i, ct in enumerate(CELL_TYPES):
        ct_display = "Pluripotent" if ct == "reference" else ct
        ax = axes[i]
        ct_data = power_grid[power_grid["cell_type"] == ct]

        pivot = ct_data.pivot_table(
            index="log2_fc", columns="baseline_mu", values="power"
        )
        pivot = pivot.sort_index(ascending=False)
        pivot = pivot[sorted(pivot.columns)]

        cells = scm.SHENDURE_BOUNDS.cells_per_cell_type.get(ct, "?")
        n_cres = N_CRES_PER_CT[ct]

        sns.heatmap(
            pivot, ax=ax, vmin=0, vmax=1,
            cmap="YlOrRd", linewidths=0.5, linecolor="white",
            cbar_kws={"label": "Power", "shrink": 0.8},
            xticklabels=[f"{x:.3f}" for x in sorted(pivot.columns)],
            yticklabels=[f"{y:.2f}" for y in pivot.index],
        )

        kircher_reference.annotate(ax, pivot.index)

        ax.set_title(f"{ct_display}\n(n_cres={n_cres}, cells={cells})",
                     fontsize=10)
        ax.set_xlabel("Baseline activity (mu)", fontsize=9)
        ax.set_ylabel("|log2 FC|", fontsize=9)
        ax.tick_params(labelsize=7)
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

    for j in range(len(CELL_TYPES), len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(
        "MWU pairwise power -- Shendure\n"
        f"(alpha=0.05, {N_LIBRARY_REPS}x{N_SIMS} reps per grid cell)",
        fontsize=13, y=1.01,
    )
    plt.tight_layout()
    svg_path = OUTPUT_DIR / "shendure_pairwise_power_mwu_all_cell_types.svg"
    fig.savefig(svg_path, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {svg_path}", flush=True)

    # --- 50% power contour ---
    power_threshold = 0.50
    fig, ax = plt.subplots(figsize=(10, 6))
    palette = sns.color_palette("tab10", n_colors=len(CELL_TYPES))

    for i, ct in enumerate(CELL_TYPES):
        ct_display = "Pluripotent" if ct == "reference" else ct
        ct_data = power_grid[power_grid["cell_type"] == ct]

        baselines = sorted(ct_data["baseline_mu"].unique())
        min_fc = []
        valid_baselines = []

        for bl in baselines:
            bl_data = ct_data[ct_data["baseline_mu"] == bl].sort_values(
                "log2_fc"
            )
            above = bl_data[bl_data["power"] >= power_threshold]
            if len(above) > 0:
                min_fc.append(above["log2_fc"].min())
                valid_baselines.append(bl)

        cells = scm.SHENDURE_BOUNDS.cells_per_cell_type.get(ct, "?")
        label = f"{ct_display} ({cells} cells)"

        if valid_baselines:
            ax.plot(valid_baselines, min_fc, "o-", color=palette[i],
                    label=label, markersize=4, lw=1.5)
        else:
            ax.plot(
                [], [], "o-", color=palette[i],
                label=f"{label} (no {int(power_threshold*100)}% power)"
            )

    ax.axhline(y=0.20, color="gray", linestyle="--", lw=1, alpha=0.7)
    ax.set_xscale("log")
    ax.set_xlabel("Baseline activity (mu)", fontsize=11)
    ax.set_ylabel(f"Min |log2 FC| for {int(power_threshold*100)}% power",
                  fontsize=11)
    ax.set_title(
        f"Minimum detectable fold change -- Shendure "
        f"(MWU, alpha=0.05, {int(power_threshold*100)}% power)",
        fontsize=12,
    )
    ax.legend(fontsize=8, loc="upper right", ncol=2)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    kircher_reference.annotate_continuous(ax)

    svg_path2 = OUTPUT_DIR / "shendure_pairwise_power_50pct_contours.svg"
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
                "--job-name=shendure_pw_worker",
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
