#!/usr/bin/env python3
"""
Cohen pairwise CRE power analysis -- Welch t-test

Companion to cohen_pairwise_power_mwu.py. Reuses the sims already on disk
(produced by the MWU script) and runs Welch's t-test on the same
'hs_pairwise' hypothesis set. No resimulation: count parquets are loaded
from disk and ttest results are written to tests/hs_pairwise/ttest/
alongside the existing tests/hs_pairwise/mwu/ results.

Usage:
    python cohen_pairwise_power_ttest.py [phase]

Phases:
    compute   -- walk existing sims and run t-test on any missing
    plot      -- aggregate t-test results and produce SVG plots
    all       -- run both phases sequentially (default)
"""

import sys
import os
import pickle
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
# Paths and parameters (mirror cohen_pairwise_power_mwu.py)
# ---------------------------------------------------------------------------

DATA_ROOT = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new")
OUTPUT_DIR = Path(__file__).resolve().parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

SIM_DATE = "2026-04-10"
SIM_DIR = DATA_ROOT / "simulated" / f"{SIM_DATE}_cohen_pw"

ORTHO_DIR = DATA_ROOT / "cohen" / "cohen_obsingle_nb_phantom"
with open(ORTHO_DIR / "by_cell_type_parameters.pkl", "rb") as f:
    _params = pickle.load(f)

CELL_TYPES = sorted(_params.nb.keys())
N_CRES = max(len(_params.nb[ct]) for ct in CELL_TYPES)
del _params

N_LIBRARY_REPS = 20
N_SIMS = 5
ALPHA = 0.05


# ---------------------------------------------------------------------------
# Sim discovery
# ---------------------------------------------------------------------------

def _sim_has_mwu(sim_dir):
    d = sim_dir / "tests" / "hs_pairwise" / "mwu"
    return d.exists() and any(d.iterdir())


def _sim_has_ttest(sim_dir):
    d = sim_dir / "tests" / "hs_pairwise" / "ttest"
    return d.exists() and any(d.iterdir())


def _discover_sims():
    """Yield (cohort_dir, sim_name) for every sim whose MWU results exist."""
    sims = []
    for cohort_dir in sorted(SIM_DIR.glob("cohort_*")):
        if not cohort_dir.is_dir():
            continue
        for d in sorted(cohort_dir.iterdir()):
            if d.is_dir() and _sim_has_mwu(d):
                sims.append((cohort_dir, d.name))
    return sims


# ---------------------------------------------------------------------------
# Phase: compute
# ---------------------------------------------------------------------------

def phase_compute(client):
    sims = _discover_sims()
    print(f"Discovered {len(sims)} complete sims under {SIM_DIR}.",
          flush=True)

    to_run = [
        (cd, sn) for (cd, sn) in sims if not _sim_has_ttest(cd / sn)
    ]
    print(f"{len(to_run)} need t-test; "
          f"{len(sims) - len(to_run)} already have it.", flush=True)

    for k, (cohort_dir, sim_name) in enumerate(to_run):
        sim = scm.de_novo_simulation(
            location=cohort_dir, name=sim_name, client=client
        )
        sim.ttest("hs_pairwise")
        sim.save()
        print(f"  [{k+1}/{len(to_run)}] {cohort_dir.name}/{sim_name}",
              flush=True)

    print("Compute phase done.", flush=True)


# ---------------------------------------------------------------------------
# Phase: plot
# ---------------------------------------------------------------------------

def phase_plot(client):
    sims = _discover_sims()
    sims = [(cd, sn) for (cd, sn) in sims if _sim_has_ttest(cd / sn)]
    print(f"Loading {len(sims)} sims with t-test results.", flush=True)

    power_rows = []
    null_rows = []

    for cohort_dir, sim_name in sims:
        sim = scm.de_novo_simulation(
            location=cohort_dir, name=sim_name, client=client
        )
        gt = sim.ground_truth

        for i in range(sim.get_state_field("n_sims")):
            m = sim._merge_in_ground_truth(
                hypothesis_set_name="hs_pairwise",
                test_type="ttest",
                index=i,
            )

            for _, row in m.iterrows():
                cre_name = row["comparison_CRE"]
                ct = row["comparison_cell_type"]

                if "_fc_" in cre_name:
                    parts = cre_name.split("_fc_")
                    anchor_name = parts[0]
                    log2_fc = float(parts[1])
                    anchor_mu = gt[
                        (gt["cre_id"] == anchor_name)
                        & (gt["cell_type"] == ct)
                    ]["mu"].iloc[0]

                    power_rows.append({
                        "cell_type": ct,
                        "baseline_mu": anchor_mu,
                        "log2_fc": log2_fc,
                        "p_value": row["p_value"],
                        "reject_null": row["p_value"] < ALPHA,
                    })
                elif cre_name.startswith("filler_"):
                    null_rows.append({
                        "cell_type": ct,
                        "p_value": row["p_value"],
                    })

    power_df = pd.DataFrame(power_rows)
    null_df = pd.DataFrame(null_rows)

    power_df.to_parquet(
        OUTPUT_DIR / "cohen_pairwise_power_ttest_df.parquet", index=False
    )
    null_df.to_parquet(
        OUTPUT_DIR / "cohen_pairwise_null_ttest_df.parquet", index=False
    )
    print(f"Saved power_df ({len(power_df)}) and null_df ({len(null_df)}).",
          flush=True)

    for ct in CELL_TYPES:
        ct_null = null_df[null_df["cell_type"] == ct]
        if len(ct_null) > 0:
            fpr = ct_null["p_value"].lt(ALPHA).mean()
            print(f"  {ct}: null FPR@{ALPHA} = {fpr:.4f} "
                  f"(n={len(ct_null)})", flush=True)

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
        ct_display = "Rod" if ct == "reference" else ct
        ax = axes[i]
        ct_data = power_grid[power_grid["cell_type"] == ct]

        pivot = ct_data.pivot_table(
            index="log2_fc", columns="baseline_mu", values="power"
        )
        pivot = pivot.sort_index(ascending=False)
        pivot = pivot[sorted(pivot.columns)]

        cells_ct = scm.COHEN_BOUNDS.cells_per_cell_type.get(ct, "?")

        sns.heatmap(
            pivot, ax=ax, vmin=0, vmax=1,
            cmap="YlOrRd", linewidths=0.5, linecolor="white",
            cbar_kws={"label": "Power", "shrink": 0.8},
            xticklabels=[f"{x:.3f}" for x in sorted(pivot.columns)],
            yticklabels=[f"{y:.2f}" for y in pivot.index],
        )

        # Draws nothing if the grid stops below the reference effects, which
        # it does for any regime not re-simulated on the extended range.
        kircher_reference.annotate(ax, pivot.index)

        ax.set_title(f"{ct_display}\n(n_cres={N_CRES}, cells={cells_ct})",
                     fontsize=10)
        ax.set_xlabel("Baseline activity (mu)", fontsize=9)
        ax.set_ylabel("|log2 FC|", fontsize=9)
        ax.tick_params(labelsize=7)
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

    for j in range(len(CELL_TYPES), len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(
        "Welch t-test pairwise power -- Cohen (episomal)\n"
        f"(alpha={ALPHA}, {N_LIBRARY_REPS}x{N_SIMS} reps per grid cell)",
        fontsize=13, y=1.01,
    )
    plt.tight_layout()
    svg_path = OUTPUT_DIR / "cohen_pairwise_power_ttest_all_cell_types.svg"
    fig.savefig(svg_path, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {svg_path}", flush=True)

    # --- 50% power contour ---
    power_threshold = 0.50
    fig, ax = plt.subplots(figsize=(10, 6))
    palette = sns.color_palette("tab10", n_colors=len(CELL_TYPES))

    for i, ct in enumerate(CELL_TYPES):
        ct_display = "Rod" if ct == "reference" else ct
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

        cells_ct = scm.COHEN_BOUNDS.cells_per_cell_type.get(ct, "?")
        label = f"{ct_display} ({cells_ct} cells)"

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
        f"Minimum detectable fold change -- Cohen "
        f"(Welch t-test, alpha={ALPHA}, "
        f"{int(power_threshold*100)}% power)",
        fontsize=12,
    )
    ax.legend(fontsize=8, loc="upper right", ncol=2)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    kircher_reference.annotate_continuous(ax)

    svg_path2 = (
        OUTPUT_DIR / "cohen_pairwise_power_ttest_50pct_contours.svg"
    )
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
                "--job-name=cohen_pw_tt_worker",
                "--time=6:00:00",
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
        print("Phase: compute", flush=True)
        print(f"{'='*60}", flush=True)
        cluster, client = _make_slurm_client()
        phase_compute(client)
        client.close()
        cluster.close()

        print(f"\n{'='*60}", flush=True)
        print("Phase: plot", flush=True)
        print(f"{'='*60}", flush=True)
        cluster, client = _make_local_client()
        phase_plot(client)
        client.close()
        cluster.close()

    elif phase == "compute":
        cluster, client = _make_slurm_client()
        phase_compute(client)
        client.close()
        cluster.close()

    elif phase == "plot":
        cluster, client = _make_local_client()
        phase_plot(client)
        client.close()
        cluster.close()

    else:
        print(f"Unknown phase: {phase!r}. Choose from: compute, plot, all")
        sys.exit(1)

    print("\nDone.", flush=True)
