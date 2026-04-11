#!/usr/bin/env python3
"""
Cohen Welch t-test null calibration -- all cell types.

Mirrors `cohen_power_ttest_all_cell_types.py` exactly in scope and infra
(N_LIBRARY_REPS=156, N_SIMS=5, episomal one-sim-covers-all-CTs, both
+reporter and -reporter conditions). Difference: mu = minP for ALL CREs, so
every comparison is a true null. Output is a per-CT p-value distribution
(histogram + QQ) plus FPR@0.05 and KS statistic vs Uniform(0,1), each with
bootstrap 95% CIs.

Usage:
    python cohen_calibration_ttest_all_cell_types.py [phase]

Phases:
    simulate  -- generate simulations and run t-test (both conditions)
    plot      -- aggregate results and produce SVG + summary parquet
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
from scipy.stats import uniform, kstest

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
SIM_DIR = DATA_ROOT / "simulated" / f"{SIM_DATE}_cohen_cal"

# Cohen is episomal: same CREs across all CTs. Hardcoded n_cres=115 to match
# the cohen_power_ttest_all_cell_types.py convention. n_cres per CT in the
# new ortho is 114-116; 115 is the canonical value used in the power runs.
N_CRES = 115
CELL_TYPES = list(scm.COHEN_BOUNDS.cells_per_cell_type.index)

N_LIBRARY_REPS = 156
N_SIMS = 5
MINP = scm.COHEN_BOUNDS.reference_activity

BATCH_SIZE = 10
N_BOOTSTRAP = 1000
HIST_BINS = 25
QQ_N_QUANTILES = 2000  # downsample QQ to evenly-spaced quantiles to keep SVG light


# ---------------------------------------------------------------------------
# Phase: simulate
# ---------------------------------------------------------------------------

def phase_simulate(client):
    print(f"minP={MINP:.6f} (all CREs at minP for null calibration)")
    print(f"n_cres={N_CRES}, cell_types={CELL_TYPES}")
    print(f"library reps: {N_LIBRARY_REPS}, sims per rep: {N_SIMS}")
    print()
    for ct in CELL_TYPES:
        cells = scm.COHEN_BOUNDS.cells_per_cell_type.get(ct, "N/A")
        print(f"  {ct}: cells={cells}")

    def _process_batch(batch, hs_all_ct):
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
            s.ttest("hs_all_ct")
            s.ttest("hs_all_ct", has_reporter=False)
        for s in batch:
            s.save()
        return hs_all_ct, [(s.location, s.name) for s in batch]

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
    print(f"Resuming: {len(sim_paths)} complete sims on disk.")

    n_remaining = N_LIBRARY_REPS - len(sim_paths)
    print(f"Will run {n_remaining} fresh sims to reach {N_LIBRARY_REPS} total.")

    hs_all_ct = None
    batch = []

    for i in range(n_remaining):
        names = [f"synthcre_{j}" for j in range(N_CRES - 1)] + ["reference"]
        # Null: every CRE at minP, every CT.
        gt_df = pd.concat([
            pd.DataFrame({"cre_id": names, "mu": MINP, "cell_type": ct})
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
            print(f"Batch done -- {len(sim_paths)}/{N_LIBRARY_REPS} complete.",
                  flush=True)

    if batch:
        hs_all_ct, paths = _process_batch(batch, hs_all_ct)
        sim_paths.extend(paths)
        print(f"Batch done -- {len(sim_paths)}/{N_LIBRARY_REPS} complete.",
              flush=True)

    print(f"All {len(sim_paths)} sims done and tested.", flush=True)
    return sim_paths


# ---------------------------------------------------------------------------
# Phase: plot
# ---------------------------------------------------------------------------

def _collect_pvals_by_ct(sims, test_type):
    """Return {cell_type: np.array of p-values}."""
    rows = []
    for sim in sims:
        for i in range(sim.get_state_field("n_sims")):
            m = sim._merge_in_ground_truth(
                hypothesis_set_name="hs_all_ct",
                test_type=test_type,
                index=i,
            )
            rows.append(m[["comparison_cell_type", "p_value"]])
    combined = pd.concat(rows, ignore_index=True)
    out = {}
    for ct in CELL_TYPES:
        v = combined[combined["comparison_cell_type"] == ct]["p_value"].values
        v = v[np.isfinite(v)]
        out[ct] = v
    return out


def _bootstrap_ci(pvals, n_boot=N_BOOTSTRAP, alpha=0.05, seed=0):
    """Return dict with FPR@0.05 and KS-D vs uniform, with bootstrap 95% CI.

    Uses numpy fixed-seed RNG so the CI is deterministic given the input.
    """
    rng = np.random.default_rng(seed)
    n = len(pvals)
    if n < 2:
        nan = float("nan")
        return dict(n=n, fpr=nan, fpr_lo=nan, fpr_hi=nan,
                    ks_d=nan, ks_d_lo=nan, ks_d_hi=nan,
                    ks_p=nan)
    # Point estimates
    fpr = float(np.mean(pvals < alpha))
    ks_d, ks_p = kstest(pvals, "uniform")
    # Bootstrap
    fpr_boot = np.empty(n_boot)
    ks_d_boot = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        sample = pvals[idx]
        fpr_boot[b] = np.mean(sample < alpha)
        ks_d_boot[b] = kstest(sample, "uniform").statistic
    return dict(
        n=n,
        fpr=fpr,
        fpr_lo=float(np.percentile(fpr_boot, 2.5)),
        fpr_hi=float(np.percentile(fpr_boot, 97.5)),
        ks_d=float(ks_d),
        ks_d_lo=float(np.percentile(ks_d_boot, 2.5)),
        ks_d_hi=float(np.percentile(ks_d_boot, 97.5)),
        ks_p=float(ks_p),
    )


def phase_plot(client):
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

    pvals_reporter = _collect_pvals_by_ct(sims, "ttest")
    pvals_deflated = _collect_pvals_by_ct(sims, "ttest_deflated")

    # Save raw null p-value parquets
    for label, d in [("reporter", pvals_reporter), ("deflated", pvals_deflated)]:
        rows = []
        for ct, v in d.items():
            rows.append(pd.DataFrame({"cell_type": ct, "p_value": v}))
        out_path = OUTPUT_DIR / f"cohen_null_pvals_{label}.parquet"
        pd.concat(rows, ignore_index=True).to_parquet(out_path)
        print(f"Saved: {out_path}", flush=True)

    # Compute summary stats with bootstrap CIs
    summary_rows = []
    for cond_label, d in [("reporter", pvals_reporter),
                          ("deflated", pvals_deflated)]:
        for ct, v in d.items():
            stats = _bootstrap_ci(v)
            stats["cell_type"] = ct
            stats["condition"] = cond_label
            summary_rows.append(stats)
    summary = pd.DataFrame(summary_rows)
    cols = ["cell_type", "condition", "n", "fpr", "fpr_lo", "fpr_hi",
            "ks_d", "ks_d_lo", "ks_d_hi", "ks_p"]
    summary = summary[cols]
    summary_path = OUTPUT_DIR / "cohen_null_summary.parquet"
    summary.to_parquet(summary_path)
    print(f"Saved: {summary_path}", flush=True)
    print(summary.to_string(index=False), flush=True)

    # --- Plot: histogram + QQ per CT, paired (+/- reporter) ---
    # Layout: nrows = n_CTs, 4 cols = (hist_rep, qq_rep, hist_def, qq_def)
    nrows = len(CELL_TYPES)
    fig, axes = plt.subplots(nrows, 4, figsize=(16, 3 * nrows))
    if nrows == 1:
        axes = axes[np.newaxis, :]

    palette = sns.color_palette("tab10", n_colors=len(CELL_TYPES))

    for i, ct in enumerate(CELL_TYPES):
        ct_display = "Rod" if ct == "reference" else ct
        color = palette[i]
        cells_ct = scm.COHEN_BOUNDS.cells_per_cell_type.get(ct, "?")

        for j, (cond_label, cond_pretty, color_cond, d) in enumerate([
            ("reporter", "+reporter", "steelblue", pvals_reporter),
            ("deflated", "-reporter (deflated)", "coral", pvals_deflated),
        ]):
            v = d[ct]
            row = summary[(summary.cell_type == ct)
                          & (summary.condition == cond_label)].iloc[0]

            ax_h = axes[i, j * 2]
            ax_h.hist(v, bins=HIST_BINS, density=True,
                      edgecolor="black", linewidth=0.3, alpha=0.7,
                      color=color_cond)
            ax_h.axhline(1.0, color="red", linestyle="--", lw=1)
            ax_h.set_title(f"{ct_display} -- {cond_pretty}", fontsize=9)
            ax_h.set_xlabel("p-value", fontsize=8)
            if j == 0:
                ax_h.set_ylabel("Density", fontsize=8)
            ax_h.tick_params(labelsize=7)
            ax_h.set_xlim(0, 1)
            ax_h.text(
                0.97, 0.93,
                f"FPR@.05={row.fpr:.3f} [{row.fpr_lo:.3f},{row.fpr_hi:.3f}]\n"
                f"KS-D={row.ks_d:.3f} [{row.ks_d_lo:.3f},{row.ks_d_hi:.3f}]\n"
                f"KS-p={row.ks_p:.2g}\n"
                f"n={row.n}, cells={cells_ct}",
                transform=ax_h.transAxes, ha="right", va="top", fontsize=6,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.85),
            )

            ax_q = axes[i, j * 2 + 1]
            n = len(v)
            if n > 0:
                # Quantile-thinned QQ: evaluate at evenly-spaced quantiles
                # instead of plotting every point. Preserves QQ shape with
                # O(N_q) points.
                n_q = min(QQ_N_QUANTILES, n)
                qs = np.linspace(0, 1, n_q + 2)[1:-1]
                obs = np.quantile(v, qs)
                ax_q.scatter(qs, obs, s=2, alpha=0.7, color=color_cond)
            ax_q.plot([0, 1], [0, 1], "r--", lw=1)
            ax_q.set_title(f"{ct_display} -- {cond_pretty} QQ", fontsize=9)
            ax_q.set_xlabel("Expected (Uniform)", fontsize=8)
            ax_q.tick_params(labelsize=7)
            ax_q.set_xlim(0, 1)
            ax_q.set_ylim(0, 1)
            ax_q.set_aspect("equal")

    fig.suptitle(
        "Cohen (episomal) -- Welch t-test null calibration\n"
        "+reporter vs -reporter (deflated)",
        fontsize=12, y=1.005,
    )
    plt.tight_layout()
    svg_path = OUTPUT_DIR / "cohen_calibration_ttest_reporter_comparison.svg"
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
            memory="64G",
            processes=1,
            env_extra=[
                f"export PYTHONPATH={_repo_root}:$PYTHONPATH",
            ],
            job_extra_directives=[
                "-p priority",
                "--account=prio_skr2",
                "--job-name=cohen_cal_worker",
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
        print(f"\n{'='*60}\nPhase: simulate\n{'='*60}", flush=True)
        cluster, client = _make_slurm_client()
        phase_simulate(client)
        client.close()
        cluster.close()

        print(f"\n{'='*60}\nPhase: plot\n{'='*60}", flush=True)
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
        print(f"Unknown phase: {phase!r}. Choose simulate, plot, or all")
        sys.exit(1)

    print("\nDone.", flush=True)
