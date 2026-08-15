#!/usr/bin/env python3
"""
Shendure Welch t-test null calibration -- all cell types.

Mirrors `shendure_power_ttest_all_cell_types.py` exactly in scope and infra
(N_LIBRARY_REPS=156, N_SIMS=5, per-CT sims, both +reporter and -reporter).
Difference: mu = minP for ALL CREs, so every comparison is a true null.
Output is a per-CT p-value distribution (histogram + QQ) plus FPR@0.05 and
KS statistic vs Uniform(0,1), each with bootstrap 95% CIs.

Usage:
    python shendure_calibration_ttest_all_cell_types.py [phase]

Phases:
    simulate  -- generate simulations and run t-test (both conditions)
    plot      -- aggregate results and produce SVG + summary parquet
    all       -- run both phases sequentially (default)
"""

import sys
import os
import pickle
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
SIM_DIR = DATA_ROOT / "simulated" / f"{SIM_DATE}_shendure_cal"

ORTHO_DIR = DATA_ROOT / "shendure" / "shendure_obs_nb"
with open(ORTHO_DIR / "by_cell_type_parameters.pkl", "rb") as f:
    _params = pickle.load(f)
CELL_TYPES = sorted(_params.nb.keys())
N_CRES_PER_CT = {ct: len(_params.nb[ct]) for ct in CELL_TYPES}
del _params

MINP = scm.SHENDURE_BOUNDS.reference_activity
N_LIBRARY_REPS = 156
N_SIMS = 5
N_BOOTSTRAP = 1000
HIST_BINS = 25
QQ_N_QUANTILES = 2000  # downsample QQ to evenly-spaced quantiles to keep SVG light


# ---------------------------------------------------------------------------
# Phase: simulate
# ---------------------------------------------------------------------------

def phase_simulate(client):
    print(f"minP={MINP:.6f} (all CREs at minP for null calibration)")
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
        n_cres = N_CRES_PER_CT[cell_type]

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

        n_remaining = N_LIBRARY_REPS - len(sims_ct)
        print(f"  Running {n_remaining} fresh sims.", flush=True)
        for i in range(n_remaining):
            names = [f"synthcre_{j}" for j in range(n_cres - 1)] + ["reference"]
            gt_df = pd.DataFrame({"cre_id": names, "mu": MINP,
                                  "cell_type": cell_type})
            libraries = [
                scm.simulate_library(
                    CREs=pd.Series(names),
                    library_model=scm.SHENDURE_BOUNDS.library_model,
                )
                for _ in range(N_SIMS)
            ]
            sim = scm.de_novo_simulation(
                location=ct_dir,
                name=f"sim_{uuid.uuid4().hex[:8]}",
                client=client,
                libraries=libraries,
                library_mapping="corresponding",
                flatten_overtransfection=True,
                n_sims=N_SIMS,
                experiment_bounds=bound_ct,
                ground_truth=gt_df,
            )
            sim.gamut()
            sims_ct.append(sim)

        all_sims[cell_type] = sims_ct
        print(f"  {len(sims_ct)} replicates ready.", flush=True)

    for sims in all_sims.values():
        for sim in sims:
            sim.save()
    print("\nAll sims saved.", flush=True)

    # --- t-test: +reporter ---
    for cell_type, sims in all_sims.items():
        print(f"t-test (+reporter): {cell_type}", flush=True)
        hs = None
        for sim in sims:
            tt_dir = sim.location / sim.name / "tests" / "hs_all_ct" / "ttest"
            if tt_dir.exists() and any(tt_dir.iterdir()):
                continue
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

    # --- t-test: deflated ---
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

def _collect_pvals(sims, test_type):
    """Return concatenated p-value array (all CREs, all sims)."""
    arrs = []
    for sim in sims:
        for i in range(sim.get_state_field("n_sims")):
            m = sim._merge_in_ground_truth(
                hypothesis_set_name="hs_all_ct",
                test_type=test_type,
                index=i,
            )
            arrs.append(m["p_value"].values)
    if not arrs:
        return np.array([])
    v = np.concatenate(arrs)
    return v[np.isfinite(v)]


def _bootstrap_ci(pvals, n_boot=N_BOOTSTRAP, alpha=0.05, seed=0):
    rng = np.random.default_rng(seed)
    n = len(pvals)
    if n < 2:
        nan = float("nan")
        return dict(n=n, fpr=nan, fpr_lo=nan, fpr_hi=nan,
                    ks_d=nan, ks_d_lo=nan, ks_d_hi=nan, ks_p=nan)
    fpr = float(np.mean(pvals < alpha))
    ks_d, ks_p = kstest(pvals, "uniform")
    fpr_boot = np.empty(n_boot)
    ks_d_boot = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        sample = pvals[idx]
        fpr_boot[b] = np.mean(sample < alpha)
        ks_d_boot[b] = kstest(sample, "uniform").statistic
    return dict(
        n=n, fpr=fpr,
        fpr_lo=float(np.percentile(fpr_boot, 2.5)),
        fpr_hi=float(np.percentile(fpr_boot, 97.5)),
        ks_d=float(ks_d),
        ks_d_lo=float(np.percentile(ks_d_boot, 2.5)),
        ks_d_hi=float(np.percentile(ks_d_boot, 97.5)),
        ks_p=float(ks_p),
    )


def phase_plot(client):
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

    pvals_reporter = {ct: _collect_pvals(sims, "ttest")
                      for ct, sims in all_sims.items()}
    pvals_deflated = {ct: _collect_pvals(sims, "ttest_deflated")
                      for ct, sims in all_sims.items()}

    for label, d in [("reporter", pvals_reporter), ("deflated", pvals_deflated)]:
        rows = []
        for ct, v in d.items():
            rows.append(pd.DataFrame({"cell_type": ct, "p_value": v}))
        out_path = OUTPUT_DIR / f"shendure_null_pvals_{label}.parquet"
        pd.concat(rows, ignore_index=True).to_parquet(out_path)
        print(f"Saved: {out_path}", flush=True)

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
    summary_path = OUTPUT_DIR / "shendure_null_summary.parquet"
    summary.to_parquet(summary_path)
    print(f"Saved: {summary_path}", flush=True)
    print(summary.to_string(index=False), flush=True)

    nrows = len(CELL_TYPES)
    fig, axes = plt.subplots(nrows, 4, figsize=(16, 3 * nrows))
    if nrows == 1:
        axes = axes[np.newaxis, :]

    for i, ct in enumerate(CELL_TYPES):
        ct_display = "Pluripotent" if ct == "reference" else ct
        cells = scm.SHENDURE_BOUNDS.cells_per_cell_type.get(ct, "?")

        for j, (cond_label, cond_pretty, color_cond, d) in enumerate([
            ("reporter", "+reporter", "steelblue", pvals_reporter),
            ("deflated", "-reporter (deflated)", "coral", pvals_deflated),
        ]):
            v = d[ct]
            row = summary[(summary.cell_type == ct)
                          & (summary.condition == cond_label)].iloc[0]

            ax_h = axes[i, j * 2]
            if len(v) > 0:
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
                f"n={row.n}, cells={cells}",
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
        "Shendure (piggyBac) -- Welch t-test null calibration\n"
        "+reporter vs -reporter (deflated)",
        fontsize=12, y=1.005,
    )
    plt.tight_layout()
    svg_path = OUTPUT_DIR / "shendure_calibration_ttest_reporter_comparison.svg"
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
                "--job-name=shendure_cal_worker",
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
