#!/usr/bin/env python3
"""
Takeshi MWU null calibration follow-up.

Reuses the existing 2026-04-11_takeshi_cal sim dirs (produced by
takeshi_calibration_ttest_all_cell_types.py) and runs Mann-Whitney U on the
same hs_all_ct hypothesis set, in both +reporter and deflated modes.

No resimulation. Each sim already has its data.parquet on disk and an
hs_all_ct hypothesis set; this script just adds mwu/ and mwu_deflated/
result subdirs alongside the existing ttest/ and ttest_deflated/.

Output: per-CT MWU p-value distributions (histogram + QQ) plus FPR@0.05
(with closed-form 95% binomial CI) and KS-D statistic as a point estimate,
structured the same way as the t-test calibration plot so the two are
visually comparable.

Motivation: the t-test cal showed FPR inflation (~2x +reporter, ~3x deflated)
on takeshi. Is the inflation t-test specific (catastrophic-cancellation /
non-normal counts) or fundamental to the data (overdispersion / many
hypotheses / per-CRE tie behavior)? MWU doesn't depend on normality, so if
its FPR is closer to nominal then the t-test is the problem; if MWU is
similarly inflated, the issue is structural.

Usage:
    python takeshi_calibration_mwu_followup.py [phase]

Phases:
    compute   -- walk existing sims and run MWU on any missing
    plot      -- aggregate MWU results and produce SVG + summary parquet
    all       -- both phases sequentially (default)
"""

import sys
import os
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import kstest

_repo_root = Path(__file__).resolve().parents[4]
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

import scMPRAforge as scm
from dask.distributed import Client, LocalCluster
from dask_jobqueue import SLURMCluster

# ---------------------------------------------------------------------------
# Paths and parameters (mirror the ttest cal script)
# ---------------------------------------------------------------------------

DATA_ROOT = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new")
OUTPUT_DIR = Path(__file__).resolve().parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

SIM_DATE = "2026-04-11"
SIM_DIR = DATA_ROOT / "simulated" / f"{SIM_DATE}_takeshi_cal"

CELL_TYPES = sorted(scm.TAKESHI_BOUNDS.cells_per_cell_type.index.tolist())

HIST_BINS = 25
QQ_N_QUANTILES = 2000
ALPHA = 0.05


# ---------------------------------------------------------------------------
# Phase: compute (walk existing sims, run MWU on any without it)
# ---------------------------------------------------------------------------

def _sim_has(sim_dir, test_subdir):
    d = sim_dir / "tests" / "hs_all_ct" / test_subdir
    return d.exists() and any(d.iterdir())


def phase_compute(client):
    print(f"Walking {SIM_DIR}", flush=True)
    sim_objs = {}
    for cell_type in CELL_TYPES:
        ct_dir = SIM_DIR / cell_type
        sims_ct = []
        if ct_dir.exists():
            for d in sorted(ct_dir.iterdir()):
                if not d.is_dir():
                    continue
                # only sims that have a ttest already are guaranteed to have
                # the hypothesis set saved (the ttest script writes the
                # hypotheses.tsv as a side effect of add_hypothesis_set).
                if not _sim_has(d, "ttest"):
                    continue
                sims_ct.append(
                    scm.de_novo_simulation(
                        location=ct_dir, name=d.name, client=client
                    )
                )
        sim_objs[cell_type] = sims_ct
        print(f"  {cell_type}: loaded {len(sims_ct)} sims", flush=True)

    # --- mwu: +reporter ---
    for cell_type, sims in sim_objs.items():
        n_to_run = sum(1 for s in sims
                       if not _sim_has(s.location / s.name, "mwu"))
        print(f"mwu (+reporter): {cell_type} -- {n_to_run} need running",
              flush=True)
        for sim in sims:
            if _sim_has(sim.location / sim.name, "mwu"):
                continue
            sim.mwu("hs_all_ct")
            sim.save()

    for sims in sim_objs.values():
        for sim in sims:
            sim.save()
    print("mwu (+reporter) done and saved.", flush=True)

    # --- mwu: deflated ---
    for cell_type, sims in sim_objs.items():
        n_to_run = sum(1 for s in sims
                       if not _sim_has(s.location / s.name, "mwu_deflated"))
        print(f"mwu (deflated): {cell_type} -- {n_to_run} need running",
              flush=True)
        for sim in sims:
            if _sim_has(sim.location / sim.name, "mwu_deflated"):
                continue
            sim.mwu("hs_all_ct", has_reporter=False)
            sim.save()

    for sims in sim_objs.values():
        for sim in sims:
            sim.save()
    print("mwu (deflated) done and saved.", flush=True)


# ---------------------------------------------------------------------------
# Phase: plot
# ---------------------------------------------------------------------------

def _collect_pvals(sims, test_type):
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


def _summary_stats(pvals, alpha=ALPHA):
    """Closed-form 95% CI for FPR (binomial normal approximation).

    Bootstrap was dropped -- at n >= ~10k the Wald interval is tight to
    ~1e-3 and matches bootstrap to 3 decimals. KS-D is reported as a
    point estimate; its distribution under H0 is well-characterized by
    the returned `ks_p`.
    """
    n = len(pvals)
    if n < 2:
        nan = float("nan")
        return dict(n=n, fpr=nan, fpr_lo=nan, fpr_hi=nan,
                    ks_d=nan, ks_p=nan)
    fpr = float(np.mean(pvals < alpha))
    se = (fpr * (1.0 - fpr) / n) ** 0.5
    fpr_lo = max(0.0, fpr - 1.96 * se)
    fpr_hi = min(1.0, fpr + 1.96 * se)
    ks_d, ks_p = kstest(pvals, "uniform")
    return dict(
        n=n, fpr=fpr, fpr_lo=float(fpr_lo), fpr_hi=float(fpr_hi),
        ks_d=float(ks_d), ks_p=float(ks_p),
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
                if not _sim_has(d, "mwu"):
                    continue
                sims_ct.append(
                    scm.de_novo_simulation(
                        location=ct_dir, name=d.name, client=client
                    )
                )
        all_sims[cell_type] = sims_ct
        print(f"  {cell_type}: loaded {len(sims_ct)} sims", flush=True)

    pvals_reporter = {ct: _collect_pvals(sims, "mwu")
                      for ct, sims in all_sims.items()}
    pvals_deflated = {ct: _collect_pvals(sims, "mwu_deflated")
                      for ct, sims in all_sims.items()}

    for label, d in [("reporter", pvals_reporter), ("deflated", pvals_deflated)]:
        rows = []
        for ct, v in d.items():
            rows.append(pd.DataFrame({"cell_type": ct, "p_value": v}))
        out_path = OUTPUT_DIR / f"takeshi_null_pvals_mwu_{label}.parquet"
        pd.concat(rows, ignore_index=True).to_parquet(out_path)
        print(f"Saved: {out_path}", flush=True)

    summary_rows = []
    for cond_label, d in [("reporter", pvals_reporter),
                          ("deflated", pvals_deflated)]:
        for ct, v in d.items():
            stats = _summary_stats(v)
            stats["cell_type"] = ct
            stats["condition"] = cond_label
            summary_rows.append(stats)
    summary = pd.DataFrame(summary_rows)
    cols = ["cell_type", "condition", "n", "fpr", "fpr_lo", "fpr_hi",
            "ks_d", "ks_p"]
    summary = summary[cols]
    summary_path = OUTPUT_DIR / "takeshi_null_summary_mwu.parquet"
    summary.to_parquet(summary_path)
    print(f"Saved: {summary_path}", flush=True)
    print(summary.to_string(index=False), flush=True)

    nrows = len(CELL_TYPES)
    fig, axes = plt.subplots(nrows, 4, figsize=(16, 3 * nrows))
    if nrows == 1:
        axes = axes[np.newaxis, :]

    for i, ct in enumerate(CELL_TYPES):
        ct_display = "HepG2" if ct == "reference" else ct
        cells = scm.TAKESHI_BOUNDS.cells_per_cell_type.get(ct, "?")

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
                f"KS-D={row.ks_d:.3f}\n"
                f"KS-p={row.ks_p:.2g}\n"
                f"n={row.n}, cells={cells}",
                transform=ax_h.transAxes, ha="right", va="top", fontsize=6,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.85),
            )

            ax_q = axes[i, j * 2 + 1]
            n = len(v)
            if n > 0:
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
        "Takeshi (piggyBac) -- MWU null calibration\n"
        "+reporter vs -reporter (deflated) -- cross-check vs t-test cal",
        fontsize=12, y=1.005,
    )
    plt.tight_layout()
    svg_path = OUTPUT_DIR / "takeshi_calibration_mwu_reporter_comparison.svg"
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
                "--job-name=takeshi_mwucal_worker",
                "--time=4:00:00",
                "--exclude=a1132u18n02,a1130u22n02",
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
        print(f"\n{'='*60}\nPhase: compute\n{'='*60}", flush=True)
        cluster, client = _make_slurm_client()
        phase_compute(client)
        client.close()
        cluster.close()

        print(f"\n{'='*60}\nPhase: plot\n{'='*60}", flush=True)
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
        print(f"Unknown phase: {phase!r}. Choose compute, plot, or all")
        sys.exit(1)

    print("\nDone.", flush=True)
