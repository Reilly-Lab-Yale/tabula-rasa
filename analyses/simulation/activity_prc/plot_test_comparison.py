#!/usr/bin/env python
"""Test comparison for the manuscript: auROC and auPRC across the three regimes.

Replaces the barcharts previously produced by all_prc_summary.ipynb. Two
changes of substance:

  - Box plots rather than bars. Each point estimate rests on 25 simulated
    experiments (5 ground-truth draws x 5 replicates), and the spread is the
    result: the Wald test's mean is unremarkable on Zhao et al. while its
    replicate spread is five times any other test's. A bar hides that.
  - No pooled column. auROC is not comparable across regimes of different
    difficulty, so a median over all three answers no question. Pooling is
    also actively misleading here: the two tests that cannot be run on Yin
    et al. would be scored only on the two easier datasets and would top the
    ranking (Wald 0.899 pooled against 0.878 for the t-test, which reverses
    on the datasets where both ran).

Two cells are deliberately absent rather than plotted:

  Yin et al. x Wald       Wald is the only test that reads the fitted model
                          rather than the counts, and those orthos were fit
                          under plain consider-missing rather than with the
                          MOI correction that dataset's canonical fit uses
                          (design fit_mode='cm_phantom'). The coefficients and
                          their standard errors therefore both omit the
                          transfection term, and the number is not comparable
                          to the other two regimes.
  Yin et al. x pseudobulk Pseudobulk collapses cells to one mean per replicate
                          and needs at least two per group. This design has one
                          biological replicate.

Reads the cached per-replicate summaries, so it needs no cluster and no dask.

    python analyses/simulation/activity_prc/plot_test_comparison.py
"""
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["svg.fonttype"] = "none"

BASE = pathlib.Path(__file__).resolve().parent
OUT = BASE / "output" / "paper_figs"

BLUE, ORANGE, GREEN, MUTED = "#0072b2", "#d55e00", "#009e73", "#6b6b6b"
INK = "#1a1a1a"

DATASETS = [("shendure", "Lalanne et al."),
            ("cohen", "Zhao et al."),
            ("seelig", "Yin et al.")]
TEST_LABEL = {"mwu": "MWU", "ttest": "t-test", "ks": "KS",
              "pseudobulk": "pseudobulk", "wald_auto": "Wald"}
TEST_ORDER = ["mwu", "ttest", "ks", "pseudobulk", "wald_auto"]

# (dataset, test) pairs that are not evaluable; see module docstring.
EXCLUDED = {("seelig", "wald_auto"), ("seelig", "pseudobulk")}


def load():
    frames = []
    for ds, _ in DATASETS:
        f = BASE / ds / "output" / f"{ds}_5x5_activity_summary.tsv"
        assert f.is_file(), f"missing summary: {f}"
        d = pd.read_csv(f, sep="\t")
        d["dataset"] = ds
        frames.append(d)
    a = pd.concat(frames, ignore_index=True)

    for col in ("auroc", "auprc", "test", "replicate", "gt_draw"):
        assert col in a.columns, f"summary is missing {col}"

    # Drop the pairs we cannot evaluate before anything is summarised, so an
    # excluded cell cannot leak into the aggregate.
    before = len(a)
    a = a[~a.set_index(["dataset", "test"]).index.isin(EXCLUDED)].copy()
    print(f"dropped {before - len(a)} rows for excluded (dataset, test) pairs")

    # Every surviving cell should rest on the full 5x5 design.
    n = a.groupby(["dataset", "test"]).size()
    assert (n >= 20).all(), f"thin cells: {n[n < 20].to_dict()}"
    return a


def panel(ax, a, metric):
    """One metric: tests on x, a box per dataset plus an aggregate."""
    groups = list(DATASETS)
    width = 0.8 / len(groups)
    colours = {"shendure": BLUE, "cohen": ORANGE, "seelig": GREEN}

    for gi, (ds, _lab) in enumerate(groups):
        for ti, test in enumerate(TEST_ORDER):
            sub = a[(a.test == test) & (a.dataset == ds)]
            if sub.empty:
                continue
            pos = ti + (gi - (len(groups) - 1) / 2) * width
            bp = ax.boxplot([sub[metric].values], positions=[pos], widths=width * 0.85,
                            patch_artist=True, showfliers=False, zorder=3,
                            medianprops=dict(color=INK, lw=1.1),
                            whiskerprops=dict(color=colours[ds], lw=0.9),
                            capprops=dict(color=colours[ds], lw=0.9),
                            boxprops=dict(facecolor=colours[ds], edgecolor=colours[ds],
                                          alpha=0.55, lw=0.9))
            del bp

    ax.set_xticks(range(len(TEST_ORDER)))
    ax.set_xticklabels([TEST_LABEL[t] for t in TEST_ORDER], fontsize=8.5)
    ax.set_ylabel(f"au{metric[2:].upper()}", fontsize=9, color=INK)
    ax.grid(True, axis="y", color="#e6e6e6", lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#cccccc")
    ax.tick_params(colors=MUTED, labelsize=8, length=0)


def main():
    a = load()

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.4))
    for ax, metric in zip(axes, ["auroc", "auprc"]):
        panel(ax, a, metric)

    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=c, edgecolor=c, alpha=0.55)
               for c in (BLUE, ORANGE, GREEN)]
    fig.legend(handles=handles,
               labels=[lab for _, lab in DATASETS],
               loc="lower center", ncol=3, frameon=False, fontsize=8,
               bbox_to_anchor=(0.5, -0.04))

    fig.tight_layout(rect=(0, 0.06, 1, 1))
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(OUT / f"test_comparison_box.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT/'test_comparison_box.svg'} and .png")

    for metric in ("auroc", "auprc"):
        print(f"\nmedian {metric} by dataset:")
        med = a.groupby(["dataset", "test"])[metric].median().unstack()
        print(med.reindex(columns=TEST_ORDER).round(3).to_string())


if __name__ == "__main__":
    main()
