#!/usr/bin/env python
"""Median-replicate ROC and PRC curves, one pair per dataset.

Scriptified from the paper-figure cell of all_prc_summary.ipynb so the curves
can be regenerated without stepping through a notebook. Two changes:

  - Wald is dropped from Yin et al. Wald is the only test that reads the
    fitted model rather than the counts, and those orthos were fit under plain
    consider-missing rather than with the MOI correction that dataset's
    canonical fit uses (design fit_mode='cm_phantom'), so its curve there is
    not comparable to the other regimes. Pseudobulk is absent from Yin et al.
    on its own: it needs at least two replicates per group and that design has
    one.
  - The curve shown is chosen the same way as before -- the ground-truth draw
    whose mean auROC is closest to the median across draws -- but the choice
    is asserted rather than assumed.

Needs the simulation objects, hence a Dask client; run it on the cluster.

    python analyses/simulation/activity_prc/plot_curves.py
"""
import pathlib
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PathCollection
import pandas as pd

sys.path.insert(0, "/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa")
import scMPRAforge as scm  # noqa: E402

plt.rcParams["svg.fonttype"] = "none"

BASE = pathlib.Path(__file__).resolve().parent
OUT = BASE / "output" / "paper_figs"
SIM_ROOT = pathlib.Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new/simulated")

N_GT_DRAWS = 5
HYPOTHESIS_SET = "hs_all_ct"

# Legends carry the internal test identifiers; the manuscript uses these.
TEST_LABEL = {"mwu": "MWU", "ttest": "t-test", "ks": "KS",
              "pseudobulk": "pseudobulk", "wald_auto": "Wald"}

# One legend serves several panels once these are assembled (Fig. 2A pairs ROC
# with PRC, Fig. S2 puts one legend over five panels), so it may carry only
# what every panel shares. AUC belongs to a ROC panel and AP to a PRC panel,
# and each panel has its own baseline, so all of that is stripped here and the
# baselines are labelled directly on the line they describe instead.
BASELINE_LABEL = re.compile(r"^(chance|baseline=([0-9.]+))$")

# Each curve type leaves a different corner empty.
LEGEND_LOC = {"ROC": "lower right", "PRC": "lower left"}

# A fixed hue per test, for the same reason the legend carries no numbers: it
# is shared. Matplotlib colours by draw order, so the Yin et al. panels, which
# run three tests rather than five, put t-test in the green the legend has
# already given to pseudobulk. These are the default cycle colours in the order
# the five-test panels get them, so only the three-test panels change.
TEST_COLOR = {"ks": "#1f77b4", "mwu": "#ff7f0e", "pseudobulk": "#2ca02c",
              "ttest": "#d62728", "wald_auto": "#9467bd"}

DATASETS = {
    "shendure": ("Lalanne et al.", "shendure_5x5_activity"),
    "cohen":    ("Zhao et al.",    "cohen_5x5_activity"),
    "seelig":   ("Yin et al.",     "seelig_5x5_activity"),
}
# Tests to draw per dataset. None means "whatever ran".
#
# Lalanne et al. and Zhao et al. both run all five. Neither is special, and the
# manuscript features Lalanne et al. in Fig. 2A only as a representative
# example -- Zhao et al. would do as well. Yin et al. is the one that draws
# three, and the two missing tests are missing for unrelated reasons:
#
#   pseudobulk  never ran. It collapses each element's cells to one value per
#               replicate and compares those, so it needs at least two
#               replicates per group, and that design has one. Absent from the
#               summary TSV entirely.
#   wald_auto   ran, and is dropped here on purpose. It is the only test that
#               reads the fitted model rather than the counts, and the orthos
#               refit for this benchmark used plain consider-missing rather
#               than the MOI correction the canonical fit uses
#               (fit_mode='cm_phantom'), so its variance estimates are not
#               comparable to the other two regimes. Its median auROC there is
#               0.633 against 0.83-0.84 for the rest, which is an artefact of
#               that mismatch and not a result about the test.
#
# plot_test_comparison.py drops the same two cells via EXCLUDED; keep the two
# in step.
TESTS_FOR = {
    "shendure": None,
    "cohen": None,
    "seelig": ["mwu", "ttest", "ks"],
}


def median_gt_draw(slug):
    """The draw whose mean auROC is nearest the median across draws."""
    f = BASE / slug / "output" / f"{slug}_5x5_activity_summary.tsv"
    assert f.is_file(), f"missing summary: {f}"
    d = pd.read_csv(f, sep="\t")
    gt_means = d.groupby("gt_draw")["auroc"].mean()
    assert len(gt_means) == N_GT_DRAWS, (
        f"{slug}: expected {N_GT_DRAWS} ground-truth draws, found {len(gt_means)}")
    pick = int((gt_means - gt_means.median()).abs().idxmin())
    print(f"  {slug}: GT draw {pick} "
          f"(mean auROC {gt_means[pick]:.3f}; draws span "
          f"{gt_means.min():.3f}-{gt_means.max():.3f})")
    return pick


def relabel(entry):
    """'wald_auto (AUC=0.929)' -> 'Wald'; other entries untouched."""
    head, _, _rest = entry.partition(" ")
    return TEST_LABEL.get(head, entry)


def recolor(ax):
    """Give each test its own hue regardless of how many ran on this panel."""
    curves = [ln for ln in ax.get_lines()
              if ln.get_label().partition(" ")[0] in TEST_COLOR]
    assert curves, f"no test curves found among {[l.get_label() for l in ax.get_lines()]}"
    # One scatter per curve, added in curve order. The PRC baseline is a
    # LineCollection and must not be counted among them.
    marks = [c for c in ax.collections if isinstance(c, PathCollection)]
    assert len(marks) in (0, len(curves)), \
        f"{len(marks)} threshold markers against {len(curves)} curves"
    for i, line in enumerate(curves):
        colour = TEST_COLOR[line.get_label().partition(" ")[0]]
        line.set_color(colour)
        if marks:
            marks[i].set_color(colour)


def label_baseline(ax, kind, value):
    """Name the baseline on the line itself, where its panel can be seen."""
    grey = "#444444"
    if kind == "ROC":
        # Match the diagonal's on-screen slope; the axes is not square, so 45
        # degrees would visibly miss it.
        (x0, y0), (x1, y1) = ax.transData.transform([[0, 0], [1, 1]])
        ax.text(0.63, 0.60, "chance", fontsize=8, color=grey,
                rotation=np.degrees(np.arctan2(y1 - y0, x1 - x0)),
                rotation_mode="anchor", ha="center", va="top")
    else:
        # Mid-line: the curves come down onto the right end and the legend
        # takes the left.
        ax.annotate(f"baseline {value:.3f}", xy=(0.62, value),
                    xytext=(0, 3), textcoords="offset points",
                    fontsize=8, color=grey, ha="center", va="bottom")


def main():
    from dask.distributed import Client, LocalCluster
    cluster = LocalCluster(n_workers=1, threads_per_worker=1, memory_limit="16GB")
    client = Client(cluster)
    OUT.mkdir(parents=True, exist_ok=True)

    for slug, (label, sim_stem) in DATASETS.items():
        gt = median_gt_draw(slug)
        sim = scm.de_novo_simulation(
            location=SIM_ROOT, name=f"{sim_stem}_gt{gt}", client=client)
        tests = TESTS_FOR[slug]
        for kind in ("ROC", "PRC"):
            sim.median_performance_curve(
                HYPOTHESIS_SET, kind, test_types=tests, include_alpha=True)
            fig = plt.gcf()
            fig.set_size_inches(4.4, 3.6)
            for ax in fig.axes:
                # One title, naming the regime; the generic "ROC Curve" the
                # plotting call sets is redundant beside it.
                ax.set_title(f"{label}, {kind}", fontsize=10)
                recolor(ax)
                if ax.get_legend() is None:
                    continue
                handles, labels = ax.get_legend_handles_labels()
                base = [BASELINE_LABEL.match(t) for t in labels]
                assert sum(m is not None for m in base) == 1, (
                    f"{slug} {kind}: expected exactly one baseline entry in "
                    f"{labels}")
                hit = next(m for m in base if m)
                label_baseline(ax, kind, float(hit.group(2) or 0.0))

                keep = [(h, relabel(t)) for h, t, m in zip(handles, labels, base)
                        if m is None]
                assert keep, f"{slug} {kind}: every legend entry was dropped"
                ax.legend(*zip(*keep), loc=LEGEND_LOC[kind], fontsize=8,
                          framealpha=0.95, borderpad=0.6, labelspacing=0.45,
                          handlelength=1.6)
            out = OUT / f"{slug}_median_{kind.lower()}.svg"
            fig.savefig(out, format="svg", bbox_inches="tight")
            plt.close(fig)
            print(f"    wrote {out.name}")

    client.close()
    print("done")


if __name__ == "__main__":
    main()
