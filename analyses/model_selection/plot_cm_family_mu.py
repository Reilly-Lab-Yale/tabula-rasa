#!/usr/bin/env python
"""Fitted means under consider missing, NB against ZINB. Fig. S1 panels C, D.

Panel B asks what the wrong reporter expansion does to Zhao et al.'s fitted
means. These ask the companion question of the count family: on the same fit
mode, where does ZINB put the mean that NB does not?

Consider missing enters a zero for every combination that could have been
observed, so both datasets are fit against a flood of zeros with no reporter to
say which of them were ever delivered. ZINB can attribute some of that flood to
its inflation term, which lifts the count component off the floor; NB has to
absorb all of it into the mean. The gap between the two distributions is the
size of that reattribution, and it is why an AIC preference for ZINB in this
regime is not evidence that ZINB has found something.

The mean plotted is mu as fitted. For ZINB that is the count-component mean,
conditional on the observation not being a structural zero, which is the
quantity that moves. There is no per-element marginal to plot instead: the
inflation parameter is estimated per replicate, not per element.

Reads cm_family_mu.tsv, which cm_family_mu.py writes on the cluster from the
four orthos. No unpickling happens here.

    python analyses/model_selection/plot_cm_family_mu.py
"""
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from plot_cohen_expansion import DISPLAY_WIDTH_IN, draw_mu_hist
from plot_nb_vs_zinb_bars import BLUE, ORANGE

plt.rcParams["svg.fonttype"] = "none"

BASE = pathlib.Path(__file__).resolve().parent
OUT = BASE / "output"
TSV = BASE / "cm_family_mu.tsv"

# Blue for NB and orange for ZINB, matching panel A of the same figure, where
# a fit is blue when dAIC favours NB and orange when it favours ZINB.
SERIES = [
    ("nb", "NB", BLUE),
    ("zinb", "ZINB", ORANGE),
]

# Panels are lettered by the manuscript, so the file name carries the dataset
# rather than the letter.
DATASETS = [("cohen", "Zhao et al."), ("shendure", "Lalanne et al.")]


def main():
    d = pd.read_csv(TSV, sep="\t")
    OUT.mkdir(exist_ok=True)

    for dataset, display in DATASETS:
        sub = d[d.dataset == dataset]
        assert len(sub), f"{dataset} missing from {TSV.name}"
        missing = {k for k, _, _ in SERIES} - set(sub.family)
        assert not missing, f"{dataset}: missing family {sorted(missing)}"

        fig, ax = plt.subplots(figsize=(DISPLAY_WIDTH_IN, 2.9))
        # Each panel is binned over its own range. The contrast the panel makes
        # is between the two families within one dataset, and the two datasets
        # sit six orders of magnitude apart, so a shared range would flatten
        # Zhao et al. into a couple of bins to no purpose.
        draw_mu_hist(ax, [(sub.loc[sub.family == fam, "mu"].values, label, color)
                          for fam, label, color in SERIES],
                     labels="stacked")
        fig.tight_layout()
        for ext in ("svg", "png"):
            fig.savefig(OUT / f"{dataset}_cm_family.{ext}", dpi=200,
                        bbox_inches="tight")
        plt.close(fig)

        meds = {}
        for fam, label, _ in SERIES:
            v = sub.loc[sub.family == fam, "mu"]
            meds[fam] = v.median()
            print(f"{dataset:9s} {label:4s} n={len(v):5d}  median={v.median():.4g}  "
                  f"range {v.min():.3g}-{v.max():.3g}")
        print(f"{dataset:9s} ZINB median is {meds['zinb'] / meds['nb']:.2f}x NB's\n")
    print(f"wrote {OUT}/{{cohen,shendure}}_cm_family.svg and .png")


if __name__ == "__main__":
    main()
