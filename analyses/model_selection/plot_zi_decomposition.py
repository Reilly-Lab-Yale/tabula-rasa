#!/usr/bin/env python
"""What zero inflation actually buys, as a supplementary panel.

Replaces the scatter emitted inline by zi_zero_decomposition.py, which drew
every cell type as a labelled point -- the labels overplotted into an
illegible block -- across a mostly empty x-range, and included takeshi, which
the manuscript does not introduce.

The argument the panel has to make: fitting an inflation component does not
change the zero mass the model predicts. Under consider-missing the fitted
inflation probability is ~0.55, two orders of magnitude above the observed
expansion, and the predicted zero rate still moves by hundredths of a
percentage point, because the negative binomial mean absorbs the difference.
That is why ZINB is not preferred.

Reads zi_zero_decomposition.tsv, which zi_zero_decomposition.py writes.

    python analyses/model_selection/plot_zi_decomposition.py
"""
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_nb_vs_zinb_bars import BLUE, ORANGE, INK, MUTED

plt.rcParams["svg.fonttype"] = "none"

BASE = pathlib.Path(__file__).resolve().parent
OUT = BASE / "output"
TSV = BASE / "zi_zero_decomposition.tsv"

# Half of \textwidth (498.66pt = 6.90in), the width this panel is placed at.
DISPLAY_WIDTH_IN = 0.49 * 498.66 / 72.27

# Only Lalanne et al. carries both expansions here, and it is the dataset the
# ZINB argument is made on. takeshi is fit upstream but is not one of the
# three datasets the manuscript introduces.
# Label anchors are per series: the two clusters sit at opposite corners, so
# one generic placement rule puts a label off the axes or on the points.
SERIES = [
    ("Shendure (obs)", "Lalanne:\nobserved, fine", BLUE, (0, 10), "center"),
    ("Shendure (CM)", "Lalanne:\nconsider missing", ORANGE, (0, 12), "center"),
]


def main():
    d = pd.read_csv(TSV, sep="\t")
    have = set(d.dataset)
    missing = {k for k, *_ in SERIES} - have
    assert not missing, f"missing from {TSV.name}: {sorted(missing)}"

    fig, ax = plt.subplots(figsize=(DISPLAY_WIDTH_IN, 2.9))

    for key, label, color, offset, ha in SERIES:
        sub = d[d.dataset == key]
        assert len(sub), f"no rows for {key}"
        ax.scatter(sub.mean_pi, sub.zinb_shift_pp, s=26, c=color,
                   alpha=0.85, linewidths=0.6, edgecolors="white", zorder=3)
        # Direct labels instead of a legend box: two clusters this well
        # separated name themselves, and a legend would cover the plot.
        ax.annotate(label, (sub.mean_pi.median(), sub.zinb_shift_pp.max()),
                    xytext=offset, textcoords="offset points",
                    ha=ha, va="bottom", fontsize=7.5, color=color)

    ax.axhline(0, color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=1)

    worst = d[d.dataset.isin([k for k, *_ in SERIES])].zinb_shift_pp.abs().max()
    # Bottom right: the cluster labels take the top, and both clusters sit at
    # or above the zero line, leaving this corner clear.
    ax.annotate(f"every fit within {worst:.2f} pp\nof its NB prediction",
                (0.97, 0.03), xycoords="axes fraction", ha="right",
                va="bottom", fontsize=7.5, color=MUTED, style="italic")

    ax.set_xlabel(r"fitted zero-inflation probability $\pi$",
                  fontsize=9, color=INK)
    ax.set_ylabel("change in predicted\nzero rate (pp)", fontsize=9, color=INK)
    ax.set_xlim(-0.04, 0.68)
    ax.set_ylim(-0.025, 0.075)
    ax.tick_params(labelsize=8, colors=MUTED, length=3)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#d9d9d9")
        ax.spines[side].set_linewidth(0.8)

    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(OUT / f"zi_decomposition.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)

    for key, label, *_ in SERIES:
        sub = d[d.dataset == key]
        print(f"{label:28s} n={len(sub):3d}  pi {sub.mean_pi.min():.3f}-{sub.mean_pi.max():.3f}"
              f"  shift {sub.zinb_shift_pp.min():+.3f} to {sub.zinb_shift_pp.max():+.3f} pp"
              f"  mu shift {sub.mu_shift.min():+.3f} to {sub.mu_shift.max():+.3f}")
    print(f"wrote {OUT/'zi_decomposition.svg'} and .png")


if __name__ == "__main__":
    main()
