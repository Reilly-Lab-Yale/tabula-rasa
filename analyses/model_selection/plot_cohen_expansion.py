#!/usr/bin/env python
"""What the coarse expansion does to Zhao et al.'s fitted means.

Replaces the 2x2 grid in cohen_obs_vs_obsingle.py, which carried its own
A./B./C./D. lettering into a manuscript panel that is itself lettered, at half
a text width where none of it could be read. Three of those four subpanels are
not referred to anywhere in the text; this is the one that is.

Zhao et al. has an element-level reporter. Expanding each detection across the
element's whole barcode set -- the coarse expansion -- drives the fitted means
three orders of magnitude below what the counts imply. Entering one zero per
delivery instead returns means that are plausible for these counts.

Reads cohen_expansion_mu.tsv, which cohen_obs_vs_obsingle.py writes from the
two orthos. No unpickling happens here.

    python analyses/model_selection/plot_cohen_expansion.py
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
TSV = BASE / "cohen_expansion_mu.tsv"

# Half of \textwidth (498.66pt = 6.90in), the width this panel is placed at.
DISPLAY_WIDTH_IN = 0.49 * 498.66 / 72.27

# Named for the expansion, matching the gloss in Fig 1 and the Results text,
# rather than for the internal fit_mode codes the orthos carry.
SERIES = [
    ("fine", "fine\n(one zero per delivery)", BLUE),
    ("coarse", "coarse\n(one zero per barcode)", ORANGE),
]


def draw_mu_hist(ax, series, nbins=46, labels="at-median"):
    """Overlaid log-x histograms of fitted means, each with its median marked.

    series: [(values, label, colour), ...]. Shared with plot_cm_family_mu.py so
    that the panels contrasting expansions and the panels contrasting count
    families get identical binning and furniture; the only thing that differs
    between them is which two populations are overlaid.

    labels="at-median" puts each label over its own median line, which reads
    directly but needs the medians to be far apart. "stacked" puts them in one
    block at the left instead, for panels whose medians are close enough that
    the two labels would other otherwise overlap.
    """
    allv = np.concatenate([v for v, _, _ in series])
    assert (allv > 0).all(), "a fitted mean of zero cannot go on a log axis"
    bins = np.logspace(np.log10(allv.min()), np.log10(allv.max()), nbins)

    for i, (v, label, color) in enumerate(series):
        med = np.median(v)
        ax.hist(v, bins=bins, color=color, alpha=0.65, linewidth=0)
        ax.axvline(med, color=color, ls=(0, (4, 3)), lw=1.1, zorder=4)
        if labels == "at-median":
            # Above the taller of the two distributions, so the two
            # annotations cannot land on each other.
            ax.annotate(f"{label}\nmedian {med:.2g}", (med, 1.0),
                        xycoords=("data", "axes fraction"),
                        xytext=(0, 4), textcoords="offset points",
                        ha="center", va="bottom", fontsize=7.5, color=color)
        else:
            # One line per series, first on top, all left-aligned above the
            # axes. The colour of the text is what keys the two histograms.
            ax.annotate(f"{label} median {med:.2g}", (0.0, 1.0),
                        xycoords="axes fraction",
                        xytext=(0, 4 + (len(series) - 1 - i) * 10),
                        textcoords="offset points",
                        ha="left", va="bottom", fontsize=7.5, color=color)

    ax.set_xscale("log")
    ax.set_xlabel("fitted mean, UMI per element", fontsize=9, color=INK)
    ax.set_ylabel("fits", fontsize=9, color=INK)
    ax.tick_params(labelsize=8, colors=MUTED, length=3)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#d9d9d9")
        ax.spines[side].set_linewidth(0.8)


def main():
    d = pd.read_csv(TSV, sep="\t")
    missing = {k for k, _, _ in SERIES} - set(d.expansion)
    assert not missing, f"missing from {TSV.name}: {sorted(missing)}"

    fig, ax = plt.subplots(figsize=(DISPLAY_WIDTH_IN, 2.9))
    draw_mu_hist(ax, [(d.loc[d.expansion == key, "mu"].values, label, color)
                      for key, label, color in SERIES])
    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(OUT / f"cohen_expansion.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)

    for key, label, _ in SERIES:
        v = d.loc[d.expansion == key, "mu"]
        print(f"{key:7s} n={len(v):5d}  median={np.median(v):.6g}  "
              f"range {v.min():.3g}-{v.max():.3g}")
    print(f"wrote {OUT/'cohen_expansion.svg'} and .png")


if __name__ == "__main__":
    main()
