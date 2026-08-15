#!/usr/bin/env python
"""NB vs ZINB across datasets and expansions, canonical fits marked.

Replaces the lrt_nb_vs_zinb_plot.ipynb barchart. Same underlying table, but
built for the manuscript rather than for reading at the bench:

  - the sign convention matches plot_nb_vs_zinb_shendure.py, so positive
    favours NB in both figures
  - canonical fits are distinguished from the expansions that were not
    selected, which is the point of the panel: NB wins where it matters, and
    the expansions favouring ZINB are ones nobody analyses. Being
    non-canonical does not by itself make an expansion a counterfactual
    regime of interest; it only means QC did not select it.
  - display names are the published first authors, and text stays editable

Reads lrt_nb_vs_zinb_results.tsv, which lrt_nb_vs_zinb.py writes. No fitting
happens here.

    python analyses/model_selection/plot_nb_vs_zinb_bars.py
"""
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# Editable <text> so the manuscript sync pipeline's transforms apply, matching
# the other figures.
plt.rcParams["svg.fonttype"] = "none"

BASE = pathlib.Path(__file__).resolve().parent
OUT = BASE / "output"
TSV = BASE / "lrt_nb_vs_zinb_results.tsv"

BLUE, ORANGE = "#0072b2", "#d55e00"     # Okabe-Ito (house palette)
INK, MUTED = "#1a1a1a", "#6b6b6b"

# Canonical fit per dataset, then that dataset's other expansions. Canonical
# is decided by fit quality at QC; the expansion name is only the handle.
# takeshi is fit upstream but is not one of the three datasets the manuscript
# introduces, so it is excluded here.
# Surnames only: the labels are rotated and seven of them have to fit across
# a single-column panel. The caption carries the full citations.
# Labels spell the rule out rather than using the internal codes. "obs"
# appearing under two datasets was the confusing part: it names a rule that
# any reporter-bearing design can be fit under, while its correctness depends
# on whether that reporter resolves barcodes. Naming the granularity makes the
# mismatch legible -- Zhao et al. has an element-level reporter, so the coarse
# expansion under it is visibly the wrong rule rather than a second option.
#
# "fine" is one zero per delivery event as the reporter records it
# (core.py reporter_expansion="single", and shendure's preexisting zeros);
# "coarse" is one zero per barcode of the element (reporter_expansion="coarse").
LAYOUT = [
    ("shendure_obs",      "Lalanne: observed, fine",      True),
    ("shendure_cm",       "Lalanne: consider missing",    False),
    ("cohen_obsingle",    "Zhao: observed, fine",         True),
    ("cohen_obs",         "Zhao: observed, coarse",       False),
    ("cohen_cm",          "Zhao: consider missing",       False),
    ("seelig_cm_moib",    "Yin: consider missing + MOI",  True),
    ("seelig_cm",         "Yin: consider missing",        False),
]
PANELS = [("by_cre", "by CRE"), ("by_cell_type", "by cell type")]


def stars(p, n):
    if n < 2 or np.isnan(p):
        return "n/a"
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def main():
    d = pd.read_csv(TSV, sep="\t")
    d[["dataset", "direction"]] = d.comparison.str.split(" / ", expand=True)
    # Published convention is AIC_NB - AIC_ZINB; flip so positive favours NB,
    # matching the by-fit figure. Recomputed from the AIC columns rather than
    # negating delta_aic, so a change in the upstream column cannot pass
    # through silently.
    d["d_nb"] = d.aic_zinb - d.aic_nb
    assert np.allclose(d.d_nb, -d.delta_aic, equal_nan=True), \
        "delta_aic is not AIC_NB - AIC_ZINB; the sign flip would be wrong"

    wanted = {k for k, _, _ in LAYOUT}
    missing = wanted - set(d.dataset)
    assert not missing, f"missing from results table: {sorted(missing)}"

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.2))

    for ax, (direction, title) in zip(axes, PANELS):
        sub = d[d.direction == direction]
        x = np.arange(len(LAYOUT))
        for i, (key, _, canonical) in enumerate(LAYOUT):
            v = sub[sub.dataset == key].d_nb.dropna().values
            assert v.size, f"no {direction} fits for {key}"
            mean, sem = v.mean(), (stats.sem(v) if v.size > 1 else 0.0)
            p = stats.ttest_1samp(v, 0).pvalue if v.size > 1 else np.nan
            colour = BLUE if mean > 0 else ORANGE
            ax.bar(i, mean, width=0.68, zorder=3,
                   color=colour if canonical else "white",
                   edgecolor=colour, linewidth=1.1,
                   hatch=None if canonical else "///")
            ax.errorbar(i, mean, yerr=sem, color=INK, lw=0.9, capsize=2.5, zorder=4)
            # Marks sit on the far side of the bar from zero. Scaling the bar
            # end outward works in both directions on a symlog axis because
            # multiplying a negative by >1 makes it more negative.
            # Adjacent bars of near-equal height put these annotations at the
            # same y, where "n=1344" and its neighbour collide into an
            # unreadable run. Push every other one further out so the pair
            # separates vertically instead.
            lift = 1.7 if i % 2 == 0 else 2.8
            ax.text(i, (mean + np.sign(mean) * sem) * lift,
                    f"{stars(p, v.size)}\nn={v.size}", ha="center",
                    va="bottom" if mean > 0 else "top",
                    fontsize=7, color=MUTED, linespacing=1.3)

        ax.axhline(0, color=MUTED, lw=1.0, ls="--", zorder=2)
        ax.set_yscale("symlog", linthresh=1)
        ax.set_xticks(x)
        # Canonical marked on the tick label, not above the bar, so it cannot
        # be confused with the significance marks.
        ax.set_xticklabels([f"{lab} *" if c else lab for _, lab, c in LAYOUT],
                           fontsize=8, rotation=40, ha="right",
                           rotation_mode="anchor")
        for tick, (_, _, c) in zip(ax.get_xticklabels(), LAYOUT):
            if c:
                tick.set_color(INK)
                tick.set_fontweight("bold")
        ax.set_title(title, fontsize=10, color=INK, pad=16)
        ax.grid(True, axis="y", color="#e6e6e6", lw=0.8, zorder=0)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color("#cccccc")
        ax.tick_params(colors=MUTED, labelsize=8, length=0)
        # Headroom for the significance marks, which sit outside the bar.
        # Wider below: the downward marks are two lines hanging from the bar
        # end, so they need more room than the upward ones.
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo * 12, hi * 3.2)

        # Which side means what, stated outright, so the panel can be read
        # without working back through the sign of an AIC difference. Placed
        # outside the right spine rather than over the plot: every column
        # holds a bar, so an in-axes watermark would sit behind one. y=0 is
        # located in axes coordinates because on a symlog scale the zero line
        # is nowhere near the middle.
        zero = ax.transAxes.inverted().transform(
            ax.transData.transform((0, 0)))[1]
        for label, frac in (("NB preferred", (zero + 1) / 2),
                            ("ZINB preferred", zero / 2)):
            ax.text(1.02, frac, label, transform=ax.transAxes,
                    ha="center", va="center", rotation=270,
                    fontsize=8.5, fontweight="bold", color=MUTED)

    axes[0].set_ylabel(r"mean $\Delta$AIC (ZINB $-$ NB)", fontsize=9.5, color=INK)

    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=MUTED, edgecolor=MUTED,
                      label="* canonical fit (solid)"),
        plt.Rectangle((0, 0), 1, 1, facecolor="white", edgecolor=MUTED, hatch="///",
                      label="non-canonical expansion (hatched)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               fontsize=8, bbox_to_anchor=(0.5, -0.02))

    for direction, _ in PANELS:
        sub = d[d.direction == direction]
        canon = [k for k, _, c in LAYOUT if c]
        won = [k for k in canon if sub[sub.dataset == k].d_nb.mean() > 0]
        print(f"{direction:13s}: canonical favouring NB on the mean: "
              f"{len(won)}/{len(canon)} {won}")

    fig.tight_layout(rect=(0, 0.06, 1, 1))
    OUT.mkdir(exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(OUT / f"nb_vs_zinb_bars.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT/'nb_vs_zinb_bars.svg'} and .png")


if __name__ == "__main__":
    main()
