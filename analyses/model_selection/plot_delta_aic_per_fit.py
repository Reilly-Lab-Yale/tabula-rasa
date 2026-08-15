#!/usr/bin/env python
"""Per-fit Delta AIC behind the NB vs ZINB summary, as a supplementary panel.

Replaces the lrt_nb_vs_zinb_plot.ipynb scatter grid, which drew one small
axes per dataset-expansion-direction combination -- eighteen of them, each
with its own y-scale, labelled with the internal dataset codes. Nothing could
be compared across panels, which is the only thing a reader wants here.

The layout, display names, canonical marking and sign convention are imported
from plot_nb_vs_zinb_bars rather than restated, so this figure and the main
one cannot drift apart. That figure gives the mean per fit set; this one
gives the spread the mean was taken over.

    python analyses/model_selection/plot_delta_aic_per_fit.py
"""
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_nb_vs_zinb_bars import (LAYOUT, PANELS, TSV, BLUE, ORANGE, INK,
                                  MUTED, load_delta_aic)

plt.rcParams["svg.fonttype"] = "none"

BASE = pathlib.Path(__file__).resolve().parent
OUT = BASE / "output"

# Drawn at the width it is placed at, so the point sizes below are the sizes a
# reader sees against a 10pt body. \textwidth is 498.66pt = 6.90in.
DISPLAY_WIDTH_IN = 498.66 / 72.27

JITTER = 0.26
RNG = np.random.default_rng(0)   # reproducible jitter


def strip(ax, i, values):
    """One category: the fits as jittered points, with median and IQR."""
    if not values.size:
        return
    # Alpha by count so a set of 1344 fits reads as density and a set of 4
    # still reads as four points.
    alpha = float(np.clip(30.0 / max(values.size, 1), 0.12, 0.85))
    x = i + RNG.uniform(-JITTER, JITTER, values.size)
    colors = np.where(values > 0, BLUE, ORANGE)
    ax.scatter(x, values, s=7, c=colors, alpha=alpha, linewidths=0, zorder=2)

    # Narrower than the jittered cloud, so the median reads as a mark on the
    # points rather than as a bar of its own.
    q1, med, q3 = np.percentile(values, [25, 50, 75])
    ax.plot([i, i], [q1, q3], color=INK, lw=1.2, alpha=0.9,
            solid_capstyle="butt", zorder=3)
    ax.plot([i - 0.17, i + 0.17], [med, med], color=INK, lw=2.0,
            solid_capstyle="butt", zorder=4)


def main():
    d = load_delta_aic(TSV)

    # Stacked, not side by side: seven labels of this length need the full
    # width to be read, and stacking means one label block instead of two.
    fig, axes = plt.subplots(
        2, 1, figsize=(DISPLAY_WIDTH_IN, 5.6), sharex=True)

    for ax, (direction, title) in zip(axes, PANELS):
        sub = d[d.direction == direction]
        for i, (key, _, _) in enumerate(LAYOUT):
            v = sub[sub.dataset == key].d_nb.dropna().values
            assert v.size, f"no {direction} fits for {key}"
            strip(ax, i, v)
            # Upright: seven counts across a half-width panel collide when set
            # horizontally, and n matters here -- two of these sets have n=2.
            ax.annotate(f"n={v.size}", (i, 1.0), xycoords=("data", "axes fraction"),
                        xytext=(0, 3), textcoords="offset points",
                        ha="center", va="bottom", fontsize=7, color=MUTED)

        ax.axhline(0, color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=1)
        ax.set_yscale("symlog", linthresh=1)
        # Left-aligned over the axes rather than centred: the counts sit along
        # the top, and a centred title lands among them.
        ax.set_title(title, fontsize=10, color=INK, pad=16, loc="left")
        ax.set_xlim(-0.6, len(LAYOUT) - 0.4)
        ax.set_xticks(range(len(LAYOUT)))
        ax.set_xticklabels([f"{lab} *" if c else lab for _, lab, c in LAYOUT],
                           fontsize=8, rotation=40, ha="right",
                           rotation_mode="anchor")
        for tick, (_, _, c) in zip(ax.get_xticklabels(), LAYOUT):
            tick.set_color(INK if c else MUTED)
        ax.tick_params(axis="y", labelsize=8, colors=MUTED, length=3)
        ax.tick_params(axis="x", length=0)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color("#d9d9d9")
            ax.spines[side].set_linewidth(0.8)

        # Placed in axes coordinates: on a symlog scale the zero line is not
        # at the middle of the data range, and it sits differently in each.
        zero = ax.transLimits.transform((0, 0))[1]
        for label, frac in (("NB preferred", (zero + 1) / 2),
                            ("ZINB preferred", zero / 2)):
            ax.text(1.012, frac, label, transform=ax.transAxes, rotation=90,
                    ha="left", va="center", fontsize=7.5, color=MUTED)

    fig.supylabel(r"$\Delta$AIC per fit (ZINB $-$ NB)", fontsize=9.5, color=INK,
                  x=0.005)

    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(OUT / f"delta_aic_per_fit.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)

    for direction, _ in PANELS:
        sub = d[d.direction == direction]
        for key, lab, canon in LAYOUT:
            v = sub[sub.dataset == key].d_nb.dropna()
            print(f"{direction:13s} {lab:32s} n={len(v):5d} "
                  f"median={v.median():10.2f}  NB-favouring={100*(v>0).mean():5.1f}%")
    print(f"wrote {OUT/'delta_aic_per_fit.svg'} and .png")


if __name__ == "__main__":
    main()
