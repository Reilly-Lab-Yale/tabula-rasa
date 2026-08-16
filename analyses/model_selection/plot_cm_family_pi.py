#!/usr/bin/env python
"""Where the zeros went: fitted inflation parameter against the zero mass.

Fig. S1 panels B-E. Each panel is one dataset under two zero regimes, ZINB
only, showing the fitted inflation parameter pi per (CRE, replicate).

Each label carries the fraction of that fit's observations that were zero, so
the reader can see how much zero mass the expansion introduced alongside what
pi did about it. The two are deliberately NOT drawn on the same axis: a plain
NB generates zeros of its own, so pi is only ever meant to cover the excess
over that, and putting a zero-fraction rule on the pi axis would invite a
comparison that overstates what pi should be.

What the panels do compare is pi against pi: the shift between the two
regimes, against the extra zero mass the flooded regime took on. Where a large
increase in zeros moves pi hardly at all, the inflation term is not
identifiable and the zeros are going into the count mean instead.

All four panels share a fixed 0 to 1 axis, pi being a probability.

Reads cm_family_pi.tsv and cm_family_zerofrac.tsv, which cm_family_pi.py writes
on the cluster. No unpickling happens here.

    python analyses/model_selection/plot_cm_family_pi.py
"""
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_cohen_expansion import DISPLAY_WIDTH_IN
from plot_nb_vs_zinb_bars import BLUE, ORANGE, INK, MUTED

plt.rcParams["svg.fonttype"] = "none"

BASE = pathlib.Path(__file__).resolve().parent
OUT = BASE / "output"
PI_TSV = BASE / "cm_family_pi.tsv"
ZERO_TSV = BASE / "cm_family_zerofrac.tsv"

NBINS = 40

# (output stem, dataset, baseline regime, flooded regime). Blue is always the
# regime with less zero mass, orange the one that adds it, as in panel A.
PANELS = [
    ("cohen_pi_coarse", "cohen", "fine", "coarse"),
    ("cohen_pi_cm", "cohen", "fine", "consider missing"),
    ("shendure_pi_cm", "shendure", "observed", "consider missing"),
    ("seelig_pi_cm", "seelig", "consider missing + MOI", "consider missing"),
]


def implied_pi(base, flood):
    """What pi would be if every zero the expansion added were structural.

    The expansion asserts that the zeros it introduces are combinations that
    were never delivered, so they belong to the inflation component by
    construction. Carrying the baseline fit's own structural zeros across and
    adding the new ones gives the value pi would take if the fit accepted that
    assertion in full:

        pi_implied = (pi_base * T_base + delta_zeros) / T_flood

    The ratio of the fitted pi to this is what the panels report. It is a
    like-for-like comparison in a way that differencing the two zero fractions
    is not, because the expansion changes the number of observations as well as
    the number of zeros -- the two fractions are not taken over the same
    universe, and differencing them can imply an uptake above 100%.
    """
    delta = flood["n_zeros"] - base["n_zeros"]
    assert delta > 0, "the flooded regime must add zeros"
    return (base["pi"] * base["n_total"] + delta) / flood["n_total"]


def draw(ax, series):
    """series: [dict(pi_values, pi, zero_fraction, n_zeros, n_total, label,
    colour), ...] with the baseline first."""
    bins = np.linspace(0, 1, NBINS + 1)
    for i, s in enumerate(series):
        ax.hist(s["pi_values"], bins=bins, color=s["colour"], alpha=0.65,
                linewidth=0)
        ax.axvline(s["pi"], color=s["colour"], ls=(0, (4, 3)), lw=1.1, zorder=4)
        # Three decimals on the percentage: several of these regimes round to
        # 100% at one, and the gap between 99.958 and 99.985 separates two of
        # the panels.
        ax.annotate(f"{s['label']}: $\\pi$ {s['pi']:.3f}, "
                    f"{s['zero_fraction']:.3%} zeros",
                    (0.0, 1.0), xycoords="axes fraction",
                    xytext=(0, 14 + (len(series) - 1 - i) * 10),
                    textcoords="offset points",
                    ha="left", va="bottom", fontsize=7.5, color=s["colour"])

    exp = implied_pi(series[0], series[1])
    # Four decimals, for the same reason the zero fractions get three: every
    # implied value here rounds to 1.000 at three, which is the rounding that
    # made these panels confusing in the first place.
    ax.annotate(f"implied $\\pi$ {exp:.4f}, fitted {series[1]['pi']:.4f} "
                f"({series[1]['pi'] / exp:.0%} uptake)",
                (0.0, 1.0), xycoords="axes fraction",
                xytext=(0, 4), textcoords="offset points",
                ha="left", va="bottom", fontsize=7.5, color=INK)

    ax.set_xlim(0, 1)
    ax.set_xlabel("zero-inflation parameter $\\pi$", fontsize=9, color=INK)
    ax.set_ylabel("fits", fontsize=9, color=INK)
    ax.tick_params(labelsize=8, colors=MUTED, length=3)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#d9d9d9")
        ax.spines[side].set_linewidth(0.8)
    return exp


def main():
    pi = pd.read_csv(PI_TSV, sep="\t")
    zero = pd.read_csv(ZERO_TSV, sep="\t").set_index(["dataset", "regime"])
    OUT.mkdir(exist_ok=True)

    for stem, dataset, base_regime, flood_regime in PANELS:
        series = []
        for regime, color in ((base_regime, BLUE), (flood_regime, ORANGE)):
            v = pi[(pi.dataset == dataset) & (pi.regime == regime)]["pi"].values
            assert len(v), f"{dataset}/{regime} missing from {PI_TSV.name}"
            row = zero.loc[(dataset, regime)]
            series.append({"pi_values": v, "pi": float(np.median(v)),
                           "zero_fraction": float(row.zero_fraction),
                           "n_zeros": int(row.n_zeros),
                           "n_total": int(row.n_total),
                           "label": regime, "colour": color})

        # The premise of every panel is that the second regime adds zeros.
        assert series[1]["zero_fraction"] > series[0]["zero_fraction"], (
            f"{stem}: {flood_regime} is not the flooded regime")
        # Only zeros are added; the real observations are the same in both.
        nz = [s["n_total"] - s["n_zeros"] for s in series]
        assert nz[0] == nz[1], (
            f"{stem}: nonzero observations differ between regimes "
            f"({nz[0]:,} vs {nz[1]:,}); the expansion changed more than zeros")

        fig, ax = plt.subplots(figsize=(DISPLAY_WIDTH_IN, 2.9))
        exp = draw(ax, series)
        fig.tight_layout()
        for ext in ("svg", "png"):
            fig.savefig(OUT / f"{stem}.{ext}", dpi=200, bbox_inches="tight")
        plt.close(fig)

        print(f"{stem}:")
        for s in series:
            v = s["pi_values"]
            print(f"    {s['label']:24s} n={len(v):5d}  pi median={s['pi']:.4f}  "
                  f"IQR {np.percentile(v, 25):.4f}-{np.percentile(v, 75):.4f}  "
                  f"zeros={s['zero_fraction']:.5f} "
                  f"({s['n_zeros']:,} of {s['n_total']:,})")
        print(f"    added {series[1]['n_zeros'] - series[0]['n_zeros']:,} zeros; "
              f"implied pi={exp:.4f}, fitted {series[1]['pi']:.4f} "
              f"-> {series[1]['pi'] / exp:.1%} uptake\n")
    print(f"wrote {len(PANELS)} panels to {OUT}")


if __name__ == "__main__":
    main()
