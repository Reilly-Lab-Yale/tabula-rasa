#!/usr/bin/env python
"""Null FPR per cell type, t-test against MWU, one figure per dataset.

Supersedes the grouped bars in fpr_bars.py for the manuscript panel. The
quantity a reader needs is the distance from nominal 0.05, and everything
happens between 0.04 and 0.08. Bars anchored at zero spend four fifths of
their length on empty range and squeeze the entire signal into the top of the
axis, and a truncated bar is not an option -- a bar's length is its value.
Dots carry position rather than length, so the axis can start where the data
does, and the segment joining each pair reads as the gap between the tests.

Cell types run down the y-axis, which also gets rid of the 45-degree tick
labels the bar version needed for names like EpiblastPrimitiveStreak.

Reads the cached null-pvals parquets written by the per-dataset calibration
scripts, so this runs anywhere the parquets are, no cluster needed.

    python analyses/simulation/activity_calibration/fpr_dumbbell.py
"""
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent

ALPHA = 0.05

# Okabe-Ito, the palette scripts/figure_styling.py remaps everything else to.
# This pair passes the categorical checks outright, including contrast against
# the page, so no secondary encoding is being leaned on to tell them apart.
TEST_COLOR = {"ttest": "#d55e00", "mwu": "#0072b2"}
TEST_LABEL = {"ttest": "t-test", "mwu": "MWU"}
INK, MUTED, HAIR = "#1a1a1a", "#666666", "#c9c9c9"

DATASETS = {
    "shendure": ("Lalanne et al.", "Pluripotent"),
    "cohen": ("Zhao et al.", "Rod"),
    "seelig": ("Yin et al.", "HepG2"),
}
CONDITIONS = [("reporter", "reporter used"), ("deflated", "reporter withheld")]

plt.rcParams.update({
    "svg.fonttype": "none",       # keep text editable downstream
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.unicode_minus": False,  # the repo is plain ASCII
})


def fpr_table(dataset):
    """(condition, test) -> {cell type: FPR}, for whichever files exist."""
    out = {}
    for cond, _ in CONDITIONS:
        for test in ("ttest", "mwu"):
            tag = "" if test == "ttest" else f"{test}_"
            p = HERE / dataset / "output" / f"{dataset}_null_pvals_{tag}{cond}.parquet"
            if not p.is_file():
                continue
            d = pd.read_parquet(p)
            g = d.groupby("cell_type", observed=True)["p_value"]
            out[(cond, test)] = g.apply(lambda s: float((s < ALPHA).mean())).to_dict()
    assert out, f"{dataset}: no null-pvals parquets found"
    # Every condition present must carry both tests, or a dumbbell has one end.
    conds = {c for c, _ in out}
    for c in conds:
        assert (c, "ttest") in out and (c, "mwu") in out, \
            f"{dataset}: condition {c} is missing one of the two tests"
    return out


def order(table):
    """Reference cell type first, then worst calibration at the top.

    Ranked on the t-test without a reporter where that condition exists, since
    that is the regime the figure is about; the reference row is pinned so the
    same row is in the same place in every dataset.
    """
    cts = {ct for d in table.values() for ct in d}
    rank_on = table.get(("deflated", "ttest")) or table[("reporter", "ttest")]
    rest = sorted((c for c in cts if c != "reference"),
                  key=lambda c: -abs(rank_on.get(c, ALPHA) - ALPHA))
    return (["reference"] if "reference" in cts else []) + rest


def draw(ax, table, cond, cts, labels, xlim, show_ylabels, marks=True):
    """cts are the keys in the data, labels are what the reader sees.

    marks=False keeps the axis furniture but plots nothing, for a condition a
    dataset never ran. The row labels then stay in the same column in every
    panel instead of jumping to whichever facet happens to carry data.
    """
    assert len(cts) == len(labels), "one label per cell type"
    y = np.arange(len(cts))[::-1]  # first listed cell type at the top
    # The dashed line is nominal alpha. It is left unlabelled: it sits on the
    # 0.05 tick, and the axis label and caption both name the threshold.
    if marks:
        ax.axvline(ALPHA, color=MUTED, ls="--", lw=0.9, zorder=1)

    if marks:
        lo = [table[(cond, "ttest")][c] for c in cts]
        hi = [table[(cond, "mwu")][c] for c in cts]
        # The connector is the gap between the two tests, so it is drawn under
        # the dots and kept lighter than either of them.
        ax.hlines(y, lo, hi, color=HAIR, lw=2, zorder=2)
        for test, vals in (("ttest", lo), ("mwu", hi)):
            ax.scatter(vals, y, s=46, color=TEST_COLOR[test],
                       label=TEST_LABEL[test], zorder=3,
                       edgecolors="white", linewidths=1.2, clip_on=False)

    ax.set_yticks(y)
    ax.set_yticklabels(labels if show_ylabels else [])
    ax.set_ylim(-0.7, len(cts) - 0.3)
    ax.set_xlim(*xlim)
    ax.grid(axis="x", color=HAIR, lw=0.5, alpha=0.6, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(HAIR)
    ax.tick_params(length=0, colors=INK)


def figure(dataset, table, xlim):
    display, ref_display = DATASETS[dataset]
    cts = order(table)
    labels = [ref_display if c == "reference" else c for c in cts]

    # Always two columns, even when a dataset has only one condition, so that
    # a dumbbell of a given length means the same thing in all three panels.
    conds = [(c, lab) for c, lab in CONDITIONS if (c, "ttest") in table]

    # Roughly a quarter inch per row on top of fixed furniture, so a two-row
    # dataset does not get the same slab of page as a ten-row one. The margins
    # are set in inches and converted, rather than left to tight_layout, so
    # that the three panels stack with their plot areas in line.
    # Sized so all three stack inside one page with the caption: three sets of
    # axis furniture is most of the cost, and the float overflows if this grows.
    top_in, bottom_in, pitch = 0.60, 0.46, 0.22
    height = top_in + bottom_in + pitch * len(cts)
    fig, axes = plt.subplots(1, len(CONDITIONS), figsize=(6.5, height),
                             squeeze=False)
    for i, (cond, lab) in enumerate(CONDITIONS):
        ax = axes[0][i]
        ran = (cond, "ttest") in table
        draw(ax, table, cond, cts, labels, xlim, show_ylabels=(i == 0),
             marks=ran)
        ax.set_title(lab, color=INK if ran else MUTED, pad=8)
        if not ran:
            # Yin et al. ran no transfection reporter, so this condition does
            # not exist for it. Say so in the slot rather than letting the one
            # condition that did run spread across the whole panel.
            ax.text(0.5, 0.5, "not applicable\nno transfection reporter",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=8, color=MUTED, linespacing=1.5,
                    bbox=dict(facecolor="white", edgecolor="none", pad=3))

    fig.subplots_adjust(left=0.25, right=0.995, wspace=0.07,
                        top=1 - top_in / height, bottom=bottom_in / height)
    fig.supxlabel(f"false-positive rate at p < {ALPHA:g}", y=0.012,
                  fontsize=9, color=INK)

    drawn = next(ax for ax in axes[0] if ax.get_legend_handles_labels()[1])
    handles = {lb: h for lb, h in zip(*reversed(drawn.get_legend_handles_labels()))}
    fig.legend([handles[TEST_LABEL[t]] for t in ("ttest", "mwu")],
               [TEST_LABEL[t] for t in ("ttest", "mwu")],
               loc="upper right", bbox_to_anchor=(0.995, 0.999), ncol=2,
               frameon=False, handletextpad=0.25, columnspacing=1.4)
    # Indented, not flush left: the manuscript overlays a bold panel letter in
    # the top-left corner via \panel{}, and a flush-left title lands under it.
    fig.suptitle(display, x=0.055, y=0.995, ha="left", va="top",
                 fontsize=10, color=INK)
    return fig


def main():
    tables = {ds: fpr_table(ds) for ds in DATASETS}

    # One x range for all three panels. On per-panel ranges Yin et al.'s
    # spread of 0.047 to 0.055 would be drawn as wide as Lalanne et al.'s
    # 0.042 to 0.076, which is the opposite of what the figure is claiming.
    vals = [v for t in tables.values() for d in t.values() for v in d.values()]
    pad = 0.08 * (max(vals) - min(vals))
    xlim = (min(vals) - pad, max(vals) + pad)
    assert xlim[0] < ALPHA < xlim[1], f"nominal {ALPHA} falls outside {xlim}"
    print(f"shared x range {xlim[0]:.4f}-{xlim[1]:.4f} over {len(vals)} cells")

    for dataset in DATASETS:
        fig = figure(dataset, tables[dataset], xlim)
        out = HERE / dataset / "output" / f"{dataset}_calibration_fpr_dumbbell.svg"
        # No bbox_inches="tight". The three panels stack in the manuscript and
        # the margins above are set so their plot areas line up; cropping each
        # one to its own contents would undo that, since the longest cell-type
        # name differs per dataset and the crop follows it.
        fig.savefig(out, format="svg")
        plt.close(fig)
        print(f"wrote {out.relative_to(HERE)}")


if __name__ == "__main__":
    main()
