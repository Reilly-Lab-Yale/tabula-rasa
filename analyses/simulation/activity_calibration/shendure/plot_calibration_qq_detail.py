#!/usr/bin/env python
"""Per-cell-type QQ detail for the Lalanne et al. null calibration, MWU.

One small panel per cell type, faceted 5x2. Companion to the single-axis
overlay in ../qq_overlay.py, which puts every cell type on one pair of axes;
this figure is the per-cell-type breakdown, so it stays faceted.

The deviation form, the Kolmogorov band and the curve downsampling are shared
with the overlay and documented in ../qq_common.py. `--form classic` draws the
undifferenced [0, 1] QQ for comparison.

    python plot_calibration_qq_detail.py                       # reporter, deviation
    python plot_calibration_qq_detail.py --form classic        # undifferenced QQ
    python plot_calibration_qq_detail.py --condition both      # both conditions
"""
import argparse
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "output"

sys.path.insert(0, str(HERE.parent))
from qq_common import (  # noqa: E402
    KS_CRIT_95, INK, MUTED, HAIR, BAND, RC_BASE,
    check_reduction, curve, expected, reduce_curve,
)

TEST = "mwu"
ALPHA = 0.05

# Code identifiers upstream, reader-facing labels here. "reference" is the
# cell type the synthetic CREs are calibrated against, pluripotent cells in
# the mouse embryoid bodies; it is a real category, not a placeholder.
# Newlines are panel-title line breaks, not part of the name.
CELL_TYPE_LABEL = {
    "reference": "Pluripotent\n(reference)",
    "Cardiomyocytes": "Cardiomyocytes",
    "EpiblastPrimitiveStreak": "Epiblast /\nprimitive streak",
    "ExEndodermParietal": "ExE endoderm\n(parietal)",
    "ExEndodermVisceral": "ExE endoderm\n(visceral)",
    "Haematoendothelial": "Haemato-\nendothelial",
    "Mesoderm": "Mesoderm",
    "NeuroectodermBrain": "Neuroectoderm\n(brain)",
    "NeuroectodermRostral": "Neuroectoderm\n(rostral)",
    "SurfaceEctoderm": "Surface\nectoderm",
}

# Okabe-Ito. One curve per panel in the default figure, so the colour is doing
# no encoding work and there is nothing to key -- the panel title names the
# cell type. The second colour is only used by --condition both, where the two
# conditions are a real two-level categorical.
CONDITION_COLOR = {"reporter": "#0072b2", "deflated": "#d55e00"}
CONDITION_LABEL = {"reporter": "reporter used", "deflated": "reporter withheld"}

NCOLS, NROWS = 5, 2

# Set here rather than left to a downstream scaler: the manuscript sync applies
# no font scaling, so a point size in this file is a point size on the page.
# The figure is drawn at \textwidth (498.66pt = 6.90in) and included unscaled.
FIG_W = 6.90
plt.rcParams.update({
    **RC_BASE,
    "font.size": 8,
    "axes.labelsize": 8.5,
    "axes.titlesize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
})


def load(condition):
    """{cell type: sorted finite p-values} for one condition."""
    path = OUT / f"shendure_null_pvals_{TEST}_{condition}.parquet"
    assert path.is_file(), f"missing null p-values: {path}"
    df = pd.read_parquet(path)
    n_rows = len(df)

    out = {}
    for ct, g in df.groupby("cell_type", observed=True):
        v = g["p_value"].to_numpy()
        finite = np.isfinite(v)
        assert finite.all(), f"{ct}: {(~finite).sum()} non-finite p-values of {len(v)}"
        assert v.min() >= 0.0 and v.max() <= 1.0, \
            f"{ct}: p-values outside [0,1], range {v.min()}-{v.max()}"
        out[ct] = np.sort(v)

    assert sum(len(v) for v in out.values()) == n_rows, \
        f"groupby lost rows: {sum(len(v) for v in out.values())} of {n_rows}"
    assert set(out) == set(CELL_TYPE_LABEL), \
        f"cell types in data {sorted(out)} != labelled {sorted(CELL_TYPE_LABEL)}"
    assert len(out) == NCOLS * NROWS, \
        f"{len(out)} cell types will not fill a {NCOLS}x{NROWS} grid"
    return out


def ks_table(condition):
    """{cell type: (D, p)} from the cached summary, keyed by name not position."""
    path = OUT / f"shendure_null_summary_{TEST}.parquet"
    assert path.is_file(), f"missing null summary: {path}"
    s = pd.read_parquet(path)
    s = s[s["condition"] == condition]
    assert len(s) == NCOLS * NROWS, \
        f"summary has {len(s)} rows for condition {condition}, expected {NCOLS * NROWS}"
    assert s["cell_type"].is_unique, "summary has duplicate cell types in one condition"
    assert ((s["ks_p"] >= 0) & (s["ks_p"] <= 1)).all(), "KS p-values outside [0,1]"
    return {r.cell_type: (r.ks_d, r.ks_p, r.n) for r in s.itertuples()}


LINEWIDTH_PT = 1.3
RDP_FRACTION = 0.30   # simplification tolerance, as a fraction of the linewidth


def fmt_p(p):
    return f"{p:.3f}" if p >= 1e-3 else f"{p:.0e}"


def order(ks):
    """Reference pinned first, then worst calibration first.

    Matches fpr_dumbbell.py, so a cell type sits in the same relative place in
    both figures, and puts the two cell types that actually fail in the top row.
    """
    rest = sorted((c for c in ks if c != "reference"), key=lambda c: -ks[c][0])
    return ["reference"] + rest


def figure(pvals_by_cond, ks_by_cond, form):
    conds = list(pvals_by_cond)
    cts = order(ks_by_cond[conds[0]])

    if form == "classic":
        ylim, ylabel = (0.0, 1.0), "observed p-value"
    else:
        span = max(
            abs(v - expected(len(v))).max()
            for p in pvals_by_cond.values() for v in p.values()
        )
        lim = 1.25 * span
        ylim = (-lim, lim)
        ylabel = "observed - expected p-value"

    # Margins in inches, then converted, so the panels come out at a size this
    # file chose rather than whatever tight_layout arrives at.
    left, right, bottom, top = 0.66, 0.05, 0.52, 0.42
    wgap, hgap = 0.13, 0.46
    panel_w = (FIG_W - left - right - wgap * (NCOLS - 1)) / NCOLS
    panel_h = 1.05
    fig_h = top + bottom + panel_h * NROWS + hgap * (NROWS - 1)

    # Simplification tolerance in data units: a fraction of the drawn
    # linewidth, given how many data units a point of panel height is worth.
    stroke_y = LINEWIDTH_PT / 72 / panel_h * (ylim[1] - ylim[0])
    tol_y = RDP_FRACTION * stroke_y
    for cond in conds:
        check_reduction(pvals_by_cond[cond], form, tol_y, 0.5 * stroke_y,
                        linewidth_pt=LINEWIDTH_PT)

    fig, axes = plt.subplots(NROWS, NCOLS, figsize=(FIG_W, fig_h), squeeze=False)
    fig.subplots_adjust(
        left=left / FIG_W, right=1 - right / FIG_W,
        bottom=bottom / fig_h, top=1 - top / fig_h,
        wspace=wgap / panel_w, hspace=hgap / panel_h,
    )

    for k, ct in enumerate(cts):
        ax = axes[k // NCOLS][k % NCOLS]
        first_col, last_row = k % NCOLS == 0, k // NCOLS == NROWS - 1

        # Kolmogorov acceptance region. n is per condition but differs by at
        # most one row between them, so the first condition sets the band.
        n = len(pvals_by_cond[conds[0]][ct])
        crit = KS_CRIT_95 / np.sqrt(n)
        if form == "classic":
            grid = np.linspace(0, 1, 200)
            ax.fill_between(grid, np.clip(grid - crit, 0, 1), np.clip(grid + crit, 0, 1),
                            color=BAND, lw=0, zorder=0)
            ax.plot([0, 1], [0, 1], color=MUTED, ls="--", lw=0.8, zorder=2)
        else:
            ax.axhspan(-crit, crit, color=BAND, lw=0, zorder=0)
            ax.axhline(0.0, color=MUTED, ls="--", lw=0.8, zorder=2)

        for cond in conds:
            xr, yr = reduce_curve(pvals_by_cond[cond][ct], form, tol_y)
            ax.plot(xr, yr, color=CONDITION_COLOR[cond], lw=LINEWIDTH_PT,
                    solid_capstyle="round", solid_joinstyle="round", zorder=3,
                    label=CONDITION_LABEL[cond] if k == 0 else None)

        ax.set_title(CELL_TYPE_LABEL[ct], color=INK, pad=3.5, linespacing=1.15)

        # KS result for the panel, on the top edge. The deviation curve is a
        # bridge pinned to zero at both ends and peaks well under half the
        # y-range, so the top strip is free in every panel. One line per curve
        # drawn, in the same order as the legend; a single line would read as
        # the result for both conditions.
        for j, cond in enumerate(conds):
            d, p, n_summary = ks_by_cond[cond][ct]
            n_cond = len(pvals_by_cond[cond][ct])
            assert n_summary == n_cond, (
                f"{ct} ({cond}): summary n={n_summary} does not match "
                f"{n_cond} p-values on disk"
            )
            if len(conds) == 1:
                color = INK if p < ALPHA else MUTED
            else:
                color = CONDITION_COLOR[cond]  # keyed by the legend
            ax.text(0.5, 0.975 - 0.085 * j, f"D={d:.4f}  p={fmt_p(p)}",
                    transform=ax.transAxes, ha="center", va="top",
                    fontsize=6.5, color=color)

        ax.set_xlim(0, 1)
        ax.set_ylim(*ylim)
        ax.set_xticks([0, 0.5, 1])
        if form != "classic":
            ax.set_yticks([-0.02, 0, 0.02])
        ax.set_xticklabels(["0", "0.5", "1"] if last_row else [])
        if not first_col:
            ax.set_yticklabels([])
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(HAIR)
        ax.tick_params(length=2.5, width=0.6, color=HAIR, colors=INK, pad=2)

    # One direct label on the band instead of a legend entry: it is furniture,
    # it appears in all ten panels, and the caption defines it. The leader ends
    # on the band's lower edge, which is below the curve in this panel, so it
    # does not cross the data on its way there. Not drawn on the classic form,
    # where the band is narrower than the line and there is nothing to point at.
    if form != "classic":
        crit0 = KS_CRIT_95 / np.sqrt(len(pvals_by_cond[conds[0]][cts[0]]))
        axes[0][0].annotate(
            "95% KS band", xy=(0.55, -crit0),
            xytext=(0.03, 0.06), textcoords="axes fraction",
            fontsize=6, color=MUTED, ha="left", va="bottom",
            arrowprops=dict(arrowstyle="-", lw=0.5, color=MUTED,
                            shrinkA=1.5, shrinkB=1.0),
        )

    fig.supxlabel("expected p-value (uniform quantile)", y=0.018, fontsize=8.5, color=INK)
    fig.supylabel(ylabel, x=0.011, fontsize=8.5, color=INK)
    if form != "classic":
        # Which way is bad. Written as one note against the zero line rather
        # than a conservative/anti-conservative pair at the two ends of the
        # axis: the rows share one y-axis, so a pair pinned to the figure edge
        # reads as belonging to whichever row it happens to sit beside, and a
        # pair repeated per row does not fit in a panel this short.
        fig.text(0.031, (bottom / fig_h + (1 - top / fig_h)) / 2,
                 "below zero: anti-conservative", rotation=90,
                 ha="left", va="center", fontsize=6.5, color=MUTED)

    if len(conds) > 1:
        h, lb = axes[0][0].get_legend_handles_labels()
        fig.legend(h, lb, loc="lower right", bbox_to_anchor=(1 - right / FIG_W, 0.0),
                   ncol=2, frameon=False, fontsize=7.5, handletextpad=0.4,
                   columnspacing=1.2, handlelength=1.6)

    return fig, fig_h


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--condition", default="reporter",
                    choices=["reporter", "deflated", "both"])
    ap.add_argument("--form", default="deviation", choices=["deviation", "classic"])
    args = ap.parse_args()

    conds = ["reporter", "deflated"] if args.condition == "both" else [args.condition]
    pvals = {c: load(c) for c in conds}
    ks = {c: ks_table(c) for c in conds}

    fig, fig_h = figure(pvals, ks, args.form)

    tag = "" if args.form == "deviation" else f"_{args.form}"
    out = OUT / f"shendure_calibration_{TEST}_qq_detail_{args.condition}{tag}.svg"
    # No bbox_inches="tight": the margins above are in inches so that the drawn
    # width is the width the manuscript includes it at, and cropping to content
    # would shrink it by however much slack the outermost text happens to leave.
    fig.savefig(out, format="svg")
    print(f"saved {out}  ({FIG_W:.2f} x {fig_h:.2f} in, "
          f"{out.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
