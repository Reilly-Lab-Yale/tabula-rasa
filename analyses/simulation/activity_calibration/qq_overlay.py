#!/usr/bin/env python
"""Null-calibration QQ overlays: every cell type on one axis, one row per dataset.

Companion to shendure/plot_calibration_qq_detail.py, which breaks a single
(dataset, test, condition) cell out into one panel per cell type. This figure
is the cross-comparison instead: how calibration differs between the two tests,
between the reporter and deflated fit conditions, and between the three
datasets. Per-cell-type detail is the other figure's job.

Form. The panels plot `observed - expected` against `expected`, with ideal
calibration the horizontal line at zero; see qq_common.py for why the
undifferenced [0, 1] QQ cannot show departures this small, and for the meaning
of the grey Kolmogorov band.

Layout. One SVG per dataset, four columns wide, drawn at \\textwidth and
included unscaled. The four columns are the same slots in every row -- reporter
used {MWU, t-test}, reporter withheld {MWU, t-test} -- so the rows stack into a
single grid and a curve is comparable to the one beside it and the one above
it. Yin et al. ran no transfection reporter, so its first two slots are marked
not applicable rather than being filled by the two conditions it does have.

Scale. All three manuscript rows share one y-scale, computed over the ten
panels they contain. Per-panel scaling would rescale each family to fill its
own box and destroy exactly the comparison the figure exists for: the t-test's
excursion in the deflated condition is three to seven times the MWU excursion
beside it, and that ratio has to survive to the page.

Cell-type identity is not encoded. Up to ten curves per panel is past any
honest categorical palette, and at 1.4 in a ten-entry key would be most of the
panel. Colour carries the test instead -- the figure's primary contrast, and
the same key fpr_dumbbell.py uses -- and the one curve a reader would ask about,
the worst-calibrated, is named in the corner of each panel.

    python qq_overlay.py                    # the three manuscript rows
    python qq_overlay.py --dataset takeshi  # one dataset, on its own scale
"""
import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from qq_common import (
    KS_CRIT_95, INK, MUTED, HAIR, BAND, RC_BASE, TEST_COLOR, TEST_LABEL,
    CONDITION_LABEL, RDP_FRACTION, check_reduction, expected, reduce_curve,
)

HERE = pathlib.Path(__file__).resolve().parent

# display name, reference-cell-type display name, expected number of cell types
DATASETS = {
    "shendure": ("Lalanne et al.", "Pluripotent", 10),
    "cohen": ("Zhao et al.", "Rod", 4),
    "seelig": ("Yin et al.", "HepG2", 2),
    "takeshi": ("Takeshi", "HepG2", 3),
}
MANUSCRIPT_ROWS = ["shendure", "cohen", "seelig"]

# Code identifiers upstream, reader-facing names here. "reference" is the cell
# type the synthetic CREs are calibrated against and is per-dataset, so it is
# resolved from DATASETS rather than listed here.
CELL_TYPE_LABEL = {
    "EpiblastPrimitiveStreak": "Epiblast / prim. streak",
    "ExEndodermParietal": "ExE endoderm (par.)",
    "ExEndodermVisceral": "ExE endoderm (visc.)",
    "Haematoendothelial": "Haematoendothelial",
    "NeuroectodermBrain": "Neuroectoderm (brain)",
    "NeuroectodermRostral": "Neuroectoderm (rostral)",
    "SurfaceEctoderm": "Surface ectoderm",
    "Mueller Glia": "Mueller glia",
    "SKNSH": "SK-N-SH",
}

# Column order. Grouped by condition first so that the two tests sit side by
# side within a condition: "does the t-test misbehave where MWU does not" is a
# comparison between neighbours, which is the one the figure is making.
COLUMNS = [("reporter", "mwu"), ("reporter", "ttest"),
           ("deflated", "mwu"), ("deflated", "ttest")]
GROUPS = [("reporter", 0, 2), ("deflated", 2, 4)]

FORM = "deviation"

# Inches, converted to figure fractions below, so the drawn width is the width
# the manuscript includes the file at and a point size here is a point size on
# the page. \textwidth is 498.66pt = 6.90in and the sync applies no scaling.
FIG_W = 6.90
LEFT, RIGHT, WGAP = 0.64, 0.34, 0.11
PANEL_W = (FIG_W - LEFT - RIGHT - WGAP * (len(COLUMNS) - 1)) / len(COLUMNS)
PANEL_H = 1.02
# The headerless rows carry a quarter inch of top margin rather than the 0.10
# the axis alone needs: the manuscript overlays a bold panel letter on the
# top-left corner of each row via \panel{}, and the y-axis label is as long as
# the panel is tall, so without the slack the letter lands on top of it.
TOP_HEADERS, TOP_PLAIN = 0.46, 0.24
BOTTOM_XLABEL, BOTTOM_PLAIN = 0.44, 0.26

LINEWIDTH_PT = 1.0
CURVE_ALPHA = 0.55

plt.rcParams.update({
    **RC_BASE,
    "font.size": 7,
    "axes.labelsize": 7.5,
    "axes.titlesize": 7.5,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
})


def ct_label(dataset, ct):
    if ct == "reference":
        return DATASETS[dataset][1]
    return CELL_TYPE_LABEL.get(ct, ct)


def pvals_path(dataset, test, condition):
    tag = "" if test == "ttest" else f"{test}_"
    return HERE / dataset / "output" / f"{dataset}_null_pvals_{tag}{condition}.parquet"


def summary_path(dataset, test):
    tag = "" if test == "ttest" else f"_{test}"
    return HERE / dataset / "output" / f"{dataset}_null_summary{tag}.parquet"


def load_panel(dataset, test, condition):
    """{cell type: sorted finite p-values}, or None where the file is absent.

    Absent is a real state, not an error: Yin et al. ran no transfection
    reporter, so it has no reporter-condition p-values to load.
    """
    path = pvals_path(dataset, test, condition)
    if not path.is_file():
        return None
    df = pd.read_parquet(path)
    n_rows = len(df)
    n_expected = DATASETS[dataset][2]

    out = {}
    for ct, g in df.groupby("cell_type", observed=True):
        v = g["p_value"].to_numpy()
        assert len(v) > 0, f"{dataset}/{test}/{condition}: empty group {ct}"
        finite = np.isfinite(v)
        assert finite.all(), \
            f"{dataset}/{test}/{condition} {ct}: {(~finite).sum()} non-finite of {len(v)}"
        assert v.min() >= 0.0 and v.max() <= 1.0, \
            f"{dataset}/{test}/{condition} {ct}: p-values outside [0,1], " \
            f"range {v.min()}-{v.max()}"
        out[ct] = np.sort(v)

    assert sum(len(v) for v in out.values()) == n_rows, \
        f"{dataset}/{test}/{condition}: groupby lost rows, " \
        f"{sum(len(v) for v in out.values())} of {n_rows}"
    assert len(out) == n_expected, \
        f"{dataset}/{test}/{condition}: {len(out)} cell types, expected {n_expected}"
    return out


def load_summary(dataset, test, condition, pvals):
    """{cell type: (D, p, n)} from the cached summary, checked against the curves.

    The summary is what the manuscript quotes, so the annotation reads from it
    rather than recomputing. It is only trustworthy if it describes the same
    p-values that are being drawn, hence the two checks: the per-cell-type row
    count, and KS D against the largest deviation in the drawn curve. Those
    two differ by at most one rank step, since D takes the sup over both sides
    of each jump and the curve is evaluated at rank midpoints.
    """
    path = summary_path(dataset, test)
    assert path.is_file(), f"missing null summary: {path}"
    s = pd.read_parquet(path)
    s = s[s["condition"] == condition]
    assert len(s) == len(pvals), \
        f"{dataset}/{test}/{condition}: summary has {len(s)} rows, " \
        f"{len(pvals)} cell types on disk"
    assert s["cell_type"].is_unique, \
        f"{dataset}/{test}/{condition}: duplicate cell types in the summary"
    assert ((s["ks_p"] >= 0) & (s["ks_p"] <= 1)).all(), \
        f"{dataset}/{test}/{condition}: KS p-values outside [0,1]"

    out = {r.cell_type: (r.ks_d, r.ks_p, r.n) for r in s.itertuples()}
    assert set(out) == set(pvals), \
        f"{dataset}/{test}/{condition}: summary cell types {sorted(out)} " \
        f"!= data {sorted(pvals)}"
    for ct, (d, _, n) in out.items():
        v = pvals[ct]
        assert n == len(v), \
            f"{dataset}/{test}/{condition} {ct}: summary n={n}, {len(v)} p-values on disk"
        dev = float(np.abs(v - expected(len(v))).max())
        assert abs(d - dev) <= 1.5 / len(v), \
            f"{dataset}/{test}/{condition} {ct}: summary D={d:.6f} disagrees with " \
            f"the drawn curve's largest deviation {dev:.6f} by more than one rank step"
    return out


def nice_step(lim):
    """Y tick step: the largest 1/2/5 x 10^k that puts two ticks on each side."""
    raw = lim / 2.0
    mag = 10.0 ** np.floor(np.log10(raw))
    step = float(max(m * mag for m in (1, 2, 5) if m * mag <= raw * 1.001))
    assert step <= lim, f"tick step {step} leaves no tick inside +/-{lim}"
    return step


def draw_panel(ax, pvals, test, tol_y, ylim, ks):
    """One (condition, test) cell: band, zero line, the family, and its worst."""
    # The Kolmogorov band is a function of n, and n differs between cell types
    # (60.7k to 160.7k in Lalanne et al., under 0.2% in the other two), so a
    # single band cannot be exact for all the curves in a panel. It is drawn at
    # the widest, i.e. from the smallest n present, which is the conservative
    # choice in the direction the figure is arguing: a curve that leaves the
    # band has a KS D above the acceptance threshold for every cell type here,
    # not just for the smallest one.
    n_min = min(len(v) for v in pvals.values())
    crit = KS_CRIT_95 / np.sqrt(n_min)
    ax.axhspan(-crit, crit, color=BAND, lw=0, zorder=0)
    ax.axhline(0.0, color=MUTED, ls="--", lw=0.7, zorder=2)

    for ct, v in pvals.items():
        xr, yr = reduce_curve(v, FORM, tol_y)
        ax.plot(xr, yr, color=TEST_COLOR[test], lw=LINEWIDTH_PT, alpha=CURVE_ALPHA,
                solid_capstyle="round", solid_joinstyle="round", zorder=3)

    assert ylim[1] > max(ks[ct][0] for ct in pvals), \
        f"a curve leaves the panel: max D {max(ks[ct][0] for ct in pvals):.4f} " \
        f"exceeds the y limit {ylim[1]:.4f}"
    return crit


def annotate_worst(ax, dataset, pvals, ks):
    """Name the worst-calibrated curve in the panel, top left.

    The panel's own extreme is the one fact a reader wants that the overlay
    cannot otherwise give, once the per-cell-type key is gone. Top left is free
    in every panel: the deviation curve is a bridge pinned to zero at both ends
    and no positive excursion anywhere in this data reaches half the y-range.
    """
    worst = max(pvals, key=lambda c: ks[c][0])
    d = ks[worst][0]
    ax.text(0.035, 0.97, f"max D = {d:.3f}\n{ct_label(dataset, worst)}",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=6, color=MUTED, linespacing=1.3)


def not_applicable(ax):
    ax.text(0.5, 0.5, "not applicable\nno transfection reporter",
            transform=ax.transAxes, ha="center", va="center",
            fontsize=6, color=MUTED, linespacing=1.5)


def figure(dataset, panels, ks_tables, ylim, show_headers, show_xlabel):
    top = TOP_HEADERS if show_headers else TOP_PLAIN
    bottom = BOTTOM_XLABEL if show_xlabel else BOTTOM_PLAIN
    fig_h = top + bottom + PANEL_H

    # Simplification tolerance in data units: a fraction of the drawn linewidth,
    # given how many data units a point of panel height is worth.
    stroke_y = LINEWIDTH_PT / 72 / PANEL_H * (ylim[1] - ylim[0])
    tol_y = RDP_FRACTION * stroke_y
    for (cond, test), p in panels.items():
        if p is not None:
            check_reduction(p, FORM, tol_y, 0.5 * stroke_y,
                            linewidth_pt=LINEWIDTH_PT,
                            tag=f"{dataset} {test} {cond}: ")

    fig, axes = plt.subplots(1, len(COLUMNS), figsize=(FIG_W, fig_h), squeeze=False)
    fig.subplots_adjust(
        left=LEFT / FIG_W, right=1 - RIGHT / FIG_W,
        bottom=bottom / fig_h, top=1 - top / fig_h,
        wspace=WGAP / PANEL_W,
    )
    axes = axes[0]

    step = nice_step(ylim[1])
    yticks = np.arange(-np.floor(ylim[1] / step) * step, ylim[1] + 1e-12, step)
    band_ax, band_crit = None, None

    for k, (cond, test) in enumerate(COLUMNS):
        ax = axes[k]
        p = panels[(cond, test)]
        if p is None:
            not_applicable(ax)
        else:
            crit = draw_panel(ax, p, test, tol_y, ylim, ks_tables[(cond, test)])
            annotate_worst(ax, dataset, p, ks_tables[(cond, test)])
            if band_ax is None:
                band_ax, band_crit = ax, crit

        ax.set_xlim(0, 1)
        ax.set_ylim(*ylim)
        ax.set_xticks([0, 0.5, 1])
        ax.set_xticklabels(["0", "0.5", "1"])
        ax.set_yticks(yticks)
        if k > 0:
            ax.set_yticklabels([])
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(HAIR)
        ax.tick_params(length=2.5, width=0.6, color=HAIR, colors=INK, pad=2)

    axes[0].set_ylabel("observed - expected", color=INK, labelpad=2)
    if show_xlabel:
        # Centred on the block of panels, not on the canvas: the left margin
        # carries the y-axis furniture and the right only a rotated row label,
        # so the two are not equal and fig.supxlabel would sit off-centre.
        fig.text((LEFT + 0.5 * (FIG_W - LEFT - RIGHT)) / FIG_W, 0.03 / fig_h,
                 "expected p-value (uniform quantile)",
                 ha="center", va="bottom", fontsize=7.5, color=INK)

    # Headers, and which way is bad, on the first row only: the rows stack into
    # one grid on a single y-scale, so repeating either beside three identical
    # axes is noise.
    if show_headers:
        draw_headers(fig, fig_h, top)
        fig.text(0.008, bottom / fig_h + 0.5 * PANEL_H / fig_h,
                 "below zero: anti-conservative", rotation=90,
                 ha="left", va="center", fontsize=6, color=MUTED)

    # Dataset identity on the right, where nothing else lives; the rows carry
    # no titles, so this is what tells them apart once they are stacked.
    fig.text(1 - 0.10 / FIG_W, bottom / fig_h + 0.5 * PANEL_H / fig_h,
             DATASETS[dataset][0], rotation=-90,
             ha="center", va="center", fontsize=7.5, color=INK)

    # One direct label on the band rather than a legend entry: it is furniture,
    # it appears in every panel, and the caption defines it. Drawn on the first
    # row only, since the rows stack into one grid.
    if show_headers and band_ax is not None:
        # Anchored on the band's upper edge and led to the top right, the one
        # corner no curve reaches: every family here runs anti-conservative or
        # barely positive, so the space above the band on the right is empty.
        band_ax.annotate(
            "95% KS band", xy=(0.91, band_crit),
            xytext=(0.44, 0.64), textcoords="axes fraction",
            fontsize=6, color=MUTED, ha="left", va="bottom",
            arrowprops=dict(arrowstyle="-", lw=0.5, color=MUTED,
                            shrinkA=1.5, shrinkB=1.0),
        )
    return fig, fig_h


def draw_headers(fig, fig_h, top):
    """Two-level column headers, in inches from the top of the panel row.

    The test name sits over a rule in that test's colour, so the key a reader
    needs is where the colour is used rather than in a legend box somewhere
    else, and the test is named in text as well as coloured.
    """
    panel_top = 1 - top / fig_h

    def cx(j):
        return (LEFT + j * (PANEL_W + WGAP) + PANEL_W / 2) / FIG_W

    def edge(j, side):
        x = LEFT + j * (PANEL_W + WGAP)
        return (x if side == "l" else x + PANEL_W) / FIG_W

    for j, (_, test) in enumerate(COLUMNS):
        y_rule = panel_top + 0.05 / fig_h
        half = 0.5 * 0.34 / FIG_W
        fig.add_artist(Line2D([cx(j) - half, cx(j) + half], [y_rule, y_rule],
                              color=TEST_COLOR[test], lw=1.8,
                              solid_capstyle="round"))
        fig.text(cx(j), panel_top + 0.10 / fig_h, TEST_LABEL[test],
                 ha="center", va="bottom", fontsize=7.5, color=INK)

    for cond, j0, j1 in GROUPS:
        y_rule = panel_top + 0.25 / fig_h
        fig.add_artist(Line2D([edge(j0, "l"), edge(j1 - 1, "r")], [y_rule, y_rule],
                              color=HAIR, lw=0.7))
        fig.text(0.5 * (edge(j0, "l") + edge(j1 - 1, "r")),
                 panel_top + 0.29 / fig_h, CONDITION_LABEL[cond],
                 ha="center", va="bottom", fontsize=7.5, color=MUTED)


def build(datasets, ylim=None):
    panels, ks = {}, {}
    for ds in datasets:
        for cond, test in COLUMNS:
            p = load_panel(ds, test, cond)
            panels[(ds, cond, test)] = p
            ks[(ds, cond, test)] = None if p is None else load_summary(ds, test, cond, p)

    drawn = [k for k, v in panels.items() if v is not None]
    assert drawn, f"no null-pvals parquets found for {datasets}"

    if ylim is None:
        span = max(float(np.abs(v - expected(len(v))).max())
                   for k in drawn for v in panels[k].values())
        lim = 1.25 * span
        ylim = (-lim, lim)
        print(f"shared y range +/-{lim:.4f} over {len(drawn)} panels "
              f"(largest deviation {span:.4f})")

    for i, ds in enumerate(datasets):
        sub_p = {(c, t): panels[(ds, c, t)] for c, t in COLUMNS}
        sub_k = {(c, t): ks[(ds, c, t)] for c, t in COLUMNS}
        fig, fig_h = figure(ds, sub_p, sub_k, ylim,
                            show_headers=(i == 0),
                            show_xlabel=(i == len(datasets) - 1))
        out = HERE / ds / "output" / f"{ds}_calibration_qq_overlay.svg"
        # No bbox_inches="tight": the margins above are in inches so that the
        # rows stack with their panel columns in line, and cropping each row to
        # its own contents would undo that.
        fig.savefig(out, format="svg")
        plt.close(fig)
        print(f"wrote {out.relative_to(HERE)}  ({FIG_W:.2f} x {fig_h:.2f} in, "
              f"{out.stat().st_size / 1024:.0f} KB)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", choices=list(DATASETS),
                    help="draw one dataset on its own y-scale; "
                         "default is the three manuscript rows on a shared one")
    args = ap.parse_args()

    if args.dataset:
        build([args.dataset])
    else:
        build(MANUSCRIPT_ROWS)


if __name__ == "__main__":
    main()
