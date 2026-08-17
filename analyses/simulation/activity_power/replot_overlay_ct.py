#!/usr/bin/env python
"""CRE-vs-minP activity power against fold change, one curve per cell type.

One SVG per (dataset, condition). Three of them are manuscript Fig 3A, placed
side by side in one row at 0.32\\textwidth each:

    shendure/output/overlay_ct/shendure_power_reporter_overlay_ct.svg
    cohen/output/overlay_ct/cohen_power_reporter_overlay_ct.svg
    seelig/output/overlay_ct/seelig_power_deflated_overlay_ct.svg

Those three are drawn as small multiples of a single panel: identical figure
size, identical margins in inches, one shared pair of axes. Axis furniture is
therefore labelled once on the outer edges -- y ticks and the y label on the
left-hand panel, the x label under the middle one -- exactly as a shared-axis
row would be if LaTeX could see all three at once. Every other output is
standalone and carries its full furniture.

Sibling figure: main_figs_power_at_fc.py::fig_b_power_curves (Fig 4B). Same
type scale, same Okabe-Ito pair for the reporter arms, same treatment of the
reference cell type and the 80% rule.

    python replot_overlay_ct.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent

# The reporter/deflated pair is MWU. Measuring the reporter's power benefit
# with the t-test confounds it: the t-test over-rejects in exactly the deflated
# condition, so part of the apparent benefit is that miscalibration rather than
# the reporter. The 2026-08-13 run added the MWU arms for this reason -- see the
# ARMS comment in shendure/shendure_power_ttest_all_cell_types.py.
DATASETS = {
    "cohen": [
        ("reporter", "output/cohen_power_df_2026-08-13_mwu.parquet"),
        ("deflated", "output/cohen_power_df_2026-08-13_mwu_deflated.parquet"),
    ],
    "shendure": [
        ("reporter", "shendure/output/power_df_mwu_2026-08-13.parquet"),
        ("deflated", "shendure/output/power_df_mwu_deflated_2026-08-13.parquet"),
    ],
    "seelig": [
        ("deflated", "seelig/output/power_df_mwu_deflated_2026-08-17.parquet"),
    ],
    "takeshi": [
        ("reporter", "takeshi/output/takeshi_power_df_reporter.parquet"),
        ("deflated", "takeshi/output/takeshi_power_df_deflated.parquet"),
    ],
}

# Author-citation form, written here rather than left to the manuscript's
# dataset_names SVG transform: that transform rewrites <text> nodes, so it is a
# silent no-op on any figure whose text is baked to paths, and the panel would
# print lab nicknames next to a panel that prints citations.
DISPLAY_NAME = {
    "shendure": "Lalanne et al.",
    "cohen": "Zhao et al.",
    "seelig": "Yin et al.",
    "takeshi": "Takeshi (unpublished)",
}

# The reference cell type is a real cell type, renamed at fit time.
REFERENCE_NAMES = {
    "cohen": "Rod",
    "shendure": "Pluripotent",
    "seelig": "HepG2",
    "takeshi": "HepG2",
}

N_CELL_TYPES = {"shendure": 10, "cohen": 4, "seelig": 2, "takeshi": 3}

# Cell-type identifiers as they come out of the fits, and what a reader sees.
# Shared with main_figs_power_at_fc.py so a cell type is spelled the same way
# in Fig 3 and Fig 4.
CT_DISPLAY = {
    "EpiblastPrimitiveStreak": "Epiblast/prim. streak",
    "ExEndodermParietal": "ExE endoderm parietal",
    "ExEndodermVisceral": "ExE endoderm visceral",
    "NeuroectodermBrain": "Neuroectoderm brain",
    "NeuroectodermRostral": "Neuroectoderm rostral",
    "SurfaceEctoderm": "Surface ectoderm",
    "Mueller Glia": "Mueller glia",
    "SKNSH": "SK-N-SH",
}

# What the panel says about the test's access to a transfection reporter. Yin
# et al. ran none, so its "deflated" arm is not a withholding experiment.
CONDITION_LABEL = {
    ("seelig", "deflated"): "no transfection reporter",
    ("*", "reporter"): "reporter used",
    ("*", "deflated"): "reporter withheld",
}

# Okabe-Ito, the palette scripts/figure_styling.py remaps everything else to,
# and the same two hues Fig 4B and activity_calibration/fpr_dumbbell.py already
# use for these two arms. The pair passes every categorical check outright
# (worst CVD dE 21.9, normal-vision 31.2, both above 3:1 on white).
CONDITION_COLOR = {"reporter": "#0072b2", "deflated": "#d55e00"}
INK, MUTED, HAIR = "#1a1a1a", "#666666", "#c9c9c9"

POWER_TARGET = 0.8

# Fold change is a ratio, so the axis is log2: a 2x activation and a 2x
# repression are the same size of effect and sit the same distance from 1x.
# The window is shared by all three manuscript panels -- a per-panel range
# would draw Zhao et al.'s 1.2x threshold at the same width as Yin et al.'s
# 2.4x one, which is the opposite of what the figure is claiming. It runs to
# 4x because that is where Yin et al.'s curves finally clear 80%; Zhao et al.'s
# simulation stops at 1.6x and the rest of its panel is marked as unsimulated
# rather than left as ambiguous white space.
LOG2_WINDOW = (-1.0, 2.0)
# Bin *width* in log2, not bin count, so a bin is the same multiplicative step
# in every dataset even though their raw fc ranges differ by orders of
# magnitude.
LOG2_BIN_WIDTH = 0.1
MIN_BIN_N = 200  # a bin thinner than this is too noisy to draw

# Sized for a 0.32\textwidth inclusion in a 6.90 in text block: 2.21 in on the
# page, and the manuscript sync applies no font scaling, so a point size here
# is the point size a reader gets.
FIG_W, FIG_H = 2.21, 1.85
LEFT_IN, RIGHT_IN, TOP_IN, BOTTOM_IN = 0.42, 0.09, 0.36, 0.34

# The three panels of manuscript Fig 3A, in the order the tabular places them.
ROW_ROLE = {
    ("shendure", "reporter"): "left",
    ("cohen", "reporter"): "middle",
    ("seelig", "deflated"): "right",
}

# Direct labels, hand-placed. A V-shaped curve family leaves only the inside of
# the V free, so there is no rule that puts a label in clear space on its own.
# Each entry is (cell type key, text x in fold change, text y in power,
# horizontal alignment); the marker end of the leader is the curve's ascending
# 80% crossing, computed from the data.
ANNOTATIONS = {
    ("shendure", "reporter"): [
        ("reference", 0.98, 0.93, "center"),
        ("Cardiomyocytes", 1.62, 0.28, "left"),
    ],
    ("cohen", "reporter"): [
        ("reference", 0.80, 0.50, "right"),
        ("Bipolar", 1.58, 0.32, "right"),
    ],
    ("seelig", "deflated"): [
        ("reference", 1.35, 0.93, "center"),
        ("K562", 3.15, 0.38, "center"),
    ],
}

plt.rcParams.update({
    "svg.fonttype": "none",       # keep text editable downstream
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.size": 7,
    "axes.labelsize": 7.5,
    "axes.titlesize": 7.5,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5,
    "axes.unicode_minus": False,  # the repo is plain ASCII
})


# ----- data -----

def display_ct(ct, dataset):
    """Reader-facing name for a cell type in a given dataset."""
    if ct == "reference":
        return REFERENCE_NAMES[dataset]
    return CT_DISPLAY.get(ct, ct)


def condition_label(dataset, condition):
    return CONDITION_LABEL.get((dataset, condition),
                               CONDITION_LABEL[("*", condition)])


def load_power(path, dataset):
    """Fold change, cell type and rejection flag for one simulation arm."""
    df = pd.read_parquet(path)
    if "fc" not in df.columns and "log2_fc" in df.columns:
        df = df.assign(fc=2.0 ** df["log2_fc"].astype(float))
    missing = {"cell_type", "reject_null", "fc"} - set(df.columns)
    assert not missing, f"{path.name}: missing columns {sorted(missing)}"
    assert (df["fc"] > 0).all(), \
        f"{path.name}: fold change must be positive, min {df['fc'].min()}"
    assert df["cell_type"].nunique() == N_CELL_TYPES[dataset], \
        (f"{dataset}: expected {N_CELL_TYPES[dataset]} cell types, got "
         f"{df['cell_type'].nunique()}: {sorted(df['cell_type'].unique())}")
    return df


def pow_curve_data(df, lo=LOG2_WINDOW[0], hi=LOG2_WINDOW[1],
                   width=LOG2_BIN_WIDTH):
    """Power per cell type in fixed-width log2 fold-change bins.

    Bins whose centre falls outside the arm's own simulated range are dropped
    rather than drawn from a handful of edge draws, so a curve ends where the
    simulation ended.
    """
    edges = np.arange(lo, hi + width / 2, width)
    assert len(edges) - 1 >= 10, \
        f"only {len(edges) - 1} bins over log2 [{lo}, {hi}] at width {width}"

    log2_fc = np.log2(df["fc"].to_numpy())
    sub = df.assign(log2_fc=log2_fc)
    sub = sub[(sub["log2_fc"] >= edges[0]) & (sub["log2_fc"] <= edges[-1])]
    assert len(sub), \
        (f"log2 window [{lo}, {hi}] is empty; data spans "
         f"[{log2_fc.min():.3f}, {log2_fc.max():.3f}]")

    binned = pd.cut(sub["log2_fc"], bins=edges, include_lowest=True)
    out = (
        sub.groupby(["cell_type", binned], observed=True)["reject_null"]
        .agg(["mean", "size"])
        .rename(columns={"mean": "power", "size": "n"})
        .reset_index()
    )
    out["log2_fc"] = out["log2_fc"].apply(lambda iv: iv.mid).astype(float)
    out["fc"] = 2.0 ** out["log2_fc"]

    # Keep only bins the arm actually simulated, then require those to be well
    # populated -- an underfilled bin is a noise spike on a smooth curve.
    covered = out[(out["log2_fc"] >= log2_fc.min())
                  & (out["log2_fc"] <= log2_fc.max())]
    assert len(covered), "every bin centre falls outside the simulated range"
    assert covered["n"].min() >= MIN_BIN_N, \
        (f"thinnest kept bin has {covered['n'].min()} draws (floor "
         f"{MIN_BIN_N}) at fold change "
         f"{covered.loc[covered['n'].idxmin(), 'fc']:.2f}")
    assert covered["power"].between(0.0, 1.0).all(), \
        (f"power outside [0, 1]: {covered['power'].min()} to "
         f"{covered['power'].max()}")

    n_cts = covered["cell_type"].nunique()
    n_bins = covered.groupby("cell_type", observed=True).size()
    assert n_bins.nunique() == 1, \
        f"cell types disagree on bin count: {n_bins.to_dict()}"
    assert len(covered) == n_cts * int(n_bins.iloc[0]), \
        (f"expected {n_cts} cell types x {int(n_bins.iloc[0])} bins = "
         f"{n_cts * int(n_bins.iloc[0])} curve points, got {len(covered)}")
    return covered.sort_values(["cell_type", "log2_fc"]).reset_index(drop=True)


def crossing_fc(curve, level=POWER_TARGET):
    """Fold change where an ascending curve first reaches `level`.

    Interpolated in log2 between the two straddling bins. None if the curve
    never gets there within the window.
    """
    up = curve[curve["log2_fc"] > 0]
    x, y = up["log2_fc"].to_numpy(), up["power"].to_numpy()
    for i in range(1, len(y)):
        if y[i - 1] < level <= y[i]:
            t = (level - y[i - 1]) / (y[i] - y[i - 1])
            return float(2.0 ** (x[i - 1] + t * (x[i] - x[i - 1])))
    return None


# ----- drawing -----

def style_axes(ax):
    """Hairline y grid, no box, ticks without tick marks."""
    ax.grid(axis="y", color=HAIR, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(HAIR)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(length=0, colors=INK, pad=2)


def draw_unsimulated(ax, curve):
    """Mark the part of the shared window this arm never simulated.

    Returns the fold change the arm reaches. Zhao et al.'s simulation stops at
    1.6x, well short of the 4x the other two datasets need. Saying so in the
    panel is the price of a shared axis; the alternative is a curve that
    appears to stop for no reason.
    """
    span = LOG2_WINDOW[1] - LOG2_WINDOW[0]
    edge = min(curve["log2_fc"].max() + LOG2_BIN_WIDTH / 2, LOG2_WINDOW[1])
    if LOG2_WINDOW[1] - edge < 0.08 * span:
        return 2.0 ** LOG2_WINDOW[1]
    ax.axvspan(2.0 ** edge, 2.0 ** LOG2_WINDOW[1],
               color=HAIR, alpha=0.22, lw=0, zorder=0)
    ax.text(2.0 ** ((edge + LOG2_WINDOW[1]) / 2), 0.5, "not\nsimulated",
            ha="center", va="center", fontsize=6, color=MUTED,
            linespacing=1.3, zorder=1)
    return 2.0 ** edge


def draw_curves(ax, curve, dataset, color):
    """One line per cell type; the reference cell type carries the panel."""
    for ct, sub in curve.groupby("cell_type", observed=True):
        is_ref = ct == "reference"
        # The reference cell type is the one every other cell type is measured
        # against, and the one Fig 3B's heatmap shows, so it is drawn solid and
        # on top while the rest of the dataset reads as a band around it.
        ax.plot(sub["fc"], sub["power"], color=color,
                lw=1.6 if is_ref else 0.8,
                alpha=1.0 if is_ref else 0.45,
                solid_capstyle="round",
                zorder=4 if is_ref else 3)


def draw_annotations(ax, curve, dataset, condition, color):
    """Name a couple of cell types by pointing at their 80% crossing."""
    for ct, text_fc, text_power, ha in ANNOTATIONS.get((dataset, condition), []):
        sub = curve[curve["cell_type"] == ct]
        assert len(sub), f"{dataset}: no curve for annotated cell type {ct!r}"
        fc = crossing_fc(sub)
        assert fc is not None, \
            f"{dataset}: {ct} never reaches {POWER_TARGET:.0%} in the window"
        ax.annotate(
            display_ct(ct, dataset),
            xy=(fc, POWER_TARGET), xytext=(text_fc, text_power),
            fontsize=6, color=INK, ha=ha, va="center",
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.5,
                            shrinkA=1.5, shrinkB=2.5),
            annotation_clip=False, zorder=6,
        )
        ax.plot([fc], [POWER_TARGET], marker="o", markersize=3.2, color=color,
                markeredgecolor="white", markeredgewidth=0.8, zorder=6)


def draw_title(fig, dataset, condition, color, role=None):
    """Dataset name, plus a condition key on the standalone panels only.

    Indented to the left edge of the axes rather than the figure: the
    manuscript overlays a bold panel letter on the top-left corner of the row
    via \\panel{}, and a flush-left title lands under it.

    The manuscript row (``role`` set) carries no condition key. Each of those
    three panels shows its dataset's only empirical regime -- Lalanne et al.
    and Zhao et al. ran a transfection reporter, Yin et al. did not -- so the
    label restates the dataset rather than distinguishing anything. On the
    standalone panels the condition is a real contrast, a withheld reporter
    against a used one, and the key stays.
    """
    x = LEFT_IN / FIG_W
    fig.text(x, 1 - 0.11 / FIG_H, DISPLAY_NAME[dataset],
             ha="left", va="center", fontsize=7.5, color=INK)
    if role is not None:
        return
    fig.add_artist(Line2D(
        [x, x + 0.085 / FIG_W], [1 - 0.245 / FIG_H] * 2,
        color=color, lw=1.6, solid_capstyle="round",
    ))
    fig.text(x + 0.115 / FIG_W, 1 - 0.245 / FIG_H,
             condition_label(dataset, condition),
             ha="left", va="center", fontsize=6.5, color=MUTED)


def plot_overlay(curve, dataset, condition, out_path):
    role = ROW_ROLE.get((dataset, condition))
    color = CONDITION_COLOR[condition]

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    right_edge = draw_unsimulated(ax, curve)

    # Fold change 1 is the null: no difference between the CRE and minP, so
    # power falls to alpha there and the curve family is V-shaped by
    # construction. The hairline says the dip is the axis, not the data.
    ax.axvline(1.0, color=HAIR, lw=0.8, zorder=1)
    ax.axhline(POWER_TARGET, color=MUTED, ls="--", lw=0.7, zorder=2)

    draw_curves(ax, curve, dataset, color)
    draw_annotations(ax, curve, dataset, condition, color)

    ax.set_xscale("log", base=2)
    ax.set_xlim(2.0 ** LOG2_WINDOW[0], 2.0 ** LOG2_WINDOW[1])
    ax.xaxis.set_major_locator(FixedLocator([0.5, 1, 2, 4]))
    ax.xaxis.set_major_formatter(FixedFormatter(["0.5x", "1x", "2x", "4x"]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_ylim(0, 1.02)
    ax.set_yticks([0, 0.5, 1.0], ["0", "0.5", "1"])
    style_axes(ax)

    # Outer-edge labelling for the manuscript row; everything else stands alone.
    if role in (None, "left"):
        ax.set_ylabel("power", color=INK, labelpad=2)
        # 80% is the conventional target and the reason the dashed rule is
        # there. Naming it once beats a legend entry for a line that is not
        # data. It sits under the right-hand end of the simulated range, the
        # one stretch of the rule that every curve has already climbed past.
        ax.annotate(f"{POWER_TARGET:.0%}", xy=(right_edge, POWER_TARGET),
                    xytext=(-2, -1), textcoords="offset points",
                    fontsize=6, color=MUTED, ha="right", va="top")
    else:
        ax.set_yticklabels([])
    if role in (None, "middle"):
        ax.set_xlabel("fold change vs minP", color=INK, labelpad=1)

    draw_title(fig, dataset, condition, color, role)
    fig.subplots_adjust(
        left=LEFT_IN / FIG_W, right=1 - RIGHT_IN / FIG_W,
        top=1 - TOP_IN / FIG_H, bottom=BOTTOM_IN / FIG_H,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    # No bbox_inches="tight". The three manuscript panels sit in one row with
    # their plot areas in line, and cropping each to its own contents would
    # undo that: the left panel carries tick labels the others do not.
    fig.savefig(out_path, format="svg")
    plt.close(fig)


def report(curve, dataset, condition):
    print(f"  {DISPLAY_NAME[dataset]} ({condition}): "
          f"{curve['cell_type'].nunique()} cell types, "
          f"{curve['log2_fc'].nunique()} bins, thinnest {curve['n'].min()} draws")
    for ct, sub in curve.groupby("cell_type", observed=True):
        fc = crossing_fc(sub)
        where = f"{fc:.2f}x" if fc is not None else "never in window"
        print(f"      {display_ct(ct, dataset):<24} {POWER_TARGET:.0%} at {where}")


def main():
    produced, skipped = [], []
    for dataset, items in DATASETS.items():
        for condition, rel in items:
            path = BASE / rel
            if not path.exists():
                skipped.append((str(path), "missing"))
                continue
            curve = pow_curve_data(load_power(path, dataset))
            report(curve, dataset, condition)
            out_path = (BASE / dataset / "output" / "overlay_ct" /
                        f"{dataset}_power_{condition}_overlay_ct.svg")
            plot_overlay(curve, dataset, condition, out_path)
            produced.append(str(out_path))

    print(f"\nPRODUCED ({FIG_W:.2f} x {FIG_H:.2f} in each):")
    for p in produced:
        print(" ", p)
    print(f"Total: {len(produced)}")
    if skipped:
        print("SKIPPED:")
        for p, why in skipped:
            print(" ", p, "->", why)


if __name__ == "__main__":
    main()
