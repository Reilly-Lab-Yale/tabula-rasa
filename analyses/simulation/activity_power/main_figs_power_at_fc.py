#!/usr/bin/env python3
"""Main-text power figures for the transfection-reporter panel (manuscript Fig 4).

Fig 4A -- fig_a_reporter_dumbbell -> output/fig_a_reporter_dumbbell.svg
    Power at a single fold change, per cell type, with the transfection
    reporter and without it. The two dots are the paired arms of the same
    simulation and the segment between them is the reporter's contribution.

Fig 4B -- fig_b_power_curves -> output/fig_b_slope.svg
    The same two arms as continuous power-vs-fold-change curves, one curve per
    cell type, over a fold-change window both datasets cover. Panel A reads one
    vertical slice of this; panel B says whether that slice is representative.

Only Zhao et al. (code name cohen) and Lalanne et al. (shendure) appear in
either manuscript panel. Yin et al. (seelig) ran no transfection reporter, so
it has no "+reporter" arm and never can have one -- a paired contrast cannot
include it. It stays in DATASETS because the two supplementary/exploratory
figures below (fig_b_dataset_bars, fig_b_telescope) do show all three.

Also written, not manuscript figures:
    output/fig_b_dataset_bars.svg     best-available power per dataset
    output/fig_b_telescope.svg        nested bars over a per-dataset FC triplet
    output/fig_b_slope_discrete.svg   the three-point slope form of Fig 4B,
                                      kept for comparison against the curves

Usage:
    python main_figs_power_at_fc.py
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

TARGET_FC = 1.5
FC_TOL = 0.05  # +/- 5% window around target

# Continuous power curves (Fig 4B) are drawn over a window both datasets cover.
# Zhao et al.'s simulation stops at fc 1.60, Lalanne et al.'s runs to 4.0. A
# shared axis wide enough for Lalanne would leave most of Zhao's panel empty,
# and a per-panel axis would make the two datasets' slopes look comparable when
# they are drawn at different fold-change-per-inch. Clipping both to the common
# informative window keeps every panel full and every slope comparable.
FC_WINDOW = (1.0, 1.6)
# Bin *width*, not bin count, so a bin means the same fold-change interval in
# both datasets even though their raw fc ranges differ by a factor of four.
FC_BIN_WIDTH = 0.03
MIN_BIN_N = 30  # a bin thinner than this is too noisy to draw

# Per-dataset FC triplets for the discrete variants. Each spans the dataset's
# informative regime. Cohen caps at 1.6; seelig needs higher FCs to show any
# power; shendure is in the middle.
FC_TRIPLET = {
    "cohen":    [1.1, 1.2, 1.3],
    "shendure": [1.2, 1.5, 2.0],
    "seelig":   [1.5, 2.0, 3.0],
}

SCRIPT_DIR = Path(__file__).resolve().parent
OUT_DIR = SCRIPT_DIR / "output"
OUT_DIR.mkdir(exist_ok=True)

# The reporter/deflated pair is MWU. Measuring the reporter's power benefit
# with the t-test confounds it: the t-test over-rejects in exactly the deflated
# condition, so part of the apparent benefit is that miscalibration rather than
# the reporter. The 2026-08-13 run added the MWU arms for this reason -- see the
# ARMS comment in shendure/shendure_power_ttest_all_cell_types.py.
# Yin et al. is still t-test: no MWU arm has been run on its power sims.
DATASETS = {
    "cohen": {
        "reporter": SCRIPT_DIR / "output/cohen_power_df_2026-08-13_mwu.parquet",
        "deflated": SCRIPT_DIR / "output/cohen_power_df_2026-08-13_mwu_deflated.parquet",
        "ref_display": "Rod",
        "title": "Zhao et al. (episomal, CRE reporter)",
        "n_cell_types": 4,
    },
    "shendure": {
        "reporter": SCRIPT_DIR / "shendure/output/power_df_mwu_2026-08-13.parquet",
        "deflated": SCRIPT_DIR / "shendure/output/power_df_mwu_deflated_2026-08-13.parquet",
        "ref_display": "Pluripotent",
        "title": "Lalanne et al. (piggyBac, oBC reporter)",
        "n_cell_types": 10,
    },
    "seelig": {
        "deflated": SCRIPT_DIR / "seelig/output/seelig_power_df_deflated.parquet",
        "ref_display": "HepG2",
        "title": "Yin et al. (episomal, no reporter)",
        "n_cell_types": 2,
    },
}

# The two datasets that ran both arms, in the order they appear in the
# manuscript panels: the four-row dataset on top of the ten-row one.
PAIRED_DSETS = ["cohen", "shendure"]

# Cell-type identifiers as they come out of the fits, and what a reader sees.
# Only the camel-case run-together names need rewriting; anything absent here
# is passed through, and "reference" is handled separately per dataset.
CT_DISPLAY = {
    "EpiblastPrimitiveStreak": "Epiblast/prim. streak",
    "ExEndodermParietal": "ExE endoderm parietal",
    "ExEndodermVisceral": "ExE endoderm visceral",
    "NeuroectodermBrain": "Neuroectoderm brain",
    "NeuroectodermRostral": "Neuroectoderm rostral",
    "SurfaceEctoderm": "Surface ectoderm",
    "Haematoendothelial": "Haematoendothelial",
    "Mueller Glia": "Mueller glia",
}

# Okabe-Ito, the palette scripts/figure_styling.py remaps everything else to,
# and the pair already used by activity_calibration/fpr_dumbbell.py so the two
# figures sit next to each other without changing colour conventions. The pair
# passes every categorical check outright (worst CVD dE 21.9, normal-vision
# 31.2, both above 3:1 on white), so nothing is leaning on shape or position to
# tell the two arms apart.
COLOR_REP = "#0072b2"   # with transfection reporter
COLOR_DEFL = "#d55e00"  # without
INK, MUTED, HAIR = "#1a1a1a", "#666666", "#c9c9c9"

LABEL_REP = "with reporter"
LABEL_DEFL = "without reporter"

# Sized for a 0.49\textwidth inclusion: 3.38 in on the page, no downstream font
# scaling, so a point size here is the point size a reader gets.
FIG_W = 3.38

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


def display_ct(ct, dataset):
    """Reader-facing name for a cell type in a given dataset."""
    if ct == "reference":
        return DATASETS[dataset]["ref_display"]
    return CT_DISPLAY.get(ct, ct)


def style_axes(ax, grid_axis):
    """Hairline grid on one axis, no box, ticks without tick marks."""
    ax.grid(axis=grid_axis, color=HAIR, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(HAIR)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(length=0, colors=INK, pad=2)


# ----- power at a single fold change -----

def power_at_fc(df, target_fc=TARGET_FC, tol=FC_TOL):
    """Mean reject_null over a window of fc around target, per cell type."""
    lo = target_fc * (1 - tol)
    hi = target_fc * (1 + tol)
    window = df[(df["fc"] >= lo) & (df["fc"] <= hi)]
    assert len(window), \
        f"fc window [{lo:.3f}, {hi:.3f}] is empty; data spans " \
        f"[{df['fc'].min():.3f}, {df['fc'].max():.3f}]"
    out = (
        window.groupby("cell_type", observed=True)["reject_null"]
        .agg(["mean", "size"])
        .rename(columns={"mean": "power", "size": "n"})
        .reset_index()
    )
    # Every cell type in the parquet must survive the window, or a dumbbell
    # silently loses a row rather than failing.
    assert len(out) == df["cell_type"].nunique(), \
        f"fc window [{lo:.3f}, {hi:.3f}] dropped cell types: " \
        f"{df['cell_type'].nunique()} in, {len(out)} out"
    assert out["power"].between(0.0, 1.0).all(), \
        f"power outside [0, 1]: {out['power'].min()} to {out['power'].max()}"
    return out


def rename_reference(df, display):
    df = df.copy()
    df.loc[df["cell_type"] == "reference", "cell_type"] = display
    return df


def load_powers():
    out = {}
    for name, d in DATASETS.items():
        out[name] = {}
        if "reporter" in d:
            rep = pd.read_parquet(d["reporter"])
            out[name]["reporter"] = rename_reference(
                power_at_fc(rep), d["ref_display"]
            )
        defl = pd.read_parquet(d["deflated"])
        out[name]["deflated"] = rename_reference(
            power_at_fc(defl), d["ref_display"]
        )
        for mode, tab in out[name].items():
            assert len(tab) == d["n_cell_types"], \
                f"{name} ({mode}): expected {d['n_cell_types']} cell types, " \
                f"got {len(tab)}: {sorted(tab['cell_type'])}"
    return out


# ----- Fig 4A: paired dumbbell at one fold change -----

def fig_a_reporter_dumbbell(powers):
    """Fig 4A: power with and without the reporter, per cell type."""
    n_rows = [DATASETS[d]["n_cell_types"] for d in PAIRED_DSETS]
    total_rows = sum(n_rows)

    # Geometry in inches rather than tight_layout, so the row pitch is the same
    # in both panels and the figure height follows the number of rows.
    pitch, top_in, gap_in, bottom_in = 0.165, 0.40, 0.30, 0.36
    plot_in = pitch * total_rows
    height = top_in + plot_in + gap_in + bottom_in

    fig, axes = plt.subplots(
        len(PAIRED_DSETS), 1, figsize=(FIG_W, height),
        gridspec_kw={"height_ratios": n_rows,
                     "hspace": gap_in / (plot_in / len(PAIRED_DSETS))},
        sharex=True,
    )

    for ax, name in zip(axes, PAIRED_DSETS):
        rep = powers[name]["reporter"].set_index("cell_type")["power"]
        defl = powers[name]["deflated"].set_index("cell_type")["power"]
        assert set(rep.index) == set(defl.index), \
            f"{name}: arms disagree on cell types, " \
            f"{sorted(set(rep.index) ^ set(defl.index))} unmatched"

        # Strongest achievable power at the top.
        cts = rep.sort_values(ascending=True).index.tolist()
        y = np.arange(len(cts))
        rep_vals = rep.loc[cts].values
        defl_vals = defl.loc[cts].values

        # The connector is the reporter's contribution. It is drawn under both
        # dots and lighter than either, so it reads as the distance between
        # them rather than as a third mark.
        ax.hlines(y, defl_vals, rep_vals, color=HAIR, lw=1.8, zorder=2)
        for vals, color, label in (
            (defl_vals, COLOR_DEFL, LABEL_DEFL),
            (rep_vals, COLOR_REP, LABEL_REP),
        ):
            ax.scatter(vals, y, s=26, color=color, label=label, zorder=3,
                       edgecolors="white", linewidths=1.0, clip_on=False)

        ax.set_yticks(y)
        ax.set_yticklabels([display_ct(c, name) for c in cts])
        ax.set_ylim(-0.5, len(cts) - 0.5)
        ax.set_xlim(0, 1.0)
        ax.set_xticks(np.arange(0, 1.01, 0.25), ["0", "0.25", "0.5", "0.75", "1"])
        ax.axvline(0.8, color=MUTED, ls="--", lw=0.7, zorder=1)
        ax.set_title(DATASETS[name]["title"], loc="left", color=INK, pad=4)
        style_axes(ax, grid_axis="x")

    # 80% power is the conventional target and the reason the dashed rule is
    # there; naming it once beats a legend entry for a line that is not data.
    # It goes on the bottom row of the lower panel, the only place right of the
    # rule that no dumbbell reaches.
    axes[-1].annotate("80%", xy=(0.8, 0.0), xytext=(3, 0),
                      textcoords="offset points",
                      fontsize=6, color=MUTED, ha="left", va="center")

    axes[-1].set_xlabel(f"power at {TARGET_FC:g}-fold change", color=INK)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.995, 0.999),
               ncol=2, frameon=False, handletextpad=0.15, columnspacing=1.0,
               borderpad=0.0, borderaxespad=0.0)
    fig.subplots_adjust(
        left=1.12 / FIG_W, right=1 - 0.04 / FIG_W,
        top=1 - top_in / height, bottom=bottom_in / height,
    )
    out = OUT_DIR / "fig_a_reporter_dumbbell.svg"
    fig.savefig(out, format="svg")
    plt.close(fig)
    print(f"Saved: {out}  ({FIG_W:.2f} x {height:.2f} in)")


# ----- Fig 4B: continuous power curves -----

def power_curve(df, lo=FC_WINDOW[0], hi=FC_WINDOW[1], width=FC_BIN_WIDTH):
    """Power per cell type in fixed-width fc bins over [lo, hi]."""
    edges = np.arange(lo, hi + width / 2, width)
    n_bins = len(edges) - 1
    assert n_bins >= 10, f"only {n_bins} bins over [{lo}, {hi}] at width {width}"

    sub = df[(df["fc"] >= edges[0]) & (df["fc"] <= edges[-1])]
    assert len(sub), \
        f"fc window [{lo}, {hi}] is empty; data spans " \
        f"[{df['fc'].min():.3f}, {df['fc'].max():.3f}]"

    binned = pd.cut(sub["fc"], bins=edges, include_lowest=True)
    out = (
        sub.groupby(["cell_type", binned], observed=True)["reject_null"]
        .agg(["mean", "size"])
        .rename(columns={"mean": "power", "size": "n"})
        .reset_index()
    )
    out["fc"] = out["fc"].apply(lambda iv: iv.mid).astype(float)

    n_cts = sub["cell_type"].nunique()
    assert len(out) == n_cts * n_bins, \
        f"expected {n_cts} cell types x {n_bins} bins = {n_cts * n_bins} " \
        f"curve points, got {len(out)} -- some bin is empty"
    assert out["n"].min() >= MIN_BIN_N, \
        f"thinnest bin has {out['n'].min()} draws (floor {MIN_BIN_N})"
    assert out["power"].between(0.0, 1.0).all(), \
        f"power outside [0, 1]: {out['power'].min()} to {out['power'].max()}"
    return out


def load_power_curves():
    """{dataset: {arm: curve_df}} for the two datasets that ran both arms."""
    out = {}
    for name in PAIRED_DSETS:
        d = DATASETS[name]
        out[name] = {}
        for arm in ("reporter", "deflated"):
            curve = power_curve(pd.read_parquet(d[arm]))
            assert curve["cell_type"].nunique() == d["n_cell_types"], \
                f"{name} ({arm}): expected {d['n_cell_types']} cell types, " \
                f"got {curve['cell_type'].nunique()}"
            out[name][arm] = curve
        print(f"  {name}: {out[name]['reporter']['fc'].nunique()} bins of "
              f"width {FC_BIN_WIDTH:g}, thinnest "
              f"{min(out[name][a]['n'].min() for a in out[name])} draws")
    return out


def fig_b_power_curves(curves):
    """Fig 4B: power vs fold change, both arms, one curve per cell type."""
    # Enough right margin to carry the direct labels; no legend box, since two
    # series that are labelled where they separate do not need one.
    left_in, right_in = 0.40, 0.70
    top_in, gap_in, bottom_in = 0.26, 0.30, 0.38
    panel_in = 1.20
    height = top_in + panel_in * len(PAIRED_DSETS) + gap_in + bottom_in

    fig, axes = plt.subplots(
        len(PAIRED_DSETS), 1, figsize=(FIG_W, height),
        gridspec_kw={"hspace": gap_in / panel_in},
        sharex=True, sharey=True,
    )

    for ax, name in zip(axes, PAIRED_DSETS):
        ref_display = DATASETS[name]["ref_display"]
        for arm, color in (("deflated", COLOR_DEFL), ("reporter", COLOR_REP)):
            df = curves[name][arm]
            for ct, sub in df.groupby("cell_type", observed=True):
                sub = sub.sort_values("fc")
                is_ref = ct == "reference"
                # The reference cell type is the one the other cell types are
                # measured against, so it is drawn solid and on top while the
                # rest of its arm reads as a band.
                ax.plot(sub["fc"], sub["power"], color=color,
                        lw=1.6 if is_ref else 0.8,
                        alpha=1.0 if is_ref else 0.45,
                        solid_capstyle="round",
                        zorder=4 if is_ref else 3)
                if is_ref and arm == "reporter":
                    ax.annotate(
                        ref_display,
                        xy=(sub["fc"].iloc[-1], sub["power"].iloc[-1]),
                        xytext=(4, 0), textcoords="offset points",
                        fontsize=6, color=INK, va="center", ha="left",
                        annotation_clip=False,
                    )

        ax.axhline(0.8, color=MUTED, ls="--", lw=0.7, zorder=1)
        ax.axvline(TARGET_FC, color=HAIR, lw=0.8, zorder=1)
        ax.set_xlim(*FC_WINDOW)
        ax.set_ylim(0, 1.02)
        ax.set_yticks([0, 0.5, 1.0], ["0", "0.5", "1"])
        ax.set_title(DATASETS[name]["title"], loc="left", color=INK, pad=4)
        style_axes(ax, grid_axis="y")

    # Panel A's fold change, so the two panels of the figure lock together.
    axes[0].annotate(f"{TARGET_FC:g}x", xy=(TARGET_FC, 1.0),
                     xycoords=("data", "axes fraction"),
                     xytext=(2, 1), textcoords="offset points",
                     fontsize=6, color=MUTED, ha="left", va="bottom")
    # Direct labels go where the two arms are furthest apart: Lalanne et al.'s
    # right edge, where the bands are separated by most of the axis. One
    # labelling serves both panels, which share colours and axes.
    lower = axes[-1]
    lower.annotate("80%", xy=(FC_WINDOW[0], 0.8),
                   xytext=(2, 2), textcoords="offset points",
                   fontsize=6, color=MUTED, ha="left", va="bottom")
    for arm, color, label in (("reporter", COLOR_REP, LABEL_REP),
                              ("deflated", COLOR_DEFL, LABEL_DEFL)):
        df = curves["shendure"][arm]
        at_edge = df[df["fc"] == df["fc"].max()]
        # Anchor on the weakest cell type in the arm, keeping clear of the
        # reference label that rides the top of the +reporter band.
        y = at_edge["power"].min()
        # A short key in the arm's colour carries the identity; the text stays
        # ink, which the vermillion is too light to be at this size.
        lower.plot([FC_WINDOW[1] + 0.012, FC_WINDOW[1] + 0.042], [y, y],
                   color=color, lw=1.6, solid_capstyle="round",
                   clip_on=False, zorder=5)
        lower.annotate(label.replace(" ", "\n"), xy=(FC_WINDOW[1] + 0.052, y),
                       xytext=(0, 0), textcoords="offset points",
                       fontsize=6, color=INK, va="center", ha="left",
                       linespacing=1.25, annotation_clip=False)

    lower.set_xlabel("fold change", color=INK)
    for ax in axes:
        ax.set_ylabel("power", color=INK)

    fig.subplots_adjust(
        left=left_in / FIG_W, right=1 - right_in / FIG_W,
        top=1 - top_in / height, bottom=bottom_in / height,
    )
    out = OUT_DIR / "fig_b_slope.svg"
    fig.savefig(out, format="svg")
    plt.close(fig)
    print(f"Saved: {out}  ({FIG_W:.2f} x {height:.2f} in)")


# ----- discrete variants, not manuscript panels -----

def fig_b_dataset_bars(powers):
    """Best-available power at TARGET_FC across all three datasets."""
    best_mode = {"cohen": "reporter", "shendure": "reporter", "seelig": "deflated"}
    mode_label = {"reporter": LABEL_REP, "deflated": LABEL_DEFL}
    mode_color = {"reporter": COLOR_REP, "deflated": COLOR_DEFL}

    dsets = ["cohen", "shendure", "seelig"]
    heights = [max(DATASETS[d]["n_cell_types"], 1) for d in dsets]

    fig, axes = plt.subplots(
        3, 1, figsize=(6, 5.5),
        gridspec_kw={"height_ratios": heights},
        sharex=True,
    )

    for ax, name in zip(axes, dsets):
        mode = best_mode[name]
        df = powers[name][mode].sort_values("power")
        cts = df["cell_type"].tolist()
        vals = df["power"].values
        y = np.arange(len(cts))
        ax.barh(y, vals, color=mode_color[mode], height=0.7)
        ax.set_yticks(y)
        ax.set_yticklabels([display_ct(c, name) for c in cts])
        ax.set_xlim(0, 1.0)
        ax.axvline(0.8, color=MUTED, ls="--", lw=0.7)
        ax.set_title(f"{DATASETS[name]['title'].split(' (')[0]}, "
                     f"{mode_label[mode]}", loc="left")
        style_axes(ax, grid_axis="x")
        ax.invert_yaxis()  # strongest at top

    axes[-1].set_xlabel(f"power at {TARGET_FC:g}-fold change")
    plt.tight_layout()
    out = OUT_DIR / "fig_b_dataset_bars.svg"
    fig.savefig(out, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


def load_multi_fc_powers():
    """Return {dataset: {mode: {fc_val: power_df}}} using FC_TRIPLET."""
    out = {}
    for name, d in DATASETS.items():
        out[name] = {}
        for parquet_key in ("reporter", "deflated"):
            if parquet_key not in d:
                continue
            raw = pd.read_parquet(d[parquet_key])
            out[name][parquet_key] = {
                fc: rename_reference(power_at_fc(raw, target_fc=fc, tol=FC_TOL),
                                     d["ref_display"])
                for fc in FC_TRIPLET[name]
            }
    return out


def _shades(base_color, n=3):
    """Build n shades from light to dark for a base matplotlib color."""
    import matplotlib.colors as mcolors
    rgb = np.array(mcolors.to_rgb(base_color))
    # fractions: larger -> closer to base (darker).
    fracs = np.linspace(0.30, 0.90, n)
    return [tuple(1.0 - f * (1.0 - rgb)) for f in fracs]


def fig_b_telescope(mpows):
    """Telescope (nested) bars: 3 overlaid bars per CT, one per FC."""
    best_mode = {"cohen": "reporter", "shendure": "reporter", "seelig": "deflated"}
    base_color = {"reporter": COLOR_REP, "deflated": COLOR_DEFL}

    dsets = ["cohen", "shendure", "seelig"]
    heights = [max(DATASETS[d]["n_cell_types"], 1) for d in dsets]

    fig, axes = plt.subplots(
        3, 1, figsize=(6.5, 6),
        gridspec_kw={"height_ratios": heights},
        sharex=True,
    )

    for ax, name in zip(axes, dsets):
        mode = best_mode[name]
        triplet = FC_TRIPLET[name]  # ascending
        light, mid, dark = _shades(base_color[mode], n=3)

        # Sort CTs by the highest-FC power (strongest at top).
        sort_df = mpows[name][mode][triplet[-1]].sort_values("power")
        cts = sort_df["cell_type"].tolist()
        y = np.arange(len(cts))

        fc_lo, fc_mi, fc_hi = triplet
        by_fc = {fc: mpows[name][mode][fc].set_index("cell_type").loc[cts, "power"].values
                 for fc in triplet}

        # Draw longest first so shorter bars paint on top.
        for fc, color in ((fc_hi, light), (fc_mi, mid), (fc_lo, dark)):
            ax.barh(y, by_fc[fc], color=color, height=0.75, label=f"FC={fc}")

        ax.set_yticks(y)
        ax.set_yticklabels([display_ct(c, name) for c in cts])
        ax.set_xlim(0, 1.0)
        ax.axvline(0.8, color=MUTED, ls="--", lw=0.7)
        ax.set_title(DATASETS[name]["title"], loc="left")
        style_axes(ax, grid_axis="x")
        ax.invert_yaxis()
        ax.legend(loc="lower right", frameon=False, ncol=3,
                  handlelength=1.2, handletextpad=0.4, columnspacing=0.8)

    axes[-1].set_xlabel("power")
    plt.tight_layout()
    out = OUT_DIR / "fig_b_telescope.svg"
    fig.savefig(out, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


def fig_b_slope_discrete(mpows):
    """Three-point slope form of Fig 4B, kept for comparison with the curves.

    Each panel's fold changes are those of FC_TRIPLET, drawn as equally spaced
    categories. Zhao et al.'s triplet spans 0.2 in fold change and Lalanne et
    al.'s spans 0.8, so a line's slope means something different in each panel
    -- which is why the manuscript uses the continuous version.
    """
    best_mode = {"cohen": "reporter", "shendure": "reporter"}
    fig, axes = plt.subplots(1, 2, figsize=(FIG_W, 1.9), sharey=True)

    for ax, name in zip(axes, PAIRED_DSETS):
        mode = best_mode[name]
        triplet = FC_TRIPLET[name]
        ref_ct = DATASETS[name]["ref_display"]
        x = np.arange(len(triplet))

        for ct in mpows[name][mode][triplet[0]]["cell_type"]:
            v = [float(mpows[name][mode][fc].set_index("cell_type").loc[ct, "power"])
                 for fc in triplet]
            is_ref = ct == ref_ct
            ax.plot(x, v, marker="o", color=COLOR_REP,
                    alpha=1.0 if is_ref else 0.45,
                    lw=1.6 if is_ref else 0.8,
                    markersize=3.5 if is_ref else 2.5,
                    zorder=4 if is_ref else 3)
            if is_ref:
                ax.annotate(ref_ct, xy=(x[-1], v[-1]), xytext=(3, 0),
                            textcoords="offset points", fontsize=6, color=INK,
                            va="center", annotation_clip=False)

        ax.set_xticks(x)
        ax.set_xticklabels([f"{fc:g}" for fc in triplet])
        ax.set_xlim(-0.25, len(triplet) - 0.55)
        ax.set_ylim(0, 1.02)
        ax.set_yticks([0, 0.5, 1.0], ["0", "0.5", "1"])
        ax.axhline(0.8, color=MUTED, ls="--", lw=0.7, zorder=1)
        ax.set_title(DATASETS[name]["title"].split(" (")[0], loc="left", color=INK, pad=4)
        ax.set_xlabel("fold change", color=INK)
        style_axes(ax, grid_axis="y")

    axes[0].set_ylabel("power", color=INK)
    fig.subplots_adjust(left=0.12, right=0.86, top=0.84, bottom=0.24, wspace=0.35)
    out = OUT_DIR / "fig_b_slope_discrete.svg"
    fig.savefig(out, format="svg")
    plt.close(fig)
    print(f"Saved: {out}")


def report(powers):
    """Per-cell-type arms and gaps at TARGET_FC -- the numbers for the caption."""
    print(f"=== power at fc {TARGET_FC:g} (+/- {FC_TOL:.0%}) ===")
    for name in PAIRED_DSETS:
        rep = powers[name]["reporter"].set_index("cell_type")["power"]
        defl = powers[name]["deflated"].set_index("cell_type")["power"]
        gap = (rep - defl).sort_values(ascending=False)
        print(f"\n{DATASETS[name]['title']}")
        print(f"{'cell type':<24}{'+rep':>8}{'-rep':>8}{'gap':>8}")
        for ct in gap.index:
            print(f"{display_ct(ct, name):<24}{rep[ct]:>8.3f}{defl[ct]:>8.3f}"
                  f"{gap[ct]:>8.3f}")
        print(f"{'mean':<24}{rep.mean():>8.3f}{defl.mean():>8.3f}{gap.mean():>8.3f}")
    print()


def main():
    powers = load_powers()
    report(powers)

    fig_a_reporter_dumbbell(powers)

    print("binning power curves:")
    curves = load_power_curves()
    fig_b_power_curves(curves)

    mpows = load_multi_fc_powers()
    fig_b_slope_discrete(mpows)
    fig_b_dataset_bars(powers)
    fig_b_telescope(mpows)


if __name__ == "__main__":
    main()
