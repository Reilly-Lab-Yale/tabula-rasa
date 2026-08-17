"""Shared drawing code for the pairwise CRE-vs-CRE power heatmaps.

One panel per cell type: power over (baseline activity, |log2 FC|), with the
two published-effect reference lines from kircher_reference.py.

The three `*_reference.svg` panels are manuscript Fig 3B, placed side by side
in one row; the full per-cell-type sets are supplementary Figs S6 (Lalanne
et al., 10 panels), S7 (Zhao et al., 4) and S8 (Yin et al., 2). Every one of
those inclusions is at 0.32\\textwidth of a 6.90 in text block, so a panel is
2.21 in wide on the page. The manuscript sync applies no font scaling, so a
point size here is the point size a reader gets. Margins are fixed in inches
rather than left to tight_layout, so the three plot areas line up across the
Fig 3B row -- and, since the left and right margins match
activity_power/replot_overlay_ct.py, with the Fig 3A panels above them too.

Sibling figures: activity_power/replot_overlay_ct.py (Fig 3A) and
activity_power/main_figs_power_at_fc.py (Fig 4). Same type scale, same
dataset naming, same ink/muted/hairline greys.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

import kircher_reference

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
# Shared with the Fig 3A and Fig 4 scripts so a cell type is spelled the same
# way wherever it appears.
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

INK, MUTED = "#1a1a1a", "#666666"

# Power is a magnitude on a fixed 0-1 scale, so the ramp is sequential and its
# lightness is monotone (checked: relative luminance falls 0.971 -> 0.047 over
# the ramp, strictly). Kept as YlOrRd so this panel matches the contour and
# all-cell-type figures built from the same aggregates.
CMAP = "YlOrRd"

# Sized for a 0.32\textwidth inclusion in a 6.90 in text block.
FIG_W = 2.21
LEFT_IN, RIGHT_IN = 0.42, 0.09   # same as Fig 3A, so the columns align
TOP_IN, BOTTOM_IN = 0.56, 0.34
PLOT_H_IN = 0.95
FIG_H = TOP_IN + PLOT_H_IN + BOTTOM_IN

# One colourbar per panel, but a slim horizontal strip in the top margin
# instead of a full-height bar down the right-hand side. Each SVG is also used
# standalone in a supplementary grid, so no panel can be left without a scale;
# what can be cut is how much a repeated 0-1 scale costs. Sitting in the
# margin it takes no plot width at all, and three identical strips across the
# Fig 3B row read as one key repeated rather than as three legends.
CBAR_W_IN, CBAR_H_IN = 0.85, 0.05
CBAR_TOP_IN = 0.345          # from the top of the figure
TITLE_TOP_IN, SUBTITLE_TOP_IN = 0.115, 0.245

# Four x ticks out of 24 log-spaced baselines: enough to read the range the
# axis spans, and the widest label (five characters at 6.5 pt) still clears
# its neighbour across 1.70 in of axis. Five ticks does not.
N_XTICKS = 4

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
        return REFERENCE_NAMES[dataset]
    return CT_DISPLAY.get(ct, ct)


def power_pivot(power_grid, cell_type):
    """Power as rows of |log2 FC| (descending) by columns of baseline mu.

    Descending rows because kircher_reference draws onto a heatmap whose y
    axis runs top-down, with row i centred at i + 0.5.
    """
    ct_data = power_grid[power_grid["cell_type"] == cell_type]
    assert len(ct_data), f"no grid rows for cell type {cell_type!r}"

    pivot = ct_data.pivot_table(
        index="log2_fc", columns="baseline_mu", values="power"
    )
    pivot = pivot.sort_index(ascending=False)
    pivot = pivot[sorted(pivot.columns)]

    n_fc = ct_data["log2_fc"].nunique()
    n_mu = ct_data["baseline_mu"].nunique()
    assert pivot.shape == (n_fc, n_mu), (
        f"{cell_type}: pivot is {pivot.shape}, expected "
        f"({n_fc}, {n_mu}) from {len(ct_data)} grid rows")
    assert not pivot.isna().any().any(), (
        f"{cell_type}: {int(pivot.isna().sum().sum())} of {pivot.size} grid "
        f"cells have no power estimate")
    assert pivot.to_numpy().min() >= 0.0 and pivot.to_numpy().max() <= 1.0, (
        f"{cell_type}: power outside [0, 1], range "
        f"{pivot.to_numpy().min()} to {pivot.to_numpy().max()}")
    return pivot


def _fmt_mu(v):
    """Two significant figures, which is all 2.21 in of axis has room for.

    Trailing zeros kept, so a tick ladder does not mix 0.80 with 3 and read
    as though the two were quoted to different precision.
    """
    return f"{v:#.2g}"


def _set_xticks(ax, mu_vals):
    """Label N_XTICKS of the baseline columns, evenly spaced by column."""
    idx = np.unique(np.linspace(0, len(mu_vals) - 1, N_XTICKS).round().astype(int))
    ax.set_xticks(idx + 0.5)
    ax.set_xticklabels([_fmt_mu(mu_vals[i]) for i in idx], rotation=0)
    # The end labels are centred over columns that sit hard against the axes,
    # so they would hang into the y-tick gutter and past the right spine.
    labels = ax.get_xticklabels()
    labels[0].set_ha("left")
    labels[-1].set_ha("right")


def _set_yticks(ax, fc_vals):
    """Label every other fold-change row, top row included.

    Ten rows over 0.95 in is a 0.095 in pitch, which a 6.5 pt label very
    nearly fills; every other row leaves the ladder readable and the grid
    lines still show where the unlabelled rows are.
    """
    idx = np.arange(0, len(fc_vals), 2)
    ax.set_yticks(idx + 0.5)
    ax.set_yticklabels([f"{fc_vals[i]:.2f}" for i in idx], rotation=0)


def _draw_colorbar(fig, mesh):
    """Slim horizontal scale in the top margin, named at its right-hand end."""
    cax = fig.add_axes([
        LEFT_IN / FIG_W,
        (FIG_H - CBAR_TOP_IN - CBAR_H_IN) / FIG_H,
        CBAR_W_IN / FIG_W,
        CBAR_H_IN / FIG_H,
    ])
    cbar = fig.colorbar(mesh, cax=cax, orientation="horizontal")
    cbar.set_ticks([0.0, 0.5, 1.0])
    cbar.set_ticklabels(["0", "0.5", "1"])
    cbar.outline.set_visible(False)
    cax.tick_params(length=0, labelsize=6, colors=MUTED, pad=1.5)
    fig.text((LEFT_IN + CBAR_W_IN + 0.05) / FIG_W,
             1 - (CBAR_TOP_IN + CBAR_H_IN / 2) / FIG_H,
             "power", ha="left", va="center", fontsize=6.5, color=MUTED)


def _draw_title(fig, dataset, ct, n_cells):
    """Dataset above, the cell type and its depth below.

    Indented to the left edge of the axes rather than the figure: the
    manuscript overlays a bold panel letter on the top-left corner of the row
    via \\panel{}, and a flush-left title lands under it.
    """
    x = LEFT_IN / FIG_W
    return [
        fig.text(x, 1 - TITLE_TOP_IN / FIG_H, DISPLAY_NAME[dataset],
                 ha="left", va="center", fontsize=7.5, color=INK),
        fig.text(x, 1 - SUBTITLE_TOP_IN / FIG_H,
                 f"{display_ct(ct, dataset)}, {n_cells:,} cells",
                 ha="left", va="center", fontsize=6.5, color=MUTED),
    ]


def _assert_fits(fig, texts):
    """No title line may run past the right edge of the figure.

    A long cell-type name overflows silently otherwise: nothing is cropped,
    because the panel is saved at a fixed size, so the text simply runs off
    the SVG and disappears when the manuscript places it.
    """
    fig.canvas.draw()
    limit = FIG_W * fig.dpi
    for t in texts:
        right = t.get_window_extent(fig.canvas.get_renderer()).x1
        assert right <= limit, (
            f"title line {t.get_text()!r} ends {right / fig.dpi:.2f} in "
            f"across a {FIG_W:.2f} in panel")


def plot_panel(pivot, dataset, ct, n_cells, out_path):
    """One cell type's power surface. Returns the reference lines drawn."""
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    sns.heatmap(
        pivot, ax=ax, vmin=0, vmax=1, cmap=CMAP,
        linewidths=0.4, linecolor="white", cbar=False,
        xticklabels=False, yticklabels=False,
    )

    # kircher_reference's own marker size is left alone: at 4 pt on a 0.095 in
    # grid row the circles and squares are the one thing that still separates
    # the two rules by shape once the dash patterns are this short.
    drawn = kircher_reference.annotate(ax, pivot.index)
    assert drawn == 2, (
        f"{dataset}/{ct}: grid reaches only {drawn} of the two reference "
        f"effects; |log2 FC| runs {min(pivot.index)} to {max(pivot.index)}")

    _set_xticks(ax, list(pivot.columns))
    _set_yticks(ax, list(pivot.index))
    ax.tick_params(length=0, colors=INK, pad=2)
    ax.set_xlabel("baseline activity (counts per cell)", color=INK, labelpad=2)
    ax.set_ylabel("|log2 FC|", color=INK, labelpad=2)

    _draw_colorbar(fig, ax.collections[0])
    texts = _draw_title(fig, dataset, ct, n_cells)

    fig.subplots_adjust(
        left=LEFT_IN / FIG_W, right=1 - RIGHT_IN / FIG_W,
        top=1 - TOP_IN / FIG_H, bottom=BOTTOM_IN / FIG_H,
    )
    _assert_fits(fig, texts)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    # No bbox_inches="tight". The three Fig 3B panels sit in one row with
    # their plot areas in line, and cropping each to its own contents would
    # undo that: a long cell-type name is wider than a short one.
    fig.savefig(out_path, format="svg")
    plt.close(fig)
    return drawn


def render_dataset(power_df, dataset, cells_per_cell_type, panel_dir):
    """Every cell type of one dataset, one SVG each. Returns the paths."""
    power_grid = (
        power_df
        .groupby(["cell_type", "baseline_mu", "log2_fc"])
        .agg(power=("reject_null", "mean"),
             n_tests=("reject_null", "count"))
        .reset_index()
    )
    cell_types = sorted(power_grid["cell_type"].unique())
    assert len(cell_types) == N_CELL_TYPES[dataset], (
        f"{dataset}: expected {N_CELL_TYPES[dataset]} cell types, got "
        f"{len(cell_types)}: {cell_types}")
    assert power_grid["n_tests"].nunique() == 1, (
        f"{dataset}: grid cells carry unequal test counts, "
        f"{sorted(power_grid['n_tests'].unique())}")
    assert power_grid["power"].between(0.0, 1.0).all(), (
        f"{dataset}: power outside [0, 1], range "
        f"{power_grid['power'].min()} to {power_grid['power'].max()}")

    n_tests = int(power_grid["n_tests"].iloc[0])
    print(f"{DISPLAY_NAME[dataset]}: {len(cell_types)} cell types, "
          f"{power_grid['baseline_mu'].nunique()} baselines x "
          f"{power_grid['log2_fc'].nunique()} fold changes, "
          f"{n_tests} tests per grid cell")

    produced = []
    for ct in cell_types:
        n_cells = int(cells_per_cell_type[ct])
        assert n_cells > 0, f"{dataset}/{ct}: {n_cells} cells"
        pivot = power_pivot(power_grid, ct)
        safe_ct = ct.replace("/", "_").replace(" ", "_")
        out_path = panel_dir / f"{dataset}_pairwise_power_mwu_{safe_ct}.svg"
        drawn = plot_panel(pivot, dataset, ct, n_cells, out_path)
        produced.append(out_path)
        print(f"  {display_ct(ct, dataset):<24} {n_cells:>7,} cells  "
              f"max power {pivot.to_numpy().max():.2f}  "
              f"{drawn} reference lines  -> {out_path.name}")

    assert len(produced) == N_CELL_TYPES[dataset], (
        f"{dataset}: wrote {len(produced)} panels, expected "
        f"{N_CELL_TYPES[dataset]}")
    print(f"  {len(produced)} panels at {FIG_W:.2f} x {FIG_H:.2f} in")
    return produced
