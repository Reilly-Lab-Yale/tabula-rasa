"""Re-generate Shendure pairwise power plots from cached parquet.

Produces both publication and presentation versions of each figure.
"""
import pickle
import scMPRAforge as scm
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

data_root = Path("/gpfs/gibbs/pi/reilly/tabula_data")
output_dir = Path("output")

n_library_reps = 20
n_sims = 5

with open(data_root / "shendure" / "shendure_ortho_20260306" / "by_cell_type_parameters.pkl", "rb") as f:
    _params = pickle.load(f)
n_cres_per_ct = {ct: len(_params.nb[ct]) for ct in _params.nb}
cell_types = sorted(n_cres_per_ct.keys())
del _params

power_df = pd.read_parquet(output_dir / "power_df.parquet")
power_grid = (
    power_df
    .groupby(["cell_type", "baseline_mu", "log2_fc"])
    .agg(power=("reject_null", "mean"), n_tests=("reject_null", "count"))
    .reset_index()
)

STYLES = {
    "pub": dict(
        fig_cell_w=7, fig_cell_h=5,
        title_fs=10, axlabel_fs=9, tick_fs=7, annot_fs=7, suptitle_fs=13,
        contour_figsize=(10, 6), contour_title_fs=12, contour_axlabel_fs=11,
        contour_legend_fs=8, contour_annot_fs=8,
        markersize=4, lw=1.5,
        suffix="",
    ),
    "pres": dict(
        fig_cell_w=9, fig_cell_h=7,
        title_fs=18, axlabel_fs=16, tick_fs=13, annot_fs=13, suptitle_fs=20,
        contour_figsize=(14, 8), contour_title_fs=20, contour_axlabel_fs=18,
        contour_legend_fs=14, contour_annot_fs=14,
        markersize=7, lw=2.5,
        suffix="_presentation",
    ),
}

for style_name, s in STYLES.items():
    # --- Heatmap ---
    ncols = 2
    nrows = int(np.ceil(len(cell_types) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(s["fig_cell_w"] * ncols, s["fig_cell_h"] * nrows))
    axes = axes.flatten()

    for i, ct in enumerate(cell_types):
        ct_display = "Pluripotent" if ct == "reference" else ct
        ax = axes[i]
        ct_data = power_grid[power_grid["cell_type"] == ct]

        pivot = ct_data.pivot_table(index="log2_fc", columns="baseline_mu", values="power")
        pivot = pivot.sort_index(ascending=False)
        pivot = pivot[sorted(pivot.columns)]

        cells = scm.SHENDURE_BOUNDS.cells_per_cell_type.get(ct, "?")
        n_cres = n_cres_per_ct[ct]

        sns.heatmap(
            pivot, ax=ax, vmin=0, vmax=1,
            cmap="YlOrRd", linewidths=0.5, linecolor="white",
            cbar_kws={"label": "Power", "shrink": 0.8},
            xticklabels=[f"{x:.2f}" for x in sorted(pivot.columns)],
            yticklabels=[f"{y:.2f}" for y in pivot.index],
        )
        ax.collections[0].colorbar.ax.tick_params(labelsize=s["tick_fs"])
        ax.collections[0].colorbar.set_label("Power", fontsize=s["axlabel_fs"])

        fc_vals = list(pivot.index)
        if 0.20 in fc_vals:
            y_pos = fc_vals.index(0.20) + 0.5
            ax.axhline(y=y_pos, color="blue", linestyle="--", lw=1.5, alpha=0.7)
            ax.text(0.02, y_pos - 0.3, "90% of variant\neffects below",
                    transform=ax.get_yaxis_transform(), fontsize=s["annot_fs"],
                    color="blue", va="top")

        ax.set_title(f"{ct_display}\n(n_cres={n_cres}, cells={cells})", fontsize=s["title_fs"])
        ax.set_xlabel("Baseline activity (mu)", fontsize=s["axlabel_fs"])
        ax.set_ylabel("log2(FC)", fontsize=s["axlabel_fs"])
        ax.tick_params(labelsize=s["tick_fs"])
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

    for j in range(len(cell_types), len(axes)):
        axes[j].set_visible(False)

    fig.suptitle("MWU pairwise power — Shendure\n"
                 f"(α=0.05, {n_library_reps}×{n_sims} reps per grid cell)",
                 fontsize=s["suptitle_fs"], y=1.01)
    plt.tight_layout()

    svg_path = output_dir / f"pairwise_power_mwu_all_cell_types{s['suffix']}.svg"
    fig.savefig(svg_path, format="svg", bbox_inches="tight")
    print(f"Saved: {svg_path}")
    plt.close()

    # --- 50% contour ---
    power_threshold = 0.50
    fig, ax = plt.subplots(figsize=s["contour_figsize"])
    palette = sns.color_palette("tab10", n_colors=len(cell_types))

    for i, ct in enumerate(cell_types):
        ct_display = "Pluripotent" if ct == "reference" else ct
        ct_data = power_grid[power_grid["cell_type"] == ct]

        baselines = sorted(ct_data["baseline_mu"].unique())
        min_fc = []
        valid_baselines = []

        for bl in baselines:
            bl_data = ct_data[ct_data["baseline_mu"] == bl].sort_values("log2_fc")
            above = bl_data[bl_data["power"] >= power_threshold]
            if len(above) > 0:
                min_fc.append(above["log2_fc"].min())
                valid_baselines.append(bl)

        cells = scm.SHENDURE_BOUNDS.cells_per_cell_type.get(ct, "?")
        label = f"{ct_display} ({cells} cells)"

        if valid_baselines:
            ax.plot(valid_baselines, min_fc, "o-", color=palette[i],
                    label=label, markersize=s["markersize"], lw=s["lw"])
        else:
            ax.plot([], [], "o-", color=palette[i], label=f"{label} (no {int(power_threshold*100)}% power)")

    ax.axhline(y=0.20, color="gray", linestyle="--", lw=1, alpha=0.7)
    ax.set_xscale("log")
    ax.set_xlabel("Baseline activity (mu)", fontsize=s["contour_axlabel_fs"])
    ax.set_ylabel(f"Min log2(FC) for {int(power_threshold*100)}% power", fontsize=s["contour_axlabel_fs"])
    ax.set_title(f"Minimum detectable fold change — Shendure (MWU, α=0.05, {int(power_threshold*100)}% power)",
                 fontsize=s["contour_title_fs"])
    ax.legend(fontsize=s["contour_legend_fs"], loc="upper right", ncol=2)
    ax.tick_params(labelsize=s["contour_legend_fs"])
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    ax.text(ax.get_xlim()[1] * 0.95, 0.205, "90% of variant effects",
            ha="right", va="bottom", fontsize=s["contour_annot_fs"], color="gray")

    svg_path2 = output_dir / f"pairwise_power_50pct_contours{s['suffix']}.svg"
    fig.savefig(svg_path2, format="svg", bbox_inches="tight")
    print(f"Saved: {svg_path2}")
    plt.close()
