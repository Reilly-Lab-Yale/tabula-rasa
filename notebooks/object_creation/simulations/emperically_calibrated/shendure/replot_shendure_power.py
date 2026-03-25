"""Re-generate Shendure pairwise power plots from cached parquet."""
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

# --- Heatmap ---
ncols = 2
nrows = int(np.ceil(len(cell_types) / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 5 * nrows))
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

    fc_vals = list(pivot.index)
    if 0.20 in fc_vals:
        y_pos = fc_vals.index(0.20) + 0.5
        ax.axhline(y=y_pos, color="blue", linestyle="--", lw=1.5, alpha=0.7)
        ax.text(0.02, y_pos - 0.3, "90% of variant\neffects below",
                transform=ax.get_yaxis_transform(), fontsize=7,
                color="blue", va="top")

    ax.set_title(f"{ct_display}\n(n_cres={n_cres}, cells={cells})", fontsize=10)
    ax.set_xlabel("Baseline activity (mu)", fontsize=9)
    ax.set_ylabel("log2(FC)", fontsize=9)
    ax.tick_params(labelsize=7)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

for j in range(len(cell_types), len(axes)):
    axes[j].set_visible(False)

fig.suptitle("MWU pairwise power — Shendure\n"
             f"(α=0.05, {n_library_reps}×{n_sims} reps per grid cell)",
             fontsize=13, y=1.01)
plt.tight_layout()

svg_path = output_dir / "pairwise_power_mwu_all_cell_types.svg"
fig.savefig(svg_path, format="svg", bbox_inches="tight")
print(f"Saved: {svg_path}")
plt.close()

# --- 50% contour ---
power_threshold = 0.50
fig, ax = plt.subplots(figsize=(10, 6))
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
                label=label, markersize=4, lw=1.5)
    else:
        ax.plot([], [], "o-", color=palette[i], label=f"{label} (no {int(power_threshold*100)}% power)")

ax.axhline(y=0.20, color="gray", linestyle="--", lw=1, alpha=0.7)
ax.set_xscale("log")
ax.set_xlabel("Baseline activity (mu)", fontsize=11)
ax.set_ylabel(f"Min log2(FC) for {int(power_threshold*100)}% power", fontsize=11)
ax.set_title(f"Minimum detectable fold change — Shendure (MWU, α=0.05, {int(power_threshold*100)}% power)", fontsize=12)
ax.legend(fontsize=8, loc="upper right", ncol=2)
ax.grid(True, alpha=0.3)
plt.tight_layout()
ax.text(ax.get_xlim()[1] * 0.95, 0.205, "90% of variant effects",
        ha="right", va="bottom", fontsize=8, color="gray")

svg_path2 = output_dir / "pairwise_power_50pct_contours.svg"
fig.savefig(svg_path2, format="svg", bbox_inches="tight")
print(f"Saved: {svg_path2}")
plt.close()
