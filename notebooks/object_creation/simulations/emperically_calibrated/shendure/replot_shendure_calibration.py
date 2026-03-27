"""Re-generate Shendure MWU calibration plots from saved sim data."""
import pickle
import scMPRAforge as scm
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import uniform

data_root = Path("/gpfs/gibbs/pi/reilly/tabula_data")
output_dir = Path("output")
sim_date = "2026-03-23"

with open(data_root / "shendure" / "shendure_ortho_20260306" / "by_cell_type_parameters.pkl", "rb") as f:
    _params = pickle.load(f)
n_cres_per_ct = {ct: len(_params.nb[ct]) for ct in _params.nb}
cell_types = sorted(n_cres_per_ct.keys())
del _params

# Reload sims from disk (no cluster needed — just reading saved test results)
from dask.distributed import LocalCluster, Client
cluster = LocalCluster(n_workers=1, threads_per_worker=1, memory_limit="8G")
client = Client(cluster)

sim_dir = data_root / f"{sim_date}_shendure_cal"
all_pvals = {ct: [] for ct in cell_types}

for ct in cell_types:
    ct_dir = sim_dir / ct
    if not ct_dir.exists():
        print(f"WARNING: no sim dir for {ct}")
        continue
    for sim_path in sorted(ct_dir.iterdir()):
        mwu_dir = sim_path / "tests" / "hs_all_ct" / "mwu"
        if not (mwu_dir.exists() and any(mwu_dir.iterdir())):
            continue
        sim = scm.de_novo_simulation(location=ct_dir, name=sim_path.name, client=client)
        for i in range(sim.get_state_field("n_sims")):
            m = sim._merge_in_ground_truth(
                hypothesis_set_name="hs_all_ct", test_type="mwu", index=i
            )
            all_pvals[ct].append(m["p_value"].values)

all_pvals = {ct: np.concatenate(v) for ct, v in all_pvals.items() if v}
for ct in cell_types:
    if ct in all_pvals:
        fpr = np.mean(all_pvals[ct] < 0.05)
        print(f"  {ct}: n={len(all_pvals[ct])}, FPR@0.05={fpr:.4f}")

client.close()
cluster.close()

STYLES = {
    "pub": dict(
        fig_col_w=5, fig_row_h=3,
        title_fs=9, axlabel_fs=8, tick_fs=7, annot_fs=7, suptitle_fs=12,
        scatter_s=1, suffix="",
    ),
    "pres": dict(
        fig_col_w=8, fig_row_h=5,
        title_fs=16, axlabel_fs=14, tick_fs=12, annot_fs=12, suptitle_fs=18,
        scatter_s=3, suffix="_presentation",
    ),
}

for style_name, s in STYLES.items():
    ncols = 2
    nrows = len(cell_types)
    fig, axes = plt.subplots(nrows, ncols, figsize=(s["fig_col_w"] * ncols, s["fig_row_h"] * nrows))
    if nrows == 1:
        axes = axes[np.newaxis, :]

    palette = sns.color_palette("tab20", n_colors=len(cell_types))

    for i, ct in enumerate(cell_types):
        ct_display = "Pluripotent" if ct == "reference" else ct
        pvals = all_pvals.get(ct, np.array([]))
        pvals = pvals[np.isfinite(pvals)]
        color = palette[i]
        fpr = np.mean(pvals < 0.05) if len(pvals) else float("nan")
        cells = scm.SHENDURE_BOUNDS.cells_per_cell_type.get(ct, "?")

        # Histogram
        ax_h = axes[i, 0]
        ax_h.hist(pvals, bins=50, density=True, edgecolor="black",
                  linewidth=0.3, alpha=0.7, color=color)
        ax_h.axhline(1.0, color="red", linestyle="--", lw=1, label="Uniform(0,1)")
        ax_h.set_title(f"{ct_display}", fontsize=s["title_fs"])
        ax_h.set_xlabel("p-value", fontsize=s["axlabel_fs"])
        ax_h.set_ylabel("Density", fontsize=s["axlabel_fs"])
        ax_h.tick_params(labelsize=s["tick_fs"])
        ax_h.set_xlim(0, 1)
        ax_h.text(
            0.97, 0.92,
            f"FPR@0.05 = {fpr:.3f}\nn={len(pvals)}, cells={cells}",
            transform=ax_h.transAxes, ha="right", va="top", fontsize=s["annot_fs"],
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8),
        )
        if i == 0:
            ax_h.legend(fontsize=s["annot_fs"], loc="upper left")

        # QQ plot
        ax_q = axes[i, 1]
        observed = np.sort(pvals)
        n = len(observed)
        expected = uniform.ppf(np.linspace(0, 1, n + 2)[1:-1])
        ax_q.scatter(expected, observed, s=s["scatter_s"], alpha=0.5, color=color)
        ax_q.plot([0, 1], [0, 1], "r--", lw=1)
        ax_q.set_title(f"{ct_display} — QQ", fontsize=s["title_fs"])
        ax_q.set_xlabel("Expected (Uniform)", fontsize=s["axlabel_fs"])
        ax_q.set_ylabel("Observed p-value", fontsize=s["axlabel_fs"])
        ax_q.tick_params(labelsize=s["tick_fs"])
        ax_q.set_xlim(0, 1)
        ax_q.set_ylim(0, 1)
        ax_q.set_aspect("equal")

    fig.suptitle("MWU p-value calibration under null — Shendure",
                 fontsize=s["suptitle_fs"], y=1.01)
    plt.tight_layout()
    svg_path = output_dir / f"calibration_mwu_all_cell_types{s['suffix']}.svg"
    fig.savefig(svg_path, format="svg", bbox_inches="tight")
    print(f"Saved: {svg_path}")
    plt.close()
