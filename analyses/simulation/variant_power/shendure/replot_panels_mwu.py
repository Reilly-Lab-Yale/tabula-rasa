#!/usr/bin/env python3
"""Replot Shendure MWU pairwise power heatmaps as separate per-CT SVGs.

Reads shendure_pairwise_power_df.parquet and emits one SVG per cell type into
output/panels_mwu/. No simulation -- pure plotting.
"""

import sys
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

_repo_root = Path(__file__).resolve().parents[4]
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import scMPRAforge as scm
import kircher_reference

OUTPUT_DIR = Path(__file__).resolve().parent / "output"
PANEL_DIR = OUTPUT_DIR / "panels_mwu"
PANEL_DIR.mkdir(exist_ok=True)

DATA_ROOT = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new")
ORTHO_DIR = DATA_ROOT / "shendure" / "shendure_obs_nb"
with open(ORTHO_DIR / "by_cell_type_parameters.pkl", "rb") as f:
    _params = pickle.load(f)
N_CRES_PER_CT = {ct: len(_params.nb[ct]) for ct in sorted(_params.nb.keys())}
del _params

N_LIBRARY_REPS = 20
N_SIMS = 5

power_df = pd.read_parquet(OUTPUT_DIR / "shendure_pairwise_power_df.parquet")
power_grid = (
    power_df
    .groupby(["cell_type", "baseline_mu", "log2_fc"])
    .agg(power=("reject_null", "mean"),
         n_tests=("reject_null", "count"))
    .reset_index()
)

CELL_TYPES = sorted(power_grid["cell_type"].unique())

for ct in CELL_TYPES:
    ct_display = "Pluripotent" if ct == "reference" else ct
    ct_data = power_grid[power_grid["cell_type"] == ct]

    pivot = ct_data.pivot_table(
        index="log2_fc", columns="baseline_mu", values="power"
    )
    pivot = pivot.sort_index(ascending=False)
    pivot = pivot[sorted(pivot.columns)]

    cells = scm.SHENDURE_BOUNDS.cells_per_cell_type.get(ct, "?")
    n_cres = N_CRES_PER_CT.get(ct, "?")

    fig, ax = plt.subplots(figsize=(7, 5))
    sns.heatmap(
        pivot, ax=ax, vmin=0, vmax=1,
        cmap="YlOrRd", linewidths=0.5, linecolor="white",
        cbar_kws={"label": "Power", "shrink": 0.8},
        xticklabels=[f"{x:.3f}" for x in sorted(pivot.columns)],
        yticklabels=[f"{y:.2f}" for y in pivot.index],
    )

    drawn = kircher_reference.annotate(ax, pivot.index)
    assert drawn == 2, (
        f"{ct}: grid reaches only {drawn} of the two reference effects; "
        f"fold-change range is {min(pivot.index)}-{max(pivot.index)}")

    ax.set_title(
        f"MWU pairwise power -- Shendure -- {ct_display}\n"
        f"(n_cres={n_cres}, cells={cells}, alpha=0.05, "
        f"{N_LIBRARY_REPS}x{N_SIMS} reps per grid cell)",
        fontsize=10,
    )
    ax.set_xlabel("Baseline activity (mu)", fontsize=9)
    ax.set_ylabel("log2(FC)", fontsize=9)
    ax.tick_params(labelsize=7)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

    plt.tight_layout()
    safe_ct = ct.replace("/", "_").replace(" ", "_")
    out_path = PANEL_DIR / f"shendure_pairwise_power_mwu_{safe_ct}.svg"
    fig.savefig(out_path, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}", flush=True)
