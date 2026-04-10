#!/usr/bin/env python3
"""
Empirical CRE effect-size distribution for Cohen (episomal) dataset.

Loads the obsingle (observed-only) training data, computes per-CRE mean UMI
per cell type, and plots the fold-change distribution relative to the median
CRE. This contextualizes the power analysis: power to detect FC=X is only
meaningful if real CREs actually exhibit FC=X.

Usage:
    python cohen_empirical_fc.py
"""

import sys
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_repo_root = Path(__file__).resolve().parents[4]
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

OUTPUT_DIR = Path(__file__).resolve().parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

ORTHO_DIR = Path(
    "/nfs/roberts/project/pi_skr2/shared/tabula_data_new"
    "/cohen/cohen_obsingle_nb_phantom"
)


def main():
    with open(ORTHO_DIR / "training_data.pkl", "rb") as f:
        td = pickle.load(f)

    df = td.data.compute()
    for col in df.columns:
        if hasattr(df[col], "sparse"):
            df[col] = df[col].sparse.to_dense()

    print(f"Total rows: {len(df):,}")
    print(f"Cell types: {sorted(df['cell_type'].unique())}")
    print(f"CREs: {df['cre_id'].nunique()}")

    # Per (cell_type, cre_id): mean UMI across all observations
    mean_raw = (
        df.groupby(["cell_type", "cre_id"])["umis_mpra_bc"]
        .mean()
        .reset_index()
    )
    mean_raw.columns = ["cell_type", "cre_id", "mean_umi"]

    # FC relative to median CRE within each cell type
    cell_types = sorted(mean_raw["cell_type"].unique())

    for ct in cell_types:
        sub = mean_raw[mean_raw["cell_type"] == ct]
        ref = sub["mean_umi"].median()
        fc = sub["mean_umi"] / ref
        ct_display = "Rod" if ct == "reference" else ct
        print(f"\n  {ct_display}: {len(sub)} CREs, median UMI={ref:.4f}")
        print(f"    FC range: [{fc.min():.3f}, {fc.max():.3f}]")
        pcts = [5, 25, 50, 75, 95]
        fc_pcts = np.percentile(fc.dropna(), pcts)
        print(
            "    FC percentiles:",
            {p: f"{v:.3f}" for p, v in zip(pcts, fc_pcts)},
        )

    # --- Plot ---
    fig, axes = plt.subplots(1, len(cell_types), figsize=(14, 3.5), sharey=True)

    for ax, ct in zip(axes, cell_types):
        sub = mean_raw[mean_raw["cell_type"] == ct].copy()
        ref = sub["mean_umi"].median()
        sub["fc"] = sub["mean_umi"] / ref
        ct_display = "Rod" if ct == "reference" else ct

        ax.hist(
            sub["fc"], bins=25,
            color="steelblue", edgecolor="white", alpha=0.8,
        )
        ax.axvline(1.0, color="black", linestyle="-", lw=0.8)
        ax.set_title(ct_display, fontsize=10)
        ax.set_xlabel("FC (vs median CRE)", fontsize=9)
        ax.set_xlim(0.8, 1.7)

    axes[0].set_ylabel("Number of CREs", fontsize=9)
    fig.suptitle(
        "Cohen -- empirical CRE effect sizes (raw, per cell type)",
        fontsize=12,
    )
    plt.tight_layout()
    out = OUTPUT_DIR / "cohen_empirical_fc_distribution.svg"
    fig.savefig(out, format="svg", bbox_inches="tight")
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
