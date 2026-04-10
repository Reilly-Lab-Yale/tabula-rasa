#!/usr/bin/env python3
"""
FC-at-80%-power bar chart for Shendure (piggyBac) dataset.

Reads the aggregated power parquet files produced by
shendure_power_ttest_all_cell_types.py and finds, for each cell type, the
fold-change at which Welch's t-test first reaches 80% power. Produces a
grouped horizontal bar chart comparing +reporter vs -reporter (deflated).

Usage:
    python fc_at_80_power_shendure.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_repo_root = Path(__file__).resolve().parents[4]
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

import scMPRAforge as scm

OUTPUT_DIR = Path(__file__).resolve().parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

# Must match the simulation script's activity range.
MINP = scm.SHENDURE_BOUNDS.reference_activity
MAX_ACTIVITY = 4.0 * MINP
MAX_FC = MAX_ACTIVITY / MINP  # 4.0


def fc_at_80(df, ct, n_bins=200):
    """Find the FC where power first crosses 0.8, on each side of FC=1."""
    sub = df[df["cell_type"] == ct].copy()
    sub["fc_bin"] = pd.cut(sub["fc"], bins=n_bins)
    binned = (
        sub.groupby("fc_bin", observed=True)["reject_null"]
        .mean()
        .reset_index(name="power")
    )
    binned["bin_center"] = binned["fc_bin"].apply(lambda x: x.mid).astype(float)
    binned = binned.sort_values("bin_center")

    # Right side (FC > 1): first bin where power >= 0.8
    right = binned[binned["bin_center"] > 1.0]
    right_cross = right[right["power"] >= 0.8]
    if len(right_cross) > 0:
        fc_right = float(right_cross.iloc[0]["bin_center"])
    else:
        fc_right = MAX_FC  # capped

    # Left side (FC < 1): first bin (from FC=1 going left) where power >= 0.8
    left = binned[binned["bin_center"] < 1.0].sort_values(
        "bin_center", ascending=False
    )
    left_cross = left[left["power"] >= 0.8]
    if len(left_cross) > 0:
        fc_left = float(left_cross.iloc[0]["bin_center"])
    else:
        fc_left = 1.0 / MAX_FC  # capped

    return fc_right, fc_left


def main():
    rep = pd.read_parquet(OUTPUT_DIR / "power_df_reporter.parquet")
    defl = pd.read_parquet(OUTPUT_DIR / "power_df_deflated.parquet")

    cell_types = sorted(rep["cell_type"].unique().tolist())

    results = []
    for ct in cell_types:
        r_right, r_left = fc_at_80(rep, ct)
        d_right, d_left = fc_at_80(defl, ct)
        ct_display = "Pluripotent" if ct == "reference" else ct
        n_cells = scm.SHENDURE_BOUNDS.cells_per_cell_type.get(ct, "?")
        results.append({
            "cell_type": ct_display,
            "n_cells": n_cells,
            "reporter_up": r_right,
            "reporter_down": r_left,
            "deflated_up": d_right,
            "deflated_down": d_left,
            "reporter_up_capped": r_right >= MAX_FC - 0.01,
            "deflated_up_capped": d_right >= MAX_FC - 0.01,
            "reporter_down_capped": r_left <= 1.0 / MAX_FC + 0.01,
            "deflated_down_capped": d_left <= 1.0 / MAX_FC + 0.01,
        })

    res = pd.DataFrame(results)
    print(res.to_string(index=False))

    # --- Plot: FC at 80% power, both directions ---
    fig, ax = plt.subplots(figsize=(7, 5))

    y_pos = np.arange(len(res))
    bar_height = 0.35

    # Upregulation (right of 1.0)
    ax.barh(
        y_pos - bar_height / 2,
        res["reporter_up"] - 1.0, bar_height,
        left=1.0, color="steelblue", label="+reporter",
    )
    ax.barh(
        y_pos + bar_height / 2,
        res["deflated_up"] - 1.0, bar_height,
        left=1.0, color="coral", label="-reporter (deflated)",
    )

    # Downregulation (left of 1.0)
    ax.barh(
        y_pos - bar_height / 2,
        res["reporter_down"] - 1.0, bar_height,
        left=1.0, color="steelblue",
    )
    ax.barh(
        y_pos + bar_height / 2,
        res["deflated_down"] - 1.0, bar_height,
        left=1.0, color="coral",
    )

    # Add > / < annotations for capped values
    for i, row in res.iterrows():
        if row["reporter_up_capped"]:
            ax.text(
                row["reporter_up"] + 0.01, i - bar_height / 2, ">",
                va="center", ha="left", fontsize=10, fontweight="bold",
                color="steelblue",
            )
        if row["deflated_up_capped"]:
            ax.text(
                row["deflated_up"] + 0.01, i + bar_height / 2, ">",
                va="center", ha="left", fontsize=10, fontweight="bold",
                color="coral",
            )
        if row["reporter_down_capped"]:
            ax.text(
                row["reporter_down"] - 0.01, i - bar_height / 2, "<",
                va="center", ha="right", fontsize=10, fontweight="bold",
                color="steelblue",
            )
        if row["deflated_down_capped"]:
            ax.text(
                row["deflated_down"] - 0.01, i + bar_height / 2, "<",
                va="center", ha="right", fontsize=10, fontweight="bold",
                color="coral",
            )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(
        [f"{row.cell_type} ({row.n_cells})" for _, row in res.iterrows()],
        fontsize=8,
    )
    ax.set_xlabel("Fold-change at 80% power", fontsize=10)
    ax.set_title(
        "Shendure (piggyBac) -- FC needed to detect activity\n"
        "+reporter vs -reporter",
        fontsize=11,
    )
    ax.axvline(1.0, color="black", linestyle="-", lw=0.8)
    ax.legend(fontsize=8, loc="lower right")
    margin = 0.05
    ax.set_xlim(
        min(res["deflated_down"].min(), res["reporter_down"].min()) - margin,
        max(res["deflated_up"].max(), res["reporter_up"].max()) + margin,
    )

    plt.tight_layout()
    out = OUTPUT_DIR / "fc_at_80_power_shendure.svg"
    fig.savefig(out, format="svg", bbox_inches="tight")
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
