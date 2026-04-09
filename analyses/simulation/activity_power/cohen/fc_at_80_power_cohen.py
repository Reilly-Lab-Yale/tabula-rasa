#!/usr/bin/env python3
"""
FC-at-80%-power bar chart for Cohen (episomal) dataset.

Reads the aggregated power parquet files produced by
cohen_power_ttest_all_cell_types.py and finds, for each cell type, the
fold-change at which Welch's t-test first reaches 80% power. Produces a
grouped horizontal bar chart comparing +reporter vs -reporter (deflated).

Usage:
    python fc_at_80_power_cohen.py
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

OUTPUT_DIR = Path(__file__).resolve().parent / "output"

MINP = 0.936273  # cohen reference_activity
MAX_FC = 1.6     # upper bound of sim range


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
    rep = pd.read_parquet(OUTPUT_DIR / "cohen_power_df_reporter.parquet")
    defl = pd.read_parquet(OUTPUT_DIR / "cohen_power_df_deflated.parquet")

    cell_types = sorted(rep["cell_type"].unique().tolist())

    results = []
    for ct in cell_types:
        r_right, r_left = fc_at_80(rep, ct)
        d_right, d_left = fc_at_80(defl, ct)
        ct_display = "Rod" if ct == "reference" else ct
        results.append({
            "cell_type": ct_display,
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

    # --- Plot: FC at 80% power, upregulation side only ---
    fig, ax = plt.subplots(figsize=(6, 4))

    y_pos = np.arange(len(res))
    bar_height = 0.35

    ax.barh(
        y_pos - bar_height / 2, res["reporter_up"], bar_height,
        color="steelblue", label="+reporter",
    )
    ax.barh(
        y_pos + bar_height / 2, res["deflated_up"], bar_height,
        color="coral", label="-reporter (deflated)",
    )

    # Add > annotations for capped values
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

    ax.set_yticks(y_pos)
    ax.set_yticklabels(res["cell_type"], fontsize=9)
    ax.set_xlabel("Fold-change at 80% power (upregulation)", fontsize=10)
    ax.set_title(
        "Cohen (episomal) -- FC needed to detect activity\n"
        "+reporter vs -reporter",
        fontsize=11,
    )
    ax.axvline(1.0, color="grey", linestyle=":", lw=0.8)
    ax.legend(fontsize=8, loc="lower right")
    ax.set_xlim(0.9, MAX_FC + 0.1)

    plt.tight_layout()
    out = OUTPUT_DIR / "fc_at_80_power_cohen.svg"
    fig.savefig(out, format="svg", bbox_inches="tight")
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
