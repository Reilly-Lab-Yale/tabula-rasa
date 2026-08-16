"""Replot activity-power 'minP vs CRE' as a single axis with cell types
overlaid as colored curves. One SVG per (dataset, condition).
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

BASE = Path(__file__).parent

REFERENCE_NAMES = {
    "cohen": "Rod",
    "shendure": "Pluripotent",
    "seelig": "HepG2",
    "takeshi": "HepG2",
}

# The reporter/deflated pair is MWU. Measuring the reporter's power benefit
# with the t-test confounds it: the t-test over-rejects in exactly the deflated
# condition, so part of the apparent benefit is that miscalibration rather than
# the reporter. The 2026-08-13 run added the MWU arms for this reason -- see the
# ARMS comment in shendure/shendure_power_ttest_all_cell_types.py.
# Yin et al. is still t-test: no MWU arm has been run on its power sims.
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
        ("deflated", "seelig/output/seelig_power_df_deflated.parquet"),
    ],
    "takeshi": [
        ("reporter", "takeshi/output/takeshi_power_df_reporter.parquet"),
        ("deflated", "takeshi/output/takeshi_power_df_deflated.parquet"),
    ],
}


def pow_curve_data(df, n_bins=100):
    df = df.copy()
    df["fc_bin"] = pd.cut(df["fc"], bins=n_bins)
    binned = (
        df.groupby("fc_bin", observed=True)["reject_null"]
        .mean()
        .reset_index(name="reject_frac")
    )
    binned["bin_center"] = binned["fc_bin"].apply(lambda x: x.mid).astype(float)
    return binned


def load_parquet(path, condition):
    df = pd.read_parquet(path)
    if "fc" not in df.columns and "log2_fc" in df.columns:
        df = df.copy()
        df["fc"] = 2.0 ** df["log2_fc"].astype(float)
    return df


def display_ct(ct, dataset):
    if ct == "reference":
        return REFERENCE_NAMES[dataset]
    return ct


def plot_overlay(df, dataset, condition, out_path):
    cell_types = sorted(df["cell_type"].unique().tolist())
    palette = sns.color_palette("tab10", n_colors=len(cell_types))

    fig, ax = plt.subplots(figsize=(7, 5))
    for color, ct in zip(palette, cell_types):
        sub = df[df["cell_type"] == ct]
        binned = pow_curve_data(sub)
        ax.plot(
            binned["bin_center"], binned["reject_frac"],
            color=color, marker="o", markersize=2, linewidth=1,
            label=display_ct(ct, dataset),
        )

    ax.axhline(0.8, color="black", linestyle="--", lw=0.8)
    ax.axvline(1.0, color="grey", linestyle=":", lw=0.8)
    ax.set_ylim(0, 1)
    ax.set_xlabel("FC (activity / minP)")
    ax.set_ylabel("Power (TPR)")
    ax.set_title(f"{dataset.capitalize()} -- {condition} (minP power)")
    ax.legend(title="Cell type", loc="best", fontsize=8)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="svg")
    plt.close(fig)


def main():
    produced = []
    skipped = []
    for dataset, items in DATASETS.items():
        for condition, rel in items:
            path = BASE / rel
            if not path.exists():
                skipped.append((str(path), "missing"))
                continue
            try:
                df = load_parquet(path, condition)
            except Exception as e:
                skipped.append((str(path), f"load error: {e}"))
                continue
            required = {"reject_null", "cell_type"}
            if not required.issubset(df.columns) or "fc" not in df.columns:
                skipped.append((str(path), f"schema: cols={list(df.columns)}"))
                continue
            out_path = (
                BASE / dataset / "output" / "overlay_ct" /
                f"{dataset}_power_{condition}_overlay_ct.svg"
            )
            plot_overlay(df, dataset, condition, out_path)
            produced.append(str(out_path))

    print("PRODUCED:")
    for p in produced:
        print(" ", p)
    print(f"Total: {len(produced)}")
    if skipped:
        print("SKIPPED:")
        for p, why in skipped:
            print(" ", p, "->", why)


if __name__ == "__main__":
    main()
