"""Big multi-panel calibration plot (rows = cell types, cols = [hist, QQ])
for any (dataset, test, condition).

Reads the cached null-pvals parquet ({ds}_null_pvals_{[mwu_]}{cond}.parquet)
so no Dask cluster or raw-sim loading is needed. QQ scatter is capped at
50k points per panel to keep SVGs under GitHub's 100MB limit.

Output filename matches the existing convention:
  calibration_{test}_all_cell_types[_deflated].svg

Usage:
    python replot_calibration_grand.py --dataset shendure --test ttest --condition reporter
    python replot_calibration_grand.py --dataset seelig --test mwu --condition deflated
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import uniform

import scMPRAforge as scm

HERE = Path(__file__).resolve().parent

DATASET_DISPLAY = {
    "shendure": "Shendure",
    "cohen":    "Cohen (episomal)",
    "seelig":   "Seelig (no reporter)",
    "takeshi":  "Takeshi",
}
REFERENCE_DISPLAY = {
    "shendure": "Pluripotent",
    "cohen":    "Rod",
    "seelig":   "HepG2",
    "takeshi":  "HepG2",
}
BOUNDS_NAME = {
    "shendure": "SHENDURE_BOUNDS",
    "cohen":    "COHEN_BOUNDS",
    "seelig":   "SEELIG_BOUNDS",
    "takeshi":  "TAKESHI_BOUNDS",
}

QQ_SCATTER_CAP = 50000  # cap QQ-scatter points per panel; FPR / shape unchanged


def _pvals_path(dataset: str, test: str, condition: str) -> Path:
    test_tag = "" if test == "ttest" else f"{test}_"
    return HERE / dataset / "output" / f"{dataset}_null_pvals_{test_tag}{condition}.parquet"


def _output_path(dataset: str, test: str, condition: str) -> Path:
    cond_tag = "" if condition == "reporter" else f"_{condition}"
    return HERE / dataset / "output" / f"calibration_{test}_all_cell_types{cond_tag}.svg"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True, choices=list(DATASET_DISPLAY))
    parser.add_argument("--test", required=True, choices=["ttest", "mwu"])
    parser.add_argument("--condition", required=True, choices=["reporter", "deflated"])
    args = parser.parse_args()

    pvals_path = _pvals_path(args.dataset, args.test, args.condition)
    if not pvals_path.exists():
        sys.exit(f"Missing pvals parquet: {pvals_path}")

    df = pd.read_parquet(pvals_path)
    bounds = getattr(scm, BOUNDS_NAME[args.dataset])
    ref_display = REFERENCE_DISPLAY[args.dataset]
    ds_display = DATASET_DISPLAY[args.dataset]
    test_label = "MWU" if args.test == "mwu" else "t-test"
    cond_label = "+reporter" if args.condition == "reporter" else "-reporter (deflated)"

    counts = df["cell_type"].value_counts()
    cell_types = ["reference"] + [c for c in counts.index if c != "reference"]
    cell_types = [c for c in cell_types if c in counts.index]

    palette = sns.color_palette("tab20", n_colors=len(cell_types))

    nrows = len(cell_types)
    ncols = 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3 * nrows))
    if nrows == 1:
        axes = axes[np.newaxis, :]

    for i, ct in enumerate(cell_types):
        ct_display = ref_display if ct == "reference" else ct
        pvals = df.loc[df["cell_type"] == ct, "p_value"].to_numpy()
        pvals = pvals[np.isfinite(pvals)]
        color = palette[i]
        fpr = float(np.mean(pvals < 0.05)) if len(pvals) else float("nan")
        cells_ct = bounds.cells_per_cell_type.get(ct, "?")

        # Histogram
        ax_h = axes[i, 0]
        ax_h.hist(pvals, bins=50, density=True, edgecolor="black",
                  linewidth=0.3, alpha=0.7, color=color)
        ax_h.axhline(1.0, color="red", linestyle="--", lw=1, label="Uniform(0,1)")
        ax_h.set_title(f"{ct_display}", fontsize=9)
        ax_h.set_xlabel("p-value", fontsize=8)
        ax_h.set_ylabel("Density", fontsize=8)
        ax_h.tick_params(labelsize=7)
        ax_h.set_xlim(0, 1)
        ax_h.text(
            0.97, 0.92,
            f"FPR@0.05 = {fpr:.3f}\nn={len(pvals)}, cells={cells_ct}",
            transform=ax_h.transAxes, ha="right", va="top", fontsize=7,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8),
        )
        if i == 0:
            ax_h.legend(fontsize=7, loc="upper left")

        # QQ plot (cap scatter at QQ_SCATTER_CAP points to keep SVG under
        # GitHub's 100MB limit; FPR / shape are unchanged).
        ax_q = axes[i, 1]
        observed = np.sort(pvals)
        n = len(observed)
        expected = uniform.ppf(np.linspace(0, 1, n + 2)[1:-1])
        if n > QQ_SCATTER_CAP:
            idx = np.linspace(0, n - 1, QQ_SCATTER_CAP).astype(int)
            observed = observed[idx]
            expected = expected[idx]
        ax_q.scatter(expected, observed, s=1, alpha=0.5, color=color)
        ax_q.plot([0, 1], [0, 1], "r--", lw=1)
        ax_q.set_title(f"{ct_display} -- QQ", fontsize=9)
        ax_q.set_xlabel("Expected (Uniform)", fontsize=8)
        ax_q.set_ylabel("Observed p-value", fontsize=8)
        ax_q.tick_params(labelsize=7)
        ax_q.set_xlim(0, 1)
        ax_q.set_ylim(0, 1)
        ax_q.set_aspect("equal")

    fig.suptitle(
        f"{test_label} p-value calibration under null -- {ds_display}, {cond_label}",
        fontsize=12, y=1.005,
    )
    plt.tight_layout()
    out = _output_path(args.dataset, args.test, args.condition)
    fig.savefig(out, format="svg", bbox_inches="tight")
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
