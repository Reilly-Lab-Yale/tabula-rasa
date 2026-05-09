"""Single-panel QQ-overlay calibration plot.

Reads a cached null pvals parquet (cell_type, p_value) and overlays one
semi-transparent line per cell type on a single QQ panel. Legend shows
FPR@0.05 per cell type.

Usage:
    python qq_overlay.py --dataset shendure --test ttest --condition reporter
    python qq_overlay.py --dataset cohen --test mwu --condition deflated
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import scMPRAforge  # noqa: F401  (sets editable-text rcParams on import)

HERE = Path(__file__).resolve().parent

# Per-dataset reference-cell-type display name.
REFERENCE_DISPLAY = {
    "shendure": "Pluripotent",
    "cohen": "Rod",
    "seelig": "HepG2",
    "takeshi": "HepG2",
}

QQ_N_QUANTILES = 2000  # cap per-CT line points so SVG stays compact


def _pvals_path(dataset: str, test: str, condition: str) -> Path:
    test_tag = "" if test == "ttest" else f"{test}_"
    return HERE / dataset / "output" / f"{dataset}_null_pvals_{test_tag}{condition}.parquet"


def _output_path(dataset: str, test: str, condition: str) -> Path:
    return HERE / dataset / "output" / f"{dataset}_calibration_{test}_qq_{condition}_overlay.svg"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True, choices=list(REFERENCE_DISPLAY))
    parser.add_argument("--test", required=True, choices=["ttest", "mwu"])
    parser.add_argument("--condition", required=True, choices=["reporter", "deflated"])
    args = parser.parse_args()

    pvals_path = _pvals_path(args.dataset, args.test, args.condition)
    if not pvals_path.exists():
        sys.exit(f"Missing pvals parquet: {pvals_path}")

    df = pd.read_parquet(pvals_path)
    ref_display = REFERENCE_DISPLAY[args.dataset]

    # Order: reference first, then by descending n.
    counts = df["cell_type"].value_counts()
    cell_types = ["reference"] + [ct for ct in counts.index if ct != "reference"]
    cell_types = [ct for ct in cell_types if ct in counts.index]

    palette = sns.color_palette("tab10", n_colors=len(cell_types))

    fig, ax = plt.subplots(figsize=(6, 5))

    legend_handles = []
    for ct, color in zip(cell_types, palette):
        v = df.loc[df["cell_type"] == ct, "p_value"].to_numpy()
        v = v[np.isfinite(v)]
        if len(v) == 0:
            continue
        n = len(v)
        n_q = min(QQ_N_QUANTILES, n)
        qs = np.linspace(0, 1, n_q + 2)[1:-1]
        obs = np.quantile(v, qs)

        ct_display = ref_display if ct == "reference" else ct
        fpr = float(np.mean(v < 0.05))
        line, = ax.plot(qs, obs, color=color, alpha=0.5, lw=1.2,
                        label=f"{ct_display}  FPR={fpr:.3f}")
        legend_handles.append(line)

    ax.plot([0, 1], [0, 1], color="black", linestyle="--", lw=0.8)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xlabel("Expected (Uniform)")
    ax.set_ylabel("Observed p-value")

    cond_pretty = "+reporter" if args.condition == "reporter" else "-reporter (deflated)"
    test_pretty = "t-test" if args.test == "ttest" else "MWU"
    ax.set_title(f"{args.dataset} -- {test_pretty}, {cond_pretty}  "
                 f"({len(legend_handles)} cell types)")

    ax.legend(handles=legend_handles, loc="center left",
              bbox_to_anchor=(1.02, 0.5), fontsize=7, frameon=False)

    out = _output_path(args.dataset, args.test, args.condition)
    fig.savefig(out, format="svg", bbox_inches="tight")
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
