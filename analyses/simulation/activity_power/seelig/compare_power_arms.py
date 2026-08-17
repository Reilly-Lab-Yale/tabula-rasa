#!/usr/bin/env python3
"""
Compare Seelig activity-power arms.

Reports, per cell type, power at a fixed set of fold changes and the fold
change at which the curve first crosses 80% power. Any number of parquets
produced by seelig_power_ttest_all_cell_types.py (or the surviving
seelig_power_df_deflated.parquet) can be passed; each is summarised on the
same grid so arms and runs are directly comparable.

Usage:
    python compare_power_arms.py <label>=<parquet> [<label>=<parquet> ...]
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Fold changes the manuscript reports power at.
REPORT_FCS = [1.2, 1.5, 2.0, 2.5, 3.0]
# Half-width of the window averaged around each reported fold change. Narrow
# enough that the curve is locally flat, wide enough to hold many draws.
FC_WINDOW = 0.05
TARGET_POWER = 0.8
# Bin width for the monotone-smoothed curve used to locate the 80% crossing.
CROSSING_BIN = 0.02


def load(path):
    df = pd.read_parquet(path)
    n_in = len(df)
    assert {"reject_null", "fc", "cell_type"} <= set(df.columns), (
        f"{path}: missing expected columns, has {sorted(df.columns)}"
    )
    df = df[["cell_type", "fc", "reject_null"]].copy()
    df["reject_null"] = df["reject_null"].astype(float)
    assert len(df) == n_in, f"{path}: column select changed rows, {n_in} -> {len(df)}"
    assert df["fc"].notna().all(), f"{path}: {df['fc'].isna().sum()} null fold changes"
    assert (df["fc"] > 0).all(), (
        f"{path}: {(df['fc'] <= 0).sum()} non-positive fold changes"
    )
    assert df["reject_null"].isin([0.0, 1.0]).all(), (
        f"{path}: reject_null is not binary, "
        f"values {sorted(df['reject_null'].unique())[:5]}"
    )
    return df


def power_at(df, fc):
    """Rejection rate in a narrow window around fc."""
    sel = df[(df["fc"] >= fc - FC_WINDOW) & (df["fc"] <= fc + FC_WINDOW)]
    if len(sel) == 0:
        return np.nan, 0
    return sel["reject_null"].mean(), len(sel)


def fc_at_power(df, target=TARGET_POWER):
    """First fold change above 1 whose binned power reaches `target`.

    The tests are two-sided, so the power curve is U-shaped: a CRE far weaker
    than the reference is as detectable as one far stronger. Only the FC>1
    branch is the quantity of interest, and searching the whole range would
    return a crossing down in the FC<1 tail instead.

    Within that branch the curve is made non-decreasing before the search:
    power is monotone in effect size by construction, so a dip is sampling
    noise, and without the correction a single noisy bin can trip the crossing
    early.
    """
    df = df[df["fc"] >= 1.0]
    if df.empty:
        return np.nan
    lo, hi = df["fc"].min(), df["fc"].max()
    edges = np.arange(lo, hi + CROSSING_BIN, CROSSING_BIN)
    binned = (
        df.groupby(pd.cut(df["fc"], bins=edges), observed=True)["reject_null"]
        .agg(["mean", "size"])
        .reset_index()
    )
    # Drop thin bins; their means are too noisy to define a crossing.
    binned = binned[binned["size"] >= 30]
    if binned.empty:
        return np.nan
    centers = binned["fc"].apply(lambda x: x.mid).to_numpy(dtype=float)
    power = np.maximum.accumulate(binned["mean"].to_numpy(dtype=float))
    above = np.flatnonzero(power >= target)
    if above.size == 0:
        return np.nan
    i = above[0]
    if i == 0:
        return float(centers[0])
    # Linear interpolation between the bracketing bin centres.
    x0, x1 = centers[i - 1], centers[i]
    y0, y1 = power[i - 1], power[i]
    if y1 == y0:
        return float(x1)
    return float(x0 + (target - y0) * (x1 - x0) / (y1 - y0))


def main(specs):
    frames = {}
    for spec in specs:
        assert "=" in spec, f"expected <label>=<path>, got {spec!r}"
        label, path = spec.split("=", 1)
        p = Path(path)
        assert p.exists(), f"{label}: no such file {p}"
        frames[label] = load(p)
        print(f"{label}: {len(frames[label]):,} rows from {p}")
    print()

    cell_types = sorted(set().union(*(set(d["cell_type"]) for d in frames.values())))

    header = ["cell_type", "arm"] + [f"FC {f}" for f in REPORT_FCS] + ["FC@80%"]
    print("  ".join(f"{h:>12}" for h in header))
    print("-" * (14 * len(header)))

    summary = {}
    for ct in cell_types:
        for label, df in frames.items():
            sub = df[df["cell_type"] == ct]
            if sub.empty:
                continue
            cells = []
            for fc in REPORT_FCS:
                pw, n = power_at(sub, fc)
                cells.append(f"{pw:.3f}" if n else "n/a")
            crossing = fc_at_power(sub)
            summary[(ct, label)] = (
                [power_at(sub, fc)[0] for fc in REPORT_FCS],
                crossing,
            )
            row = [ct, label] + cells + [
                f"{crossing:.2f}" if np.isfinite(crossing) else "never"
            ]
            print("  ".join(f"{c:>12}" for c in row))
        print()

    # Pairwise deltas against the first arm listed, which is the reference the
    # others are being checked against.
    labels = list(frames)
    if len(labels) > 1:
        ref = labels[0]
        print(f"Deltas vs {ref} (positive = higher power than {ref}):")
        for ct in cell_types:
            for label in labels[1:]:
                if (ct, ref) not in summary or (ct, label) not in summary:
                    continue
                a_fcs, a_cross = summary[(ct, ref)]
                b_fcs, b_cross = summary[(ct, label)]
                d = [f"{(b - a):+.3f}" for a, b in zip(a_fcs, b_fcs)]
                dc = (
                    f"{(b_cross - a_cross):+.2f}"
                    if np.isfinite(a_cross) and np.isfinite(b_cross)
                    else "n/a"
                )
                print(
                    "  ".join(
                        f"{c:>12}" for c in [ct, label] + d + [dc]
                    )
                )
        print()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1:])
