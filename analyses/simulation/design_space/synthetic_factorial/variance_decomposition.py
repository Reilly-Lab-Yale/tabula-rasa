#!/usr/bin/env python3
"""How much of the power surface is additive, and how much is interaction?

The attribution bars in the design-space figure decompose the empirical power
gap one axis at a time, which is only honest if the axes act close to
independently. This quantifies that: nested least-squares fits of activity
power against the seven design axes, with and without cross terms.

Every axis is drawn log-uniformly, so all fits are on log10 axes,
standardised. Reported models:

    A   seven linear main effects
    B   A plus all 21 pairwise products
    C   cubic main effects, no cross terms
    D   cubic main effects plus pairwise products

A vs B is the question the attribution rests on. C is included because the
sweep's dominant departure from linearity is curvature within an axis --
dynamic range saturates above roughly 3x -- and fitting cross terms on top of
linear main effects lets that curvature leak into the interaction estimates.
Cubic is a coarse probe of that curvature, adequate for showing it is the
larger term but not a publication-grade smoother; natural splines would be
the honest choice if a curvature number is ever quoted directly.

    python variance_decomposition.py

Writes output/variance_decomposition.tsv.
"""
import itertools
import pathlib

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "output"
SAMPLES = OUT / "samples_power_union.parquet"
METRIC = "power_auc_1to3"
AXES = ["n_cells", "n_cres", "bcs_per_cre", "moi", "lib_alpha_nb", "minP",
        "activity_max_mult"]


def r2(X, y):
    A = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ beta
    return 1.0 - (resid ** 2).sum() / ((y - y.mean()) ** 2).sum(), beta


def main():
    assert SAMPLES.is_file(), f"missing {SAMPLES}; run the plot phase first"
    d = pd.read_parquet(SAMPLES)
    assert METRIC in d.columns, f"{METRIC} not in {list(d.columns)}"
    y = d[METRIC].values.astype(float)
    assert np.isfinite(y).all(), "non-finite power values"

    L = np.column_stack([np.log10(d[a].values.astype(float)) for a in AXES])
    assert np.isfinite(L).all(), "non-positive axis value; all axes are log-drawn"
    L = (L - L.mean(0)) / L.std(0)

    pairs = list(itertools.combinations(range(len(AXES)), 2))
    P = np.column_stack([L[:, i] * L[:, j] for i, j in pairs])
    S = np.column_stack([L, L ** 2, L ** 3])

    rA, _ = r2(L, y)
    rB, bB = r2(np.column_stack([L, P]), y)
    rC, _ = r2(S, y)
    rD, _ = r2(np.column_stack([S, P]), y)

    rows = [
        ("A", "linear main effects", len(AXES), rA, np.nan),
        ("B", "A + pairwise interactions", len(AXES) + len(pairs), rB, rB - rA),
        ("C", "cubic main effects", 3 * len(AXES), rC, rC - rA),
        ("D", "cubic main effects + interactions", 3 * len(AXES) + len(pairs),
         rD, rD - rC),
    ]
    tab = pd.DataFrame(rows, columns=["model", "terms", "n_params", "r2",
                                      "gain_over_reference"])
    print(f"{len(d)} samples, metric {METRIC}\n")
    print(tab.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    # Interactions are what the additive decomposition assumes away, so name
    # the largest few rather than only reporting that the total is small.
    inter = bB[1 + len(AXES):]
    order = np.argsort(-np.abs(inter))[:8]
    print("\nlargest pairwise interaction coefficients (model B, standardised):")
    for k in order:
        i, j = pairs[k]
        print(f"  {AXES[i]:18s} x {AXES[j]:18s} {inter[k]:+.4f}")

    assert rB - rA < rC - rA, (
        "interactions now explain more than main-effect curvature; the "
        "additive reading of the attribution bars needs revisiting")

    OUT.mkdir(exist_ok=True)
    out = OUT / "variance_decomposition.tsv"
    tab.to_csv(out, sep="\t", index=False)
    inter_tab = pd.DataFrame(
        [(AXES[i], AXES[j], inter[k]) for k, (i, j) in enumerate(pairs)],
        columns=["axis_a", "axis_b", "beta"],
    ).sort_values("beta", key=np.abs, ascending=False)
    inter_out = OUT / "pairwise_interactions.tsv"
    inter_tab.to_csv(inter_out, sep="\t", index=False)
    print(f"\nwrote {out} and {inter_out}")


if __name__ == "__main__":
    main()
