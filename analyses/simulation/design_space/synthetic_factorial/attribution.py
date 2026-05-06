#!/usr/bin/env python3
"""Cohen-vs-shendure per-axis power attribution.

Reads output/samples_power_combined.parquet (1099 LHS samples covering the
synthetic factorial axis space, including the cohen-corner top-up). For each
of the 7 axes, fits a 1-D LOESS marginal of power_auc_1to3 vs that axis,
then evaluates the smoother at the cohen-Rod and shendure-Pluripotent
empirical values. Per-axis delta = cohen_pred - shendure_pred.

Outputs:
- output/attribution_bar.svg     bar chart of per-axis deltas
- output/attribution_marginals.svg  marginals layout with cohen/shendure
                                      points + delta arrows per panel
- stdout: per-axis delta, sum, KNN-validated power at each empirical point

Notes:
- activity_max_mult is the only axis where empirical values fall outside the
  LHS bracket [2, 8] (cohen=1.11, shendure=99). Both are clamped to LHS
  edge before evaluating LOESS. The bar is hatched and labelled accordingly.
- Marginal-sum is approximate -- it overcounts to the extent that axes are
  correlated in the LHS (LHS is space-filling so correlation is small).
- KNN check: 20 nearest LHS samples (in standardized log-feature space)
  give a non-parametric power prediction at each anchor for sanity checking.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Pull shared definitions from synthetic_factorial.py.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from synthetic_factorial import (
    AXIS_NAMES, AXIS_BOUNDS, TOPUP_AXIS_BOUNDS, EMPIRICAL, _loess_band, OUT,
)

# Combined bounds = union of original LHS and top-up brackets, since the
# combined parquet contains samples from both. Without this, cohen's
# bcs_per_cre (17244) and moi (149) get flagged as out-of-range when they
# are actually well inside the topup coverage.
COMBINED_BOUNDS = {}
for axis in AXIS_NAMES:
    a_lo, a_hi, log_scale = AXIS_BOUNDS[axis]
    t_lo, t_hi, _ = TOPUP_AXIS_BOUNDS[axis]
    COMBINED_BOUNDS[axis] = (min(a_lo, t_lo), max(a_hi, t_hi), log_scale)

METRIC = "power_auc_1to3"
EMP_NAMES = ["shendure-Pluripotent", "cohen-Rod"]


def _x_to_plot(axis: str, x):
    return np.log10(x) if AXIS_BOUNDS[axis][2] else x


def _clamp(axis: str, x):
    """Clamp x to the COMBINED bracket; return (clamped, was_clamped)."""
    lo, hi, _ = COMBINED_BOUNDS[axis]
    c = max(lo, min(hi, x))
    return c, (c != x)


def _loess_predict_at(x_data, y_data, x_query, frac=0.4):
    """1-D LOESS evaluated at one x. Caller provides x already log-transformed
    if axis is log-scale."""
    from statsmodels.nonparametric.smoothers_lowess import lowess
    valid = np.isfinite(x_data) & np.isfinite(y_data)
    return float(lowess(y_data[valid], x_data[valid],
                        frac=frac, xvals=np.array([x_query]),
                        return_sorted=False)[0])


def _knn_predict(df: pd.DataFrame, anchor: dict, k: int = 20) -> float:
    """Power_auc_1to3 averaged over k nearest LHS samples to `anchor`,
    in standardized log-feature space. Out-of-range anchor coords are
    clamped first."""
    feats = []
    qf = []
    for axis in AXIS_NAMES:
        v = anchor.get(axis)
        if v is None:
            continue
        v_c, _ = _clamp(axis, v)
        x = df[axis].values.astype(float)
        if AXIS_BOUNDS[axis][2]:
            x = np.log10(x)
            vp = np.log10(v_c)
        else:
            vp = v_c
        mu, sd = x.mean(), x.std()
        if sd == 0:
            continue
        feats.append((x - mu) / sd)
        qf.append((vp - mu) / sd)
    X = np.column_stack(feats)
    q = np.array(qf)
    d = np.sqrt(((X - q) ** 2).sum(axis=1))
    nearest = np.argsort(d)[:k]
    return float(df[METRIC].values[nearest].mean())


def main():
    df = pd.read_parquet(OUT / "samples_power_combined.parquet")
    print(f"loaded {len(df)} samples; metric={METRIC}\n")

    deltas = {}
    predictions = {name: {} for name in EMP_NAMES}
    clamp_flags = {name: {} for name in EMP_NAMES}

    for axis in AXIS_NAMES:
        x = df[axis].values.astype(float)
        y = df[METRIC].values.astype(float)
        x_plot = _x_to_plot(axis, x)
        for name in EMP_NAMES:
            v = EMPIRICAL[name].get(axis)
            if v is None:
                predictions[name][axis] = np.nan
                clamp_flags[name][axis] = False
                continue
            v_c, was_clamped = _clamp(axis, v)
            x_query = _x_to_plot(axis, v_c)
            predictions[name][axis] = _loess_predict_at(x_plot, y, x_query)
            clamp_flags[name][axis] = was_clamped
        deltas[axis] = predictions["cohen-Rod"][axis] - predictions["shendure-Pluripotent"][axis]

    sorted_axes = sorted(deltas, key=lambda a: -abs(deltas[a]))

    # -- text summary --
    print("Per-axis predicted power (LOESS, marginal):")
    print(f"  {'axis':<20s} {'shendure':>10s} {'cohen':>10s} {'delta':>10s}  flag")
    for axis in sorted_axes:
        ps = predictions["shendure-Pluripotent"][axis]
        pc = predictions["cohen-Rod"][axis]
        flag = "CLAMPED" if (clamp_flags["shendure-Pluripotent"][axis] or
                              clamp_flags["cohen-Rod"][axis]) else ""
        print(f"  {axis:<20s} {ps:>10.3f} {pc:>10.3f} {deltas[axis]:>+10.3f}  {flag}")
    print(f"\n  sum of deltas:  {sum(deltas.values()):+.3f}")
    print(f"  median LHS power_auc_1to3:  {df[METRIC].median():.3f}")

    knn_s = _knn_predict(df, EMPIRICAL["shendure-Pluripotent"])
    knn_c = _knn_predict(df, EMPIRICAL["cohen-Rod"])
    print(f"\nKNN (k=20 in standardized log space) sanity check:")
    print(f"  shendure neighborhood: {knn_s:.3f}")
    print(f"  cohen neighborhood:    {knn_c:.3f}")
    print(f"  KNN delta:             {knn_c - knn_s:+.3f}")
    print(f"  (compare to sum-of-LOESS-deltas {sum(deltas.values()):+.3f}; "
          f"big mismatch implies marginal decomposition is missing interactions)")

    clamped = [a for a in AXIS_NAMES
               if clamp_flags["cohen-Rod"][a] or clamp_flags["shendure-Pluripotent"][a]]
    if clamped:
        print(f"\nClamped axes (empirical outside LHS bracket; attribution approximate):")
        for axis in clamped:
            for name in EMP_NAMES:
                if not clamp_flags[name][axis]:
                    continue
                raw = EMPIRICAL[name][axis]
                clamped_val, _ = _clamp(axis, raw)
                print(f"  {name} {axis}: empirical={raw} -> LHS edge={clamped_val}")

    # -- bar chart --
    values = [deltas[a] for a in sorted_axes]
    colors = ["seagreen" if v > 0 else "indianred" for v in values]
    is_clamped = [clamp_flags["shendure-Pluripotent"][a] or clamp_flags["cohen-Rod"][a]
                  for a in sorted_axes]

    fig, ax = plt.subplots(figsize=(8, 5.5))
    bars = ax.barh(range(len(sorted_axes)), values, color=colors, edgecolor="black")
    for bar, hatched in zip(bars, is_clamped):
        if hatched:
            bar.set_hatch("///")
    ax.set_yticks(range(len(sorted_axes)))
    ax.set_yticklabels(sorted_axes)
    ax.invert_yaxis()
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel(f"delta {METRIC} attributed to axis (cohen - shendure)")
    ax.set_title("Per-axis power attribution: cohen vs shendure\n"
                 "positive = cohen advantage; hatched = clamped to LHS bound")
    # annotate each bar with its numeric value
    for i, (v, hatched) in enumerate(zip(values, is_clamped)):
        ann = f"{v:+.3f}" + (" *" if hatched else "")
        ha = "left" if v >= 0 else "right"
        offset = 0.005 if v >= 0 else -0.005
        ax.text(v + offset, i, ann, va="center", ha=ha, fontsize=9)
    plt.tight_layout()
    out_bar = OUT / "attribution_bar.svg"
    fig.savefig(out_bar, format="svg", bbox_inches="tight")
    print(f"\nSaved: {out_bar}")
    plt.close(fig)

    # -- annotated marginals --
    fig, axes = plt.subplots(3, 3, figsize=(13, 11))
    axes_flat = axes.ravel()
    ec = {"shendure-Pluripotent": "darkorange", "cohen-Rod": "purple"}
    for ax_i, axis in enumerate(sorted_axes):
        ax = axes_flat[ax_i]
        x = df[axis].values.astype(float)
        y = df[METRIC].values.astype(float)
        log_scale = AXIS_BOUNDS[axis][2]
        x_plot = np.log10(x) if log_scale else x
        valid = np.isfinite(x_plot) & np.isfinite(y)
        ax.scatter(x_plot[valid], y[valid], s=8, alpha=0.3, color="steelblue")
        try:
            xg, yhat, lo, hi = _loess_band(x_plot[valid], y[valid])
            ax.plot(xg, yhat, color="firebrick", lw=2)
            ax.fill_between(xg, lo, hi, color="firebrick", alpha=0.15)
        except Exception:
            pass

        ax_xs = {}
        ax_ys = {}
        cur_lo, cur_hi = ax.get_xlim()
        for name in EMP_NAMES:
            v = EMPIRICAL[name].get(axis)
            if v is None:
                continue
            v_c, was_c = _clamp(axis, v)
            xv = _x_to_plot(axis, v_c)
            cur_lo = min(cur_lo, xv); cur_hi = max(cur_hi, xv)
            yv = predictions[name][axis]
            color = ec[name]
            ax.scatter([xv], [yv], color=color, s=80, zorder=10,
                       edgecolor="black", linewidth=0.8)
            ax_xs[name] = xv
            ax_ys[name] = yv
            label = name.split("-")[0][:4] + ("*" if was_c else "")
            ax.text(xv, 1.02, label, fontsize=7, ha="center",
                    color=color, transform=ax.get_xaxis_transform())
        if {"shendure-Pluripotent", "cohen-Rod"} <= ax_xs.keys():
            ax.annotate(
                "",
                xy=(ax_xs["cohen-Rod"], ax_ys["cohen-Rod"]),
                xytext=(ax_xs["shendure-Pluripotent"], ax_ys["shendure-Pluripotent"]),
                arrowprops=dict(arrowstyle="->", color="black", lw=1.3))
            mid_x = 0.5 * (ax_xs["shendure-Pluripotent"] + ax_xs["cohen-Rod"])
            mid_y = 0.5 * (ax_ys["shendure-Pluripotent"] + ax_ys["cohen-Rod"])
            ax.text(mid_x, mid_y + 0.07, f"Δ={deltas[axis]:+.2f}",
                    fontsize=8, ha="center",
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                              edgecolor="black", lw=0.5))
        span = cur_hi - cur_lo
        if span > 0:
            ax.set_xlim(cur_lo - 0.03 * span, cur_hi + 0.03 * span)
        ax.set_xlabel(("log10 " if log_scale else "") + axis)
        ax.set_ylabel(METRIC)
        ax.set_ylim(0, 1)
    for j in range(len(sorted_axes), len(axes_flat)):
        axes_flat[j].axis("off")
    fig.suptitle(
        f"cohen vs shendure power attribution ({METRIC})\n"
        "panels ordered by |delta| descending; * = empirical clamped to LHS bracket"
    )
    plt.tight_layout()
    out_marg = OUT / "attribution_marginals.svg"
    fig.savefig(out_marg, format="svg", bbox_inches="tight")
    print(f"Saved: {out_marg}")
    plt.close(fig)


if __name__ == "__main__":
    main()
