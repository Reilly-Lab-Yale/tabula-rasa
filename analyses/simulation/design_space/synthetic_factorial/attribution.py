#!/usr/bin/env python3
"""Pairwise per-axis power attribution between empirical anchors.

Reads output/samples_power_pooled.parquet (6220 samples: full + topup +
topup3 + union, pooled). The union LHS dominates (5000/6220) so the
pooled set retains near-LHS-grade axis independence (cross-axis |r| <=
0.094, vs 0.44 in the pre-union 'combined' patchwork), while the prior
sweeps contribute extra density near the empirical anchors.

For each (A, B) anchor pair, fits a 1-D LOESS marginal of
power_auc_1to3 vs each axis, evaluates the smoother at the A and B
empirical coordinates, and decomposes the A-B power difference into
per-axis contributions: delta = A_pred - B_pred.

Outputs (per pair, tagged <a_short>_vs_<b_short>):
- attribution_bar_<tag>.svg          per-axis delta bar chart
- attribution_marginals_<tag>.svg    marginals layout with anchor points
                                      + delta arrows per panel
- stdout: per-axis delta, sum, KNN-validated power at each anchor

Notes:
- Empirical anchors are clamped to the union LHS bracket before
  evaluating LOESS; out-of-range values are hatched in the bar chart.
  The union box is wide enough that most empirical values fall inside.
- Marginal-sum is approximate -- it overcounts to the extent that axes
  are correlated. Pooled |r| <= 0.094 keeps this small.
- KNN check: 20 nearest LHS samples (in standardized log-feature space)
  give a non-parametric power prediction at each anchor for sanity check.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from synthetic_factorial import (
    AXIS_NAMES, AXIS_BOUNDS, UNION_AXIS_BOUNDS, EMPIRICAL, EMPIRICAL_COLORS,
    _loess_band, OUT,
)

# Union LHS already spans the full box; use it directly as the clamp
# bracket so anchors that fall inside the union range don't get flagged
# as out-of-range.
COMBINED_BOUNDS = UNION_AXIS_BOUNDS

METRIC = "power_auc_1to3"

# Short-name tags for filenames. One per EMPIRICAL key.
SHORT = {
    "shendure-Pluripotent": "shendure",
    "cohen-Rod":            "cohen",
    "takeshi-HepG2":        "takeshi",
}

# Pairs to compute. (a, b) gives delta = a_pred - b_pred; positive means
# anchor a has a power advantage on that axis.
PAIRS = [
    ("cohen-Rod",       "shendure-Pluripotent"),
    ("takeshi-HepG2",   "shendure-Pluripotent"),
    ("takeshi-HepG2",   "cohen-Rod"),
]


def _x_to_plot(axis: str, x):
    return np.log10(x) if AXIS_BOUNDS[axis][2] else x


def _clamp(axis: str, x):
    lo, hi, _ = COMBINED_BOUNDS[axis]
    c = max(lo, min(hi, x))
    return c, (c != x)


def _loess_predict_at(x_data, y_data, x_query, frac=0.4):
    from statsmodels.nonparametric.smoothers_lowess import lowess
    valid = np.isfinite(x_data) & np.isfinite(y_data)
    return float(lowess(y_data[valid], x_data[valid],
                        frac=frac, xvals=np.array([x_query]),
                        return_sorted=False)[0])


def _knn_predict(df: pd.DataFrame, anchor: dict, k: int = 20) -> float:
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


def run_pair(df: pd.DataFrame, anchor_a: str, anchor_b: str):
    """Compute and plot attribution for one (a, b) anchor pair.

    Sign convention: delta = a_pred - b_pred; positive bars favour `a`.
    """
    a_short = SHORT[anchor_a]
    b_short = SHORT[anchor_b]
    tag = f"{a_short}_vs_{b_short}"
    anchors = [anchor_a, anchor_b]
    print(f"\n{'='*60}\nPair: {anchor_a}  vs  {anchor_b}  (tag={tag})\n{'='*60}")

    deltas = {}
    predictions = {name: {} for name in anchors}
    clamp_flags = {name: {} for name in anchors}

    for axis in AXIS_NAMES:
        x = df[axis].values.astype(float)
        y = df[METRIC].values.astype(float)
        x_plot = _x_to_plot(axis, x)
        for name in anchors:
            v = EMPIRICAL[name].get(axis)
            if v is None:
                predictions[name][axis] = np.nan
                clamp_flags[name][axis] = False
                continue
            v_c, was_clamped = _clamp(axis, v)
            x_query = _x_to_plot(axis, v_c)
            predictions[name][axis] = _loess_predict_at(x_plot, y, x_query)
            clamp_flags[name][axis] = was_clamped
        deltas[axis] = predictions[anchor_a][axis] - predictions[anchor_b][axis]

    sorted_axes = sorted(deltas, key=lambda a: -abs(deltas[a]))

    # -- text summary --
    print("Per-axis predicted power (LOESS, marginal):")
    print(f"  {'axis':<20s} {a_short:>10s} {b_short:>10s} {'delta':>10s}  flag")
    for axis in sorted_axes:
        pa = predictions[anchor_a][axis]
        pb = predictions[anchor_b][axis]
        flag = "CLAMPED" if (clamp_flags[anchor_a][axis] or
                              clamp_flags[anchor_b][axis]) else ""
        print(f"  {axis:<20s} {pa:>10.3f} {pb:>10.3f} {deltas[axis]:>+10.3f}  {flag}")
    print(f"\n  sum of deltas:  {sum(deltas.values()):+.3f}")
    print(f"  median LHS power_auc_1to3:  {df[METRIC].median():.3f}")

    knn_a = _knn_predict(df, EMPIRICAL[anchor_a])
    knn_b = _knn_predict(df, EMPIRICAL[anchor_b])
    print(f"\nKNN (k=20 in standardized log space) sanity check:")
    print(f"  {a_short} neighborhood: {knn_a:.3f}")
    print(f"  {b_short} neighborhood: {knn_b:.3f}")
    print(f"  KNN delta:             {knn_a - knn_b:+.3f}")
    print(f"  (compare to sum-of-LOESS-deltas {sum(deltas.values()):+.3f}; "
          f"big mismatch implies marginal decomposition is missing interactions)")

    clamped_axes = [a for a in AXIS_NAMES
                    if clamp_flags[anchor_a][a] or clamp_flags[anchor_b][a]]
    if clamped_axes:
        print(f"\nClamped axes (empirical outside LHS bracket; attribution approximate):")
        for axis in clamped_axes:
            for name in anchors:
                if not clamp_flags[name][axis]:
                    continue
                raw = EMPIRICAL[name][axis]
                clamped_val, _ = _clamp(axis, raw)
                print(f"  {name} {axis}: empirical={raw} -> LHS edge={clamped_val}")

    # -- bar chart --
    values = [deltas[a] for a in sorted_axes]
    color_a = EMPIRICAL_COLORS[anchor_a]
    color_b = EMPIRICAL_COLORS[anchor_b]
    colors = [color_a if v > 0 else color_b for v in values]
    is_clamped = [clamp_flags[anchor_a][a] or clamp_flags[anchor_b][a]
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
    ax.set_xlabel(f"delta {METRIC} attributed to axis ({a_short} - {b_short})")
    ax.set_title(f"Per-axis power attribution: {a_short} vs {b_short}\n"
                 f"positive = {a_short} advantage; hatched = clamped to LHS bound")
    for i, (v, hatched) in enumerate(zip(values, is_clamped)):
        ann = f"{v:+.3f}" + (" *" if hatched else "")
        ha = "left" if v >= 0 else "right"
        offset = 0.005 if v >= 0 else -0.005
        ax.text(v + offset, i, ann, va="center", ha=ha, fontsize=9)
    plt.tight_layout()
    out_bar = OUT / f"attribution_bar_{tag}.svg"
    fig.savefig(out_bar, format="svg", bbox_inches="tight")
    print(f"\nSaved: {out_bar}")
    plt.close(fig)

    # -- annotated marginals --
    fig, axes = plt.subplots(3, 3, figsize=(13, 11))
    axes_flat = axes.ravel()
    ec = {anchor_a: color_a, anchor_b: color_b}
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
        for name in anchors:
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
            label = SHORT[name] + ("*" if was_c else "")
            ax.text(xv, 1.02, label, fontsize=7, ha="center",
                    color=color, transform=ax.get_xaxis_transform())
        if set(anchors) <= ax_xs.keys():
            ax.annotate(
                "",
                xy=(ax_xs[anchor_a], ax_ys[anchor_a]),
                xytext=(ax_xs[anchor_b], ax_ys[anchor_b]),
                arrowprops=dict(arrowstyle="->", color="black", lw=1.3))
            mid_x = 0.5 * (ax_xs[anchor_a] + ax_xs[anchor_b])
            mid_y = 0.5 * (ax_ys[anchor_a] + ax_ys[anchor_b])
            ax.text(mid_x, mid_y + 0.07, f"D={deltas[axis]:+.2f}",
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
        f"{a_short} vs {b_short} power attribution ({METRIC})\n"
        f"panels ordered by |delta| descending; * = empirical clamped to LHS bracket"
    )
    plt.tight_layout()
    out_marg = OUT / f"attribution_marginals_{tag}.svg"
    fig.savefig(out_marg, format="svg", bbox_inches="tight")
    print(f"Saved: {out_marg}")
    plt.close(fig)


def main():
    df = pd.read_parquet(OUT / "samples_power_pooled.parquet")
    print(f"loaded {len(df)} samples; metric={METRIC}; pairs={len(PAIRS)}")
    for a, b in PAIRS:
        run_pair(df, a, b)


if __name__ == "__main__":
    main()
