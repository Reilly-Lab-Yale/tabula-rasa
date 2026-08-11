#!/usr/bin/env python
"""Poisson is ruled out: observed dispersion vs its simulated null.

Draws the result of overdispersion.py, which is the step that motivates an
overdispersed count family at all. For each cell type the Pearson dispersion
of a Poisson fit conditioning on CRE identity, phi = X2/(n-p), is plotted
against the envelope of phi obtained by refitting data simulated from that
same Poisson fit. Under Poisson phi is 1; the simulated null pins down how
far it can stray by chance at these sample sizes and mean levels, where the
asymptotic chi-square reference is not calibrated.

    python analyses/model_selection/plot_overdispersion.py
"""
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Editable <text> so the manuscript sync pipeline's transforms apply, matching
# the other figures.
plt.rcParams["svg.fonttype"] = "none"

BASE = pathlib.Path(__file__).resolve().parent
OUT = BASE / "output"
TSV = BASE / "overdispersion.tsv"

BLUE, ORANGE = "#0072b2", "#d55e00"     # Okabe-Ito (house palette)
INK, MUTED = "#1a1a1a", "#6b6b6b"


def main():
    d = pd.read_csv(TSV, sep="\t")
    assert d.dataset.nunique() == 1, f"expected one dataset, got {sorted(d.dataset.unique())}"
    assert d.cell_type.is_unique, "cell types are not unique; one row per cell type expected"
    assert d.phi.notna().all(), "missing phi"
    # The whole point of the figure: observed dispersion sits outside the null.
    assert (d.phi > d.phi_null_max).all(), (
        "some cell type's observed phi is within its simulated null -- "
        "the figure would be claiming something the data does not support")
    assert (d.aic_pois > d.aic_nb).all(), "Poisson beats NB by AIC somewhere"

    d = d.sort_values("phi")
    y = np.arange(len(d))

    fig, ax = plt.subplots(figsize=(6.2, 3.0))

    # Null envelope: the full spread of phi across all simulations, pooled
    # over cell types. It is narrow enough that per-cell-type bands would
    # overplot into a single stripe anyway.
    lo = float((d.phi_null_mean - 3 * d.phi_null_sd).min())
    hi = float(d.phi_null_max.max())
    ax.axvspan(lo, hi, color=MUTED, alpha=0.18, zorder=1)
    ax.axvline(1.0, color=MUTED, lw=1.0, ls="--", zorder=2)

    ax.hlines(y, hi, d.phi, color=BLUE, lw=1.2, alpha=0.55, zorder=3)
    ax.scatter(d.phi, y, s=34, color=BLUE, zorder=4)

    ax.set_yticks(y)
    ax.set_yticklabels(d.cell_type, fontsize=8)
    ax.set_xscale("log")
    ax.set_xlim(0.6, float(d.phi.max()) * 1.9)
    ax.set_xlabel("Pearson dispersion $\\phi$ of the Poisson fit "
                  "(1 = Poisson)", fontsize=9, color=INK)

    for n, yi, xi in zip(d.n, y, d.phi):
        ax.text(xi * 1.12, yi, f"n={n:,}", va="center", fontsize=7, color=MUTED)

    ax.grid(True, axis="x", color="#e6e6e6", lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#cccccc")
    ax.tick_params(colors=MUTED, labelsize=8, length=0)
    # Annotated rather than put in a legend: with ten rows there is no corner
    # of the axes free of data.
    ax.annotate(f"Poisson null\n(max $\\phi$ = {hi:.2f})",
                xy=(hi, y[-1]), xytext=(1.9, y[-1] + 0.55),
                fontsize=7.5, color=MUTED, va="center",
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8,
                                shrinkA=0, shrinkB=2))
    ax.set_ylim(-0.8, len(d) - 0.1)

    print(f"phi: {d.phi.min():.1f} to {d.phi.max():.1f} (median {d.phi.median():.1f})")
    print(f"null: mean {d.phi_null_mean.mean():.2f}, max {hi:.2f}")

    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(OUT / f"overdispersion_shendure.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT/'overdispersion_shendure.svg'} and .png")


if __name__ == "__main__":
    main()
