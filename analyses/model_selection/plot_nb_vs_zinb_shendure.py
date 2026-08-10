#!/usr/bin/env python
"""NB vs ZINB for Lalanne et al., in both stratification directions.

The worked example behind the manuscript's statement that the counts are
overdispersed but not zero-inflated. lrt_nb_vs_zinb.py compares the two
families across datasets; this shows one dataset in detail, fit by fit, in
both directions:

  by cell type -- one model per cell type, CRE identity as predictor (10 fits)
  by CRE       -- one model per CRE, cell type as predictor      (208 fits)

Plots dAIC = AIC_ZINB - AIC_NB, so positive means NB is preferred. The two
panels differ in form because the counts differ: ten fits are shown
individually, two hundred as a distribution.

    python analyses/model_selection/plot_nb_vs_zinb_shendure.py
"""
import pickle
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Editable <text> so the manuscript sync pipeline's transforms apply, matching
# the other figures.
plt.rcParams["svg.fonttype"] = "none"

BASE = Path(__file__).parent
OUT = BASE / "output"
DATA = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new/shendure")
NB, ZINB = DATA / "shendure_obs_nb_phantom", DATA / "shendure_obs_zinb_phantom"

BLUE, ORANGE = "#0072b2", "#d55e00"     # Okabe-Ito (house palette)
INK, MUTED = "#1a1a1a", "#6b6b6b"


class _Placeholder:
    def __setstate__(self, state):
        self.__dict__.update(state if isinstance(state, dict) else {})


class TolerantUnpickler(pickle.Unpickler):
    """Read saved fits without importing the full scMPRAforge stack.

    Only the recorded scalars are used, so the classes need not be
    reconstructed faithfully.
    """

    def find_class(self, module, name):
        if module.startswith(("scMPRAforge", "tensorzinb")):
            return type(name, (_Placeholder,), {})
        return super().find_class(module, name)


def load(path):
    with open(path, "rb") as f:
        return TolerantUnpickler(f).load()


def delta_aic(family):
    """dAIC per fit, positive where NB beats ZINB."""
    nb, zi = load(NB / family), load(ZINB / family)
    keys = sorted(set(nb.model) & set(zi.model))
    assert keys, f"no shared fits in {family}"
    out = {}
    for k in keys:
        a_nb, a_zi = nb.model[k], zi.model[k]
        # df must differ by the zero-inflation parameters, or the two fits are
        # not the comparison we think they are.
        assert a_zi["df_model_total"] > a_nb["df_model_total"], (
            f"{family}/{k}: ZINB has no extra parameters over NB")
        out[k] = float(a_zi["aic_total"]) - float(a_nb["aic_total"])
    dropped = (set(nb.model) | set(zi.model)) - set(keys)
    if dropped:
        print(f"  {family}: {len(dropped)} fits present in only one family, skipped")
    return out


def main():
    by_ct = delta_aic("by_cell_type.pkl")
    by_cre = delta_aic("by_cre.pkl")
    print(f"by_cell_type: {len(by_ct)} fits, NB preferred in "
          f"{sum(v > 0 for v in by_ct.values())}")
    print(f"by_cre      : {len(by_cre)} fits, NB preferred in "
          f"{sum(v > 0 for v in by_cre.values())}")

    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(7.4, 3.0))

    # Ten fits: show each one.
    names = sorted(by_ct, key=lambda k: by_ct[k])
    vals = [by_ct[k] for k in names]
    colours = [BLUE if v > 0 else ORANGE for v in vals]
    ax_l.barh(range(len(names)), vals, color=colours, height=0.68, zorder=3)
    ax_l.set_yticks(range(len(names)))
    ax_l.set_yticklabels([n[:22] for n in names], fontsize=8)
    ax_l.set_xlabel(r"$\Delta$AIC (ZINB $-$ NB)", fontsize=9, color=INK)
    ax_l.set_title("by cell type", fontsize=10, color=INK)

    # Two hundred fits: show the distribution. A handful of CREs favour ZINB
    # by a very large margin, so a linear axis collapses the bulk into a spike;
    # symlog keeps both the mass near zero and the tail readable, and the tail
    # is the interesting part -- those are the CREs that really are inflated.
    v = np.array(list(by_cre.values()))
    lim = float(np.abs(v).max()) * 1.3
    bins = np.concatenate([
        -np.geomspace(lim, 10, 18), np.linspace(-10, 10, 9), np.geomspace(10, lim, 12)])
    ax_r.hist(v, bins=bins, color=BLUE, zorder=3)
    ax_r.set_xscale("symlog", linthresh=10)
    ax_r.set_xlabel(r"$\Delta$AIC (ZINB $-$ NB), symlog", fontsize=9, color=INK)
    ax_r.set_ylabel("CREs", fontsize=9, color=INK)
    ax_r.set_title(f"by CRE ({len(v)} fits)", fontsize=10, color=INK)

    for ax in (ax_l, ax_r):
        ax.axvline(0, color=MUTED, lw=1.0, ls="--", zorder=4)
        ax.grid(True, axis="x", color="#e6e6e6", lw=0.8, zorder=0)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color("#cccccc")
        ax.tick_params(colors=MUTED, labelsize=8, length=0)

    frac = sum(x > 0 for x in v) / len(v)
    n_strong_zinb = int((v < -50).sum())
    ax_r.text(0.5, 1.14, f"NB preferred in {frac:.0%}; {n_strong_zinb} CREs favour "
              r"ZINB by $>$50", transform=ax_r.transAxes, ha="center", va="bottom",
              fontsize=8, color=INK)
    fig.text(0.005, 0.01, "positive favours NB", fontsize=7.5, color=MUTED)

    fig.tight_layout(rect=(0, 0.03, 1, 1))
    OUT.mkdir(exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(OUT / f"nb_vs_zinb_shendure.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT/'nb_vs_zinb_shendure.svg'} and .png")


if __name__ == "__main__":
    main()
