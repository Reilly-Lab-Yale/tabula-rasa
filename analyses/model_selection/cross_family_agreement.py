#!/usr/bin/env python
"""Do the two stratification directions agree on the means they share?

The by-cell-type family (one model per cell type, CRE as predictor) and the
by-CRE family (one model per CRE, cell type as predictor) are both saturated
in the mean, so each (cell type, CRE) combination has a free mean in both.
Fitting them separately is a computational convenience; if that convenience
is harmless, the two families must return the same mu for the same
combination. This checks that directly.

The comparison needs no likelihood and no distributional assumption, which
is what makes it usable across fits whose zero expansions differ and whose
AICs are therefore not comparable.

Report relative deviation, not correlation: mu spans several orders of
magnitude, so Pearson r sits at 0.999999 even for a fit whose families
disagree by a factor of two. The percentiles of |mu_cre - mu_ct| / mu_ct are
what separate them.

Reads the saved parameter pickles directly -- no scMPRAforge import, no
cluster, no dask.

    python analyses/model_selection/cross_family_agreement.py
"""
import pickle
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["svg.fonttype"] = "none"

BASE = Path(__file__).resolve().parent
OUT = BASE / "output"
ROOT = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new")

BLUE, ORANGE, GREEN, MUTED = "#0072b2", "#d55e00", "#009e73", "#6b6b6b"
INK = "#1a1a1a"

CANONICAL = {
    "shendure_obs_nb_phantom",
    "cohen_obsingle_nb_phantom",
    "seelig_cm_moib_nb_phantom",
}
DATASET_COLOUR = {"shendure": BLUE, "cohen": ORANGE, "seelig": GREEN}
DATASET_LABEL = {"shendure": "Lalanne et al.", "cohen": "Zhao et al.",
                 "seelig": "Yin et al."}
# Fits shown in the left-hand scatter: the three canonical fits, plus the
# reporter-free Lalanne fit as the contrast. That one is not a broken fit --
# it is what the same experiment looks like when no transfection reporter was
# included -- so it shows how much the two directions drift apart when the
# zeros are unconditional rather than reporter-informed.
SCATTER = CANONICAL | {"shendure_cm_nb_phantom"}


def discover():
    """Every ortho under ROOT carrying both families' saved parameters."""
    found = []
    for ds_dir in sorted(ROOT.iterdir()):
        if not ds_dir.is_dir():
            continue
        for fit_dir in sorted(ds_dir.iterdir()):
            if ((fit_dir / "by_cell_type_parameters.pkl").is_file()
                    and (fit_dir / "by_cre_parameters.pkl").is_file()):
                found.append((ds_dir.name, fit_dir))
    assert found, f"no orthos with saved parameters under {ROOT}"
    return found


class _Placeholder:
    def __setstate__(self, state):
        self.__dict__.update(state if isinstance(state, dict) else {})


class TolerantUnpickler(pickle.Unpickler):
    """Read saved parameters without importing the scMPRAforge stack."""

    def find_class(self, module, name):
        if module.startswith(("scMPRAforge", "tensorzinb")):
            return type(name, (_Placeholder,), {})
        return super().find_class(module, name)


def _load(path):
    with open(path, "rb") as f:
        return TolerantUnpickler(f).load()


def _unwrap(v):
    return v.result() if hasattr(v, "result") else v


def _family_mean(params, names):
    """Marginal mean per (level, anti-level), as a long Series.

    For a ZINB fit, params.nb holds the mean of the NB component, not the
    mean of the distribution. The two families can split the same fitted
    mean differently between pi and mu -- on shendure_obs they differ by a
    factor of forty in pi -- so comparing the NB component alone reports a
    disagreement that the marginal means do not have. Zero inflation is
    regressed on replicate, so it carries no anti-level index and enters as
    a single per-model factor.
    """
    out = {}
    for level, v in params.nb.items():
        mu = _unwrap(v)["mu"]
        # NB-only fits still carry a zi dict, but its entries are None.
        zi = getattr(params, "zi", None) or {}
        zi_level = _unwrap(zi.get(level)) if zi.get(level) is not None else None
        if zi_level is not None:
            pi = float(zi_level["zi"].mean())
            assert 0.0 <= pi < 1.0, f"zero-inflation out of range: {pi}"
            mu = mu * (1.0 - pi)
        out[level] = mu
    return pd.concat(out, names=names)


def paired_mu(fit_dir):
    """Long table of the marginal mean from each family, by (cell_type, cre_id)."""
    ct = _load(fit_dir / "by_cell_type_parameters.pkl")
    cr = _load(fit_dir / "by_cre_parameters.pkl")
    assert ct.broken_on == "cell_type", f"unexpected split {ct.broken_on}"
    assert cr.broken_on == "cre_id", f"unexpected split {cr.broken_on}"

    by_ct = _family_mean(ct, ["cell_type", "cre_id"]).rename("mu_by_cell_type")
    by_cre = _family_mean(cr, ["cre_id", "cell_type"]).rename("mu_by_cre")
    by_cre = by_cre.reorder_levels(["cell_type", "cre_id"])

    outer = pd.concat([by_ct, by_cre], axis=1)
    merged = outer.dropna()
    # The two families are fit independently, so the minimum-observation
    # filter can drop a level in one direction and not the other. That is
    # expected raggedness rather than an error, but it must stay marginal:
    # if the families are comparing substantially different grids, the
    # agreement statistic below is not measuring what it claims to.
    n_unmatched = len(outer) - len(merged)
    frac = n_unmatched / max(len(outer), 1)
    assert frac < 0.05, (
        f"{fit_dir.name}: {n_unmatched} of {len(outer)} (cell type, CRE) "
        f"pairs ({frac:.1%}) appear in only one family")
    if n_unmatched:
        print(f"    {fit_dir.name}: {n_unmatched} pairs in one family only "
              f"({len(by_ct)} by-cell-type vs {len(by_cre)} by-CRE), "
              f"compared on the {len(merged)} shared")
    assert (merged > 0).all().all(), f"{fit_dir.name}: non-positive mu"
    return merged


def short(name):
    """Fit directory name as it reads in the manuscript's terms."""
    ds = name.split("_")[0]
    rest = name[len(ds) + 1:].replace("_phantom", "").replace("_", "+")
    return f"{DATASET_LABEL.get(ds, ds)} {rest}"


def main():
    rows = []
    for dataset, fit_dir in discover():
        m = paired_mu(fit_dir)
        rel = ((m.mu_by_cre - m.mu_by_cell_type).abs()
               / m.mu_by_cell_type).sort_values()
        rows.append(dict(name=fit_dir.name, dataset=dataset, pairs=m, rel=rel,
                         canonical=fit_dir.name in CANONICAL,
                         colour=DATASET_COLOUR.get(dataset, MUTED)))

    rows.sort(key=lambda r: r["rel"].median())
    print(f"{'fit':32s} {'n':>6s} {'median':>9s} {'p95':>9s} {'max':>9s}  canonical")
    for r in rows:
        print(f"{r['name']:32s} {len(r['rel']):6d} {r['rel'].median():9.2e} "
              f"{r['rel'].quantile(0.95):9.2e} {r['rel'].max():9.2e}"
              f"  {'*' if r['canonical'] else ''}")

    canon = [r for r in rows if r["canonical"]]
    assert canon, "no canonical fits loaded"
    worst = max(r["rel"].median() for r in canon)
    assert worst < 0.01, (
        f"a canonical fit's families disagree by a median of {worst:.1%}; "
        f"the stratified split is not behaving as a reparameterisation")

    fig, (ax_l, ax_r) = plt.subplots(
        1, 2, figsize=(8.6, 4.0), gridspec_kw={"width_ratios": [1, 1.25]})

    shown = [r for r in rows if r["name"] in SCATTER]
    for r in shown:
        colour = r["colour"] if r["canonical"] else MUTED
        ax_l.scatter(r["pairs"].mu_by_cell_type, r["pairs"].mu_by_cre, s=5,
                     alpha=0.35, color=colour, linewidths=0, zorder=3,
                     label=short(r["name"]))
    # Limits from the data: the axes are log-scaled, so the autoscaled
    # bounds can include zero before the scale is applied.
    allv = np.concatenate([r["pairs"].values.ravel() for r in shown])
    lims = [allv.min() * 0.6, allv.max() * 1.6]
    ax_l.set_xscale("log"); ax_l.set_yscale("log")
    ax_l.plot(lims, lims, ls="--", lw=1.0, color=INK, zorder=4)
    ax_l.set_xlim(lims); ax_l.set_ylim(lims)
    ax_l.set_xlabel(r"$\mu$, by-cell-type family", fontsize=9, color=INK)
    ax_l.set_ylabel(r"$\mu$, by-CRE family", fontsize=9, color=INK)
    ax_l.set_title("shared means, canonical fits", fontsize=10, color=INK)

    # One row per fit: median disagreement with the 5th-95th percentile span.
    y = np.arange(len(rows))
    for i, r in enumerate(rows):
        lo, hi = r["rel"].quantile(0.05), r["rel"].quantile(0.95)
        ax_r.plot([max(lo, 1e-9), hi], [i, i], color=r["colour"], lw=1.4,
                  alpha=0.55, zorder=3)
        ax_r.scatter([r["rel"].median()], [i], s=34, color=r["colour"],
                     zorder=4, marker="o" if r["canonical"] else "x")
    ax_r.set_yticks(y)
    ax_r.set_yticklabels([short(r["name"]) + (" *" if r["canonical"] else "")
                          for r in rows], fontsize=7.5)
    ax_r.set_xscale("log")
    ax_r.set_xlim(left=1e-5)
    ax_r.set_xlabel("relative disagreement between families\n"
                    "(median, 5th-95th percentile)", fontsize=9, color=INK)
    ax_r.set_title("all fits with saved parameters", fontsize=10, color=INK)
    ax_r.text(0.99, 0.02, "* canonical", transform=ax_r.transAxes, ha="right",
              va="bottom", fontsize=7.5, color=MUTED)

    for ax in (ax_l, ax_r):
        ax.grid(True, color="#e6e6e6", lw=0.8, zorder=0)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color("#cccccc")
        ax.tick_params(colors=MUTED, labelsize=8, length=0)

    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(OUT / f"cross_family_agreement.{ext}", dpi=200,
                    bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT/'cross_family_agreement.svg'} and .png")


if __name__ == "__main__":
    main()
