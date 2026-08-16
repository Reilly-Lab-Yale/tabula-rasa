#!/usr/bin/env python
"""Fitted means from the consider-missing fits, NB against ZINB.

Panels C and D of Fig. S1 contrast the two count families on one fit mode, the
way panel B contrasts the two reporter expansions on one family. This reads the
four orthos and caches what those panels need, so plotting needs neither the
cluster nor an unpickle. plot_cm_family_mu.py reads the result.

The mean cached here is mu as fitted. For ZINB that is the count-component
mean, conditional on the observation not being a structural zero, which is the
quantity that moves when the inflation parameter takes zeros out of the count
process. There is no per-element marginal to cache instead: zi is estimated per
replicate, not per element, so (1 - zi) * mu is not defined at the granularity
of these rows.

Runs on the cluster, where the orthos live.

    sbatch analyses/model_selection/run_cm_family_mu.sh
"""
import pickle
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa")
import scMPRAforge.core as scm  # noqa: F401,E402  (needed to unpickle)

DATA_ROOT = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new")
OUT_TSV = Path(__file__).resolve().parent / "cm_family_mu.tsv"

# Both datasets whose consider-missing fits exist in both families. Yin et al.
# is absent on purpose: consider-missing is its canonical mode, so the contrast
# these panels draw is not a counterfactual there.
ORTHOS = {
    ("cohen", "nb"): "cohen/cohen_cm_nb_phantom_20260401",
    ("cohen", "zinb"): "cohen/cohen_cm_zinb_phantom_20260401",
    ("shendure", "nb"): "shendure/shendure_cm_nb_phantom",
    ("shendure", "zinb"): "shendure/shendure_cm_zinb_phantom",
}


def mu_table(rel, family):
    """Per (cell type, element) fitted mean, from the by-cell-type models."""
    path = DATA_ROOT / rel / "by_cell_type_parameters.pkl"
    assert path.is_file(), f"missing ortho artifact: {path}"
    with open(path, "rb") as fh:
        params = pickle.load(fh)

    # The family is a property of the fit, not of the file name, so check it
    # rather than trust the directory: a ZINB fit carries an inflation table
    # and an NB fit carries None.
    has_zi = params.zi[params.keys[0]] is not None
    assert has_zi == (family == "zinb"), (
        f"{rel} is named {family} but its zi table is "
        f"{'present' if has_zi else 'absent'}")
    if has_zi:
        zi = pd.concat(params.zi[k]["zi"] for k in params.keys)
        print(f"    zi per replicate: {zi.min():.3f}-{zi.max():.3f} "
              f"over {len(zi)} (level, replicate) values")

    rows = []
    for level in params.keys:
        for regressor, mu in params.nb[level]["mu"].items():
            rows.append({"level": level, "regressor": regressor, "mu": float(mu)})
    d = pd.DataFrame(rows)
    assert len(d), f"{rel}: no fitted means"
    assert (d.mu >= 0).all(), f"{rel}: negative fitted mean"
    n_zero = int((d.mu == 0).sum())
    print(f"    {family:4s} {rel.split('/')[-1]:32s} {len(d):5d} fits, "
          f"mu {d.mu[d.mu > 0].min():.3g}-{d.mu.max():.3g}"
          + (f", {n_zero} exactly zero" if n_zero else ""))
    return d


def main():
    frames = []
    for (dataset, family), rel in ORTHOS.items():
        d = mu_table(rel, family)
        frames.append(d.assign(dataset=dataset, family=family))
    out = pd.concat(frames, ignore_index=True)

    # The two families are fit on the same design, so they must describe the
    # same set of (cell type, element) combinations. If they do not, the
    # panels would be overlaying histograms of different populations.
    for dataset in {d for d, _ in ORTHOS}:
        sets = {}
        for family in ("nb", "zinb"):
            sub = out[(out.dataset == dataset) & (out.family == family)]
            sets[family] = set(zip(sub.level, sub.regressor))
        assert sets["nb"] == sets["zinb"], (
            f"{dataset}: families cover different combinations "
            f"(nb only {len(sets['nb'] - sets['zinb'])}, "
            f"zinb only {len(sets['zinb'] - sets['nb'])})")
        print(f"  {dataset}: {len(sets['nb'])} combinations, both families")

    out = out[["dataset", "family", "level", "regressor", "mu"]]
    out.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"wrote {OUT_TSV} ({len(out)} rows)")


if __name__ == "__main__":
    main()
