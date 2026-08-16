#!/usr/bin/env python
"""Fitted inflation parameters, and the zero mass each expansion introduced.

Fig. S1 panels B-E ask where the zeros went. A ZINB handed a pile of zeros it
cannot attribute has two places to put them: the inflation parameter pi, or the
count mean. Which one it uses is the identifiability question, and pi is the
direct read on it.

For each pair of expansions of one dataset this caches:

  - pi per (CRE, replicate), from the by-CRE models. The by-cell-type models
    carry the same parameter at a granularity of 8 values for Zhao et al.,
    which is not a distribution.
  - the fraction of the fit's observations that were zero, which is the mass
    pi would have to reach to account for the expansion on its own. Panels draw
    it as a reference line, so the gap between it and pi is what got pushed
    into the mean instead.

The zero fraction is reconstructed the way `_nb_versus_means` reconstructs its
QC denominator, by the same helpers the fits used, because most of these zeros
are phantom-compressed and never exist as rows:

  preexisting          zeros are rows in the table; count them
  per_delivery         _reporter_zero_counts(reporter_expansion="single")
  per_barcode          _reporter_zero_counts(reporter_expansion="coarse")
  all_combinations     _cm_group_totals
  all_combinations_moi _cm_group_totals with the MOI correction

Runs on the cluster, where the orthos live.

    sbatch analyses/model_selection/run_cm_family_pi.sh
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa")
import scMPRAforge.core as scm  # noqa: E402

DATA_ROOT = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new")
HERE = Path(__file__).resolve().parent
PI_TSV = HERE / "cm_family_pi.tsv"
ZERO_TSV = HERE / "cm_family_zerofrac.tsv"

SPLIT = "cell_type"   # denominator is a whole-fit property; fewer, larger groups

# (dataset, regime, relative path). regime names match the manuscript's gloss.
ORTHOS = [
    ("cohen", "fine", "cohen/cohen_obsingle_zinb_phantom"),
    ("cohen", "coarse", "cohen/cohen_obs_zinb_phantom_20260401"),
    ("cohen", "consider missing", "cohen/cohen_cm_zinb_phantom_20260401"),
    ("shendure", "observed", "shendure/shendure_obs_zinb"),
    ("shendure", "consider missing", "shendure/shendure_cm_zinb_phantom"),
    ("seelig", "consider missing + MOI", "seelig/seelig_cm_moib_zinb_phantom"),
    ("seelig", "consider missing", "seelig/seelig_cm_zinb_phantom"),
]


def _resolve(x):
    """A loaded ortho holds its parameters as futures; a saved one does not."""
    return x.result() if hasattr(x, "result") else x


def pi_rows(ortho, name):
    """One row per (CRE, replicate). by_cre models, so this is a distribution."""
    params = ortho.by_cre_parameters
    assert params is not None, f"{name}: no by_cre_parameters"
    rows = []
    for level in params.keys:
        zi = _resolve(params.zi[level])
        assert zi is not None, f"{name}: level {level} has no zi table -- not a ZINB fit"
        for rep_id, pi in zi["zi"].items():
            rows.append({"level": str(level), "rep_id": str(rep_id), "pi": float(pi)})
    d = pd.DataFrame(rows)
    assert len(d), f"{name}: no inflation parameters"
    assert d.pi.between(0, 1).all(), (
        f"{name}: pi outside [0,1] ({d.pi.min():.4g}-{d.pi.max():.4g}); "
        "params.zi is meant to be probability-scale, not logit")
    return d


def zero_fraction(ortho, name, expansion):
    """(n_zeros, n_total) over the whole fit, on the denominator it was fit on."""
    dat = ortho.training_data
    data = dat.get_data(include_missing=False)
    if hasattr(data["umis_mpra_bc"].dtype, "fill_value"):
        data["umis_mpra_bc"] = data["umis_mpra_bc"].astype("int64")
    obs = scm._to_pandas(data)
    n_nonzero = int((obs["umis_mpra_bc"] > 0).sum())

    if expansion == "preexisting":
        # The table already holds the zeros the fit saw; nothing was added.
        n_zeros = int((obs["umis_mpra_bc"] == 0).sum())
        return n_zeros, n_zeros + n_nonzero

    cell_map = obs[["rep_id", "cell_bc", "cell_type"]].drop_duplicates()
    mpra_map = obs[["rep_id", "mpra_bc", "cre_id"]].drop_duplicates()

    if expansion in ("per_delivery", "per_barcode"):
        reporter = getattr(dat, "_coarse_reporter", None)
        assert reporter is not None, f"{name}: reporter-keyed fit with no reporter table"
        anti = scm.anti_split(SPLIT)
        nz = obs[obs["umis_mpra_bc"] > 0][[SPLIT, anti, "rep_id", "umis_mpra_bc", "cell_bc"]].copy()
        for frame in (nz, cell_map, mpra_map):
            for col in frame.columns:
                if col != "umis_mpra_bc":
                    frame[col] = frame[col].astype(str)
        levels = sorted(obs[SPLIT].astype(str).unique())
        mode = "single" if expansion == "per_delivery" else "coarse"
        totals = scm._reporter_zero_counts(
            nz, reporter, mpra_map, cell_map, SPLIT, levels,
            reporter_expansion=mode)
        counted = int(totals["n_total"].sum())
        # The two modes report different things, per _build_obs_phantom_inputs:
        # in single mode the detection count IS the zero count and the total is
        # built up from it, while in coarse mode the count already includes the
        # nonzero observations.
        if mode == "single":
            return counted, n_nonzero + counted
        assert counted >= n_nonzero, (
            f"{name}: coarse total {counted:,} is below the {n_nonzero:,} "
            "nonzero observations it is meant to contain")
        return counted - n_nonzero, counted

    assert expansion in ("all_combinations", "all_combinations_moi"), \
        f"{name}: unhandled expansion {expansion}"
    obs_for_maps = (obs.groupby(["rep_id", "cell_bc", "mpra_bc"])["umis_mpra_bc"]
                    .sum().reset_index())
    cm_maps = {"cell_map": cell_map, "mpra_map": mpra_map, "observed": obs_for_maps}
    moi_correction = None
    if expansion == "all_combinations_moi":
        moi = float(dat.describe_transfection().mu_nb)
        n_lib = int(mpra_map["mpra_bc"].nunique())
        moi_correction = {"moi": moi, "n_library_barcodes": n_lib}
        p_t = 1.0 - (1.0 - 1.0 / n_lib) ** moi
        print(f"    MOI {moi:.4f}, {n_lib} library barcodes, "
              f"P(transfected)={p_t:.6f}")

    n_zeros = n_total = 0
    for level in sorted(obs[SPLIT].astype(str).unique()):
        gt = scm._cm_group_totals(cm_maps, SPLIT, level, moi_correction=moi_correction)
        n_zeros += int(gt["n_zeros"].sum())
        n_total += int(gt["n_zeros"].sum() + gt["n_nonzero"].sum())
    return n_zeros, n_total


def main():
    from dask.distributed import Client, LocalCluster
    cluster = LocalCluster(n_workers=1, threads_per_worker=1, processes=False,
                           memory_limit="60GB")
    client = Client(cluster)

    pi_frames, zero_rows = [], []
    for dataset, regime, rel in ORTHOS:
        name = Path(rel).name
        print(f"\n=== {dataset} / {regime}  ({name})", flush=True)
        ortho = scm.ortho.load(client, str(DATA_ROOT / Path(rel).parent), name)
        expansion = ortho.meta["zero_expansion"]
        print(f"    zero_expansion={expansion}")

        d = pi_rows(ortho, name)
        print(f"    pi: {len(d)} values over {d.level.nunique()} CREs x "
              f"{d.rep_id.nunique()} reps, median {d.pi.median():.4f}, "
              f"range {d.pi.min():.4f}-{d.pi.max():.4f}")
        pi_frames.append(d.assign(dataset=dataset, regime=regime, ortho=name))

        n_zeros, n_total = zero_fraction(ortho, name, expansion)
        assert 0 <= n_zeros <= n_total, f"{name}: {n_zeros} zeros of {n_total}"
        frac = n_zeros / n_total
        print(f"    zeros: {n_zeros:,} of {n_total:,} observations = {frac:.4f}")
        zero_rows.append({"dataset": dataset, "regime": regime, "ortho": name,
                          "zero_expansion": expansion, "n_zeros": n_zeros,
                          "n_total": n_total, "zero_fraction": frac})

    pi = pd.concat(pi_frames, ignore_index=True)[
        ["dataset", "regime", "ortho", "level", "rep_id", "pi"]]
    zero = pd.DataFrame(zero_rows)

    # Within a dataset the flooded regime must add zeros relative to the
    # baseline, or the panel's premise is wrong.
    for dataset, g in zero.groupby("dataset"):
        if len(g) > 1:
            print(f"\n{dataset}: zero fraction "
                  + ", ".join(f"{r.regime} {r.zero_fraction:.4f}"
                              for r in g.itertuples()))

    pi.to_csv(PI_TSV, sep="\t", index=False)
    zero.to_csv(ZERO_TSV, sep="\t", index=False)
    print(f"\nwrote {PI_TSV} ({len(pi)} rows) and {ZERO_TSV} ({len(zero)} rows)")
    client.close()


if __name__ == "__main__":
    main()
