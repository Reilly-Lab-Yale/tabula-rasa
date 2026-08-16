#!/usr/bin/env python
"""Per-element expression among nonzero observations, for every dataset.

Written to test whether the identifiability of the zero-inflation parameter
tracks how large the nonzero counts are. The mechanism would be that a
negative binomial with a small mean already puts most of its mass at zero, so
extra zeros are indistinguishable from ordinary sampling and the inflation
term buys nothing; with a large mean the NB's own P(0) is small, so excess
zeros can only be structural and pi is determined.

Joined against lambda = 2(lnL_ZINB - lnL_NB) from lrt_nb_vs_zinb_results.tsv,
this gives the test both between datasets and, more usefully, within them.

Statistics are taken over nonzero rows of the source table, which is a
property of the dataset rather than of any expansion, so one ortho per dataset
is enough and the numbers apply to all of that dataset's fits.

    sbatch analyses/model_selection/run_nonzero_expression.sh
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa")
import scMPRAforge.core as scm  # noqa: E402

DATA_ROOT = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new")
OUT_TSV = Path(__file__).resolve().parent / "nonzero_expression.tsv"

# One ortho per dataset, chosen only as a handle on that dataset's counts.
# takeshi is not one of the manuscript's three, but it is a fourth independent
# observation of the effect and the test is stronger with it.
SOURCES = {
    "shendure": "shendure/shendure_obs_nb",
    "cohen": "cohen/cohen_obsingle_nb_phantom",
    "seelig": "seelig/seelig_cm_nb_phantom",
    "takeshi": "takeshi/takeshi_obs_nb",
}


def main():
    from dask.distributed import Client, LocalCluster
    cluster = LocalCluster(n_workers=1, threads_per_worker=1, processes=False,
                           memory_limit="60GB")
    client = Client(cluster)

    frames = []
    for dataset, rel in SOURCES.items():
        name = Path(rel).name
        print(f"\n=== {dataset} ({name})", flush=True)
        ortho = scm.ortho.load(client, str(DATA_ROOT / Path(rel).parent), name)
        data = ortho.training_data.get_data(include_missing=False)
        if hasattr(data["umis_mpra_bc"].dtype, "fill_value"):
            data["umis_mpra_bc"] = data["umis_mpra_bc"].astype("int64")
        obs = scm._to_pandas(data[["cre_id", "umis_mpra_bc"]])

        n_rows = len(obs)
        nz = obs[obs["umis_mpra_bc"] > 0]
        assert len(nz), f"{dataset}: no nonzero observations"
        print(f"    {n_rows:,} rows in the source table, {len(nz):,} nonzero "
              f"({len(nz) / n_rows:.1%})")

        d = (nz.groupby("cre_id")["umis_mpra_bc"]
             .agg(n_nonzero="size", mean_nonzero="mean", median_nonzero="median")
             .reset_index().rename(columns={"cre_id": "model"}))
        d["model"] = d["model"].astype(str)
        # Total observations per element in the source table, so the zero
        # fraction of the raw data is available too.
        tot = obs.groupby("cre_id").size().rename("n_source").reset_index()
        tot["cre_id"] = tot["cre_id"].astype(str)
        d = d.merge(tot, left_on="model", right_on="cre_id").drop(columns="cre_id")
        d["frac_zero_source"] = 1 - d["n_nonzero"] / d["n_source"]

        print(f"    {len(d)} elements; mean nonzero UMI median "
              f"{d.mean_nonzero.median():.3f}, "
              f"range {d.mean_nonzero.min():.3f}-{d.mean_nonzero.max():.3f}")
        frames.append(d.assign(dataset=dataset))

    out = pd.concat(frames, ignore_index=True)[
        ["dataset", "model", "n_nonzero", "n_source", "frac_zero_source",
         "mean_nonzero", "median_nonzero"]]
    assert (out.mean_nonzero >= 1).all(), "a nonzero count below 1 is impossible"
    out.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"\nwrote {OUT_TSV} ({len(out)} rows)")
    client.close()


if __name__ == "__main__":
    main()
