#!/usr/bin/env python
"""
Compute seelig collision rate (birthday-paradox estimate) directly from the
raw UMI-wise TSV, replicating the calculation in Bounds.from_ortho().

Does not require Dask or a cluster.
"""
import sys
import numpy as np
import pandas as pd

tsv = "/nfs/roberts/project/pi_skr2/shared/tabula_data/seelig/seelig_scmpra_umiwise.tsv.gz"

print("Reading TSV...", flush=True)
df = pd.read_csv(tsv, sep="\t")
print(f"  {len(df)} rows, columns: {list(df.columns)}", flush=True)

total_uniq_mpra_bc = df["mpra_bc"].nunique()
print(f"  Total unique mpra_bcs: {total_uniq_mpra_bc}", flush=True)

tfection = (
    df.groupby(["rep_id", "cell_bc"])["mpra_bc"]
    .nunique()
    .reset_index()
    .rename(columns={"mpra_bc": "unique_mpra_bc"})
)

observed = tfection["unique_mpra_bc"]
tfection["tot_plasmid"] = (
    np.log(1 - observed / total_uniq_mpra_bc)
    / np.log((total_uniq_mpra_bc - 1) / total_uniq_mpra_bc)
)

total_tfection = tfection["tot_plasmid"].sum()
excess = total_tfection - tfection["unique_mpra_bc"].sum()
pct = excess / total_tfection * 100

print(
    f"scMPRAforge: INFO: Computed a total of {excess} estimated collision events, "
    f"out of a total of {total_tfection}, or {pct}%"
)
