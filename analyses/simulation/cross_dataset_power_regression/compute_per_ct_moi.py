#!/usr/bin/env python3
"""
Per-cell-type MOI from raw scMPRA td for shendure/cohen/seelig.

Mirrors the corrected `describe_transfection` logic (groupby on
`(rep_id, cell_bc)`) but aggregates per cell_type instead of fitting one
global NB. The "MOI" reported here is the per-cell mean of unique mpra_bc
detected (no coupon-collector correction -- matches what the bounds
`get_effective_moi` would yield, but per CT). Reads each dataset's raw tsv
once and writes `per_ct_moi.parquet`.
"""

from pathlib import Path
import pandas as pd

OUT = Path(__file__).resolve().parent / "output" / "per_ct_moi.parquet"

DATASETS = {
    "shendure": {
        "path": "/nfs/roberts/project/pi_skr2/shared/tabula_data_new/shendure/shendure_processed.tsv",
        "compression": None,
        "ref_ct": "Pluripotent",
    },
    "cohen": {
        "path": "/nfs/roberts/project/pi_skr2/shared/tabula_data_new/cohen/retina_single_counting_u6.tsv",
        "compression": None,
        "ref_ct": "Rod",
    },
    "seelig": {
        "path": "/nfs/roberts/project/pi_skr2/shared/tabula_data_new/seelig/seelig_scmpra_umiwise.tsv.gz",
        "compression": "gzip",
        "ref_ct": "HepG2",
    },
}

rows = []
for ds, cfg in DATASETS.items():
    print(f"=== {ds} ===", flush=True)
    df = pd.read_csv(
        cfg["path"], sep="\t", compression=cfg["compression"],
        usecols=["rep_id", "cell_bc", "cell_type", "mpra_bc"],
        dtype=str,
    )
    print(f"  rows: {len(df)}", flush=True)

    # per-cell unique barcodes (correct cell identity = (rep_id, cell_bc))
    per_cell = (
        df.groupby(["rep_id", "cell_bc"], observed=True)
        .agg(n_mpra_bc=("mpra_bc", "nunique"), cell_type=("cell_type", "first"))
        .reset_index()
    )
    print(f"  unique cells: {len(per_cell)}", flush=True)

    # rename the dataset's reference cell type to "reference" to match
    # how the bounds objects label it
    per_cell["cell_type"] = per_cell["cell_type"].where(
        per_cell["cell_type"] != cfg["ref_ct"], "reference"
    )

    # per-CT mean and median of per-cell unique mpra_bc
    summ = (
        per_cell.groupby("cell_type", observed=True)["n_mpra_bc"]
        .agg(["mean", "median", "size"])
        .reset_index()
        .rename(columns={"mean": "moi_mean", "median": "moi_median", "size": "n_cells"})
    )
    summ["dataset"] = ds
    print(summ.to_string(index=False), flush=True)
    rows.append(summ)

out = pd.concat(rows, ignore_index=True)
OUT.parent.mkdir(exist_ok=True)
out.to_parquet(OUT)
print(f"\nSaved: {OUT}", flush=True)
