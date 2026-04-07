"""
Takeshi dataset preprocessing.

Input:  all.integrate.qced.xls (tab-separated, already QC'd to oBC UMI > 5)
Output: takeshi_scmpra_umiwise.scmpra/ (parquet directory)
        takeshi_scmpra_umiwise.tsv     (for inspection)

Column mapping:
    mBC         -> mpra_bc
    CELL        -> cell_bc
    mBC_count   -> umis_mpra_bc
    Celltype    -> cell_type
    oBC_count   -> umis_transfection_bc
    CRE_name    -> cre_id
    (added)     -> rep_id = "R1"
    error_corrected_mBC -> dropped
    nCount_RNA          -> dropped
"""

import sys
sys.path.insert(0, "/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa-takeshi")

import pandas as pd
import numpy as np
import scMPRAforge as scm

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DATA_ROOT = "/nfs/roberts/project/pi_skr2/shared/tabula_data_new/takeshi"
INPUT = f"{DATA_ROOT}/all.integrate.qced.xls"
OUT_SCMPRA = f"{DATA_ROOT}/takeshi_scmpra_umiwise.scmpra"
OUT_TSV = f"{DATA_ROOT}/takeshi_scmpra_umiwise.tsv"

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
print("Loading raw data...")
df = pd.read_csv(INPUT, sep="\t")
print(f"  Raw rows: {len(df):,}")
print(f"  Columns: {list(df.columns)}")

# ---------------------------------------------------------------------------
# Rename columns to standard form
# ---------------------------------------------------------------------------
df = df.rename(columns={
    "mBC": "mpra_bc",
    "CELL": "cell_bc",
    "mBC_count": "umis_mpra_bc",
    "Celltype": "cell_type",
    "oBC_count": "umis_transfection_bc",
    "CRE_name": "cre_id",
})

# Drop unused columns
df = df.drop(columns=["error_corrected_mBC", "nCount_RNA"])

# Add replicate ID (single replicate)
df["rep_id"] = "R1"

# Cast count columns to int
df["umis_mpra_bc"] = df["umis_mpra_bc"].astype(int)
df["umis_transfection_bc"] = df["umis_transfection_bc"].astype(int)

# Standardize cell type names (remove hyphens for consistency)
df["cell_type"] = df["cell_type"].replace({"SK-N-SH": "SKNSH"})

print(f"\n  Renamed columns: {list(df.columns)}")
print(f"  Cell types: {df['cell_type'].value_counts().to_dict()}")
print(f"  Unique CREs: {df['cre_id'].nunique()}")
print(f"  Unique cells: {df['cell_bc'].nunique()}")
print(f"  Unique mBCs: {df['mpra_bc'].nunique()}")

# ---------------------------------------------------------------------------
# Diagnostic: are the NA:NA CREs negative controls?
# Compare mean UMI counts for NA:NA vs variant CREs
# ---------------------------------------------------------------------------
print("\n--- CRE group comparison (potential negative controls) ---")

is_na_na = df["cre_id"].str.contains(":NA:NA$", regex=True)
groups = {
    "NA:NA CREs (n={})".format(df.loc[is_na_na, "cre_id"].nunique()):
        df.loc[is_na_na],
    "Variant CREs (n={})".format(df.loc[~is_na_na, "cre_id"].nunique()):
        df.loc[~is_na_na],
}

for label, sub in groups.items():
    n_rows = len(sub)
    mean_umi = sub["umis_mpra_bc"].mean()
    median_umi = sub["umis_mpra_bc"].median()
    frac_zero = (sub["umis_mpra_bc"] == 0).mean()
    print(f"  {label}:")
    print(f"    rows={n_rows:,}  mean_umi={mean_umi:.3f}  "
          f"median_umi={median_umi:.1f}  frac_zero={frac_zero:.3f}")

# Per-CRE breakdown for NA:NA
print("\n  Per-CRE mean UMI (NA:NA):")
na_na_stats = (df.loc[is_na_na]
               .groupby("cre_id")["umis_mpra_bc"]
               .agg(["mean", "median", "size"])
               .sort_values("mean"))
print(na_na_stats.to_string(float_format="%.3f"))

# Per-CRE breakdown for variant CREs (just summary stats)
variant_means = (df.loc[~is_na_na]
                 .groupby("cre_id")["umis_mpra_bc"]
                 .mean())
print(f"\n  Variant CRE mean UMI: "
      f"min={variant_means.min():.3f}  "
      f"median={variant_means.median():.3f}  "
      f"max={variant_means.max():.3f}")

# Per cell-type breakdown
print("\n  Per cell-type mean UMI by CRE group:")
for ct in sorted(df["cell_type"].unique()):
    sub_ct = df[df["cell_type"] == ct]
    na_mean = sub_ct.loc[is_na_na[sub_ct.index], "umis_mpra_bc"].mean()
    var_mean = sub_ct.loc[~is_na_na[sub_ct.index], "umis_mpra_bc"].mean()
    print(f"    {ct:8s}  NA:NA={na_mean:.3f}  variant={var_mean:.3f}  "
          f"ratio={na_mean/var_mean:.3f}")

# ---------------------------------------------------------------------------
# Diagnostic: ortho_filter impact preview
# ortho_filter drops (cell_type, cre_id) combos with < 3 nonzero observations
# ---------------------------------------------------------------------------
MIN_PTS = 3
print(f"\n--- ortho_filter preview (MIN_PTS={MIN_PTS}) ---")

nonzero = df[df["umis_mpra_bc"] > 0]
nz_counts = nonzero.groupby(["cell_type", "cre_id"]).size().reset_index(name="nonzero_count")

# All combos present in data
all_combos = df[["cell_type", "cre_id"]].drop_duplicates()
n_total_combos = len(all_combos)

# Merge to get counts for all combos (some may have 0 nonzero)
combo_stats = all_combos.merge(nz_counts, on=["cell_type", "cre_id"], how="left")
combo_stats["nonzero_count"] = combo_stats["nonzero_count"].fillna(0).astype(int)

would_drop = combo_stats[combo_stats["nonzero_count"] < MIN_PTS]
would_keep = combo_stats[combo_stats["nonzero_count"] >= MIN_PTS]

print(f"  Total (cell_type, cre_id) combos: {n_total_combos}")
print(f"  Would keep: {len(would_keep)}")
print(f"  Would drop: {len(would_drop)}")

if len(would_drop) > 0:
    print(f"\n  Dropped combos by cell type:")
    drop_by_ct = would_drop.groupby("cell_type").size()
    keep_by_ct = would_keep.groupby("cell_type").size()
    for ct in sorted(df["cell_type"].unique()):
        n_drop = drop_by_ct.get(ct, 0)
        n_keep = keep_by_ct.get(ct, 0)
        print(f"    {ct:8s}  keep={n_keep:4d}  drop={n_drop:4d}")

    # Are NA:NA CREs disproportionately dropped?
    drop_na = would_drop[would_drop["cre_id"].str.contains(":NA:NA$", regex=True)]
    drop_var = would_drop[~would_drop["cre_id"].str.contains(":NA:NA$", regex=True)]
    print(f"\n  Dropped NA:NA combos: {len(drop_na)} / {len(would_drop[would_drop['cre_id'].str.contains(':NA:NA$', regex=True)])} total NA:NA")
    print(f"  Dropped variant combos: {len(drop_var)} / {len(would_drop[~would_drop['cre_id'].str.contains(':NA:NA$', regex=True)])} total variant")

    # Show the actual dropped combos
    print(f"\n  All dropped combos:")
    for _, row in would_drop.sort_values(["cell_type", "cre_id"]).iterrows():
        print(f"    {row['cell_type']:8s} x {row['cre_id']:40s}  nonzero={row['nonzero_count']}")

# Row-level impact
rows_in_dropped = df.merge(would_drop[["cell_type", "cre_id"]], on=["cell_type", "cre_id"], how="inner")
print(f"\n  Rows lost to ortho_filter: {len(rows_in_dropped):,} / {len(df):,} "
      f"({100*len(rows_in_dropped)/len(df):.2f}%)")

# ---------------------------------------------------------------------------
# Save TSV
# ---------------------------------------------------------------------------
print(f"\nSaving TSV to {OUT_TSV} ...")
df.to_csv(OUT_TSV, sep="\t", index=False)

# ---------------------------------------------------------------------------
# Save as scMPRA_data object
# ---------------------------------------------------------------------------
print(f"Saving scMPRA_data to {OUT_SCMPRA} ...")
obj = scm.scMPRA_data()
obj.data = df
obj.table_type = "mpra_umiwise"
obj.source = OUT_SCMPRA
obj.metadata = {
    "dataset": "takeshi",
    "description": "Takeshi scMPRA: 3 cell lines (K562, HepG2, SKNSH), "
                   "180 CREs (80 emVar loci x 2 alleles + 18 NA:NA + 2 orphan), "
                   "piggyBac integration, oBC transfection reporter",
    "qc": "Pre-filtered to oBC UMI > 5 by Takeshi",
    "preprocessing": "Column rename only, no aggregation or filtering applied",
}
obj.to_parquet(OUT_SCMPRA)

print("Done.")
