"""
Regenerate the Cohen scMPRA data object with two fixes applied:

1. Filter ambiguous rBCs: remove reads where the same rBC maps to more than one
   CRE (pBC) within a replicate (sequencing artifacts — see Cohen et al. 2023
   Methods; DNA library has ~99.5% unique cBC-rBC pairs).

2. U6-informed zero expansion: for (cell, CRE) pairs where U6 detects transfection
   but no MPRA barcodes are observed, add zeros for ALL barcodes of that CRE
   observed elsewhere in that replicate. Rationale: U6 confirms the CRE arrived;
   every unobserved barcode of that CRE in that cell is therefore a true zero.
   Processed one CRE at a time to avoid OOM.

Inputs:
  .../cohen/unjoined/read_wise_mpra_retina.tsv   — raw read-wise MPRA
  .../cohen/unjoined/u6.tsv                      — U6 reporter data

Outputs:
  .../cohen/unjoined/read_wise_mpra_retina_filtered.tsv  — ambiguity-filtered reads
  .../cohen/retina_single_counting_u6.scmpra/            — scMPRA_data parquet dir
  .../cohen/retina_single_counting_u6.tsv                — TSV (legacy compat)
"""
import os, sys, shutil
import pandas as pd
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path

sys.path.insert(0, '/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa')
import scMPRAforge as scm

DATA_ROOT = Path('/nfs/roberts/project/pi_skr2/shared/tabula_data/cohen')
UNJOINED  = DATA_ROOT / 'unjoined'
PARTS_DIR = DATA_ROOT / '_zero_parts'

# ── 1. Load and filter ambiguous rBCs ─────────────────────────────────────────
filtered_tsv = UNJOINED / 'read_wise_mpra_retina_filtered.tsv'
if filtered_tsv.exists():
    print("[1] Loading pre-filtered reads...", flush=True)
    mpra_filt = pd.read_csv(filtered_tsv, sep='\t')
    print(f"    {len(mpra_filt):,} reads", flush=True)
else:
    print("[1] Loading and filtering ambiguous rBCs...", flush=True)
    mpra_raw = pd.read_csv(UNJOINED / 'read_wise_mpra_retina.tsv', sep='\t')
    bc_cre = mpra_raw.groupby(['rep_id', 'mpra_bc'])['cre_id'].nunique().reset_index()
    bc_cre.columns = ['rep_id', 'mpra_bc', 'n_cre']
    unambig = bc_cre[bc_cre['n_cre'] == 1][['rep_id', 'mpra_bc']]
    mpra_filt = mpra_raw.merge(unambig, on=['rep_id', 'mpra_bc'], how='inner')
    n_removed = len(mpra_raw) - len(mpra_filt)
    print(f"    Removed {n_removed:,} reads ({n_removed/len(mpra_raw)*100:.2f}%) ambiguous rBCs", flush=True)
    print(f"    Remaining: {len(mpra_filt):,} reads", flush=True)
    mpra_filt.to_csv(filtered_tsv, sep='\t', index=False)

# ── 2. Convert to UMI-wise ─────────────────────────────────────────────────────
print("[2] Converting to UMI-wise...", flush=True)
umi_cols = ['cell_bc', 'rep_id', 'cre_id', 'cell_type', 'mpra_bc']
mpra_umi = (
    mpra_filt
    .groupby(umi_cols)['umi']
    .nunique()
    .reset_index()
    .rename(columns={'umi': 'umis_mpra_bc'})
)
print(f"    {len(mpra_umi):,} UMI-wise rows", flush=True)

# ── 3. Load and aggregate U6 ───────────────────────────────────────────────────
print("[3] Loading U6 data...", flush=True)
u6_raw = pd.read_csv(UNJOINED / 'u6.tsv', sep='\t')
u6_umi = (
    u6_raw
    .groupby(['rep_id', 'cell_bc', 'cell_type', 'cre_id'])['umi']
    .nunique()
    .reset_index()
    .rename(columns={'umi': 'umis_transfection_bc'})
)
print(f"    {len(u6_umi):,} U6 (cell, CRE) entries", flush=True)

# ── 4. Find right_only rows (U6 present, no MPRA) ─────────────────────────────
print("[4] Identifying U6-only (cell, CRE) pairs...", flush=True)
mpra_pairs = mpra_umi[['rep_id', 'cell_bc', 'cre_id']].drop_duplicates()
u6_joined = u6_umi.merge(mpra_pairs, on=['rep_id', 'cell_bc', 'cre_id'],
                          how='left', indicator=True)
right_only = u6_joined[u6_joined['_merge'] == 'left_only'].drop(columns='_merge')
print(f"    {len(right_only):,} right-only (U6 no MPRA) rows", flush=True)

# ── 5. Build bc_map: unique barcodes per (rep, CRE) ───────────────────────────
print("[5] Building barcode map...", flush=True)
bc_map = mpra_umi[['rep_id', 'cre_id', 'mpra_bc']].drop_duplicates()
print(f"    {len(bc_map):,} unique (rep, CRE, barcode) entries", flush=True)

# ── 6. Expand right_only rows to barcode level, streaming to parquet ──────────
PARTS_DIR.mkdir(exist_ok=True)

# pyarrow schema for zero rows
schema = pa.schema([
    ('cell_bc', pa.string()),
    ('rep_id', pa.string()),
    ('cre_id', pa.string()),
    ('cell_type', pa.string()),
    ('mpra_bc', pa.string()),
    ('umis_mpra_bc', pa.int64()),
    ('umis_transfection_bc', pa.float64()),
])

existing_parts = list(PARTS_DIR.glob('zeros_*.parquet'))
if existing_parts:
    total_written = sum(pq.read_metadata(p).num_rows for p in existing_parts)
    print(f"[6] Skipping expansion — {len(existing_parts)} zero parts already exist "
          f"({total_written:,} rows). Proceeding to step 7.", flush=True)
else:
    print("[6] Expanding U6-only zeros to barcode level (chunked)...", flush=True)
    total_written = 0
    CELL_CHUNK = 500  # process this many cells at a time per CRE to control memory

    groups = list(right_only.groupby(['rep_id', 'cre_id']))
    for i, ((rep_id, cre_id), group) in enumerate(groups):
        barcodes = bc_map[
            (bc_map['rep_id'] == rep_id) & (bc_map['cre_id'] == cre_id)
        ]['mpra_bc'].values

        if len(barcodes) == 0:
            continue

        cells = group[['cell_bc', 'cell_type', 'umis_transfection_bc']].values
        n_barcodes = len(barcodes)

        # Process cells in chunks to cap per-iteration memory
        part_tables = []
        for start in range(0, len(cells), CELL_CHUNK):
            chunk_cells = cells[start:start + CELL_CHUNK]
            n_chunk = len(chunk_cells)

            cell_bc_arr   = np.repeat(chunk_cells[:, 0], n_barcodes)
            cell_type_arr = np.repeat(chunk_cells[:, 1], n_barcodes)
            u6_arr        = np.repeat(chunk_cells[:, 2].astype(float), n_barcodes)
            mpra_bc_arr   = np.tile(barcodes, n_chunk)

            tbl = pa.table({
                'cell_bc': pa.array(cell_bc_arr, type=pa.string()),
                'rep_id': pa.array([str(rep_id)] * len(cell_bc_arr), type=pa.string()),
                'cre_id': pa.array([str(cre_id)] * len(cell_bc_arr), type=pa.string()),
                'cell_type': pa.array(cell_type_arr, type=pa.string()),
                'mpra_bc': pa.array(mpra_bc_arr, type=pa.string()),
                'umis_mpra_bc': pa.array(np.zeros(len(cell_bc_arr), dtype=np.int64)),
                'umis_transfection_bc': pa.array(u6_arr, type=pa.float64()),
            }, schema=schema)
            part_tables.append(tbl)

        combined = pa.concat_tables(part_tables)
        part_file = PARTS_DIR / f'zeros_{rep_id}_{cre_id.replace("/","_")}.parquet'
        pq.write_table(combined, part_file, compression='snappy')
        total_written += len(combined)

        if i % 20 == 0 or i == len(groups) - 1:
            print(f"    [{i+1}/{len(groups)}] {cre_id} rep{rep_id}: "
                  f"{len(combined):,} zeros written (total so far: {total_written:,})", flush=True)

    print(f"    Total zero rows written: {total_written:,}", flush=True)

# ── 7. Assemble final table ────────────────────────────────────────────────────
print("[7] Assembling final table...", flush=True)

# Non-right_only (both + left_only) from original UMI data joined with U6
non_right = mpra_umi.copy()
non_right['umis_mpra_bc'] = non_right['umis_mpra_bc'].astype(int)
both_u6 = u6_umi[['rep_id', 'cell_bc', 'cre_id', 'umis_transfection_bc']]
non_right = non_right.merge(both_u6, on=['rep_id', 'cell_bc', 'cre_id'], how='left')
non_right['umis_transfection_bc'] = non_right['umis_transfection_bc'].fillna(0)
non_right['rep_id'] = non_right['rep_id'].astype(str)

print(f"    Non-zero rows: {len(non_right):,}", flush=True)

# Write non-zero rows as parquet too then combine with dask
non_right_tbl = pa.Table.from_pandas(non_right[list(schema.names)], schema=schema)
pq.write_table(non_right_tbl, PARTS_DIR / 'non_zero.parquet', compression='snappy')

# Copy parts directly into scMPRA_data parquet directory — no in-memory re-read
print("[8] Assembling scMPRA_data parquet directory (file copy, no re-read)...", flush=True)
import json as _json

out_scmpra = DATA_ROOT / 'retina_single_counting_u6.scmpra'
data_parquet_dir = out_scmpra / 'data.parquet'
data_parquet_dir.mkdir(parents=True, exist_ok=True)

all_parts = sorted(PARTS_DIR.glob('*.parquet'))
print(f"    Copying {len(all_parts)} part files...", flush=True)
for i, part_file in enumerate(all_parts):
    dest = data_parquet_dir / f'part.{i}.parquet'
    shutil.copy2(part_file, dest)
    if i % 20 == 0 or i == len(all_parts) - 1:
        print(f"    [{i+1}/{len(all_parts)}] {part_file.name} → {dest.name}", flush=True)

with open(out_scmpra / 'members.json', 'w') as _f:
    _json.dump({"table_type": "mpra_umiwise", "source": str(out_scmpra)}, _f)
print(f"    Parquet directory saved: {out_scmpra}", flush=True)

# Cleanup temp parts
shutil.rmtree(PARTS_DIR)
print("[DONE]", flush=True)
