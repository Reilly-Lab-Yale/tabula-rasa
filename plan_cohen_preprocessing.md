# Plan: Cohen preprocessing fix + cohen_cm pipeline

## Problems discovered (March 2026)

1. **Ambiguous rBCs**: 49,646 (rep, mpra_bc) pairs mapping to multiple CREs in the raw
   MPRA data. These are sequencing artifacts — homopolymer/error rBCs that co-occur
   with multiple cBCs. Cohen et al. 2023 library has ~99.5% unique cBC-rBC pairs;
   cross-CRE collisions in scRNA-seq data are errors. Must filter before any downstream
   use.

2. **Single-counting U6 zeros underestimates**: preprocessing inserts 1 dummy row per
   (cell, CRE) when U6 detects CRE but no MPRA signal. Biologically, each unobserved
   barcode of that CRE in that cell is a true zero. Fix: expand to all barcodes of
   that (rep, CRE) observed elsewhere in the dataset.

3. **cohen_cm incompatibility**: `mpra_bc="dummy"` (from the U6 zero-padding) breaks
   `_get_missing_maps()` uniqueness check. Even with unique dummies per CRE, the
   dummy barcodes would be treated as real barcodes by consider_missing. Solution:
   cohen_cm scripts must load pre-U6-join filtered MPRA data and let consider_missing
   handle zero expansion.

## Implementation

### A. Preprocessing regen (job 6818181 running)

`regen_scmpra_object.py` + `wrap_regen.sh` — single standalone script replaces
notebook re-run (raw data files on palmer_scratch no longer available). Steps:

1. ✅ Filter ambiguous rBCs from `read_wise_mpra_retina.tsv`
   → `unjoined/read_wise_mpra_retina_filtered.tsv` (3,597,480 reads, -59.6%)
2. ⏳ Convert to UMI-wise, join with U6, expand right_only rows to barcode level
3. ⏳ Save `retina_single_counting_u6.tsv` + `retina_single_counting_u6.scmpra/`

### B. cohen_obs fit scripts

After A completes, update `cohen_obs_nb/fit.py` and `cohen_cm_zinb/` / `cohen_cm_nb/`
obs scripts to use `from_parquet("retina_single_counting_u6.scmpra")` instead of
`from_tsv(...)`.

### C. cohen_cm fit scripts

Must load filtered pre-join data (not the U6-padded file):
```python
cohen = scm.scMPRA_data.from_tsv(str(path / "unjoined/read_wise_mpra_retina_filtered.tsv"))
cohen.read_wise_to_umi_wise()
cohen.set_consider_missing(True)
```
Files to update: `cohen_cm_zinb/fit_by_cre.py`, `fit_by_cell_type.py`,
`cohen_cm_nb/fit_by_cre.py`, `fit_by_cell_type.py`.

**Not yet done** — waiting for A to complete and confirm data is clean.

## Order of operations

1. ⏳ Wait for job 6818181 (A) to complete
2. □ Verify new TSV/parquet (check row counts, no dummy barcodes)
3. □ Update cohen_obs scripts to use `from_parquet` (B)
4. □ Update cohen_cm scripts to use filtered pre-join data (C)
5. □ Resubmit cohen_cm_zinb and cohen_cm_nb fits
6. □ Re-run cohen_obs ortho fits (data changed — old fits are stale)
