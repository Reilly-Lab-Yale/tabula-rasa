# Plan: Cohen data refresh + ortho refit

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

## A. Preprocessing regen — DONE

`regen_scmpra_object.py` + `wrap_regen.sh` — single standalone script replaces
notebook re-run (raw data files on palmer_scratch no longer available).

**Outputs:**
- `unjoined/read_wise_mpra_retina_filtered.tsv` — ambiguity-filtered reads
- `retina_single_counting_u6.scmpra/` — 229 parquet parts, UMI-wise with barcode-level U6 zeros

### Verification (2026-03-30)

Old = `unjoined/retina_single_counting_u6.tsv` (original preprocessing, Oct 2025)
New = `retina_single_counting_u6.scmpra/` (regen, Mar 2026)

| Metric | Old | New |
|--------|----:|----:|
| Total rows | 7,351,469 | 1,493,337,547 |
| Zero rows | 75,954 (1.0%) | 1,490,148,989 (99.8%) |
| Nonzero rows | 7,275,515 | 3,188,558 |
| Dummy barcode rows | 75,954 | 0 |
| Ambiguous (rep, mpra_bc) pairs | 49,646 | 0 |
| Unique cells | 22,127 | 22,101 |
| Unique barcodes | 2,013,178 | 1,965,339 |
| Unique CREs | 115 | 115 |

- **3,188,558 shared nonzero rows** — all UMI counts agree perfectly (0 mismatches)
- **4,086,957 lost nonzero rows** — all from ambiguous rBC filtering (confirmed: all 49,646 lost barcode identities were ambiguous in old data)
- **0 gained nonzero rows** — no new signal introduced
- New object is 200× larger due to barcode-level zero expansion replacing single-count dummies

### Zero expansion: U6 zeros are a subset of consider_missing

The U6-informed zeros only cover (cell, CRE) pairs where U6 detected transfection.
`consider_missing` produces the full Cartesian product (cell × barcode) per rep:

| | Rows |
|---|---:|
| Current new object (U6 zeros) | 1.49 billion |
| consider_missing would produce | 21.8 billion |
| CM is this many × larger | 14.6× |

Per rep: Rep 1 = 11,329 cells × 1,041,845 barcodes; Rep 2 = 10,802 cells × 924,902 barcodes.

For comparison, Seelig CM = 68.7M rows (10,336 cells × 6,648 barcodes, 1 rep).
Cohen CM would be **317× larger** than Seelig CM. The difference is library depth:
Cohen has ~1M barcodes vs Seelig's ~6.6K.

### Feasibility of cohen_cm fits

**Sparse encoding**: `.compute()` preserves `pd.SparseDtype` columns — the 1.49B-row
object (99.8% zeros) stays sparse through Dask collection. Design matrices from
`formulaic` use `scipy.sparse` CSR format. TensorZINB accepts sparse matrices natively.

**int32 index overflow**: Not an issue. scipy 1.17.1 (installed in `tz` env) auto-promotes
`indptr`/`indices` to int64 when matrix dimensions exceed 2^31. By-cell-type slices at
~11B rows exceed int32 max (2.1B) but are handled by auto-promotion. TensorFlow sparse
tensors use int64 indices by default. No code changes needed.

**Remaining concern — `_prepare_subset_for_modeling()`**: This function (core.py:332-338)
calls `pd.to_numeric(...).astype("int64")` on `umis_mpra_bc`, which densifies the sparse
column. This is a one-line fix (check if already sparse, skip conversion) but should be
addressed before attempting cohen_cm fits.

**Bottom line**: cohen_cm is architecturally feasible if the sparse path is kept end-to-end.
The by-CRE fits are smaller (11K cells × ~8.5K barcodes per CRE ≈ 93M rows each — very
manageable). The by-cell-type fits (~11B rows each) are the challenge — comparable in
spirit to Shendure CM but with wider design matrices. GPU (H200) strongly recommended.

## B. cohen_obs ortho fit

Fit using the new `retina_single_counting_u6.scmpra/` object with `consider_missing=False`.
This is the "observed + U6-informed zeros" condition — the standard Cohen analysis with
corrected preprocessing.

- by_cre: CPU-only, should be straightforward (3.2M nonzero + 1.49B sparse-zero rows,
  but with consider_missing=False the fit only sees observed data)
- by_cell_type: same — moderate resource needs

## C. cohen_cm ortho fit

Must load filtered pre-join data (not the U6-padded file) and let consider_missing
handle all zero expansion:
```python
cohen = scm.scMPRA_data.from_tsv(str(path / "unjoined/read_wise_mpra_retina_filtered.tsv"))
cohen.read_wise_to_umi_wise()
cohen.set_consider_missing(True)
```

- by_cre: ~93M rows per CRE slice, CPU feasible
- by_cell_type: ~11B rows per slice, GPU strongly recommended, may need
  `_prepare_subset_for_modeling` sparse-densify fix first

## Order of operations

1. ✅ Preprocessing regen (regen_scmpra_object.py)
2. ✅ Verify new object (old-vs-new comparison, 2026-03-30)
3. □ Write + submit cohen_obs ortho fit scripts (B)
4. □ Fix `_prepare_subset_for_modeling` sparse densification (one-line fix in core.py)
5. □ Write + submit cohen_cm ortho fit scripts (C)
6. □ Merge results into ortho objects
