# Cohen downstream re-runs required

The Cohen preprocessing was fixed in March 2026 (two bugs, see below). All downstream
analyses that used Cohen data need to be re-run.

## What changed

1. **Ambiguous rBC filter** (`break_by_barcode.ipynb` step 11): reads where the same
   rBC (MPRA random barcode) maps to more than one cBC/CRE within a replicate are now
   removed. These are sequencing artifacts — homopolymer/error reads producing a
   barcode sequence that happens to co-occur with multiple CREs. The Cohen et al. 2023
   paper reports ~99.5% unique cBC-rBC pairs in their DNA library; any cross-CRE
   collision in scRNA-seq data is error.

2. **Barcode-level U6 zero expansion** (`make_scmpra_object.ipynb` step 12): when U6
   detects a CRE in a cell but no MPRA barcodes are observed, the old code inserted a
   single dummy row (one zero per CRE). The fix expands to one zero per known barcode
   for that (rep, CRE), matching the granularity of non-zero observations.

## What needs re-running

- [ ] Cohen `COHEN_BOUNDS` preset (bounds extraction)
- [ ] Cohen power analyses
- [ ] Cohen p-value calibrations
- [ ] Any figures / results notebooks using Cohen data
- [ ] cohen_obs ortho fits (the underlying data changed)
- [ ] cohen_cm ortho fits (separate issue — see plan_cohen_preprocessing.md)
