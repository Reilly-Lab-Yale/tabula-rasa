# scMPRAforge — CLAUDE.md

## What this project is

scMPRAforge is a Python package for analysis of single-cell MPRA (Massively Parallel Reporter Assay) data. It fits stratified ZINB (zero-inflated negative binomial) or NB (negative binomial) models to scMPRA count data and provides hypothesis testing (Wald, MWU, bootstrap), simulation, and power analysis.

## Key concepts

- **ortho**: the core fitted object. Contains two model sets for the same data: `by_cre` (one model per CRE, with cell-type as treatment) and `by_cell_type` (one model per cell type, with CRE as treatment). Produced by `criss_cross()` or the individual `fit_by_cre_models()` / `fit_by_cell_type_models()` methods.
- **scMPRA_data**: wrapper around a Dask DataFrame of MPRA count data (UMI-wise or read-wise). Tracks operations applied (filtering, consider_missing, etc.).
- **Bounds**: empirical summary of an scMPRA experiment (MOI distribution, cells per cell type, library composition). Used to parameterize simulations. Presets: `SHENDURE_BOUNDS`, `COHEN_BOUNDS`.
- **consider_missing**: inflates the count matrix to include unobserved (cell, CRE) pairs as zeros. Essential for datasets without a transfection reporter (seelig), optional for datasets with one (shendure, cohen). Dramatically increases dataset size and compute cost.
- **TensorZINB**: external package (`tensorzinb-plusplus`) that performs the actual ZINB/NB fitting via TensorFlow. The `tz` conda environment has the version with sparse matrix support and GPU device placement.

## Three empirical datasets

| Dataset | Integration | Transfection reporter | consider_missing default | flatten_overtransfection |
|---------|------------|----------------------|--------------------------|--------------------------|
| Shendure (Lalanne et al. 2024) | piggyBac | oBC (barcode-level) | Off (reporter removes untransfected zeros) | True (measurement cannot resolve duplicate integrations) |
| Cohen | Non-integrating | CRE-level reporter | Off (reporter removes untransfected zeros) | True (CRE-level reporter cannot resolve barcode-level duplicates) |
| Seelig | Non-integrating | None | On (no reporter = no way to introduce zeros without it) | True (no reporter = cannot resolve duplicates) |

## Compute environment

- **Cluster**: Bouchet (Yale YCRC), formerly McCleary
- **Account**: `prio_skr2` (Priority Tier — no hard job/memory limits)
- **Partitions**: `priority` (CPU), `priority_gpu` (GPU, H200 or RTX 5000 Ada)
- **Conda env**: `tz` (TensorFlow 2.20, tensorzinb-plusplus v0.0.6, sparse support, GPU)
- **Dask**: distributed scheduler with SLURM workers. By-CRE fits are CPU-only. By-cell-type fits with consider_missing use GPU (H200, ~2h wall time) via `pre_fit_hook` pattern for just-in-time GPU worker submission.

## Key resource patterns (from run graveyards)

- **By-CRE fits**: CPU-only, 1 worker, ~1-2h, 15-60 GB RAM
- **By-cell-type fits (no consider_missing)**: CPU-only, ~1-2h, moderate RAM
- **By-cell-type fits (with consider_missing)**: GPU strongly preferred. Seelig: ~2h on 2x H200. Shendure: ~43h on 16x96GB CPU workers (only option due to design matrix size)
- **Driver jobs**: 1 CPU, 16-24 GB RAM, coordinator only
- **GPU idle preemption**: priority_gpu kills jobs not using GPU after ~60 min. Use `pre_fit_hook` to submit GPU workers just-in-time.

## Code layout

```
scMPRAforge/
  core.py           -- all main classes and functions (~7000 lines)
  utils.py          -- helper functions
  presets/           -- serialized Bounds objects
analyses/
  preprocessing/    -- raw data processing
  model_fitting/    -- ortho fitting scripts, QC, bounds extraction
    fits/           -- per-ortho fit scripts and wrappers
    qc/             -- QC plots and summary notebook
    bounds/         -- bounds extraction scripts and output
  model_selection/  -- NB vs ZINB (LRT, AIC, ZI decomposition, bias-variance)
  simulation/       -- empirically calibrated + synthetic simulations
    activity_prc/   -- 5x5 activity PRC (per dataset)
    activity_calibration/  -- null p-value calibration
    activity_power/ -- power curves
    variant_power/  -- pairwise variant effect power
    fully_synthetic/
  empirical_testing/ -- real-data hypothesis testing
  supplementary/    -- coupon collector, collision rates
```

## Running fits

All fit scripts follow the pattern: sbatch wrapper (.sh) -> ipython script (.py) -> Dask cluster -> TensorZINB. The wrapper activates `tz` conda env and `cd`s to the script directory (never use `dirname "$0"` in sbatch scripts). Scripts are in `analyses/model_fitting/fits/`.

## After a successful fit

1. Run `runstats` (`python3 ~/.slurm_run_stats.py`) in the job directory — produces a timestamped `run_stats_*.txt`.
2. Write `run_summary.md` from it — must include: start/end times, elapsed, node names, CPU model (from `scontrol show node` → `AvailableFeatures`), allocated vs peak RSS, CPU/GPU utilization (max and avg). Keep the raw `run_stats_*.txt` alongside it.
3. Delete `slurm-*.out`, `slurm-*.err`, `worker_*.out` logs only — do NOT delete `run_stats_*.txt`.
4. Commit `run_summary.md` + `run_stats_*.txt` with the fit scripts.

## Cohen data and consider_missing

Cohen's transfection reporter (U6) operates at CRE-level, not barcode-level. The original preprocessing inserted a single `mpra_bc="dummy"` zero row per (cell, CRE) when U6 detected transfection but no MPRA signal — this broke `consider_missing` because `"dummy"` mapped to all CREs, failing the `(rep_id, mpra_bc) -> cre_id` uniqueness check.

**Regenerated preprocessing** (`regen_scmpra_object.py`, March 2026) fixes both issues:
1. Filters 49,646 ambiguous rBC pairs (sequencing artifacts)
2. Performs CRE-coarse reporter-informed expansion: expands U6-only (cell, CRE) detections to all barcodes of that CRE (replacing a single dummy row with real barcode-level zeros)

The output `retina_single_counting_u6.scmpra/` has no dummy barcodes and passes `_get_missing_maps()` validation. It contains the full CRE-coarse reporter-informed zeros, which are a strict subset of what `consider_missing` would produce.

**Memory cap note**: `_inflate_missing_split_level` estimates ~8,750 GB for Cohen CM cell-type slices, but this assumes dense string-per-row representation. The actual data is 99.8% sparse zeros. Cohen CM fits must set `consider_missing_max_memory_gb = None` to bypass the cap. The estimator should eventually be updated to account for sparse encoding.

## Terminology: zero expansion modes

| Term | Verb form | Meaning |
|------|-----------|---------|
| **reporter-informed zeros** | reporter-informed expansion | Any zeros derived from a transfection reporter signal — the reporter confirms the CRE was present in the cell, so unobserved MPRA signal is a true zero not a missing observation. Encompasses all reporter-based zero imputation. |
| **CRE-coarse reporter-informed zeros** | CRE-coarse reporter-informed expansion | The specific case where the transfection reporter operates at CRE level (not barcode level), requiring imputation across all barcodes of that CRE. Because the reporter cannot distinguish which barcodes were present, every unobserved barcode of a detected CRE is imputed as zero. Cohen U6 is the canonical example. Produces many more zeros than a barcode-level reporter would (where only the specific unobserved barcode gets a zero). |
| **consider_missing zeros** | consider_missing expansion | Full Cartesian-product expansion: every (cell × barcode) combination in a replicate, regardless of reporter signal. Appropriate when there is no reporter (Seelig). |

**Relationship**: CRE-coarse reporter-informed zeros ⊂ reporter-informed zeros ⊂ consider_missing zeros.
The obs condition uses only reporter-informed zeros (or none); CM uses all possible zeros.

## Reporter-informed zero logic (Cohen U6 orphan handling)

When computing phantom zero weights via `_reporter_zero_counts`, there is an
ambiguous case: cells where MPRA signal is observed but U6 (transfection
reporter) is not detected (-U6 +MPRA, "orphan" observations).

Two failure modes: (A) spurious MPRA (index hopping / ambient RNA, U6 absence
is correct), or (B) false-negative U6 (CRE genuinely transfected, U6
undersequenced). Empirical evidence from Cohen data strongly favors (B): 64.5%
of U6-confirmed (cell, CRE) pairs have zero MPRA obs (even confirmed
transfections are frequently missed), and U6 is undersequenced.

The compromise: treat orphan cells as confirmed transfections.
`_reporter_zero_counts` augments the reporter-confirmed set with orphan cells
before computing phantom zero weights. This is conservative -- if (A) dominated,
the correct treatment would be to drop orphan obs entirely. Under this treatment
the phantom zero weight for a group cannot go negative (n_total >= n_nonzero by
construction).

## Important implementation notes

- `nb_only=True` parameter propagates through `standard_fit` → `_tensorzinb_fit` → `TensorZINB`. Disables MoM init automatically. Result dict has `nb_only=True` stamp and no `x_pi` weights.
- MoM initialization: used for by-cell-type ZINB fits. Falls back to NB init when P(X=0|NB) > 0.05 (ZI poorly identified). See seelig attempt 12-15 in run graveyard.
- `_extract_zi` and `_label_tensorzinb_regressors` guard against `nb_only=True` (no `x_pi`).
- AIC/BIC are not currently computed or stored. `llf_total` is available from every fit.
- `overtransfected()` reports per-event collision rate (not per-cell — birthday paradox makes per-cell misleading). Threshold: `WARN_MULTI_TRANSFECTION_PERCENT = 2.0%`.

## Data locations

- **Primary (clean)**: `/nfs/roberts/project/pi_skr2/shared/tabula_data_new/{shendure,cohen,seelig}/`
  Preprocessing outputs + phantom-zero orthos only. See `README.md` there for
  full provenance tracing each file back to GEO downloads.
- **Legacy**: `/nfs/roberts/project/pi_skr2/shared/tabula_data/` -- old orthos,
  simulations, intermediate files. Do not use for new fits.
- Simulations: `/nfs/roberts/project/pi_skr2/shared/tabula_data/simulated/`
  (will migrate to tabula_data_new when re-run)
