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
  core.py           — all main classes and functions (~7000 lines)
  utils.py          — helper functions
  presets/           — serialized Bounds objects
notebooks/
  preprocessing/    — raw data processing
  object_creation/  — ortho fitting scripts, simulation setup
    orthos/         — fit scripts, wrappers, run graveyards
    simulations/    — empirically calibrated simulations
  results/          — analysis notebooks, figures
```

## Running fits

All fit scripts follow the pattern: sbatch wrapper (.sh) → ipython script (.py) → Dask cluster → TensorZINB. The wrapper activates `tz` conda env and `cd`s to the script directory (never use `dirname "$0"` in sbatch scripts). Scripts are in `notebooks/object_creation/orthos/`.

## After a successful fit

1. Run `runstats` (`python3 ~/.slurm_run_stats.py`) in the job directory — produces a timestamped `run_stats_*.txt`.
2. Write `run_summary.md` from it — must include: start/end times, elapsed, node names, CPU model (from `scontrol show node` → `AvailableFeatures`), allocated vs peak RSS, CPU/GPU utilization (max and avg). Keep the raw `run_stats_*.txt` alongside it.
3. Delete `slurm-*.out`, `slurm-*.err`, `worker_*.out` logs only — do NOT delete `run_stats_*.txt`.
4. Commit `run_summary.md` + `run_stats_*.txt` with the fit scripts.

## Cohen data and consider_missing

Cohen's transfection reporter (U6) operates at CRE-level, not barcode-level. The preprocessing notebook (`cohen_retina_denovo_process/12_make_scmpra_object/`) joins MPRA data with U6 data: when U6 detects a CRE in a cell but no MPRA barcodes are observed, it inserts a single row with `mpra_bc="dummy"` and `umis_mpra_bc=0` ("single counting U6 approach" from the original paper). This creates 75,954 zero rows.

**Problem**: `_get_missing_maps()` (called by `consider_missing`) requires unique `(rep_id, mpra_bc) -> cre_id` mapping. The sentinel `"dummy"` maps to all CREs, failing validation with 49,647 ambiguous keys. This makes `consider_missing=True` incompatible with the U6-padded TSV.

**Solution for cohen_cm fits**: load the pre-U6-join MPRA data (`unjoined/read_wise_mpra_retina.tsv`), convert to UMI-wise, and let `consider_missing` handle all zero expansion. This bypasses U6 entirely — which is the point of the "ignore reporter" condition.

**Separate issue — single-counting underestimates zeros**: the current preprocessing adds 1 zero per (cell, CRE) when U6 says the CRE is present but no MPRA signal is observed. Biologically, all barcodes for that CRE produced zero UMIs, so the zero should be expanded to all barcodes. The preprocessing should be updated to perform this barcode-level expansion for the standard (obs) condition.

## Important implementation notes

- `nb_only=True` parameter propagates through `standard_fit` → `_tensorzinb_fit` → `TensorZINB`. Disables MoM init automatically. Result dict has `nb_only=True` stamp and no `x_pi` weights.
- MoM initialization: used for by-cell-type ZINB fits. Falls back to NB init when P(X=0|NB) > 0.05 (ZI poorly identified). See seelig attempt 12-15 in run graveyard.
- `_extract_zi` and `_label_tensorzinb_regressors` guard against `nb_only=True` (no `x_pi`).
- AIC/BIC are not currently computed or stored. `llf_total` is available from every fit.
- `overtransfected()` reports per-event collision rate (not per-cell — birthday paradox makes per-cell misleading). Threshold: `WARN_MULTI_TRANSFECTION_PERCENT = 2.0%`.

## Data locations

- Empirical data: `/nfs/roberts/project/pi_skr2/shared/tabula_data/{shendure,cohen,seelig}/`
- Ortho objects: same directories, named `*_ortho_*`
- Simulations: `/nfs/roberts/project/pi_skr2/shared/tabula_data/simulated/`
