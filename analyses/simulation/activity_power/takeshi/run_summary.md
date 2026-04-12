# takeshi activity power -- run summary (rerun on scratch with flatten=True)

Welch's t-test power analysis for the Takeshi dataset (3 cell types: HepG2,
K562, SKNSH), both `+reporter` (obs condition) and deflated (drop-zeros)
test conditions. 156 library replicates x 5 sims per replicate per cell type.

This is a **rerun** of the original 2026-04-11 takeshi power. The first run
used `flatten_overtransfection=False` (a bug inherited from the broken
shendure power script). The corrected rerun uses `flatten_overtransfection=True`
and writes to scratch instead of project storage. The earlier (broken)
parquets and SVG in `output/` have been overwritten.

## Job

| Field        | Value                                              |
|--------------|----------------------------------------------------|
| Driver job   | 8009056 (`takeshi_pow`)                            |
| Workers      | 7960245-7960264, 7966245-..., 8009057-... (across 3 attempts) |
| Partition    | priority / prio_skr2                               |
| Script       | `takeshi_power_ttest_all_cell_types.py all`        |
| Outcome      | COMPLETED, exit 0                                  |
| Sim location | `/nfs/roberts/scratch/pi_skr2/mcn26/simulated/2026-04-11_takeshi_pow/` (scratch) |

## Timing

| Field      | Value                                                |
|------------|------------------------------------------------------|
| Start      | 2026-04-11 17:05:32                                  |
| End        | 2026-04-11 18:06:10                                  |
| Elapsed    | 01:00:38                                             |
| Time limit | 14:00:00 (used 7%)                                   |

Net wall clock matches the original (broken) run almost exactly -- the
flatten=True fix doesn't materially change compute cost.

## Hardware

Driver and workers ran on the priority partition. Driver node:
`a1130u05n03`. Workers spread across many nodes.

CPU: Intel Xeon **8562Y+** (Emerald Rapids), `cpugen=emeraldrapids`.

## Resource utilization

| Process       | Alloc        | Peak RSS  | CPU util         |
|---------------|--------------|-----------|------------------|
| Driver (1x)   | 1 CPU, 64 GB | 20.1 GB   | 50.5%            |
| Workers (20x) | 1 CPU, 23 GB | typically 4-6 GB | typically 90% |

- Driver memory: 20.1/64 GB (31%) -- 64 GB request was conservative;
  ~32 GB would have been ample.
- Driver CPU: 50.5% -- in line with the original run; the script is
  fundamentally a single-CPU coordinator that blocks on Dask futures.
- Worker CPU: ~90% across the 20 workers -- end-to-end saturated.

## Outputs

In `output/`:
- `takeshi_power_df_reporter.parquet` -- aggregated `+reporter` test results
- `takeshi_power_df_deflated.parquet` -- aggregated deflated test results
- `takeshi_power_ttest_reporter_comparison.svg` -- per-cell-type power curves

Simulations live at
`/nfs/roberts/scratch/pi_skr2/mcn26/simulated/2026-04-11_takeshi_pow/{K562,SKNSH,reference}/`
(156 sims per cell type, all with both ttest and ttest_deflated results).

## Bounds

Used `scm.TAKESHI_BOUNDS` (canonical = `takeshi_obs_nb_phantom`).
- minP (`reference_activity`) = 0.011637
- min activity = 0.001582 (`min_mpra_umi`)
- max activity = 4 * minP = 0.046547
- 3 cell types, reference = HepG2 (encoded as `"reference"` in bounds)
- 20 negative controls
- nb_only = True, consider_missing = False

## Caveats

- **t-test FPR is inflated on this dataset** -- the takeshi t-test
  calibration showed FPR ~0.10-0.11 in `+reporter` and ~0.14-0.15 in
  deflated, vs nominal 0.05. The MWU calibration cross-check showed
  nominal FPR (~0.05) on the same sims, so the inflation is t-test
  specific, not a property of the data. The power-curve numbers in
  this run are therefore over-optimistic; consider re-running with
  `sim.mwu("hs_all_ct")` on the same on-disk sims to get an unbiased
  surface (no resimulation needed).
- Five worker subjobs landed on the bad node `a1130u22n02` during the
  first two attempts and got silent-killed; that node is now in the
  `--exclude` list in both the wrapper and the Dask SLURMCluster config.
- The first attempt also crashed on disk-quota in `tabula_data_new`;
  switching to scratch fixed it (and tabula_data legacy was cleaned up
  separately, freeing 211 GB).

## Run history

| Attempt | Job ID | Outcome | Note |
|---------|--------|---------|------|
| 1 (flatten=False) | 7960231 | COMPLETED but stale | original buggy run |
| 2 (flatten=True, project) | 7966268 | killed (rate limit) | 0 workers got submitted |
| 3 (flatten=True, project) | 7978455 | silent kill on a1130u22n02 | bad node |
| 4 (flatten=True, project) | 7999601 | crashed (disk quota) | tabula_data_new full |
| **5 (flatten=True, scratch)** | **8009056** | **COMPLETED** | **canonical run** |
