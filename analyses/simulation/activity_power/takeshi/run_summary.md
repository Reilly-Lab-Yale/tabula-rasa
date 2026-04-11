# takeshi activity power -- run summary

Welch's t-test power analysis for the Takeshi dataset (3 cell types: HepG2,
K562, SKNSH), both `+reporter` (obs condition) and deflated (drop-zeros)
test conditions. 156 library replicates x 5 sims per replicate per cell type.

## Job

| Field        | Value                                              |
|--------------|----------------------------------------------------|
| Driver job   | 7960231 (`takeshi_pow`)                            |
| Workers      | 7960245-7960264 (20x `takeshi_pow_worker`)         |
| Partition    | priority / prio_skr2                               |
| Script       | `takeshi_power_ttest_all_cell_types.py all`        |
| Outcome      | COMPLETED, exit 0                                  |

## Timing

| Field      | Value                                                |
|------------|------------------------------------------------------|
| Start      | 2026-04-11 12:24:57                                  |
| End        | 2026-04-11 13:25:15                                  |
| Elapsed    | 01:00:18                                             |
| Time limit | 14:00:00 (used 7%)                                   |

Workers ran ~59:06 each (started ~20s after the driver, cancelled cleanly
when SLURMCluster shut down).

## Hardware

Driver and all workers ran on the priority partition. Nodes used:
`a1132u22n01`, `a1132u22n03`, `a1132u24n01`, `a1132u33n01`, `a1132u33n04`.

CPU: Intel Xeon **8562Y+** (Emerald Rapids), `cpugen=emeraldrapids`.

## Resource utilization

| Process       | Alloc        | Peak RSS  | CPU util         |
|---------------|--------------|-----------|------------------|
| Driver (1x)   | 1 CPU, 64 GB | 16.4 GB   | 49.7%            |
| Workers (20x) | 1 CPU, 23 GB | 4.7-5.3 GB| 90.5-91.3% (avg 90.9%) |

- Driver memory: 16.4/64 GB (26%) -- 64 GB request was conservative; ~32 GB
  would have been ample. Worker memory: ~5/23 GB (~22%); 12-16 GB request
  would suffice next time.
- Driver CPU: 49.7% -- expected for a single-CPU coordinator (it spends time
  blocked on Dask futures).
- Worker CPU: ~91% across all 20 -- workers were saturated end-to-end.

## Outputs

In `output/`:
- `takeshi_power_df_reporter.parquet` -- aggregated `+reporter` test results
- `takeshi_power_df_deflated.parquet` -- aggregated deflated test results
- `takeshi_power_ttest_reporter_comparison.svg` -- per-cell-type power curves

Simulations live at
`/nfs/roberts/project/pi_skr2/shared/tabula_data_new/simulated/2026-04-11_takeshi_pow/{K562,SKNSH,reference}/`
(156 sims per cell type, all with both ttest and ttest_deflated results).

## Bounds

Used `scm.TAKESHI_BOUNDS` (canonical = `takeshi_obs_nb_phantom`, added in
this branch). Key params:

- minP (`reference_activity`) = 0.011637
- min activity = 0.001582 (`min_mpra_umi`)
- max activity = 4 * minP = 0.046547
- 3 cell types, reference = HepG2 (encoded as `"reference"` in bounds)
- 20 negative controls
- nb_only = True, consider_missing = False
