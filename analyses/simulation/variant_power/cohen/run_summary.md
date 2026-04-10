# Cohen variant power -- run summary

Two jobs in this directory, both COMPLETED with exit 0.

Hardware for all nodes (driver + workers): Intel Xeon 8562Y+ (Emerald Rapids),
2 sockets x 32 cores = 64 cores/node, 1014 GB RealMemory.

## Job 7850789 -- cohen_pw (MWU pairwise power)

- Script: `cohen_pairwise_power_mwu.py all`
- Start:   2026-04-10 11:56:12
- End:     2026-04-10 15:54:09
- Elapsed: 03:57:57
- Node:    a1130u09n04
- Alloc:   1 CPU, 32 GB
- Peak RSS (driver): 20.0 GB / 32 GB
- Driver CPU util: 28.7% (1:08:27 CPU / 3:57:57 wall)
- State:   COMPLETED, exit 0:0
- Dask workers: 20 (cohen_pw_worker, 23 GB each), priority partition.
  Spot-check on 7850795: 40.9% CPU util, 8.0 GB peak RSS / 23 GB.
- Sims: 40/40 (2 cohorts x 20 library reps x 5 sims, 4 cell types in each sim)
- Outputs: `output/cohen_pairwise_power_mwu_all_cell_types.svg`,
  `output/cohen_pairwise_power_50pct_contours.svg`,
  `output/cohen_pairwise_power_df.parquet`,
  `output/cohen_pairwise_null_df.parquet`,
  `output/cohen_pairwise_power_sim_registry.json`
- Sims on disk: `2026-04-10_cohen_pw/cohort_{0,1}/sim_*`
- Null FPR@0.05 (filler-vs-filler calibration check inside power sims):
  Bipolar 0.0500, Interneuron 0.0600, Mueller Glia 0.0483, reference 0.0550.

## Job 7858624 -- cohen_null (t-test null calibration)

- Script: `cohen_pairwise_null_calibration_ttest.py all`
- Start:   2026-04-10 14:41:32
- End:     2026-04-10 16:37:42
- Elapsed: 01:56:10
- Node:    a1132u07n02
- Alloc:   1 CPU, 32 GB
- Peak RSS (driver): 15.8 GB / 32 GB
- Driver CPU util: 29.2% (0:34:01 CPU / 1:56:10 wall)
- State:   COMPLETED, exit 0:0
- Dask workers: 20 (cohen_null_worker), priority partition.
- Sims: 20/20 (1 cohort x 20 reps x 5 sims, 4 cell types, 24 anchor null pairs each)
- Outputs: `output/cohen_null_ttest_fpr_heatmap.svg`,
  `output/cohen_null_ttest_fpr_lines.svg`,
  `output/cohen_null_ttest_df.parquet`,
  `output/cohen_null_ttest_fpr_grid.parquet`,
  `output/cohen_null_ttest_sim_registry.json`
- Sims on disk: `2026-04-10_cohen_null_ttest/cohort_0/sim_*`
- Overall null FPR@0.05: Bipolar 0.0475, Interneuron 0.0483,
  Mueller Glia 0.0529, reference 0.0533 (all within sampling noise of 0.05).

## Notes

- Driver CPU efficiency looks low because the driver coordinates Dask workers
  and waits on their results; the real compute happens in the workers. That
  is expected for this pipeline.
- Both jobs wildly under time limit (14 h requested). 6 h limit would be
  plenty for a rerun.
- Peak driver RSS is close to half the 32 GB alloc for cohen_pw and under
  half for cohen_null; 24 GB could probably hold both.

Raw data: `run_stats_20260410_115611.txt` (alongside this file).
