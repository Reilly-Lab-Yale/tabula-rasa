# Shendure variant power -- run summary

Two jobs in this directory, both COMPLETED with exit 0.

Hardware for all nodes (driver + workers): Intel Xeon 8562Y+ (Emerald Rapids),
2 sockets x 32 cores = 64 cores/node, 1014 GB RealMemory.

## Job 7850790 -- shend_pw (MWU pairwise power)

- Script: `shendure_pairwise_power_mwu.py all`
- Start:   2026-04-10 11:56:14
- End:     2026-04-10 13:52:04
- Elapsed: 01:55:50
- Node:    a1130u18n01
- Alloc:   1 CPU, 32 GB
- Peak RSS (driver): 3.9 GB / 32 GB
- Driver CPU util: 12.2% (0:14:11 CPU / 1:55:50 wall)
- State:   COMPLETED, exit 0:0
- Dask workers: 20 (shendure_pw_worker, 23 GB each), priority partition.
- Sims: 460/460 (10 cell types, 23 total cohorts x 20 library reps x 5 sims)
- Outputs: `output/shendure_pairwise_power_mwu_all_cell_types.svg`,
  `output/shendure_pairwise_power_50pct_contours.svg`,
  `output/shendure_pairwise_power_df.parquet`,
  `output/shendure_pairwise_null_df.parquet`,
  `output/shendure_pairwise_power_sim_registry.json`
- Sims on disk: `2026-04-10_shendure_pw/<cell_type>/cohort_*/sim_*`
- Null FPR@0.05 (filler-vs-filler calibration check inside power sims):
  all 10 cell types 0.046-0.062, consistent with target 0.05.

## Job 7858625 -- shend_null (t-test null calibration)

- Script: `shendure_pairwise_null_calibration_ttest.py all`
- Start:   2026-04-10 14:41:35
- End:     2026-04-10 15:32:49
- Elapsed: 00:51:14
- Node:    a1132u07n02
- Alloc:   1 CPU, 32 GB
- Peak RSS (driver): 3.1 GB / 32 GB
- Driver CPU util: 11.9% (0:06:26 CPU / 0:51:14 wall)
- State:   COMPLETED, exit 0:0
- Dask workers: 20 (shendure_null_worker), priority partition.
- Sims: 200/200 (10 cell types x 1 cohort x 20 reps x 5 sims, 24 null pairs each)
- Outputs: `output/shendure_null_ttest_fpr_heatmap.svg`,
  `output/shendure_null_ttest_fpr_lines.svg`,
  `output/shendure_null_ttest_df.parquet`,
  `output/shendure_null_ttest_fpr_grid.parquet`,
  `output/shendure_null_ttest_sim_registry.json`
- Sims on disk: `2026-04-10_shendure_null_ttest/<cell_type>/cohort_0/sim_*`
- Overall null FPR@0.05 per cell type: 0.0396 - 0.0533 (all within sampling
  noise of 0.05; lowest is `reference` = Pluripotent at 0.040).

## Notes

- Both drivers are light CPU and light memory (peak 3-4 GB). 8 GB alloc would
  have been plenty. The entire driver workload is coordinating Dask tasks and
  aggregating results in the plot phase.
- Shendure power ran much faster than Cohen power (1h56m vs 3h58m) even
  though it has 11.5x more sims (460 vs 40): Shendure sims are per-cell-type
  (smaller libraries, smaller per-sim cost), while Cohen sims carry all 4
  cell types together with the full episomal library.
- Wall clock very comfortably under the 14 h time limit; 4 h would be plenty
  for reruns.

Raw data: `run_stats_20260410_115612.txt` (alongside this file).
