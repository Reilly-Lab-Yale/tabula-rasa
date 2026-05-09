# Seelig variant power -- run summary

One job, COMPLETED with exit 0.

Hardware: Intel Xeon 8562Y+ (Emerald Rapids), 2 sockets x 32 cores =
64 cores/node, 1014 GB RealMemory.

## Job 11219363 -- seelig_pw (MWU pairwise power, has_reporter=False)

- Script: `seelig_pairwise_power_mwu.py all`
- Start:   2026-05-09 13:39:48
- End:     2026-05-09 13:51:47
- Elapsed: 00:11:59
- Node:    a1130u05n03 (driver and many workers co-located)
- Alloc:   1 CPU, 96 GB
- Peak RSS (driver): 1.9 GB / 96 GB
- Driver CPU util: 6.6% (00:00:47 CPU / 00:11:59 wall)
- State:   COMPLETED, exit 0:0
- Dask workers: 20 (seelig_pw_worker, 45 GB each, 12h time limit, priority).
  Workers were CANCELLED at job end as the driver tore down the cluster.
  Per-worker CPU util range 5%-59% (median ~20%); peak RSS ~1.1 GB / 45 GB.
- Sims: 20/20 (1 cohort x 20 library reps x 5 sims, 2 cell types in each sim).
  Library size n_cres=1344 (full SEELIG MOIB-NB library).
- Outputs:
    - `output/seelig_pairwise_power_mwu_all_cell_types.svg` (combined panel)
    - `output/panels_mwu/seelig_pairwise_power_mwu_K562.svg`
    - `output/panels_mwu/seelig_pairwise_power_mwu_reference.svg` (HepG2)
    - `output/seelig_pairwise_power_50pct_contours.svg`
    - `output/seelig_pairwise_power_df.parquet` (38400 rows)
    - `output/seelig_pairwise_null_df.parquet` (1000 rows)
    - `output/seelig_pairwise_power_sim_registry.json`
- Sims on disk: `2026-05-09_seelig_pw/cohort_0/sim_*`
- Null FPR@0.05 (filler-vs-filler calibration, MWU deflated):
  K562 0.0560 (n=500), reference (HepG2) 0.0460 (n=500). Both well-calibrated.

## Notes

- Vastly over-provisioned. Pre-run estimate was 24-36h based on cohen scaling;
  actual run was 12 minutes. CPU and memory both used <10% of allocation.
  Future seelig pairwise MWU runs can comfortably target 1h walltime, 16 GB
  driver, and 16-24 GB workers.
- One cohort sufficed because (n_cres - 1) // 9 = 149 anchors/cohort vs 24
  total anchors (vs cohen which needed 2 cohorts at n_cres=115).
- Workers all co-located on the driver node (a1130u05n03) plus a handful of
  spillover nodes (a1130u05n04, a1130u11n02, etc.); since the run is short
  this is fine.
