# takeshi pairwise variant power -- run summary

Welch t-test pairwise variant power for the Takeshi dataset (3 cell types:
HepG2, K562, SKNSH). 24 baseline-mu anchors x 8 log2-FC offsets, split into
2 cohorts per cell type to fit n_cres budget. 20 library replicates per
cohort, 5 sims per replicate.

The run was rebooted multiple times today (rate-limit, bad-node silent kill,
disk-quota crash). The successful final run resumed from 80 partial sims
already on disk (K562 cohort 0+1, SKNSH cohort 0+1) and only had to run the
remaining ~40 (reference cohort 0+1).

## Job

| Field        | Value                                              |
|--------------|----------------------------------------------------|
| Driver job   | 8009077 (`takeshi_pw`)                             |
| Workers      | 60 worker jobs across the 4 attempts               |
| Partition    | priority / prio_skr2                               |
| Script       | `takeshi_pairwise_power_ttest.py all`              |
| Outcome      | COMPLETED, exit 0                                  |
| Sim location | `/nfs/roberts/scratch/pi_skr2/mcn26/simulated/2026-04-11_takeshi_pw/` (scratch) |

## Timing (final attempt only)

| Field      | Value                                                |
|------------|------------------------------------------------------|
| Start      | 2026-04-11 17:05:48                                  |
| End        | 2026-04-11 17:25:50                                  |
| Elapsed    | 00:20:02                                             |
| Time limit | 14:00:00 (used 2%)                                   |

The 20-min wall time is misleadingly fast: it only ran ~40 fresh sims;
the other 80 were resumed from earlier (failed) runs that completed K562
and SKNSH cohorts before crashing. End-to-end the work took ~3-4h across
the 4 attempts.

## Hardware

Driver and workers ran on the priority partition. Final-attempt nodes
included `a1132u33n02` (driver), and many worker nodes across the
a1130u/a1132u pools.

CPU: Intel Xeon **8562Y+** (Emerald Rapids), `cpugen=emeraldrapids`.

## Resource utilization (final attempt only)

| Process       | Alloc        | Peak RSS  | CPU util         |
|---------------|--------------|-----------|------------------|
| Driver (1x)   | 1 CPU, 64 GB | 2.3 GB    | 14.6%            |
| Workers (20x) | 1 CPU, 23 GB | varies    | typically 15-30% |

- Driver memory: 2.3/64 GB (4%) -- request was very oversized; 8 GB
  would suffice next time. Note: failed-run drivers used different
  numbers (the disk-crash driver hit higher RSS), so this is final-run
  only.
- Driver CPU: 14.6% -- low because the script is fundamentally a
  sequential coordinator iterating per-rep through the cohort loop and
  blocking on each sim's gamut+ttest.
- Worker CPU: highly variable (5-58%) -- per-rep work distribution is
  uneven, since the driver submits one rep at a time rather than
  fanning out the full batch.

## Outputs

In `output/`:
- `takeshi_pairwise_power_ttest_df.parquet` -- per-(cell_type, baseline_mu, log2_fc) test results
- `takeshi_pairwise_null_ttest_df.parquet` -- null-check FPR (filler_* CRE pairs)
- `takeshi_pairwise_power_ttest_all_cell_types.svg` -- per-CT power heatmap grid
- `takeshi_pairwise_power_ttest_50pct_contours.svg` -- min log2(FC) for 50% power vs baseline mu
- `takeshi_pairwise_power_sim_registry.json` -- list of contributing sims

Simulations live at
`/nfs/roberts/scratch/pi_skr2/mcn26/simulated/2026-04-11_takeshi_pw/{K562,SKNSH,reference}/cohort_{0,1}/`
(120 sims total: 6 cohorts x 20 reps).

## Anchor grid

- Baseline mu range: `[MU_P5, MU_P95]` of takeshi by-CT NB mu values, in 24 geometric steps
- log2(FC) offsets: 0.02, 0.03, 0.05, 0.07, 0.10, 0.15, 0.20, 0.30
- 1 reference CRE at minP per cohort
- N_FILLER_NULL=5 pairs of filler CREs at minP for FPR sanity check
- ALPHA = 0.05

## Caveats

- The takeshi t-test calibration showed FPR inflation (~2-3x nominal alpha
  on this dataset; see `analyses/simulation/activity_calibration/takeshi/`).
  The pairwise power numbers in this run reflect that bias -- powers are
  somewhat over-optimistic. The MWU follow-up cal showed nominal calibration,
  so a future MWU rerun on the same hs_pairwise sims would give a less-biased
  power surface. Sims are saved on scratch and can be re-tested with
  sim.mwu("hs_pairwise") without resimulating.
