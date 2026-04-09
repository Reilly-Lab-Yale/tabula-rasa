# Seelig 5x5 activity simulation -- run summary

Re-run with MOIB NB canonical ortho (seelig_cm_moib_nb_phantom) and
MWU has_reporter=False fix (zeros dropped to match real no-reporter data).

## Jobs

| Phase | JobID | Elapsed | Peak RSS | ReqMem | Partition | Node | CPU | GPU | Exit |
|-------|-------|---------|----------|--------|-----------|------|-----|-----|------|
| create (failed, no PYTHONPATH) | 7682360 | 00:00:00 | 0.03 GB | 64G | priority | a1132u22n03 | -- | -- | 1 |
| create (failed, old data present) | 7682390 | 00:00:00 | 0.03 GB | 64G | priority | a1132u22n03 | -- | -- | 1 |
| create | 7682422 | 00:55:14 | 5.3 GB | 64G | priority | a1132u22n03 | 48% | -- | 0 |
| fit | 7685070 | 00:45:24 | 8.0 GB | 256G | priority_gpu | a1122u11n01 | 14% | 8% util, 99% mem | 0 |
| wald+test+metrics | 7688130 | 00:20:39 | 16.6 GB | 256G | priority_gpu | a1122u11n01 | 73% | 1% util, 99% mem | 0 |

Total wall time (successful jobs): ~2h01m (sequential phases).

## Hardware

- a1132u22n03: 2x Xeon 8562Y+ (Emerald Rapids), 64 cores, 1015 GB RAM
- a1122u11n01: 2x Xeon 6542Y (Emerald Rapids), 48 cores, 2044 GB RAM, H200 140 GB

## Configuration

- Canonical ortho: seelig_cm_moib_nb_phantom (MOIB NB, new canonical)
- Bounds: SEELIG_BOUNDS = SEELIG_CM_MOIB_NB_BOUNDS
- 5 GT draws x 5 sim reps = 25 orthos
- Model: NB (nb_only=True), by_cell_type only, consider_missing=True
- MWU: has_reporter=False (zeros dropped before testing)
- flatten_overtransfection=True
- Dask: LocalCluster, 1 worker, 2 threads, auto memory limit

## Results

All phases completed successfully (create, fit, wald_precomp, test, metrics).

| Test | AUROC (mean +/- std) | AUPRC (mean +/- std) |
|------|---------------------|---------------------|
| MWU | 0.829 +/- 0.019 | 0.799 +/- 0.022 |
| Wald | 0.629 +/- 0.021 | 0.603 +/- 0.017 |

MWU outperforms Wald on both metrics. MWU zeros were dropped (has_reporter=False)
to match what real Seelig data looks like on disc (no transfection reporter means
only non-zero MPRA signal is observed). Wald operates on MOIB-expanded data.

Saved: analyses/simulation/activity_prc/seelig/output/seelig_5x5_activity_summary.tsv

Date: 2026-04-08
