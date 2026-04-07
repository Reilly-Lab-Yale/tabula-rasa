# Seelig 5x5 activity simulation -- run summary

Reconstructed from sacct (logs deleted before runstats).

## Jobs

| Phase | JobID | Elapsed | MaxRSS | ReqMem | Partition | Node | Exit |
|-------|-------|---------|--------|--------|-----------|------|------|
| create | 7427287 | 00:50:05 | 7.9 GB | 64G | priority | a1132u11n03 | 0 |
| fit | 7438059 | 01:14:19 | 8.7 GB | 256G | priority_gpu | a1126u11n01 | 0 |
| wald+test | 7455867 | 00:05:25 | 11.3 GB | 256G | priority_gpu | a1124u11n01 | 1 |

Total wall time: ~2h10m (sequential phases).

## Configuration

- 5 GT draws x 5 sim reps = 25 orthos
- Model: NB (nb_only=True), by_cell_type only, consider_missing=True
- flatten_overtransfection=True
- reference_pooling from SEELIG_BOUNDS.n_negative_controls
- Dask: LocalCluster, 1 worker, 2 threads, auto memory limit

## Results

Wald precomp + wald test + MWU all completed for all 5 GT draws.
Metrics phase failed: _classifier_summary got empty array after NA filtering
(known issue -- 55% of seelig CM hypotheses have <3 nonzero obs).
Test results (p-values) are saved to disk.

Date: 2026-04-06
