# Shendure 5x5 activity simulation -- run summary

Reconstructed from sacct (logs deleted before runstats).

## Jobs

| Phase | JobID | Elapsed | MaxRSS | ReqMem | Partition | Node | Exit |
|-------|-------|---------|--------|--------|-----------|------|------|
| create | 7425331 | 00:31:31 | 29.6 GB | 64G | priority | a1132u09n01 | 0 |
| fit | 7438057 | 02:01:25 | 14.5 GB | 64G | priority_gpu | a1122u11n01 | 0 |
| wald+test | 7455865 | 00:17:35 | 29.8 GB | 128G | priority_gpu | a1122u11n01 | 0 |

Total wall time: ~2h50m (sequential phases).

## Configuration

- 5 GT draws x 5 sim reps = 25 orthos
- Model: NB (nb_only=True), by_cell_type only
- flatten_overtransfection=True
- reference_pooling from SHENDURE_BOUNDS.n_negative_controls
- Dask: LocalCluster, 1 worker, 2 threads, auto memory limit

## Results

MWU AUROC=0.922, AUPRC=0.914
Wald AUROC=0.923, AUPRC=0.914

Date: 2026-04-06
