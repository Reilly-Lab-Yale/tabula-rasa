# Cohen 5x5 activity simulation -- run summary

Reconstructed from sacct (logs deleted before runstats).

## Jobs

| Phase | JobID | Elapsed | MaxRSS | ReqMem | Partition | Node | Exit |
|-------|-------|---------|--------|--------|-----------|------|------|
| create | 7425377 | 01:01:28 | 45.7 GB | 128G | priority | a1132u07n03 | 0 |
| fit | 7438058 | 01:34:10 | 21.8 GB | 64G | priority_gpu | a1126u02n01 | 0 |
| wald+test | 7455866 | 00:14:08 | 38.9 GB | 128G | priority_gpu | a1122u11n01 | 0 |

Total wall time: ~2h50m (sequential phases).

## Configuration

- 5 GT draws x 5 sim reps = 25 orthos
- Model: ZINB (nb_only=False), by_cell_type only
- flatten_overtransfection=True
- reference_pooling from COHEN_BOUNDS.n_negative_controls
- Dask: LocalCluster, 1 worker, 2 threads, auto memory limit

## Results

MWU AUROC=0.592, AUPRC=0.469
Wald AUROC=0.739, AUPRC=0.721

Date: 2026-04-06
