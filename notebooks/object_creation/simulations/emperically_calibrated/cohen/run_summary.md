# Cohen 5x5 activity simulation -- run summary

## Run 2: obsingle NB canonical (2026-04-07)

| Phase | JobID | Elapsed | MaxRSS | ReqMem | Partition | Node | CPU% | GPU% | GPU mem% | Exit |
|-------|-------|---------|--------|--------|-----------|------|------|------|----------|------|
| create | 7502046 | 00:57:57 | 39.0 GB | 128G | priority | a1132u07n02 | 72% | -- | -- | 0 |
| fit | 7511917 | 00:26:48 | 25.0 GB | 64G | priority_gpu | a1124u11n01 | 32% | 32% | 99% | 0 |
| wald+test | 7513838 | 00:17:28 | 41.6 GB | 128G | priority_gpu | a1124u19n01 | 67% | 17% | 99% | 0 |

Total wall time: ~1h42m (sequential phases).

### Configuration

- 5 GT draws x 5 sim reps = 25 orthos
- Model: NB (nb_only=True), by_cell_type only
- Ortho: cohen_obsingle_nb_phantom (reporter_expansion="single")
- flatten_overtransfection=True
- reference_pooling from COHEN_BOUNDS.n_negative_controls
- Dask: LocalCluster, 1 worker, 2 threads, auto memory limit
- GPU: H200 (140 GB)

### Results

MWU AUROC=0.880, AUPRC=0.863
Wald AUROC=0.751, AUPRC=0.649 (high variance: GT0=0.89 vs GT1-4=0.58-0.65)

## Run 1: coarse ZINB (2026-04-06, superseded)

Reconstructed from sacct (logs deleted before runstats).

| Phase | JobID | Elapsed | MaxRSS | ReqMem | Partition | Node | Exit |
|-------|-------|---------|--------|--------|-----------|------|------|
| create | 7425377 | 01:01:28 | 45.7 GB | 128G | priority | a1132u07n03 | 0 |
| fit | 7438058 | 01:34:10 | 21.8 GB | 64G | priority_gpu | a1126u02n01 | 0 |
| wald+test | 7455866 | 00:14:08 | 38.9 GB | 128G | priority_gpu | a1122u11n01 | 0 |

Total wall time: ~2h50m (sequential phases).

### Configuration

- 5 GT draws x 5 sim reps = 25 orthos
- Model: ZINB (nb_only=False), by_cell_type only
- Ortho: cohen_obs_zinb_phantom (reporter_expansion="coarse")
- flatten_overtransfection=True
- reference_pooling from COHEN_BOUNDS.n_negative_controls
- Dask: LocalCluster, 1 worker, 2 threads, auto memory limit

### Results

MWU AUROC=0.592, AUPRC=0.469
Wald AUROC=0.739, AUPRC=0.721
