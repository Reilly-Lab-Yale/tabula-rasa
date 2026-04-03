# shendure_obs_nb run summary

**Fit:** Shendure, observed-only (no consider_missing), NB-only, by_cre + by_cell_type
**Result:** SUCCESS
**Date:** 2026-03-29
**Output:** `/nfs/roberts/project/pi_skr2/shared/tabula_data/shendure/shendure_obs_nb_20260329/`

## Successful run

### Driver (job 6819698)
| Field | Value |
|---|---|
| Node | a1130u05n04 |
| Hardware | 2× Intel Xeon 8562Y+ (Emerald Rapids), 64 cores, 1 TB RAM, no GPU |
| Partition | priority (CPU) |
| Allocated | 1 core, 16 GB |
| Start | 2026-03-29 21:33:29 |
| End | 2026-03-29 21:46:29 |
| Elapsed | 00:13:00 |
| Peak RSS | 1.22 GB |
| CPU utilization | ~17% (coordinator overhead; workers do the work) |

### Workers (jobs 6819708–6819711, 4× CPU)
| Field | Value |
|---|---|
| Nodes | a1132u24n01, a1132u24n04 |
| Hardware | 2× Intel Xeon 8562Y+ (Emerald Rapids), 64 cores, 1 TB RAM, no GPU |
| Partition | priority (CPU) |
| Allocated per worker | 1 core, 56 GB |
| Active time | ~10:25 (cancelled cleanly on driver exit) |
| Peak RSS (max across workers) | 1.99 GB |
| Peak RSS (avg across workers) | 1.37 GB |
| CPU utilization (avg across workers) | ~46% |

## Failed attempts (IndexingError — dask-expr boolean index bug)
| Job | Elapsed | Failure |
|---|---|---|
| 6818681 | 00:08:36 | IndexingError: Unalignable boolean Series |
| 6819421 | 00:02:39 | IndexingError: Unalignable boolean Series |
| 6819575 | 00:08:14 | IndexingError: Unalignable boolean Series |

Fix: `raw.compute().reset_index(drop=True)` before worker dispatch loop in `standard_fit` and `_compute_mats_futures` (see commit on `nb-vs-zinb-fits`).
