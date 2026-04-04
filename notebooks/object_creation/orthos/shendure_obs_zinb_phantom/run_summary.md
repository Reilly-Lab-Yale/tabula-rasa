# Run Summary: Shendure obs ZINB phantom-zero ortho

## What was fitted

Shendure dataset, obs (no consider_missing) condition, ZINB model, phantom-zero
weighting. Fits both by_cre and by_cell_type models.

## Timing

- Job submitted: 2026-04-03 14:49:05
- Job started: 2026-04-03 14:51:04
- Job ended: 2026-04-03 15:22:59
- Elapsed (wall): 00:31:55
- by_cre fit: 6.9s
- by_cell_type fit: 1425.6s (23.8 min)
- Total fit time: 1437.4s (24.0 min)

## Resources

- Job ID: 7131601
- Node: a1132u31n01
- CPU model: Intel Xeon 8562Y+ (Emerald Rapids)
- Allocated: 4 CPUs, 64 GB RAM
- Peak RSS: 3.4 GB (MaxRSS 3583868K)
- CPU utilization: 37.2% (00:47:32 used / 02:07:40 allocated)
- Partition: priority (prio_skr2)
- Cluster: bouchet

## Notes

- All by_cell_type MoM initializations fell back to NB init (P(X=0|NB) > 0.05
  for all cell types), indicating ZI component poorly identified in obs data.
