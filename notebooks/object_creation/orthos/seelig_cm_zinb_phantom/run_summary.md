# Run Summary: Seelig CM ZINB phantom-zero ortho

## What was fitted

Seelig dataset, consider_missing (CM) condition, ZINB model, phantom-zero
weighting. Fits both by_cre and by_cell_type models with full Cartesian-product
zero expansion.

## Timing

- Job submitted: 2026-04-03 14:01:54
- Job started: 2026-04-03 14:07:36
- Job ended: 2026-04-03 15:41:30
- Elapsed (wall): 01:33:54
- by_cre fit: 74.1s
- by_cell_type fit: 4037.3s (67.3 min)
- Total fit time: 4117.1s (68.6 min)

## Resources

- Job ID: 7129290
- Node: a1132u11n03
- CPU model: Intel Xeon 8562Y+ (Emerald Rapids)
- Allocated: 4 CPUs, 256 GB RAM
- Peak RSS: 49.9 GB (MaxRSS 55467296K)
- CPU utilization: 34.0% (02:07:43 used / 06:15:36 allocated)
- Partition: priority (prio_skr2)
- Cluster: bouchet

## Notes

- MoM init from CM maps not yet implemented for K562 and reference; fell back
  to NB init.
