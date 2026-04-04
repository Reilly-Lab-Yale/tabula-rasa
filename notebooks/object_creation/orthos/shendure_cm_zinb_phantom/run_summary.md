# Run Summary: Shendure CM ZINB phantom-zero ortho

## What was fitted

Shendure dataset, consider_missing (CM) condition, ZINB model, phantom-zero
weighting. Fits both by_cre and by_cell_type models with full Cartesian-product
zero expansion.

## Timing

- Job submitted: 2026-04-03 14:01:51
- Job started: 2026-04-03 14:07:36
- Job ended: 2026-04-03 16:15:58
- Elapsed (wall): 02:08:22
- by_cre fit: 28.2s
- by_cell_type fit: 5653.6s (94.2 min)
- Total fit time: 5686.8s (94.8 min)

## Resources

- Job ID: 7129286
- Node: a1132u11n01
- CPU model: Intel Xeon 8562Y+ (Emerald Rapids)
- Allocated: 4 CPUs, 256 GB RAM
- Peak RSS: 14.3 GB (MaxRSS 15810248K)
- CPU utilization: 31.1% (02:39:36 used / 08:33:28 allocated)
- Partition: priority (prio_skr2)
- Cluster: bouchet

## Notes

- MoM init from CM maps not yet implemented; all cell types fell back to NB init.
