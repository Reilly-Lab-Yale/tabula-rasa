# Cohen CM ZINB phantom -- run summary

**Job:** 7000006
**Node:** a1132u09n01
**CPU:** Intel Xeon 8562Y+ (Emerald Rapids), 64 cores/node
**Allocated:** 4 CPUs, 256 GB RAM

## Timing

| Phase | Time |
|-------|------|
| by_cre | 66.1s |
| by_cell_type | 1861.9s (31.0 min) |
| Total fit | 1931.9s (32.2 min) |
| Wall time | 38:10 |

Start: 2026-04-01 20:54
End:   2026-04-01 21:33

## Resource usage

| Metric | Value |
|--------|-------|
| Peak RSS | 31.6 GB |
| Allocated RAM | 256 GB |
| RAM utilization | 12% |
| CPU utilization | 39% |

## Notes

- Phantom-compressed CM path: TSV (3.2M rows) + consider_missing.
- NB init (MoM not implemented for CM maps; NB init is correct for CM).
- consider_missing_max_memory_gb = None (bypass dense-string estimator).
- 256 GB over-allocated; 64 GB would suffice.
