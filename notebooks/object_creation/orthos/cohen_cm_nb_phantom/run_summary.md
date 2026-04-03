# Cohen CM NB phantom -- run summary

**Job:** 7000644
**Node:** a1132u18n03
**CPU:** Intel Xeon 8562Y+ (Emerald Rapids), 64 cores/node
**Allocated:** 4 CPUs, 256 GB RAM

## Timing

| Phase | Time |
|-------|------|
| by_cre | 71.2s |
| by_cell_type | 1310.2s (21.8 min) |
| Total fit | 1389.6s (23.2 min) |
| Wall time | 28:44 |

Start: 2026-04-01 21:10
End:   2026-04-01 21:39

## Resource usage

| Metric | Value |
|--------|-------|
| Peak RSS | 30.8 GB |
| Allocated RAM | 256 GB |
| RAM utilization | 12% |
| CPU utilization | 38% |

## Notes

- NB-only counterpart to cohen_cm_zinb_phantom.
- Phantom-compressed CM path: TSV (3.2M rows) + consider_missing.
- ~28% faster than ZINB (23.2 vs 32.2 min fit time).
- consider_missing_max_memory_gb = None (bypass dense-string estimator).
- 256 GB over-allocated; 64 GB would suffice.
