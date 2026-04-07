# Cohen obs NB phantom -- run summary

**Job:** 7000645
**Node:** a1132u22n03
**CPU:** Intel Xeon 8562Y+ (Emerald Rapids), 64 cores/node
**Allocated:** 4 CPUs, 64 GB RAM

## Timing

| Phase | Time |
|-------|------|
| by_cre | 24.1s |
| by_cell_type | 1004.3s (16.7 min) |
| Total fit | 1033.0s (17.2 min) |
| Wall time | 33:25 |

Start: 2026-04-01 21:10
End:   2026-04-01 21:44

## Resource usage

| Metric | Value |
|--------|-------|
| Peak RSS | 3.1 GB |
| Allocated RAM | 64 GB |
| RAM utilization | 5% |
| CPU utilization | 55% |

## Notes

- NB-only counterpart to cohen_obs_zinb_phantom.
- Phantom-compressed obs path with coarse U6 reporter table.
- Orphan cell fix applied: -U6+MPRA cells treated as confirmed transfections.
- ~41% faster than ZINB (17.2 vs 29.3 min fit time).
- 64 GB over-allocated; 8 GB would suffice.
