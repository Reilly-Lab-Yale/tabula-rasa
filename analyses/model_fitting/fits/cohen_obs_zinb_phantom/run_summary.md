# Cohen obs ZINB phantom -- run summary

**Job:** 7000007
**Node:** a1130u37n02
**CPU:** Intel Xeon 8562Y+ (Emerald Rapids), 64 cores/node
**Allocated:** 4 CPUs, 64 GB RAM

## Timing

| Phase | Time |
|-------|------|
| by_cre | 23.3s |
| by_cell_type | 1727.7s (28.8 min) |
| Total fit | 1756.3s (29.3 min) |
| Wall time | 46:12 |

Start: 2026-04-01 20:54
End:   2026-04-01 21:41

## Resource usage

| Metric | Value |
|--------|-------|
| Peak RSS | 6.1 GB (6220548K) |
| Allocated RAM | 64 GB |
| RAM utilization | 7% |
| CPU utilization | 49% |

## Notes

- Phantom-compressed obs path with coarse U6 reporter table.
- Orphan cell fix applied: -U6+MPRA cells treated as confirmed transfections
  (see reporter_zero_logic.md in project root).
- NB init for by_cell_type (MoM not implemented for phantom path).
- 64 GB over-allocated; 16 GB would suffice.
