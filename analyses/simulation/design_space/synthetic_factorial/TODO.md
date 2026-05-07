# synthetic_factorial: TODOs

## Adapt-race FutureCancelledError on first sample of each slice

**Symptom**: When a slice driver starts, `cluster.adapt(min=1, max=N)`
begins spawning workers via sbatch, but they take ~30s+ to register with
the scheduler. If the driver immediately submits the first sample's
transfection / transcription tasks, those futures may target workers
that are still in the middle of connecting -- the scheduler then reaps
them mid-spawn and the helper task surfaces as
`FutureCancelledError: _simulate_transcription_helper-...
cancelled for reason: already forgotten`.

The slice driver's per-sample try/except catches this so the run
continues, but it costs us one sample per affected slice. Empirically:

- topup B (2026-05-05): 1 of 100 samples lost (t0033)
- topup3 (2026-05-06): 3 of 120 samples lost in the first ~15 min
  (u0000, u0003, u0008 -- the first sample of slices 0, 3, 8)

All of these are "first sample" failures, all happen during adapt()'s
ramp-up window.

## Fix options

### Option A (recommended): `client.wait_for_workers(N, timeout=...)`

In `_make_slurm_client` after `cluster.adapt(...)`, block until at least
some workers have actually registered:

```python
print("waiting for >=8 workers to register before starting work...", flush=True)
client.wait_for_workers(min(8, N_WORKERS), timeout="600s")
```

`wait_for_workers` is the canonical dask pattern for this. Picking 8 (or
~15% of N_WORKERS) gives the scheduler a small, stable pool to start
work on; adapt() continues to scale up for the actual workload as more
tasks arrive.

### Option B: retry once on FutureCancelledError in `_run_one_sample`

Catch `FutureCancelledError` in the rep loop and retry the same rep
once before giving up. Less elegant; adds latency on real errors. A's
preferred unless A turns out flaky.

### Option C: explicit `cluster.scale(N_initial)` then switch to adapt

`cluster.scale(8)` followed by `cluster.adapt(1, N_WORKERS)` ensures a
synchronous initial pool. The original `scale()` is what triggered the
2026-05-05 sbatch rate-limit storm, so the count must stay small (~8).

## Other small TODOs

- `disk per sample` in topup3 is ~18 GB (vs 9 GB in topup) due to
  high-moi samples generating more transfection events. Worth
  documenting in resource projections.
- After topup3 completes, move full + topup + topup3 sim trees to
  cold storage via globus to free scratch.
- 5-metric replot regenerates the marginals SVGs. Note that the legacy
  `marginals.svg` (P@FC=2) is the saturated-at-edges metric; the
  primary metric is `power_auc_1to3`. Document this somewhere
  user-facing if anyone other than us looks.
