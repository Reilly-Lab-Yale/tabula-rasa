# synthetic_factorial: TODOs

## Heavy-tail OOM during transcription (resolved 2026-05-10)

**Symptom**: `FutureCancelledError: _simulate_transcription_helper-...
cancelled for reason: already forgotten` raised from `sim.save()` (core.py
~line 7872, the `[f.result() for f in fut]` await).

**Original misdiagnosis**: thought to be an adapt() ramp-up race ("first
sample of each slice fails before workers register"). Counter-evidence
from the 2026-05-09 retry: failures recurred on samples started 1+ hours
after workers were stable, and slice 0 alone lost 4 samples spread across
9 hours, not just the first.

**Real cause**: heavy samples (ncells > 35k AND moi > 200) need >40 GiB
of unmanaged worker memory during transcription. With workers allocated
at `worker_mem="48G"` (~44.7 GiB after dask overhead), the dask nanny
hit its 95% memory budget mid-task, SIGKILL'd the worker process, and
restarted it. Any future running on the killed worker became "already
forgotten" from the scheduler's perspective.

Worker logs make this visible: workers cycled through "Unmanaged memory:
32+ GiB" warnings, "Worker is at 80% memory usage. Pausing worker.",
then eventually "exceeded 95% memory budget. Restarting...". Across the
2026-05-05 topup3 run, 10 of 17 workers experienced these kills with
MaxRSS pinned at ~47 GiB (the alloc cap).

The 14 originally-lost samples (1 in topup, 13 in topup3) all sat in
the heavy-tail corner: ncells in [31771, 49382] and moi in [171, 340].

**Fix**: bumped `worker_mem` for `topup` and `topup3` modes from `48G`
to `128G` (synthetic_factorial.py:179-180). The 2026-05-10 resim of all
14 lost samples completed cleanly: 0 OOM kills, 0 FutureCancelledError,
zero worker memory warnings. Real peak unmanaged memory was 50-60 GiB
on the heaviest samples; 96G would likely also suffice but 128G is safe
headroom. Full sweep stays at 48G because the LHS box doesn't reach the
takeshi corner.

## Union-sweep rerun playbook (lessons from 2026-05-11 attempt)

**Binding constraint: YCRC per-hour sbatch rate (~200/hr).** Not compute,
not memory, not wall time. Every parallelism decision must pencil out
against this budget.

Per-driver steady-state sbatch cost:
- 5 initial worker sbatches at startup (`adapt(min=1, max=5)`)
- ~5 replacement sbatches every 8h (worker walltime expiry)
- ~1-2 transient replacements per hour (worker preempt, OOM, etc.)
- ~ 1.5-2 sbatches/hour/driver steady state, but BURSTY

Fatal failure mode: `dask_jobqueue.SLURMCluster.adapt()` calls
`_correct_state_internal` which raises `RuntimeError` on any single failed
sbatch and does not retry. Once raised, the cluster's worker pool stops
replenishing. Subsequent samples in that driver all fail with
"TimeoutError: No valid workers found" until the driver is killed.

Three mitigations to apply on rerun (most-to-least impactful):

1. **Stagger driver startup.** Submit the slice array with `--array=...%10`
   (or %5) to cap concurrent driver starts. Eliminates the initial-burst
   cascade where 50 drivers all submit 5 workers simultaneously
   (~250 sbatches in seconds, ~15,000/hr instantaneous rate).

2. **Lower max workers per cluster.** Set `n_workers=2` in MODES[mode]
   instead of 5. `n_sims=5` parallelism within a rep is already enough
   to keep two workers busy. Halves per-driver sbatch pressure across
   startup, replacement, and walltime-expiry bursts.

3. **Retry wrapper on sbatch failures.** Monkeypatch
   `dask_jobqueue.SLURMCluster._submit_job` (or whatever the current
   sbatch call is in the installed version) with exponential backoff on
   rate-limit errors. Single transient failure no longer kills the
   cluster's adapt() permanently. ~30 lines of code in core.py or a
   small monkeypatch module.

## Tiny tim top-up after union sweep

The union sweep (array 11454632, 5000 samples) drops occasional samples to
rare FutureCancelledErrors when dask's adapt() tries to top up a slice's
worker pool and hits the YCRC sbatch rate limit. Observed rate: ~0.1-2%
depending on how often heavy-tail samples trigger worker replacement.

After the main sweep finishes, identify any `v*` samples without a
`cached_results.parquet` and submit a small fixup array over just those
slice indices. Sample-grain is fine for the fixup since the count is
small (~5-100 expected) and worker-startup overhead is negligible at
that scale.

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
