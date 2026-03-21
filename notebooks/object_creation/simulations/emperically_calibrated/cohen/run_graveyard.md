# Cohen Power MWU — Run Graveyard

Records of failed/aborted runs and why they died.

---

## Run 1 — Job 6202496 (FAILED: OOM)

| Resource | Value |
|---|---|
| Main job mem | 24G |
| Worker mem | 10G / 3 processes = **3.3G per worker** |
| Worker jobs | 10 |
| Total workers | 30 |

**Failure:** Main job kernel OOM-killed → `DeadKernelError`. Workers hit signal 15 at 95% memory budget (~2.2 GiB unmanaged per process).

---

## Run 2 — Job 6212048 (ABORTED: OOM imminent)

| Resource | Value |
|---|---|
| Main job mem | 32G |
| Worker mem | 15G / 3 processes = **5G per worker** (4.66 GiB effective) |
| Worker jobs | 10 |
| Total workers | 30 |

**Failure:** Unmanaged memory per worker already at ~3.5 GiB / 4.66 GiB effective limit within minutes of startup. `processes=3` splits the 15G allocation three ways — not enough headroom for actual computation on top of the ~3.5 GiB baseline overhead (imports, numpy, etc.). Aborted before OOM kill.

**Lesson:** `processes=3` is the root cause. The per-process limit is too low for the Cohen workload regardless of total job memory. Must use `processes=1` to give each worker the full allocation.

---

## Run 3 — Job 6212578 (FAILED: OOM)

| Resource | Value |
|---|---|
| Main job mem | 32G |
| Worker mem | 15G / **1 process = 15G per worker** |
| Worker jobs | 4 |
| Total workers | 4 |

**Failure:** Main job kernel OOM-killed → `DeadKernelError` + `slurmstepd: Detected 1 oom_kill event in StepId=6212578.batch`. Workers closed gracefully (not the problem). The main process accumulates memory across 100 library rep iterations — 32G is not enough for the orchestration loop itself.

**Lesson:** Main job needs more memory too. The loop building `gt_df` and calling `de_novo_simulation` 100× is the bottleneck in the main process, not the workers.

---

## Run 4 — Job 6215670 (ABORTED: still too slow/heavy)

| Resource | Value |
|---|---|
| Main job mem | 64G |
| Worker mem | 32G / **1 process = 32G per worker** |
| Worker jobs | 2 |
| Total workers | 2 |

**Issue:** Still slow and memory-heavy. Root cause: all 100 sims submitted via `gamut()` before any `save()` — 500 futures in flight simultaneously, main process holding all library DataFrames and sim state at once. Memory and scheduler pressure even with large per-worker allocation.

**Lesson:** Need to batch gamut+save so futures are resolved and freed before submitting the next batch.

---

## Run 5 — Job 6218044 (FAILED: workers at 81% memory)

| Resource | Value |
|---|---|
| Main job mem | 64G |
| Worker mem | 32G / 1 process = **29.8 GiB effective** |
| Worker jobs | 2 |
| Total workers | 2 |
| Batch size | 5 |

**Failure:** Each Cohen sim (all 4 cell types, flatten_overtransfection=True) uses ~24 GiB per worker — 81% of 29.8 GiB effective limit. Workers oscillating at pause/resume threshold; event loop stalls. Dense_simu never hit this because it only simulated the "reference" cell type, not all 4 simultaneously.

**Lesson:** Need 64G workers to give ~40% headroom. Scale to 4 workers + batch_size=10 to restore throughput.

---

## Run 6 — Job 6220049 (FAILED: main OOM at batch 4)

| Resource | Value |
|---|---|
| Main job mem | 64G |
| Worker mem | 64G / **1 process = 64G per worker** |
| Worker jobs | 4 |
| Total workers | 4 |
| Batch size | 10 |

**Failure:** Main OOM-killed after 30 sims (3 batches). Workers closed gracefully — not the problem. Main process accumulated all completed sim objects in `sims` list; each holds ground truth, state, and rehydrated futures. 30+ of these exhausted 64G.

**Lesson:** Don't accumulate sim objects in the main process at all. Record `(location, name)` tuples only; reload from disk for aggregation.

---

## Run 7 — Job 6239474 (ABORTED: too slow, ~7h estimated vs 4h10m limit)

| Resource | Value |
|---|---|
| Main job mem | 64G |
| Worker mem | 64G / **1 process = 64G per worker** |
| Worker jobs | 4 |
| Total workers | 4 |
| Batch size | 10 |

**Issue:** Main process memory fixed (sim objects discarded per batch). Workers stable but slow — 20/100 sims in 1h46m, estimated ~7h total. Also unmanaged worker memory growing to ~48 GiB due to glibc heap fragmentation from large DataFrame allocations in `_simulate_transcription_helper` not being returned to OS after each task.

**Lesson:** glibc holds freed pages in the heap instead of returning to OS. Fix: call `malloc_trim(0)` at end of `_simulate_transcription_helper` to force page return after each task. Also increase `--time` limit.

---

## Run 8 — Job 6242197 (ABORTED: too slow, ~10h50m estimated vs 10h limit)

| Resource | Value |
|---|---|
| Main job mem | 64G |
| Worker mem | 64G / **1 process = 64G per worker** |
| Worker jobs | 4 |
| Total workers | 4 |
| Batch size | 10 |
| Time limit | 10:00:00 |

**Issue:** malloc_trim added. Memory stable (41/50 transfections in batch 2 without warnings). But pace unchanged — ~65 min/batch × 10 batches = ~10h50m, just over the limit. Not a memory problem anymore, just a throughput problem.

**Lesson:** Double the workers to halve the time.

---

## Run 9 — Job 6242773 (FAILED: main OOM at 50/100)

| Resource | Value |
|---|---|
| Main job mem | 64G |
| Worker mem | 64G / **1 process = 64G per worker** |
| Worker jobs | 8 |
| Total workers | 8 |
| Batch size | 10 |
| Time limit | 10:00:00 |

**Failure:** Main OOM-killed at exactly 50/100 sims (halfway through batch 5). Workers stable — malloc_trim working. Main process accumulated enough state across 50 sims to exhaust 64G despite discarding sim objects per batch. Also: workers were using `--nthreads 2` which allocated 2 Dask threads per process (wasting cores and causing CPU contention with the single-core SLURM allocation).

**Lesson:** Main needs 128G. Workers should run with 1 thread (no `worker_extra_args`/`--nthreads`) to match the `cores=1` SLURM allocation.

---

## Run 10 — Job 6253772 (FAILED: TIMEOUT at 10h, 90/100 complete)

| Resource | Value |
|---|---|
| Main job mem | 128G |
| Main job CPUs | 1 |
| Worker mem | 64G / **1 process = 64G per worker** |
| Worker jobs | 12 |
| Total workers | 12 |
| Batch size | 10 |
| Time limit | 10:00:00 |
| Key changes | Main 64→128G; drop `--nthreads 2`; 8→12 workers for throughput. |

**Failure:** TIMEOUT at exactly 10h. Memory was fine — main peaked at ~44GB (well under 128G), malloc_trim working on workers. 90/100 sims fully complete with MWU. 10 stubs had state.parquet but no simulated_scmpra (timeout mid-save) — deleted. Resume run (Run 11) will pick up the 90 complete sims and run only 10 fresh.

---

## Run 11 — Job 6276686 (SUCCEEDED)

| Resource | Value |
|---|---|
| Main job mem | 128G |
| Main job CPUs | 1 |
| Worker mem | 64G / **1 process = 64G per worker** |
| Worker jobs | 12 |
| Total workers | 12 |
| Batch size | 10 |
| Time limit | 14:00:00 |
| Elapsed | 00:34:59 |
| Key changes | Resume from 90 complete sims; ran only 10 fresh; extend to 14h for headroom. |

**Result:** All 100 sims complete with MWU results. Resume logic worked — 90 picked up from disk, 10 fresh sims run in ~35 min.
