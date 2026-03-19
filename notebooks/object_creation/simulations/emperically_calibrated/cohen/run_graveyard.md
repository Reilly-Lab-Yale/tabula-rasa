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

## Run 5 — (current)

| Resource | Value |
|---|---|
| Main job mem | 64G |
| Worker mem | 32G / **1 process = 32G per worker** |
| Worker jobs | 2 |
| Total workers | 2 |
| Batch size | 5 (gamut+save per batch, caps futures at 25 in flight) |
