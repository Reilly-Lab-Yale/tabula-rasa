# Attempt 3 HTML Performance Report Analysis

## Analysis: Attempt 3 HTML Performance Report vs. Graveyard

**Short answer: OOM was caused by a combination of unequal distribution AND a design flaw in the 4-thread configuration. It was NOT uniform cluster saturation — but neither was it purely one rogue worker. The heaviest-loaded nodes were the ones that died most often, and one node escaped entirely.**

---

### 1. Task distribution was highly skewed

From the 96,042 task-worker assignments recorded in the task stream:

| Node IP | Tasks | Share | Restarts (new worker addresses) |
|---|---|---|---|
| 10.178.138.46 | 21,159 | **22.0%** | **4** |
| 10.178.138.45 | 19,617 | **20.4%** | 1 |
| 10.178.138.8  | 15,707 | 16.4% | 1 |
| 10.178.138.5  | 13,829 | 14.4% | 3 |
| 10.178.138.48 | 7,761  | 8.1%  | 2 |
| 10.178.138.30 | 7,616  | 7.9%  | **4** |
| 10.178.138.7  | 5,236  | 5.5%  | **4** |
| 10.178.138.23 | 5,117  | **5.3%** | **0** ← never crashed |

The top node handled **4.1× more tasks** than the bottom node. The one node that never crashed was the one with the lightest load. Clear correlation.

Note: `.7` is the scheduler node (it also ran a worker, but scheduler overhead competed with it — explains its low task count *and* high restart count).

---

### 2. The 4-thread concurrent execution was the amplifying mechanism

The key data from the task stream:
- **Max single task duration: 195 minutes** (not 195 ms — 195 *minutes*, i.e. ~11.7M ms)
- Many tasks ran **147–156 minutes each**
- These are the `_smart_matrix` / `_inflate_missing` tasks materializing cartesian-expanded DataFrames under `consider_missing=True`

With `cores=4`, each 128 GiB node ran **4 such tasks concurrently**. If each task consumed 30–40 GiB (which the memory guard doesn't prevent — it only estimates the raw DataFrame, not downstream design matrix construction), you get 120–160 GiB peak → SLURM hard OOM kill.

The heaviest node (.46, 22% of tasks) ran 4 concurrent threads of these monsters throughout the run → 4 restarts. The lightest node (.23, 5.3%) got assigned fewer and/or shorter tasks → 0 restarts.

---

### 3. The driver process shows the final failure

The driver process memory trace (1,000 samples over 15.7 hours):
- Flat at **~82 MiB** for almost the entire run
- Then spikes: **302 → 398 → 741 MiB** in the final 3 samples

This matches the graveyard: the final `KilledWorker` exception happened during `.save()` → `flattened_copy()` → `future.result()`, which tried to pull failed `_smart_matrix` futures into the driver and caused driver-side memory growth before crashing.

---

### 4. Implications for Attempt 4 (currently running)

The fix (`cores=1, processes=1`, 16 workers) should directly resolve the concurrent-task OOM mechanism. With 1 thread per 128 GiB node, only 1 task runs at a time — even a 40 GiB task can't cause OOM alone on 119 GiB headroom.

The unequal distribution you can see in Attempt 4's Dask poll right now (`processing: 212, waiting: 314`) should be less dangerous now — an overloaded worker just queues tasks, it doesn't run them concurrently. The 15th worker being stuck pending (`QOSMaxMemoryPerUser`) is a mild loss of throughput, not a correctness/memory risk.
