# Graveyard — shendure power simulation attempts

---

## Attempt 1 — 2026-03-13 (jobs 5901175, 5901675–5901785)

**Setup:**
- `processes=1, cores=3, threads_per_worker=3` (default)
- 10 SLURM worker jobs → 10 Dask workers
- `partition=day`

**What went wrong:**
1. Wrong ortho loaded (`ortho_primordial_v4`) — predates the dask DataFrame backend refactor for `scMPRA_data`. Job started but workers failed immediately.
2. Root error: `ValueError('Unsupported table_type for parquet save: None')` in `_simulate_transcription_helper` — `scMPRA_data.table_type` was never set after construction during simulation, so `to_parquet` raised on every task.

**Fix applied:**
- Changed ortho to `shendure_ortho_20260306`
- Added `scd.table_type = "mpra_umiwise"` in `_simulate_transcription_helper` (`core.py`)

---

## Attempt 2 — 2026-03-18 (jobs 6154490 + workers, then 6154774 + workers)

**Setup:**
- `processes=1, cores=3, threads_per_worker=3` (default)
- 10 SLURM worker jobs → 10 Dask workers
- `partition=priority, account=prio_skr2`
- `shendure_ortho_20260306`

**What went wrong:**
- Job ran for ~1.5 hours before stalling completely.
- Dask scheduler pegged at ~100% CPU; all 10 workers at ~1–2% CPU (essentially idle).
- `json/counts.json` polled every 10s: `processing`, `waiting`, `waiting_data`, and `memory` counts were completely static — confirmed deadlock, not just slowness.
- Root cause: `processes=1` per SLURM job means a single Python process handles all tasks. With only 3 threads (default `threads_per_worker = cores // processes = 3`) and a large inter-dependent task graph (~10,857 tasks, `desired_workers=178`), intra-worker task dependencies caused a deadlock. The Dask scheduler spun at 100% CPU trying to resolve a frozen graph.
- Secondary inefficiency: with `cores=3, processes=1`, only ~1 core was ever active (GIL), giving ~33% CPU utilisation per worker SLURM job.

**Fix being tried (Attempt 3):**
- `processes=3, cores=3, threads_per_worker=2`
- 10 jobs × 3 workers per job = 30 Dask workers
- Each worker is an independent Python process on its own core → no GIL contention, full core utilisation
- `threads_per_worker=2` gives each worker headroom to avoid intra-worker deadlock (same class of fix that resolved a prior deadlock in this codebase)

---

## Attempt 3 — 2026-03-18 (failed immediately)

**Setup:**
- `processes=3, cores=3, threads_per_worker=2`
- 10 SLURM jobs × 3 workers = 30 Dask workers

**What went wrong:**
- `threads_per_worker` is not a valid `SLURMCluster` parameter (it's `LocalCluster`-only).
- `ValueError: Got unexpected keyword argument 'threads_per_worker'` on cluster init.

**Fix:** pass via `worker_extra_args=["--nthreads", "2"]` instead.

---

## Attempt 4 — 2026-03-18 (job 6160665, workers 6160719–6160728)

**Setup:**
- `processes=3, cores=3, worker_extra_args=["--nthreads", "2"]`
- 10 SLURM jobs × 3 workers = 30 Dask workers
- `partition=priority, account=prio_skr2`

**What went wrong:**
- Workers ran (30 connected, ~110–143% CPU per node), tasks started, 635 tasks released — then deadlock again.
- `json/counts.json` frozen: `processing=385, released=635, waiting=349` static across multiple polls.
- Root cause identified by querying `client.processing()`: **354 of 385 processing tasks were `_mwu_helper`**, with the remaining 31 being Dask sub-tasks (`read_parquet`, `astype`, `fillna`, `repartitiontofewer`).
- Deadlock mechanism: `_mwu_helper` runs on a worker thread and calls `scMPRA_data.from_parquet()` → `dd.read_parquet()` (lazy Dask DataFrame). Then `_mwu_make_bundle()` calls `_to_pandas_df()` which calls `.compute()`, submitting sub-tasks back to the scheduler. Those sub-tasks need a free worker thread, but **all 60 threads are occupied by blocked `_mwu_helper` calls waiting for their own `.compute()` to finish**. Classic thread-starvation deadlock — nothing to do with nthreads count.
- Increasing `--nthreads` would not fix this; the problem is structural: Dask tasks must not call `.compute()` inside a worker.

**Fix being tried (Attempt 5):**
- Patch `_mwu_helper` to use `pd.read_parquet()` directly instead of `scMPRA_data.from_parquet()`, bypassing Dask I/O entirely inside the worker. The data is small enough per sim that pandas is fine.

---

## Attempt 5 — (pending)
