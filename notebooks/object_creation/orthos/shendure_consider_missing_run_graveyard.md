# Shendure Consider-Missing Run Graveyard

## 2026-03-12 Attempt 1

Driver job:
- Wrapper: [wrap_shend_consider_missing.sh](/home/mcn26/project/tabula_rasa/notebooks/object_creation/orthos/wrap_shend_consider_missing.sh)
- Slurm resources: `2` CPU cores, `24G` RAM, `1-12:10:00` walltime
- Main job id: `1804829`

Dask worker configuration used:
- Defined in [shendure_consider_missing.py](/home/mcn26/project/tabula_rasa/notebooks/object_creation/orthos/shendure_consider_missing.py)
- `SLURMCluster(cores=4, memory="128G", processes=4)`
- `cluster.scale(jobs=4)`
- Effective layout: `16` worker processes total, about `29.8 GiB` memory per worker process

Observed result:
- `ortho_filter` completed
- `consider_missing` modeling ran for roughly `17` hours
- Repeated Dask nanny restarts due to workers exceeding `95%` memory budget
- Final notebook error was `P2PConsistencyError: No active shuffle ... found`
- Root cause was worker loss during shuffle after repeated memory-triggered restarts, not an application-level model exception

Evidence:
- Driver stderr: [slurm-1804829.err](/home/mcn26/project/tabula_rasa/notebooks/object_creation/orthos/slurm-1804829.err)
- Worker logs:
  - [worker_1804850.out](/home/mcn26/project/tabula_rasa/notebooks/object_creation/orthos/worker_1804850.out)
  - [worker_1804851.out](/home/mcn26/project/tabula_rasa/notebooks/object_creation/orthos/worker_1804851.out)
  - [worker_1804852.out](/home/mcn26/project/tabula_rasa/notebooks/object_creation/orthos/worker_1804852.out)
  - [worker_1804853.out](/home/mcn26/project/tabula_rasa/notebooks/object_creation/orthos/worker_1804853.out)

Next configuration to try:
- Reduce per-node Dask process count from `4` to `2`
- Increase cluster size from `4` jobs to `8`
- Keep worker job memory at `128G`
- Add a Dask HTML performance report for memory and shuffle diagnostics

---

## 2026-03-12 Attempt 2

Driver job:
- Wrapper: [wrap_shend_consider_missing.sh](/home/mcn26/project/tabula_rasa/notebooks/object_creation/orthos/wrap_shend_consider_missing.sh)
- Slurm resources: `2` CPU cores, `24G` RAM, `36:10:00` walltime
- Main job id: `1870568`, workers: `1870576`–`1870583`
- Started: `2026-03-12 10:29`, failed: `~15:25–15:29` (~5 hours in)

Dask worker configuration used:
- `SLURMCluster(cores=4, memory="128G", processes=2)`
- `cluster.scale(jobs=8)`
- Effective layout: `16` worker processes total, `~59.6 GiB` memory per worker process

Observed result:
- `ortho_filter` removed 4 combinations involving `'reference'`; dropped `641 of 2103` (cell_type, cre_id) combos with fewer than 3 nonzero entries
- Fitting ran for ~5 hours before failure
- Workers showed severe event loop unresponsiveness (up to 35+ seconds), indicating large data movement blocking the GIL
- One unconverged model warning: `Txndc12_chr4_7971`
- `_smart_matrix` task marked as **failed because 4 workers died** during shuffle
- Cascade of `P2PConsistencyError: No active shuffle ... found` errors in `_subset_to_pandas` (`scMPRAforge/core.py:2030`)
- Root cause: workers likely died from OOM or timeout mid P2P shuffle; the shuffle state became inconsistent, causing `RuntimeError: Set changed size during iteration` and `P2PIllegalStateError` in the Dask scheduler shuffle plugin
- HTML performance report was **not written** — process was stuck in a C-level Dask future wait and could not be gracefully interrupted via SIGINT; job was hard-cancelled

Fix applied after this run:
- Added SIGTERM handler to `shendure_consider_missing.py` so that `scancel` triggers a graceful Python shutdown (performance report + client/cluster close)

Next configuration to try:
- Memory pressure on workers is the most likely culprit based on dashboard observation
- Reduce `processes` from `2` to `1` so each worker gets the full `128G` instead of ~60G
- Keep `cluster_jobs=8` so still 8 workers total (down from 16), each with 128G headroom

---

## 2026-03-12 Attempt 3

Driver job:
- Wrapper: [wrap_shend_consider_missing.sh](/home/mcn26/project/tabula_rasa/notebooks/object_creation/orthos/wrap_shend_consider_missing.sh)
- Slurm resources: `2` CPU cores, `24G` RAM, `36:10:00` walltime
- Main job id: `1881657`, workers: `1881683`–`1881690`
- Started: `2026-03-12 16:53`, failed: `2026-03-13 ~08:41` (~15.7 hours in)

Dask worker configuration used:
- `SLURMCluster(cores=4, memory="128G", processes=1)`
- `cluster.scale(jobs=8)`
- Effective layout: `8` worker processes total, `~119.21 GiB` memory per worker process, `4` threads per worker

Observed result:
- `ortho_filter` removed 4 combinations involving `'reference'`; dropped `641 of 2103` (cell_type, cre_id) combos with fewer than 3 nonzero entries
- Fitting ran for ~15 hours before failure
- Workers periodically hit `80%` memory budget (paused) then recovered, but eventually exceeded `95%` → nanny killed and restarted
- "Unmanaged memory" warnings at 86–89 GiB on multiple workers — memory not released back to OS between tasks (likely Python/numpy/TF heap retention)
- At least two workers hit the SLURM OOM killer directly: `slurmstepd: error: Detected 2 oom_kill events` (process exceeded SLURM's 128G hard limit, not just the Dask nanny threshold)
- Multiple `_smart_matrix` tasks marked as **failed because 4 workers died** (at ~02:00, ~04:32, ~04:33, ~06:16)
- Final `KilledWorker` raised during `.save()` → `flattened_copy()` → `future.result()` attempting to materialize a failed `_smart_matrix` future
- SIGTERM handler worked: graceful shutdown printed `! Done, shutting down`; performance report was written

Root cause:
- With `4` threads per worker, up to 4 tasks run concurrently on a single 128 G node
- Each task with `consider_missing=True` materializes a cartesian-product-expanded DataFrame, then builds dense design matrices in `_smart_matrix` (notably the dense `Z = C(rep_id)-1` indicator matrix)
- The memory guard in `_inflate_missing_split_level` only estimates raw DataFrame size and caps at 100 GB; it does not account for downstream design matrix construction
- Accumulated per-task memory from 4 concurrent tasks routinely exceeded 128 G → SLURM OOM kills → after 4 attempts, task permanently failed

Fix applied after this run:
- `worker_cores` reduced from `4` to `1` (single-threaded workers: exactly 1 task at a time per worker)
- `cluster_jobs` increased from `8` to `16` to compensate for reduced per-worker parallelism

Next configuration to try:
- `SLURMCluster(cores=1, memory="128G", processes=1)`, `cluster.scale(jobs=16)`
- Effective layout: 16 single-threaded workers, each with 128 G exclusive to one task at a time
