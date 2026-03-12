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
- Investigate whether workers are truly OOM or hitting wall-time limits (check nanny logs more carefully)
- Consider reducing `cluster_jobs` or further splitting work to avoid large shuffles
- Consider using `repartition` more conservatively or avoiding P2P merges in `_subset_to_pandas`
