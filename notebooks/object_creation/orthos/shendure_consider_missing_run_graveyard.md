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

---

## 2026-03-13 Attempt 4

Driver job:
- Wrapper: [wrap_shend_consider_missing.sh](/home/mcn26/project/tabula_rasa/notebooks/object_creation/orthos/wrap_shend_consider_missing.sh)
- Slurm resources: `2` CPU cores, `24G` RAM, `36:10:00` walltime
- Main job id: `1910351`, workers: `1910362`–`1910376` (15 running; `1910377` pending due to `QOSMaxMemoryPerUser`)
- Started: `2026-03-13 10:26`, failed: `2026-03-13 ~19:58` (~9.5 hours in)

Dask worker configuration used:
- `SLURMCluster(cores=1, memory="128G", processes=1)`
- `cluster.scale(jobs=16)` (only 15 launched; 1 blocked by QOS memory limit)
- Effective layout: `15` single-threaded workers, `~119.21 GiB` memory per worker

Observed result:
- `ortho_filter` removed 4 combinations involving `'reference'`; dropped `641 of 2103` (cell_type, cre_id) combos
- Fitting ran for ~9.5 hours before failure
- Workers again accumulated high unmanaged memory (~92 GiB), hit 80–95% thresholds, and were restarted by the Dask nanny (SIGKILL via signal 15)
- Restarted workers disrupted an in-progress P2P shuffle → cascade of `P2PConsistencyError: No active shuffle ... found`
- Final error: `MemoryError: Unable to allocate 560 MiB for an array with shape (73424955,) and dtype int64` in `_smart_matrix` → `formulaic` → `_get_columns_for_term`
- Root cause: `formulaic` allocates a fully dense `int64` matrix before the `.astype(pd.SparseDtype(...))` conversion in `_smart_matrix`; for large cell types (73M rows) this intermediate allocation OOM'd workers already holding ~92 GiB of unmanaged data
- SIGTERM handler worked: graceful shutdown printed `! Done, shutting down`; performance report written to `shendure_ortho_consider_missing_20260310_dask_performance_report.html`

Fix applied after this run:
- `_smart_matrix`: replaced `get_model_matrix(..., output='pandas')` + `.astype(SparseDtype)` with native sparse construction via `get_model_matrix(..., output='sparse')` — formulaic builds the CSC matrix directly from category indices with no dense intermediate
- `y` (regressand) extracted directly as `data[["umis_mpra_bc"]]` (single column, no formula call needed)
- `_tensorzinb_fit`: densify `nb_regressors` and `zi_regressors` once at top of function via `.toarray()` (scipy sparse API), since TensorZINB requires dense numpy throughout

Next configuration to try:
- Same worker layout (`cores=1, memory="128G", processes=1`, `scale(jobs=16)`)
- Sparse design matrix construction should eliminate the dense intermediate OOM in `_smart_matrix`

---

## 2026-03-15 Attempt 5

Driver job:
- Wrapper: [wrap_shend_consider_missing.sh](/home/mcn26/project/tabula_rasa/notebooks/object_creation/orthos/wrap_shend_consider_missing.sh)
- Slurm resources: `2` CPU cores, `24G` RAM, `36:10:00` walltime
- Main job id: `2155413`, workers: `2155415`–`2155429`
- Started: `2026-03-15 18:02`, failed: `2026-03-16 ~04:51` (~10.8 hours in)

Dask worker configuration used:
- `SLURMCluster(cores=1, memory="128G", processes=1)`
- `cluster.scale(jobs=16)` (15 running)
- Effective layout: `15` single-threaded workers, `~119.21 GiB` memory per worker

Observed result:
- `ortho_filter` removed 4 combinations involving `'reference'`; dropped `641 of 2103` (cell_type, cre_id) combos
- **Sparse `_smart_matrix` fix worked** — no more OOM during design matrix construction; fitting ran successfully for hours
- First `_label_tensorzinb_regressors` failure at `18:58` (worker 2155415), then many more across all workers
- All failures identical: `AttributeError: 'csc_matrix' object has no attribute 'columns'` at `core.py:2110`
- Root cause: `_label_tensorzinb_regressors` and `_matricies_to_order` accessed `.columns` on the design matrices, which is a pandas DataFrame attribute; the sparse `csc_matrix` returned by formulaic's `output='sparse'` does not have `.columns`
- 1 worker also died from unrelated cause (`read-csv` task failed because 4 workers died)
- **No memory warnings, no nanny restarts, no OOM kills** — major improvement over previous runs
- SIGTERM handler worked; performance report written

Resource usage (note: run failed at labeling stage, so later stages like `extract_params`/`save` materialization did not fully execute — actual peak for a successful run may be higher):
- All workers requested `120 GiB`; actual MaxRSS ranged `8.0–31.8 GiB` (4–25% utilization)
- Highest: worker `2155425` at `31.8 GiB` (30.3 GB per jobstats); most workers under `18 GiB`
- CPU utilization `71–79%` across workers
- Memory could potentially be trimmed to `64G` (~2x observed max) if a successful run confirms similar usage — but keeping `128G` for now given >4 hr run times and risk of spiky allocations in stages that didn't run

Fix applied after this run:
- `_smart_matrix` now stores column names as plain lists in the return dict: `nb_regressor_names` and `zi_regressor_names`, extracted from `ModelMatrix.model_spec.column_names` before the formulaic wrapper is stripped by Dask serialization
- `_label_tensorzinb_regressors` and `_matricies_to_order` updated to read names from the dict keys, with fallback to `.columns` for backwards compatibility

Next configuration to try:
- Same worker layout
- Column name fix should resolve the `AttributeError`

---

## 2026-03-16 Attempt 6

Driver job:
- Main job id: `2246714`, workers: `2246723`–`2246737`
- Started: `2026-03-16 11:46`, still running at `2026-03-17 09:09` when last `_tensorzinb_fit` failure logged

Dask worker configuration used:
- `SLURMCluster(cores=1, memory="128G", processes=1)`
- `cluster.scale(jobs=16)` (15 running)
- Effective layout: `15` single-threaded workers, `~119.21 GiB` memory per worker

Observed result:
- `ortho_filter` removed 4 combinations involving `'reference'`; dropped `641 of 2103` (cell_type, cre_id) combos
- Column name fix worked — no `AttributeError` in `_label_tensorzinb_regressors`
- Fitting ran for many hours before failure; `_tensorzinb_fit` tasks permanently failed at `04:04`, `05:25`, `05:47`, `07:16`, `08:16`, `09:09` (each after 4 worker deaths)
- Root cause: TF/Keras heap memory accumulates across tasks and is not released back to the OS between fits; after many tasks, workers reach ~87–90 GiB of unmanaged memory; `nb_X.toarray()` then tries to allocate the dense design matrix on top of that, tipping the worker over the 95% nanny budget → SIGKILL → restart → shuffle/task failure cascade
- Exact error in `_tensorzinb_fit` at line 1968: `nb_X.toarray()` → `np.zeros(self.shape, ...)` → OOM

Fix applied after this run:
- Added `K.clear_session()` at the end of `_tensorzinb_fit` (after results extracted, before return) to explicitly release the TF/Keras session memory after each model fit
- Identified fundamental problem: `consider_missing=True` inflates datasets to 12M–128M rows/cell type; NB design matrix (×211 CRE columns, float64) = 20–216 GB — 7 of 10 cell types exceed the 128 GB worker limit even with a single task; `nb_X.toarray()` in `_tensorzinb_fit` is therefore infeasible for most cell types regardless of session clearing

Next configuration to try:
- Upgrade TensorZINB to accept sparse matrices natively (eliminate `nb_X.toarray()` altogether)
- User upgraded TensorZINB in new conda env `tz` (TF 2.20, tf_keras 2.20.1); `.toarray()` calls already commented out in `core.py` in anticipation

---

## 2026-03-17 Pre-run notes / fixes applied

Upgraded TensorZINB (`tz` conda env) assessed and two compatibility fixes applied before Attempt 7:

**TensorZINB upgrade (sparse support):**
- Added `SparseDense` custom layer: replaces `Dense` when input is sparse; uses `tf.sparse.sparse_dense_matmul` to avoid densification in the Keras graph
- Sparse input detection flags (`_exog_is_sparse` etc.) in `__init__`; `_matrix_rank` / `_sparse_std` / `_sparse_col_mean` / `_sparse_dot` helpers replace numpy dense equivalents
- `fit()`: `Input(sparse=True)` for sparse inputs; training via `tf.data.Dataset.from_tensors()` + `model.fit(dataset)`; LL retrieval via `SparseTensorValue` feeds to `K.function`
- `_poisson_init_each`: skips statsmodels (can't handle sparse); uses intercept-only fallback
- Migrated from `keras` (TF1-bundled) to `tf_keras` (standalone Keras 2 on TF2); keras 3.x present in env but not used by TensorZINB

**Fix 1 — `K.clear_session()` backend mismatch (`core.py`):**
- `tz` env has keras 3.13.2 alongside tf_keras 2.20.1; `import keras.backend as K` would call Keras 3's session clear, which does not clear the TF1 graph/session
- Changed to `import tf_keras.backend as K` so `clear_session()` actually releases TF memory

**Fix 2 — conda env (`wrap_shend_consider_missing.sh`):**
- Changed `conda activate env_tensorzinb` → `conda activate tz`
- Dask workers inherit the driver's Python binary; driver must be in `tz` for workers to use the sparse TensorZINB

