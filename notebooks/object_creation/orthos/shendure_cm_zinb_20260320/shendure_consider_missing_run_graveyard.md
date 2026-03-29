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

---

## 2026-03-17 Attempt 7 — failed (2GB TF graph-constant limit)

Root cause: `core.py` calls `tf.compat.v1.disable_eager_execution()` at module import time. In TF1 graph mode, `tf.SparseTensor(numpy_indices, numpy_values, shape)` tries to serialize the indices array as a protobuf graph constant (limit: 2GB). For a 128M-row sparse matrix with ~5 nonzeros/row, the COO indices array alone is ~10GB → `ValueError: Cannot create a tensor proto whose content is larger than 2GB`.

Full error path:
```
_tensorzinb_fit → TensorZINB.fit() → _scipy_to_tf_sparse()
→ tf.SparseTensor(indices=...) → ops.convert_to_tensor(numpy_indices)
→ _create_graph_constant() → make_tensor_proto()
→ ValueError: Cannot create a tensor proto whose content is larger than 2GB.
```

Fixes applied:
1. **`core.py` — defer `disable_eager_execution()`**: Removed from module level; added at the top of `_estimate_cov_se()` and `_hessian_se_graph()` (the only two functions requiring TF1 graph mode). TensorZINB sparse fits now run in TF2 eager mode; Wald functions enable graph mode lazily when called (after all fitting is done, so no interleaving risk).
2. **TensorZINB `fit()` — remove `disable_eager_execution()`**: Was being called inside `fit()`, which conflicted with the deferred approach in core.py.
3. **TensorZINB `fit()` — replace `K.function` LL computation**: `K.function([..., zinb.y], [zinb.llf])` requires graph-mode symbolic tensors; in eager mode `zinb.y`/`zinb.llf` are concrete values from the last training step, not traceable nodes. Replaced with: `preds = model(ll_inputs, training=False); _ = zinb.loss(endog_t, preds); llft = zinb.llf.numpy()`.
4. **`core.py` `_extract_mu` — scipy sparse API**: `X.sparse.to_coo()` is a pandas sparse accessor; `X` is now `csc_matrix`. Fixed: detect sparse, do `X.tocsr() @ w` directly, then reconstruct pandas sparse DataFrame via `pd.DataFrame.sparse.from_spmatrix(X, columns=design_matrix["nb_regressor_names"])` for `undo_one_hot_encoding`.
5. **`core.py` `_extract_zi` — scipy sparse API**: `Z.to_numpy()` is a pandas method; `Z` is now `csc_matrix`. Fixed: detect sparse, do `Z @ x_pi` directly, reconstruct pandas sparse DataFrame for label extraction.
6. **`tz` env — bokeh/distributed `TemplateNotFound`**: `distributed` `performance_report.html` extends `"file.html"` but bokeh 3.x renamed it to `file.html.jinja`. The broken `__exit__` masked the real error (`ValueError`) in the driver logs. Fixed: symlink `file.html → file.html.jinja` in bokeh templates dir.

Note: fixes 2 and 3 above were applied directly to the installed TensorZINB package (`site-packages/tensorzinb/tensorzinb.py`) rather than the source repo — this was a mistake. They have since been reverted to source, and the correct changes will be released as a new version of TensorZINB (see Attempt 8 notes).

---

## 2026-03-18 Attempt 8 — failed (`AttributeError: 'KerasTensor' object has no attribute 'numpy'`)

Root cause: The replacement LL computation introduced in Attempt 7 fix #3 was itself broken. After training, `zinb.loss(endog_t, preds)` was called expecting it to set `zinb.llf` to a concrete eager tensor. But `zinb.pi` and `zinb.log_theta` are `KerasTensor`s — they are set during Keras functional-API model construction (when Keras traces `zinb.loss()` with symbolic inputs to build the computation graph). Calling `zinb.loss()` with concrete eager inputs after training mixes those stored `KerasTensor`s into the computation, so `zinb.llf` remains symbolic and `.numpy()` raises `AttributeError`.

Full error path:
```
primordial.save() → experiment_model.save() → flattened_copy()
→ future.result() → _tensorzinb_fit (worker)
→ TensorZINB.fit() line 521: llft = zinb.llf.numpy()
→ AttributeError: 'KerasTensor' object has no attribute 'numpy'
```

Workers also showed a related error during model construction:
```
TypeError: You are passing KerasTensor(...), an intermediate TF-Keras symbolic
input/output, to a TF API that does not allow registering custom dispatchers,
such as `tf.cond`, `tf.function`, gradient tapes, or `tf.map_fn`.
```

Root cause analysis of the full LL chain:
- Original code: `disable_eager_execution()` + `K.function([..., zinb.y], [zinb.llf])` — works in TF1 graph mode
- We removed `disable_eager_execution()` (required for sparse SparseTensor inputs to avoid 2GB limit) → `K.function` breaks because `zinb.y`/`zinb.llf` are no longer symbolic graph nodes
- Replacement used `zinb.loss()` side-effect → breaks because `zinb.log_theta`/`zinb.pi` are KerasTensors from functional API construction

Correct fix (to be released in new TensorZINB version): compute LL entirely in numpy from `weights_dict` after training, without invoking `zinb.loss()` at all. Move `model.get_weights()` before the LL block and mirror the `ZINBLogLik.loss()` formula in numpy using `scipy.special.gammaln` and `np.logaddexp`.

All TensorZINB changes (sparse support + numpy LL) reverted from `site-packages` back to source repo. A new TensorZINB release incorporating all changes will be published and reinstalled into the `tz` env before the next attempt.

Changes committed in this session (all in `scMPRAforge/core.py`):
1. **Defer `disable_eager_execution()`**: Removed from module-level TF import; added at the top of `_estimate_cov_se()` and `_hessian_se_graph()`. These are the only two callers that require TF1 graph mode (they use `tf.compat.v1.placeholder` / `Session`), and they are only called after all fitting is complete.
2. **`_extract_mu` — scipy sparse API**: Detect `scipy.sparse.issparse(X)`; if sparse, compute `X.tocsr() @ w` directly then reconstruct pandas sparse DataFrame via `pd.DataFrame.sparse.from_spmatrix()` for `undo_one_hot_encoding`.
3. **`_extract_zi` — scipy sparse API**: Detect `scipy.sparse.issparse(Z)`; if sparse, compute `Z @ x_pi` directly then reconstruct pandas sparse DataFrame for label extraction.

---

## 2026-03-20 Attempt 9 — failed (36-hour wall-limit timeout, no crash)

No errors. The job ran cleanly from 2026-03-18 15:08 to 2026-03-20 03:17 (36 hours) and was cancelled by SLURM at the wall-limit. All 16 workers were still alive at cancellation — no OOM, no Python exception, no TF error.

Worker logs show workers entered `_mom_from_training_data` within minutes of starting (visible FutureWarning about `working_nb["beta"].loc["reference"]`), then no further logged activity until the cancellation message. No model summaries, no loss values, no fitting completion messages appeared in any worker log.

Based on observation before the job was cancelled: by-CRE models completed, but by-cell-type models were still running and taking much longer than expected. With `consider_missing=True`, cell-type model datasets are 12M–128M rows each, making both design matrix construction and TensorZINB fitting substantially slower than by-CRE models.

Contributing factor: this run predated the `standard_fit` sort-by-size change — cell-type models were submitted in arbitrary order rather than largest-first. Large cell types were not necessarily the first to start fitting.

Fixes applied before next attempt:
1. **`standard_fit` — sort levels by descending size**: Levels are now sorted by row count (descending) before futures are submitted, so the largest (slowest) models start first. Failures surface early (fail-fast), and the scheduler picks up high-priority work sooner.
2. **TensorZINB — new release installed**: The `tz` env now has the updated TensorZINB with `TensorZINBTrainingModel` (`tf.GradientTape`-based `train_step`), `ZINBLogLik._loss_components()` as a `@staticmethod` called with concrete tensors post-training (fixes the KerasTensor LL bug), `SparseDense`, and `run_eagerly=True` compilation.

---

## 2026-03-20 Attempt 10 — failed (Dask re-entrant deadlock in `_subset_to_pandas`)

Driver job:
- Main job id: `2559822`, workers: `2559825`–`2559840`
- Started: `2026-03-20`, killed manually after deadlock confirmed

Root cause: Dask re-entrant deadlock. All 16 workers were executing `_subset_to_pandas` tasks. Inside `_subset_to_pandas`, `ddf.compute()` was called to materialise the Dask DataFrame. `ddf.compute()` attempts to dispatch DDF partition sub-tasks back to the distributed scheduler — but all worker slots were already occupied by the outer `_subset_to_pandas` tasks, so the sub-tasks could never be scheduled. Result: all workers blocked indefinitely waiting for sub-tasks that could never run.

Evidence:
- Dashboard monitoring showed 138 active `_subset_to_pandas` tasks and 0 task completions over 2+ minutes
- No model summaries, no loss values, no fitting completion messages in any worker log

Fix applied:
- In `standard_fit` (`core.py`), replaced the `subset_inputs` dict + `{t: client.submit(_subset_to_pandas, subset_inputs[t])}` comprehension with a per-level loop
- For `use_missing=True` path: call `client.compute(ddf)` on the driver to dispatch DDF partition tasks directly into the scheduler's task graph, then `client.submit(_prepare_subset_for_modeling, pandas_future)` to prepare the result — no worker holds a slot while waiting for sub-tasks
- For `use_missing=False` path: unchanged (`client.submit(_subset_to_pandas, raw[raw[split] == t])` — subset is already a pandas slice, no internal compute)

---

## 2026-03-20 Attempt 11 — SUCCESS (driver OOM during save; patched)

Driver job:
- Wrapper: [wrap_shend_consider_missing.sh](/home/mcn26/project/tabula_rasa/notebooks/object_creation/orthos/wrap_shend_consider_missing.sh)
- Slurm resources: `2` CPU cores, `24G` RAM, `2-23:50:00` walltime
- Main job id: `2561074`, workers: `2561081`–`2561096`
- Started: `2026-03-20 10:05`, fitting completed: `2026-03-22 ~04:46`, driver OOM: `2026-03-22 04:51` (~42.7 hours total)

Dask worker configuration used:
- `SLURMCluster(cores=1, memory="96G", processes=1)`
- `cluster.scale(jobs=16)`
- Effective layout: `16` single-threaded workers, `~89.4 GiB` memory per worker

Observed result:
- `ortho_filter` removed 4 combinations involving `'reference'`; dropped `641 of 2103` (cell_type, cre_id) combos
- **All 16 workers completed cleanly** — no OOM kills, no nanny restarts, no Python exceptions, no shuffle failures
- `by_cre.pkl` saved at `2026-03-21 05:25` (~19 hours in); `by_cell_type.pkl` saved at `2026-03-22 04:46` (~42.7 hours in)
- Driver OOM'd at `04:51` with exit code `0:125` and MaxRSS `~24 GiB` (hit the 24G driver memory limit) during design matrix serialization in `save()`/`extract_params()`
- Design matrices were not written by the driver job

Recovery:
- `patch_design_matrices.py` submitted twice (`2636537`, `2636583`) on `2026-03-22`
  - First attempt (`10:51–10:53`, COMPLETED, 2 min): wrote most design matrices
  - Second attempt (`10:57–11:11`, COMPLETED, 14 min): completed remaining
- All 208 CRE design matrices (`by_cre_design/0.pkl`–`207.pkl`) and all 10 cell-type design matrices (`by_cell_type_design/0.pkl`–`9.pkl`) present and timestamped `2026-03-22 11:08–11:11`

Model verification (2026-03-26):
- All 5 core pickle files load cleanly: `by_cre.pkl` → `experiment_model`, `by_cell_type.pkl` → `experiment_model`, `by_cre_parameters.pkl` → `parameters`, `by_cell_type_parameters.pkl` → `parameters`, `training_data.pkl` → `scMPRA_data`
- `by_cre_design/_keys.json`: 208 entries (reference + 207 CREs); `by_cell_type_design/_keys.json`: 10 entries (reference + 9 cell types)
- **Model is intact and complete**

Model saved at:
- `/vast/palmer/pi/reilly/tabula_data/shendure/shendure_ortho_consider_missing_20260320/`

