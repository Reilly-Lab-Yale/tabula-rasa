# Cohen CM ZINB — Run Graveyard

Cohen consider_missing ZINB fits. 21.8B total CM rows (317× Seelig CM).
By-cell-type slices: ~11B rows each. By-CRE slices: ~93M rows each.

---

## Attempt 1 — 2026-03-30 (by_cell_type only, exploratory)

**Goal**: See how far the current pipeline gets before hitting memory/time limits.
No core.py optimizations applied — `_prepare_subset_for_modeling` still densifies
sparse columns. Want empirical crash data before optimizing.

**Config**:
- Data: `unjoined/read_wise_mpra_retina_filtered.tsv` → read_wise_to_umi_wise → consider_missing=True
- 4 CPU workers × 64G + 2 GPU workers × 64G (H200)
- Driver: 16G, priority partition, 12h timeout
- Worker timeout: 12h CPU, 24h GPU

**Scripts**: `fit_by_cell_type.py`, `wrap_by_cell_type.sh`
(PYTHONPATH and cd paths updated to tabula-rasa-cohen-regen worktree)

**Predictions / things to watch**:
- `_prepare_subset_for_modeling()` densifies `umis_mpra_bc` sparse column — may OOM
  on by-cell-type slices (~11B rows × 8 bytes = ~82 GB dense, exceeds 64G worker)
- Dask `.compute()` of the inflated DDF may itself OOM before even reaching the fit
- Design matrix construction via formulaic should stay sparse (output='sparse')
- scipy 1.17.1 handles int64 indices automatically, so no int32 overflow expected
- If it somehow gets to TensorZINB, sparse matrix support is in place

**Status**: FAILED

**Job**: 6839359 (driver), workers 6839486-6839489

**Result**: Immediate failure — never reached actual fitting. The `_inflate_missing_split_level`
memory estimator (core.py:1355-1387) predicted 8,749 GB for a single cell-type slice and
raised `ValueError` at the `consider_missing_max_memory_gb = 100.0` cap.

The estimator assumes dense string-per-row representation and has no concept of sparse
encoding. Real memory usage would be far lower since `umis_mpra_bc` is 99.8% zeros
stored as `pd.SparseDtype`. The cap needs to be bypassed for this dataset.

Also: script was loading `unjoined/read_wise_mpra_retina_filtered.tsv` (read-wise) and
converting to UMI-wise at runtime. This is unnecessary — the new
`retina_single_counting_u6.scmpra/` parquet object is already UMI-wise and has no dummy
barcodes, so `consider_missing` can work directly on it.

**Driver log** (slurm-6839359.out):
```
[+] CPU workers connected.
[+] Creating...
ValueError: consider_missing split-level expansion estimated peak memory 8749.66 GB exceeds cap 100.0 GB.
```

**Stderr**: ortho_filter dropped 4 of 456 (cell_type, cre_id) combos with fewer than 3 nonzero entries. Normal.

**Fixes for Attempt 2**:
1. Load from `retina_single_counting_u6.scmpra/` (UMI-wise parquet, no dummy barcodes)
2. Bypass memory cap: `cohen.consider_missing_max_memory_gb = None`
3. Add comment explaining cap estimator doesn't account for sparse encoding

---

## Attempt 2 — 2026-03-30 (by_cell_type, cap bypassed)

**Goal**: Same as Attempt 1, but bypass the memory cap to see where the pipeline
actually fails (or succeeds). No core.py changes — want empirical data first.

**Changes from Attempt 1**:
- Data loading: `from_parquet("retina_single_counting_u6.scmpra")` (UMI-wise, no dummy barcodes)
- `consider_missing_max_memory_gb = None` (bypass dense-assumption cap)
- CLAUDE.md updated: U6 zeros are a harmless subset of CM, no dummy barcode conflict

**Config**: Same resources as Attempt 1 (4 CPU × 64G, 2 GPU × 64G H200, 12h driver)

**Predictions / things to watch**:
- Will get past the cap check this time
- `_get_missing_maps()` should pass (no ambiguous barcodes)
- The Dask `.compute()` call materializing ~11B rows per cell-type slice is the next
  likely failure point — even with sparse `umis_mpra_bc`, the string columns
  (`cell_bc`, `mpra_bc`, `cre_id`) are dense per row
- If it gets past `.compute()`, `_prepare_subset_for_modeling()` densifies the sparse
  column — another potential OOM
- If it reaches TensorZINB, sparse design matrices should work (scipy 1.17.1 + int64)

**Status**: FAILED

**Jobs**: 6839730 (driver), workers 6839779-6839782

**Result**: OOM — all 4 workers killed before reaching consider_missing expansion.
Crash occurred at core.py:2154 in `_series_unique_str`, which calls `.compute()` on
the `cell_type` column of the full dataset to enumerate levels and sort by size.
The 1.49B-row parquet (U6 zeros baked in) is too large for 64GB workers even for a
simple column read — Dask was repartitioning 229 partitions (`repartitiontofewer`)
before the compute and that alone OOMed every worker.

This is upstream of consider_missing expansion entirely. The new parquet object
includes the U6 zeros as "observed" data, so `get_data(include_missing=False)` still
returns 1.49B rows. Workers died with event-loop-unresponsive warnings (GIL contention
or large data movement) then SLURM cancelled them.

**Driver error**:
```
KilledWorker: Attempted to run task ('repartitiontofewer-e55f808398db48eaca0fcfc1496b132b', 0)
on 4 different workers, but all those workers died while running it.
```

**Worker logs**: All 4 workers show event-loop-unresponsive warnings (~3-4s) during
the repartition, then `CANCELLED` by SLURM at 18:20:39.

**Root cause**: `_series_unique_str` (core.py:291-292) materializes the entire dataset
to find unique cell_type values. With 1.49B rows this OOMs 64GB workers. The U6-expanded
parquet is just too big to read naively — even before any CM expansion.

**Discussion**: Two viable fixes:
1. Don't load from the U6-expanded parquet for CM fits. Instead load the filtered
   read-wise TSV → convert to UMI-wise (small: ~3.6M reads → ~3.2M UMI rows). Then
   let consider_missing handle all zero expansion. The U6 zeros would be regenerated
   by CM rather than pre-baked.
2. Fix `_series_unique_str` to read only a single small partition or use metadata to
   get unique cell_type values without scanning all 1.49B rows.

Option 1 is simpler and sidesteps the problem entirely. The U6-expanded parquet is
the right input for obs fits; the filtered read-wise file is the right input for CM.

---

## Attempt 3 — 2026-03-30 (by_cell_type, pre-CRE-coarse-expansion input)

**Goal**: Same as Attempt 2, but feed pre-CRE-coarse-expansion data into CM so
`_series_unique_str` and other pre-fit operations work on the compact ~3.2M-row
dataset rather than the 1.49B-row U6-expanded parquet.

**Changes from Attempt 2**:
- Data loading: `from_tsv("unjoined/read_wise_mpra_retina_filtered.tsv")` → `read_wise_to_umi_wise()`
- Comment explaining why the U6-expanded parquet can't be used for CM fits

**Config**: Same resources (4 CPU × 64G, 2 GPU × 64G H200, 12h driver)

**Predictions / things to watch**:
- Pre-fit ops (`_series_unique_str`, level size sorting) should be fast on 3.2M rows
- `_get_missing_maps()` should pass (no ambiguous barcodes after filtering)
- `_inflate_missing_split_level` cap bypassed — will proceed to Dask expansion
- First real test of whether workers can handle the CM expansion (~11B rows / cell type)
- If Dask `.compute()` of the expanded DDF OOMs, we'll need to fix `_prepare_subset_for_modeling`
  sparse densification or increase worker memory

**Status**: FAILED

**Jobs**: 6843038 (driver), CPU workers 6843108-6843111, GPU workers 6843321-6843322

**Result**: OOM during `hashjoinp2p` shuffle for CM expansion — never reached model fitting.

Pre-fit ops were fast on the compact 3.2M-row input (as expected). `_get_missing_maps()` passed
cleanly. `fit_by_cell_type_models()` and `extract_params()` submitted futures lazily and returned
immediately. The actual compute triggered only when `save()` → `flattened_copy()` → `future.result()`.
At that point the Dask scheduler began the CM expansion, which requires a `hashjoinp2p` hash-shuffle
(peer-to-peer merge) to build the full cell × barcode Cartesian product per cell-type slice (~11B rows).

Each CPU worker hit ~48 GB RSS within 3-4 minutes of startup, triggering the Dask nanny kill-restart
cycle (kill at 95% of 59.6 GiB limit = 56.6 GiB; MaxRSS per sacct: 57.4 GB). Workers restarted and
OOMed again immediately. GPU workers (submitted JIT via pre_fit_hook) had been queued for ~8 minutes
and only connected at 18:46:31 — the driver saw 6 workers connected and printed "Submitting fits",
but the CPU workers died again 2 seconds later and SLURM cancelled the whole cluster (sacct: GPU
workers MaxRSS 0.25 GB, 3s wall time). Nothing was ever fit.

`by_cell_type.pkl` on disk: 0 bytes (save started, wrote nothing). Driver backtrace ends at
`future.result()` → `KilledWorker: hashjoinp2p-f27ee7c297600453756e591309bd8ca2`.

**SLURM memory (sacct)**:
- CPU workers: MaxRSS 57.4 GB / 60G allocated — at ceiling
- GPU workers: MaxRSS 0.25 GB / 64G — never ran
- Driver: MaxRSS 0.88 GB / 16G — fine

**Root cause**: 64 GB per CPU worker is insufficient for the `hashjoinp2p` shuffle that underlies CM
expansion. The by-cell-type slice is ~11B rows; with 4 workers the per-worker shuffle bucket fills
to ~48 GB before the merge completes.

**Candidate fixes for Attempt 4**:
1. More CPU workers (8 × 64G): more hash buckets → less per-worker memory during shuffle. No code
   changes. Doesn't help with `_prepare_subset_for_modeling` densification at the fit stage.
2. Fix `_prepare_subset_for_modeling` sparse densification (core.py:332-338): this happens AFTER
   the CM shuffle, so won't help with this crash — but would reduce memory during actual fitting.
3. Increase per-worker RAM (128G): brute force but expensive if other jobs are competing.

---

## Attempt 4 — 2026-03-30 (by_cell_type, 4 × 128G workers)

**Goal**: Test whether doubling per-worker RAM (64G → 128G) allows the CM expansion
shuffle to complete. Same total worker count, same code.

**Changes from Attempt 3**:
- CPU worker memory: 64G → 128G (SLURMCluster `memory="128G"`)
- Also deleted stale 0-byte `by_cell_type.pkl` from Attempt 3 (first submission
  hit `EOFError` trying to load it before the fix — job 6844157, trivial, not counted)

**Config**: 4 CPU × 128G, 2 GPU × 64G H200, 12h driver

**Predictions / things to watch**:
- Will the `hashjoinp2p` shuffle fit in 128G workers?
- Or will the p2p shuffle scale with available memory (proportional consumption)?

**Status**: FAILED

**Jobs**: 6844542 (driver), CPU workers 6844603-6844606

**Result**: Same OOM pattern as Attempt 3, proportionally scaled. Workers hit 95% of
the 119.2 GiB Dask memory limit (95.4 GiB) within 8 minutes of startup — exactly the
same percentage as Attempt 3 (95% of 59.6 GiB = 56.6 GiB). The p2p shuffle consumes
memory proportional to available worker RAM, not a fixed amount.

Crash point was `_inflate_missing_split_level` line 1357:
`rep_sizes["target_rows"].sum().compute()` triggered a `dropduplicates` task from the
persisted maps (deferred from `_get_missing_maps` lines 1277-1278). All 4 workers died
attempting it. Driver backtrace: `KilledWorker: dropduplicates-8c582a08...`.

**SLURM memory (sacct)**:
- CPU workers: MaxRSS 114.2–114.5 GB / 120G allocated — at ceiling (same 95% as Attempt 3)
- Driver: MaxRSS 1.63 GB / 16G — fine
- GPU workers: never submitted (crash before pre_fit_hook)

**Root cause**: The Dask p2p shuffle for CM expansion is fundamentally incompatible with
Cohen's data dimensions. The by-cell-type Cartesian product is ~50K cells × ~230K MPRAs
= ~11B rows per slice. Compare to Shendure CM which works fine: ~4K cells × ~3K MPRAs
= ~128M rows per slice (86× smaller). The 77× larger barcode library is the dominant factor.

More RAM per worker does not help — the shuffle fills whatever memory is available.
`client.compute(ddf)` (standard_fit line 2181) ultimately materializes the 11B-row DDF
as a single pandas DataFrame on one worker (~880 GB at ~80 bytes/row for 5 dense string
identity columns), which is impossible regardless of worker size.

**Conclusion**: The current `_inflate_missing_split_level` → `client.compute(ddf)` →
single-pandas-on-worker architecture cannot handle Cohen CM. Need to build sparse design
matrices and response vectors directly from the maps (cell_map, mpra_map, observed)
without materializing the full 11B-row expansion.
