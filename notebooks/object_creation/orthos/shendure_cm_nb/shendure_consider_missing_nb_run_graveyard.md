# Shendure Consider-Missing NB-only Run Graveyard

Fitting a plain negative binomial (NB) model on the Shendure `consider_missing=True` dataset,
for head-to-head comparison with the ZINB model
(`shendure_ortho_consider_missing_20260320`).

Key difference from the ZINB runs: `criss_cross(..., nb_only=True)` is passed through the
call chain to `TensorZINB(..., nb_only=True)`, and method-of-moments initialization is
automatically disabled (MoM produces `x_pi` which has no meaning in a plain NB model).

Scripts (Bouchet, split approach):
- By-CRE: `fit_by_cre.py` + `wrap_by_cre.sh` (CPU, 16 workers × 96GB)
- By-cell-type: `fit_by_cell_type.py` + `wrap_by_cell_type.sh` (CPU, 16 workers × 96GB)
- By-cell-type GPU race: `fit_by_cell_type_gpu.py` + `wrap_by_cell_type_gpu.sh` (4 CPU setup workers + 2 H200 via pre_fit_hook)
- Merge: `merge.py`

Output path: `/nfs/roberts/project/pi_skr2/shared/tabula_data/shendure/shendure_cm_nb_20260329/`

Previous scripts (McCleary, dead cluster — attempts 1–3):
- Driver: `shendure_consider_missing_nb.py` (used criss_cross, single monolithic run)
- Wrapper: `wrap_shend_consider_missing_nb.sh`
- Output: `/vast/palmer/pi/reilly/tabula_data/shendure/shendure_ortho_consider_missing_nb_20260326/`

---

## Attempt 4 — SUCCESS (2026-03-30): by_cre

**Driver job:** 6832721 (node a1132u07n04, Intel Xeon 8562Y+ / Emerald Rapids)
**Workers:** 6832732–6832747 (16 workers × 96GB, CPU-only, priority partition)
**Start:** 2026-03-30 12:16:52 UTC → **End:** 2026-03-30 18:56:58 UTC
**Elapsed:** 06:40:06

**Resource usage (driver):**
- Allocated: 2 CPUs, 24 GB RAM
- Peak RSS: 2.7 GB (11% of allocated)
- CPU time: 00:25:58 total (3.2% efficiency — expected; driver is a Dask coordinator, not a compute node)

**Output:** `shendure_cm_nb_20260329_by_cre`

**Notes:** Clean run, 0 erred tasks. Run stats: `run_stats_cre_6832721.txt`.

---

## Attempt 5 — SUCCESS (2026-03-30): by_cell_type (GPU)

**Driver job:** 6834191 (node a1130u35n04, Intel Xeon 8562Y+ / Emerald Rapids)
**CPU setup workers:** 6834213–6834216 (4 workers × 96GB, nodes a1132u05n04, a1132u07n04)
**GPU workers:** 6834307, 6834308 (2 × H200, 64GB each, nodes a1124u11n01, a1124u28n01)
**Start:** 2026-03-30 13:34:58 → **End:** 2026-03-31 00:43:45
**Elapsed:** 11:08:47

**Resource usage:**
- Driver: 1.6 GB peak RSS / 24 GB allocated, CPU time 00:29:55
- CPU workers peak RSS: 28–42 GB / 96 GB allocated (design matrix construction)
- GPU workers peak RSS: 39–41 GB / 64 GB allocated (TensorZINB fitting)
- GPU: H200 (1 per worker), ~141/144 GB VRAM utilized during fitting

**Timeline:**
- Setup phase (CPU workers build design matrices): ~5 min
- pre_fit_hook → GPU workers submitted and connected: ~50 min queue wait
- Fitting phase (10 cell types on 2 H200s): ~10 hours
- Observed rate: ~48 min/type for small types, Pluripotent (~128M rows) much longer

**Output:** `shendure_cm_nb_20260329_by_cell_type_gpu`

**Notes:** Clean run, 0 erred tasks. Head-to-head race vs CPU (Attempt 6). GPU finished in 11:09; CPU still running at 25+ hours. Run stats: `run_stats_gpu_6834191.txt`.

---

## Attempt 6 — IN PROGRESS (2026-03-30): by_cell_type (CPU benchmark)

**Driver job:** 6832722 (node a1132u20n03, Intel Xeon 8562Y+ / Emerald Rapids)
**Workers:** 6832748–6832763 (16 workers × 96GB, CPU-only)
**Start:** 2026-03-30 12:16:52
**Elapsed:** 1d+ (still running)

**Notes:** CPU benchmark counterpart to Attempt 5 (GPU). Same data, same config, no GPU.
Running to establish GPU vs CPU speedup factor. Will be completed or cancelled when
benchmark data is sufficient.

---

## Attempt 1 — FAILED (2026-03-26)

**Driver job:** 2877106
**Workers:** 2877135–2877150 (16 workers)

Ran for ~6.5 hours. CRE fits completed; failed during cell-type fitting phase with 26 erred tasks.

**Root cause:** `_label_tensorzinb_regressors` unconditionally accessed `model["weights"]["x_pi"]`, which does not exist in NB-only models. Job cancelled.

**Fix:** Added `if not model.get("nb_only", False):` guard in `_label_tensorzinb_regressors`. Also stamped `result["nb_only"] = nb_only` in `_tensorzinb_fit` so intent is explicit rather than inferred from weight key absence.

---

## Attempt 3 — CANCELLED (2026-03-27 to 2026-03-30)

**Driver job:** 2909824
**Workers:** 2909836–2909851 (16 workers, CPU-only)

Fitting phase completed cleanly (0 erred tasks, ~7h). `by_cre.pkl` (2MB) saved successfully at ~21:36 Mar 27. Driver then blocked in `by_cell_type.flattened_copy()` waiting on cell-type model futures. Workers were actively computing (6–7 at ~99% CPU with 17–39GB memory each) but far too slow on CPU — after ~3 days total runtime and ~325 monitor checks, `by_cell_type.pkl` was still 0 bytes. Tasks ticked only 3 times (14820→14823) over the entire save phase.

**Outcome:** Manually cancelled. CPU-only fitting of the largest cell-type models (consider_missing=True, 12M–128M rows each) is impractically slow. Retrying on GPU cluster.

---

## Attempt 2 — FAILED (2026-03-27)

**Driver job:** 2882624
**Workers:** 2882633–2882648 (16 workers)

Ran for ~8 hours. CRE fits completed (~20 min), cell-type fits completed (~8 hours). Failed at `extract_params` stage with 208 erred tasks.

**Root cause:** `_extract_zi` (nested inside `extract_parameters`) unconditionally accessed `model["weights"]["x_pi"]` — same bug pattern as Attempt 1, in a different function. The driver then blocked indefinitely waiting on errored futures and could not save. Output directory has only a 0-byte `by_cre.pkl`. Job manually cancelled after ~3.5 hours stuck.

**Fix:** Added `if model.get("nb_only", False): return None` guard in `_extract_zi`. Also patched `flatten_param_representation` to skip the nb×zi cross-join when `zi` is `None` (for downstream `describe_parameters` compatibility).

---
