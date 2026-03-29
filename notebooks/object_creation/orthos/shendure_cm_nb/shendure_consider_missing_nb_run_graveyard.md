# Shendure Consider-Missing NB-only Run Graveyard

Fitting a plain negative binomial (NB) model on the Shendure `consider_missing=True` dataset,
for head-to-head comparison with the ZINB model
(`shendure_ortho_consider_missing_20260320`).

Key difference from the ZINB runs: `criss_cross(..., nb_only=True)` is passed through the
call chain to `TensorZINB(..., nb_only=True)`, and method-of-moments initialization is
automatically disabled (MoM produces `x_pi` which has no meaning in a plain NB model).

Scripts:
- Driver: [shendure_consider_missing_nb.py](/gpfs/gibbs/project/reilly/mcn26/tabula_rasa/notebooks/object_creation/orthos/shendure_consider_missing_nb.py)
- Wrapper: [wrap_shend_consider_missing_nb.sh](/gpfs/gibbs/project/reilly/mcn26/tabula_rasa/notebooks/object_creation/orthos/wrap_shend_consider_missing_nb.sh)

Output path: `/vast/palmer/pi/reilly/tabula_data/shendure/shendure_ortho_consider_missing_nb_20260326/`

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
