# Plan: Fix flatten_overtransfection and Replot

## Goal

The shendure empirically-calibrated simulations used `flatten_overtransfection=False`, which is incorrect. All three datasets should use `flatten=True` because the sequencing measurement cannot resolve duplicate integrations of the same MPRA barcode in the same cell. Fix this and regenerate plots.

## Impact assessment

The empirical collision rate for shendure is **0.053%** of transfection events. This means:
- ~243 collision events expected across ~780,000 transfection events per replicate
- The stored `.scmpra` files on disk show **zero** collisions in practice
- Re-flattening would be a **no-op** — there is nothing to flatten

Despite this, the fix should be applied for correctness, and to prevent confusion if parameters change in future simulations (higher MOI or smaller library would increase collision rate).

## Affected files

All shendure simulation notebooks that set `flatten_overtransfection=False`:

1. `notebooks/object_creation/simulations/emperically_calibrated/shendure/shendure_calibration_mwu_all_cell_types.ipynb`
2. `notebooks/object_creation/simulations/emperically_calibrated/shendure/shendure_power_mwu_all_cell_types.ipynb`
3. `notebooks/object_creation/simulations/emperically_calibrated/shendure/shendure_pairwise_power_mwu.ipynb`
4. `notebooks/object_creation/simulations/emperically_calibrated/shendure/two_thirds_inactive.ipynb`

Cohen notebooks already use `flatten_overtransfection=True` — no changes needed.

Seelig has no simulations yet.

## Approach: replot without re-simulating

The simulation data is stored on disk. The pipeline is:

```
simulate transfection → simulate transcription → save .scmpra → fit orthos → precompute Wald → run tests → plot
```

The `.scmpra` files are saved **before** flattening (for shendure, `flatten=False` means the data was written as-is). Since there are zero actual collisions in the stored data, re-flattening and re-fitting would produce identical results. The existing orthos and test results are valid.

### What to do

1. **Fix the notebooks**: change `flatten_overtransfection=False` → `flatten_overtransfection=True` in all four shendure simulation notebooks
2. **Fix the replot scripts**: check `replot_shendure_calibration.py` and `replot_shendure_power.py` for the same issue
3. **Do NOT re-run simulations**: the results are identical (zero collisions in practice)
4. **Regenerate plots** by running the replot scripts (fast — just loads results and plots)
5. **Verify**: spot-check that the stored `.scmpra` files indeed have zero collisions (already confirmed for `shendure_calibrated_sim_with_orthos_20251118/scMPRA/0.scmpra`)

### Also fix: one_library_replicate default

In `core.py` line 6933, `one_library_replicate` has `flatten_overtransfection=False` as default. Change to `True`.

## Already fixed

- `overtransfected()` method in `core.py`: now reports per-event collision rate instead of per-cell probability (birthday paradox fix)
- `estimated_percent_conflict.ipynb`: removed spurious `*100` in bar labels

## Steps

1. Edit the 4 shendure simulation notebooks: `False` → `True`
2. Edit `core.py:6933`: default `flatten_overtransfection=False` → `True`
3. Check and fix replot scripts if needed
4. Run replot scripts to regenerate SVGs/PNGs
5. Verify plots are unchanged (expected, since zero collisions)

## Estimated effort

~30 minutes. All changes are one-line edits. No compute needed.
