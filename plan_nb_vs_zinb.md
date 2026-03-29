# Plan: NB vs ZINB Model Comparison

## Goal

Evaluate whether ZINB is justified over plain NB for scMPRA data, across all three datasets and +/- consider_missing conditions.

## Rationale

Preliminary results suggest AIC/BIC favor NB. ZI estimates are near-zero in most conditions. The one exception is high-signal data with consider_missing enabled (lots of structural zeros), where ZI picks up signal. ZINB is "wider" — it includes NB as a special case — and may be needed for future experiments without transfection reporters. This evaluation determines if that extra complexity is earning its keep.

## Evaluation matrix

| Dataset | consider_missing | ZINB | NB | Status |
|---------|-----------------|------|-----|--------|
| shendure | No | exists | **need fit** | NB fit: CPU-only, fast (~2h) |
| shendure | Yes | exists | **need fit** | NB fit: was running on McCleary (dead), redo on Bouchet. GPU for cell-type |
| cohen | No | exists | **need fit** | NB fit: CPU-only, fast |
| cohen | Yes | **need fit** | **need fit** | Both new. Tests "ignore reporter" scenario |
| seelig | Yes (default) | exists | **need fit** | NB fit: GPU for cell-type |
| seelig | No | skip | skip | Meaningless — raw data has no zeros without consider_missing |

**3 existing orthos** can compute AIC immediately. **7 new fits** needed.

## Step 1: Post-hoc AIC/BIC from existing orthos

AIC = 2k - 2·ln(L), BIC = k·ln(n) - 2·ln(L)

- `llf_total` is stored in every fit result
- `k` = number of parameters: count NB coefficients (x_mu) + ZI coefficients (x_pi, if ZINB) + 1 (theta)
- `n` = number of observations (rows in design matrix)
- Write a utility that loads each ortho, iterates over models, extracts llf/k/n, computes AIC/BIC
- Aggregate per-model AIC/BIC into a summary table

Existing orthos to compute from:
1. `shendure_ortho_20260306` (ZINB, no consider_missing)
2. `shendure_ortho_consider_missing_20260320` (ZINB, +consider_missing)
3. `seelig_ortho_20260320` (ZINB, +consider_missing)
4. Cohen ortho in `/nfs/roberts/project/pi_skr2/shared/tabula_data/cohen_ortho/` (ZINB, no consider_missing)

## Step 2: New NB fits

For each needed NB fit, create a script based on the existing ZINB fit scripts with `nb_only=True` added to `criss_cross()` or the individual fit methods.

### 2a. Shendure NB (no consider_missing)
- Template: `shendure_ortho_creation_and_validation.ipynb` (original ZINB fit)
- Add `nb_only=True`
- Resources: 1 CPU worker, 60 GB, ~2h
- Output: `shendure_ortho_nb_20260329/`

### 2b. Shendure NB (+consider_missing)
- Template: `shendure_consider_missing.py` / `shendure_consider_missing_nb.py` (already exists but needs redo on Bouchet)
- Resources: 16 CPU workers × 96 GB, ~43h. Consider GPU for cell-type.
- Output: `shendure_ortho_consider_missing_nb_20260329/`
- **Note**: check if GPU works for NB cell-type fits (should, since `use_gpu` is independent of `nb_only`)

### 2c. Cohen NB (no consider_missing)
- Template: `cohen_retina.ipynb`
- Add `nb_only=True`
- Resources: CPU-only, moderate. Cohen is smaller than shendure.
- Output: `cohen_ortho_nb_20260329/`

### 2d. Cohen ZINB (+consider_missing)
- New script needed. Template: cohen ortho creation + `set_consider_missing(True)`
- Resources: TBD — depends on cohen data size with consider_missing
- Output: `cohen_ortho_consider_missing_20260329/`

### 2e. Cohen NB (+consider_missing)
- Same as 2d but with `nb_only=True`
- Output: `cohen_ortho_consider_missing_nb_20260329/`

### 2f. Seelig NB (+consider_missing)
- Template: `fit_seelig_ortho_by_cre.py` / `fit_seelig_ortho_by_cell_type.py`
- Add `nb_only=True`
- Resources: by_cre CPU (~2h), by_cell_type GPU (~2h on 2x H200)
- Output: `seelig_ortho_nb_20260329/`

## Step 3: Post-hoc AIC/BIC from new fits

Same utility as step 1, applied to all new orthos.

## Step 4: Likelihood Ratio Test

For each (dataset, consider_missing) condition where we have both NB and ZINB:

- Λ = 2 × (LL_ZINB - LL_NB)
- Under H₀ (NB is true): Λ ~ 50:50 mixture of χ²(0) and χ²(q), where q = number of ZI parameters
- Conservative test: use plain χ²(q) — if significant even with this, ZI is definitely needed
- Compute per-model (per CRE or per cell-type) and aggregate

## Step 5: Summary

- Table: dataset × consider_missing × {AIC_NB, AIC_ZINB, ΔAIC, LRT_p}
- Aggregate across models (median ΔAIC, fraction of models where ZINB wins)
- Expected result: ZINB wins for +consider_missing/high-signal, NB wins otherwise
- Decision: if ZINB wins in the "no reporter" scenario, keep ZINB as default

## Parallelism plan

Many fits are independent and can run simultaneously on Bouchet:
- All NB fits (2a-2f) can run in parallel
- Cohen +consider_missing ZINB (2d) can run in parallel with NB fits
- AIC computation is post-hoc and fast (minutes)
- Each fit gets its own Sonnet babysitter agent

## GPU verification

Add to `_tensorzinb_fit` (or verify already present):
```python
if use_gpu:
    devices = tf.config.list_physical_devices('GPU')
    if not devices:
        raise RuntimeError("use_gpu=True but no GPU found")
    logger.info(f"GPU device: {devices[0]}")
```
