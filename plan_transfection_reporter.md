# Plan: Transfection Reporter Necessity Evaluation

==WE NEED TO RE-EVAL THE CONCLUSION THAT THERE WAS NO ACTUAL COLLISIONS IN EXISTING SIMS: I THINK THERE WAS==~MCN

## Goal

Use simulation + power analysis to determine whether transfection reporters are worth including in future scMPRA experiments, or whether correct statistical modeling (ZINB with consider_missing) can compensate for their absence.

## Background

Transfection reporters add cloning complexity but provide two benefits:
1. **Zero classification**: distinguish "not transfected" (structural zero) from "transfected but not expressed" (biological zero)
2. **Reduced zero-inflation**: filtering removes structural zeros before modeling

Without a reporter, `consider_missing` re-introduces structural zeros, and the ZINB ZI component must absorb them. The question is whether this hurts statistical power enough to justify the experimental complexity.

## Approach

Simulate matched datasets under +/- reporter conditions using empirical parameters from existing orthos. Fit with the appropriate model (NB or ZINB, as determined by plan_nb_vs_zinb.md). Compare power and calibration.

## Prerequisites

- Complete the NB vs ZINB evaluation (plan_nb_vs_zinb.md) first — we need to know which model to use in each condition
- Existing `de_novo_simulation` framework handles the heavy lifting

## Design

### Experimental conditions

For each empirical dataset (shendure, cohen):

| Condition | Reporter | consider_missing | Model | flatten_overtransfection |
|-----------|----------|-----------------|-------|--------------------------|
| A: with reporter | Yes | No | TBD by NB/ZINB eval | True |
| B: without reporter | No | Yes | TBD by NB/ZINB eval | True |

Seelig already *is* the "no reporter" case, so it serves as empirical validation rather than a simulation target.

### Key simulation parameters to match

- Ground truth effect sizes: drawn from empirical ortho parameters
- MOI distribution: from empirical Bounds
- Library composition: from empirical library tables
- Cell counts per cell type: from empirical Bounds
- Dispersion (theta): from empirical ortho
- ZI rate: from empirical ortho (for condition B, higher ZI expected)
- Number of replicates: match empirical (6 for shendure, varies for cohen)

### What "no reporter" means in simulation

1. Simulate transfection normally (same MOI, same library)
2. Simulate transcription with ZINB parameters (some cells get zero from ZI)
3. **Do NOT filter** by reporter — keep all cells, including untransfected ones
4. Apply `consider_missing=True` to inflate the matrix
5. `flatten_overtransfection=True` always (measurement cannot resolve duplicates)

### What "with reporter" means in simulation

1. Simulate transfection normally
2. Simulate transcription with ZINB parameters
3. **Filter**: remove (cell, CRE) pairs that were never transfected (simulating reporter-based filtering)
4. `consider_missing=False` (reporter already removed structural zeros)
5. `flatten_overtransfection=True` always

### Role of flatten_overtransfection

Always True in both conditions. The measurement process (sequencing) cannot distinguish two copies of the same MPRA barcode in the same cell regardless of reporter presence. The reporter tells you *what* is present, not *how many copies*.

## Metrics

For each condition, across N simulation replicates (N=5-10):

1. **P-value calibration**: under null (no true effect), are p-values uniform? Plot QQ.
2. **Power at α=0.05**: fraction of true positives detected, as function of effect size
3. **AUROC / AUPRC**: discrimination between true effects and nulls
4. **Precision-recall curves**: at the median-AUPRC replicate

## Open questions to resolve before implementation

1. **What model for each condition?** Depends on NB/ZINB eval results. If ZINB wins for +consider_missing, use ZINB for condition B and NB for condition A. If NB wins everywhere, use NB for both.

2. **How to handle shendure's clonotype bottlenecking?** Different cell types get different CRE sets due to clonal expansion. The `two_thirds_inactive` notebook already handles this — review and reuse.

3. **Effect size range**: use empirical fold-changes from real data, or sweep a controlled range? Probably both — empirical for calibration, controlled for power curves.

4. **Number of simulation replicates**: 5 is fast but noisy for power estimates. 10-20 would be better. Compute cost depends on model fitting time per replicate.

## Timeline dependency

This plan cannot start until the NB/ZINB evaluation is complete, since the model choice for each condition depends on those results. However, the simulation infrastructure (ground truth construction, library tables, bounds objects) can be prepared in parallel.
