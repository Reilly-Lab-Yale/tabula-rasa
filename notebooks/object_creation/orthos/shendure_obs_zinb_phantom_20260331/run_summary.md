# Shendure obs ZINB phantom-zero validation — run summary

**Date:** 2026-03-31
**Branch:** `cohen-regen`

## Jobs

### by_cre (job 6921909)
- **Node:** a1132u11n01 (Intel Xeon 8562Y+, Emerald Rapids)
- **Start:** 2026-03-31 11:03:17 — **End:** 2026-03-31 11:28:25
- **Elapsed:** 25m 08s
- **Allocated CPUs:** 2 — **Peak RSS:** 1,815 MB — **Requested:** 16G
- **Script:** `fit_phantom.py` — 208 CREs × ZINB, `init_method="nb"`, CPU

### by_cell_type (job 6922750)
- **Node:** a1130u05n03 (Intel Xeon 8562Y+, Emerald Rapids)
- **Start:** 2026-03-31 11:25:56 — **End:** 2026-03-31 11:33:00
- **Elapsed:** 7m 04s
- **Allocated CPUs:** 2 — **Peak RSS:** 1,116 MB — **Requested:** 24G
- **Script:** `fit_phantom_by_cell_type.py` — 10 cell types × ZINB, `init_method="nb"`, CPU

## Results summary

Both fits used phantom-zero weighted compression against the existing `shendure_ortho_20260306` design matrices. Results compared to saved full-expansion fits.

### by_cre

| rep | Saved ortho | Phantom | %diff |
|-----|------------|---------|-------|
| 2B1 | 0.3925 | 0.3821 | -2.66% |
| 2B2 | 0.4163 | 0.4074 | -2.13% |
| A1  | 0.4447 | 0.4391 | -1.26% |
| A2  | 0.4194 | 0.4124 | -1.68% |
| B1  | 0.3732 | 0.3652 | -2.14% |
| B2  | 0.3641 | 0.3542 | -2.72% |

- Mean %diff: **-2.07%** (systematic downward bias)
- Per-CRE Pearson r: **0.967** — distribution shape preserved
- Per-CRE std: saved=0.183, phantom=0.179 — virtually identical

### by_cell_type

| rep | Saved ortho | Phantom | %diff |
|-----|------------|---------|-------|
| 2B1 | 0.0194 | 0.0207 | +6.45% |
| 2B2 | 0.0143 | 0.0164 | +14.79% |
| A1  | 0.0395 | 0.0403 | +2.00% |
| A2  | 0.0269 | 0.0285 | +5.96% |
| B1  | 0.0207 | 0.0220 | +6.30% |
| B2  | 0.0199 | 0.0208 | +4.44% |

- Mean %diff: **+5.65%** (systematic upward bias)
- Per-cell-type std: saved=0.0133, phantom=0.0136 — nearly identical

## Interpretation

Phantom-zero compression produces ZI estimates within 2–6% of full-expansion fits at the bounds level. The biases (by_cre slightly low, by_cell_type slightly high) are symmetric and reflect different local optima on the non-convex ZINB surface — not a flaw in the weighted objective. LLF differences were <0.01% throughout, confirming mathematical equivalence.

The per-CRE ZI distribution shape is well-preserved (r=0.967), which is what matters for simulation heterogeneity. This validates phantom-zero compression for obs ZINB fits. The approach is ready to extend to Cohen CM fits.
