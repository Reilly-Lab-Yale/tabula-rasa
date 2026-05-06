# Plotting ideas: "what determines the cohen-vs-shendure power gap?"

Stashed 2026-05-06 mid-iteration. Once the bounds-expansion top-up run lands
and a saturation-resistant primary metric is chosen (probably power@FC=1.5 or
power-curve AUC), revisit these as figure deliverables.

## Context

Synthetic factorial sweep over 7 design axes (1000 LHS samples). Cohen empirical
power was consistently higher than shendure. The marginals show how power
depends on each axis, but they don't directly answer *which axis differences
between cohen and shendure are doing the work*.

Cohen vs shendure axis values (see `synthetic_factorial.py::EMPIRICAL`):

| Axis | Shendure | Cohen | Ratio |
|---|---|---|---|
| n_cells | 8,201 | 18,633 | 2.3x |
| n_cres | 207 | 116 | 0.56x |
| bcs_per_cre | 136 | 17,244 | 127x |
| moi | 18 | 149 | 8.3x |
| lib_alpha_nb | 0.22 | 1.39 | 6.3x |
| minP | 0.041 | 0.94 | 23x |

Caveat: bcs_per_cre and moi at cohen sit far outside the original LHS bracket.
Any attribution involving those axes needs the bounds-expansion top-up run
that's in flight (or honest "lower bound" labelling if we stick with current
data).

## Display ideas

### 1. Per-axis attribution bar chart -- "the headline"

For each axis, compute the LOESS-predicted power at the shendure value vs
the cohen value. Bar = `power(cohen) - power(shendure)`, sorted by magnitude.

Bars approximately sum to the observed cohen-shendure gap. Clearly answers
"which axis difference is doing the work."

Implementation: standalone script `attribution.py` reading
`output/samples_power.parquet` + `output/samples.parquet`. Single PNG/SVG.

Caveat: holding-others-at-median assumption may understate interactions.

### 2. Annotated marginals -- the supplement

Same 7-panel marginal layout as `output/marginals.svg`, but on each panel:
- shendure point at `(shendure_x, predicted_power_shendure)` -- orange dot
- cohen point at `(cohen_x, predicted_power_cohen)` -- purple dot
- arrow connecting them, labeled with delta power
- panels sorted by absolute delta

Lets the reader see the slope that produces the bar chart number, and which
axes have steep vs gentle gradients in the relevant region.

### 3. Counterfactual ladder -- waterfall, most explanatory

Start from the synthetic median sample's predicted power. Step through axes
one at a time, replacing the median value with cohen's value, re-predict at
each step. Each step's height = that axis's marginal contribution conditional
on prior swaps. Final value approximates predicted cohen power. Visually a
waterfall chart, ordered by descending magnitude.

This makes interactions visible: if swapping bcs_per_cre first contributes
+0.4 power but moi contributes only +0.05 *after* that swap, the order matters.
Reorder by magnitude to highlight the largest contributors.

Caveat: same as #1 -- meta-model is just LOESS marginals, so true
multidimensional interactions are partially hidden.

### Out-of-range handling

For axes where cohen sits outside the LHS bracket:
- Option 1 (conservative): clamp cohen value to the LHS edge before
  predicting power. Resulting delta is a "this much *or more*" lower bound.
  Mark with a hatched bar / "(>=)" annotation.
- Option 2 (expansion): rely on the bounds-expansion top-up run to provide
  real data in the cohen-extreme corner. If those samples land successfully,
  the LOESS smoother will extend naturally into that regime.

Going with option 2 once top-up is done.

## Recommended order

1. Build #1 first (single horizontal bar chart). ~50 lines, runs locally.
2. Layer in #2 by augmenting the existing marginals SVG.
3. Skip #3 for v1 unless reviewer asks "why this order."

## Saturation note

P@FC=2 saturates (~40% of samples at 0 or 1 ceiling). For attribution use a
non-saturating summary -- power@FC=1.5 or power-curve AUC over [1, 3], whichever
better discriminates cohen-shendure. Decision pending after the multi-metric
re-aggregate finishes.

Decision (2026-05-06): primary metric is **`power_auc_1to3`**. 0.1% saturated
at max, 1.9% at min, defined for all 1099 combined samples.

## Confounder analysis: library composition and compression

Two confounders considered when claiming "x y z assay-design factors drive
the cohen-shendure power gap":

### (1) Library composition (CRE activity range)

The original 5x5 power sims at `analyses/simulation/activity_prc/{shendure,cohen}/`
used each dataset's *full empirical mu distribution*, not a common
activity_max_mult. So the original cohen-vs-shendure power gap mixes
assay-design effects with library-composition effects.

Empirical activity ranges (p95(mu) / minP) differ dramatically:

| dataset  | minP    | p95 mu | p95/minP | library type |
|----------|---------|--------|----------|--------------|
| shendure | 0.041   | 4.10   | 99x      | strong promoters in heavy tail (eef1aP=272x, pgk1P, ubcP) |
| cohen    | 0.936   | 1.04   | 1.11x    | variants of one regulatory element (knockout/mutation) |

**The confound goes against cohen.** Cohen's narrow library means most CRE
pairs have small biological FC, which is *harder* to detect. Yet cohen
empirically has *higher* power. So if we removed the library-composition
confound (gave cohen a wider activity range), the assay-design contribution
to cohen's advantage would be even larger, not smaller. The empirical gap is
a lower bound on the true assay-design advantage.

This is borne out in the LOESS attribution: at LHS [2, 8] coverage of
activity_max_mult, the predicted power barely differs between cohen-clamped
(2.0) and shendure-clamped (8.0). When the heavy-lifters (bcs_per_cre, moi,
n_cells, n_cres) are held at empirical, varying activity_max_mult between
those clamped extremes barely moves the prediction (delta +/- 0.01-ish).

### (2) Compression of measured fold changes

Concern: maybe cohen's measured FCs are systematically smaller than the true
biological FCs (some saturation / quantification artifact), making the
power test look easier than it really is.

This is *not* a confound for simulated power. The synthetic factorial
controls ground-truth `mu` per CRE explicitly; what we call "FC" in the
power-vs-FC plots is the ratio of true mus, not the measured one. Power
emerges from the count noise around those known truths -- and that count
noise is parameterized by `lib_alpha_nb` (NB dispersion), which IS one of
our axes. So whatever measurement-noise effect is operating in the
empirical data, it is captured by the dispersion axis, which the
attribution analysis already includes.

For empirical 5x5 results specifically, compression *would* matter: it
could make a real biological FC of 2 look like a measured FC of 1.5,
which is harder to detect. But cohen's empirical mu distribution is
narrow because the *library design* is narrow (variants of one element),
not because of compression. Confirmed by inspecting the fitted mu series
in `tabula_data_new/cohen/cohen_obsingle_nb_phantom/by_cell_type_parameters.pkl`:
cohen-Rod CREs span 0.026 to 1.08 with a tight cluster near 1.0 = unit FC,
which is exactly what you would expect from saturation-mutagenesis variants
of a single CRE.

### Summary statement for methods/results

"The cohen-vs-shendure power gap arises from cohen's superior assay-design
parameters: 127x more barcodes per CRE, 8x higher MOI, ~2x more cells, and
modestly fewer CREs. These advantages are partially offset by cohen's
higher count dispersion and roughly cancel on activity_max_mult. The
narrow CRE activity range of the cohen library would, if anything, hurt
cohen's empirical power -- the confound goes against cohen, so the
observed empirical advantage is a lower bound on the true assay-design
advantage."
