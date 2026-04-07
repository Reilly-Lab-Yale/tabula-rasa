# Paper Plan

## Style note

Prefer scripts (.py) over notebooks (.ipynb) for all new analysis work --
scripts are easier to automate, version-control, and review. Notebooks are fine
for exploratory work but final analysis should be scriptified.

All re-runs and new fits must use `flatten_overtransfection=True`.

---

## NB vs ZINB model selection results

Restricting to the condition each dataset would actually use in practice
(obs for datasets with a transfection reporter, CM for seelig which has none):

| Dataset   | Condition       | LRT Bonferroni sig | AIC prefers ZINB | mu shift NB->ZINB |
|-----------|-----------------|--------------------|------------------|--------------------|
| Seelig    | CM (mandatory)  | 0%                 | 0%               | 0.22%              |
| Cohen     | obs             | 54%                | 82.5%            | 2.87%              |
| Shendure  | obs             | 6.7%               | 25%              | 12.3%              |
| Takeshi   | obs             | 0%                 | 0%               | ~0%                |

### Conclusions and canonical model choices

- **Seelig (CM)**: NB is the right model. Zero models significant, zero AIC
  wins, mu parameters identical to 4 decimal places.
  **Canonical: `SEELIG_BOUNDS` = CM NB** (`seelig_cm_nb_phantom`)

- **Cohen (obs)**: ZINB is warranted. Majority significant under Bonferroni,
  mu stable (~3% shift). By-cell-type: 100% LRT significant, 100% AIC.
  **Canonical: `COHEN_BOUNDS` = obs ZINB** (`cohen_obs_zinb_phantom_20260401`)

- **Shendure (obs)**: NB is sufficient. Only 7% survive Bonferroni,
  BIC agrees (3% prefer ZINB by_cre, 0% by_cell_type).
  **Canonical: `SHENDURE_BOUNDS` = obs NB** (`shendure_obs_nb_phantom`)

- **Takeshi (obs)**: NB is the right model. Zero significant under any test
  (LRT, AIC, BIC), both by_cre (150 models) and by_cell_type (3 models).
  ZINB lambda negative for by_cell_type -- NB fits strictly better.
  **Canonical: `TAKESHI_BOUNDS` = obs NB** (`takeshi_obs_nb_phantom`)

These canonical choices are set as default aliases in `scMPRAforge.core` and
will be used to parameterize all downstream simulations.

### Collision rates

| Dataset   | Estimated collision % |
|-----------|-----------------------|
| Cohen     | 0.117%                |
| Shendure  | 0.053%                |
| Seelig    | 7.630%                |

---

## 3x sets of hypotheses

Activity hypotheses : "CRE is significantly different from minP, within a cell type"
Variant effect hypotheses : "two arbitrary CREs are different from each other, within a cell type"
Cell-type specificity hypotheses : "the same CRE is expressed differently in two different cell types"

---

## Tasks

Remember to collect runtime statistics for all tasks with nontrivial compute.

- [x] Consolidate md files & revise figure plan to match analysis goals (this task)
- [x] Map analysis flow, make sure everything is organized
- [x] Create new data directory
- [x] Organize pre-processing
- [x] Re-fit models with phantom zeroes (all datasets, obs + cm)
- [x] Delete non-phantom fit code and orthos (superseded by phantom versions)
- [x] Regen LRT figures (NB vs ZINB comparison) with properly phantomed models
- [x] Pick 1x canonical model (NB or ZINB) per dataset from phantomed LRT results
- [x] Re-extract bounds (from_ortho bug fixed, all 10 presets extracted)
- [x] QC plots for all 10 phantom orthos (convergence, theta, r-values, mu vs mean, ZI)
- [x] Fix Wald for phantom and Wald for NB
  - NB-only loglik, weighted covariance, design matrix reconstruction from recipe
  - Adaptive ridge + sandwich->hessian fallback for extreme phantom weights
  - 15 unit tests + 3 simulation tests passing
- [x] Delete pre-phantom run_stats/run_summary files (obsolete after phantom refit)
- [x] Re-extract bounds presets with rep_ids, zi Series, CM flag, nb_only, n_negative_controls
  - All 10 presets re-extracted (shendure + cohen done, seelig job 7407765)
  - NB bounds now have zero-valued zi Series (not None)
  - n_negative_controls enables reference_pooling in simulate_library
  - Simulation scripts no longer need manual ZINB zi splicing
- [ ] Fix _classifier_summary empty-array crash on seelig CM metrics
  - roc_auc_score throws on empty input when all hypotheses are NA after MIN_PTS filtering
  - Need guard upstream of sklearn call; subset-based evaluation (only score hypotheses valid in all tests)
  - Wald+MWU test results save fine, only the summary/metrics phase crashes
- [ ] Strip out consider_missing memory heuristics in core.py (no longer needed with phantom efficiency gains)
- [ ] Review simulation notebooks (shendure + cohen calibration, power, pairwise)
- [ ] Create all figures
- [ ] Takeshi data analysis (4th dataset)
  - Re-fitting all 4 orthos (obs/CM x NB/ZINB) with corrected negative controls
    (was 18, now 20 -- added 10:57784083:A:G:R:wC and 4:108065576:T:C:R:wC)
  - Also fixed sys.path from stale worktree back to main repo
- [ ] Code style corrections
- [ ] (v2) Port Wald SE computation from TF1 graph mode to TF2 eager + GradientTape for GPU acceleration
  - Current: TF1 tf.hessians + tf.map_fn per-obs scores, CPU-only (~1s/model cohen, ~0.04s shendure)
  - Target: tf.GradientTape.jacobian for vectorized per-obs scores on GPU, est. 10-50x speedup
  - Enables dense simulation power sweeps at scale

---

## Known bugs blocking progress

### ~~Bounds.from_ortho() fails on seelig_ortho_20260320~~ FIXED

Fixed in this session. Multiple issues in `from_ortho()`:
`float(current.max())` on DataFrame, missing `ret.` prefix, None ZI for NB-only.
All resolved with `_resolve_scalar`/`_resolve_df` helpers. Also fixed
`simple_count` NB convergence failure on near-constant data (seelig library:
95% of CREs have exactly 5 barcodes) -- falls back to `fixed_count` mode.

---

## Transfection reporter necessity evaluation

### Goal

Use simulation + power analysis to determine whether transfection reporters are
worth including in future scMPRA experiments, or whether correct statistical
modeling (ZINB with consider_missing) can compensate for their absence.

### Background

Transfection reporters add cloning complexity but provide two benefits:
1. Zero classification: distinguish "not transfected" from "transfected but not
   expressed"
2. Reduced zero-inflation: filtering removes structural zeros before modeling

Without a reporter, `consider_missing` re-introduces structural zeros, and the
ZINB ZI component must absorb them.

### Design

For each empirical dataset (shendure, cohen), simulate matched conditions:

| Condition | Reporter | consider_missing | Model | flatten_overtransfection |
|-----------|----------|-----------------|-------|--------------------------|
| A: with reporter | Yes | No | per NB/ZINB results | True |
| B: without reporter | No | Yes | per NB/ZINB results | True |

Seelig is already the "no reporter" case -- serves as empirical validation.

### Prerequisites

- NB vs ZINB evaluation complete (done -- see table above)
- Existing `de_novo_simulation` framework

### Metrics

- P-value calibration (QQ under null)
- Power at alpha=0.05 vs effect size
- AUROC / AUPRC
- Precision-recall curves at median-AUPRC replicate

### Open questions

- How to handle shendure's clonotype bottlenecking? `two_thirds_inactive`
  notebook handles this -- review and reuse.
- Number of simulation replicates: 5 is fast but noisy, 10-20 better.

### Note

This may be deferred to v2 / companion paper depending on timeline. If included
in v1, could serve as Fig 5 (see below).

---

## Figures

Note: 1a, 1b, etc are temporary labels. Give each figure a catchy name and
refer to it that way in directories and code.

### Fig 1: modeling approach and justification

- 1a: cartoon of stratified modeling approach (hand-drawn, not from code)
- 1b: table of run statistics -- CPU/GPU usage for different tasks (at the end)
- 1c: goodness of fit: NB vs ZINB AIC comparison
  - Data ready: `lrt_nb_vs_zinb_results.tsv`, `mu_stability_nb_vs_zinb.tsv`
  - Plotting: `lrt_nb_vs_zinb_plot.ipynb`
- 1d: bias-variance / ground truth recovery
  - Show that as cell number increases, parameter estimates converge to ground
    truth. Existing: `notebooks/results/bias/vsground_subsampling.ipynb`
    (synthetic data, 100/1000/10000 cells). MSE drops ~50x, bias near-zero at
    all sizes -- the interesting finding is how few cells suffice.
  - Needs: updated paths (old McCleary /gpfs/gibbs/ paths), re-run, saved SVGs.
  - Source data: `/nfs/roberts/project/pi_skr2/shared/tabula_data/simulated/fake_cres_v3_*.tsv`
- ~~1e: minimum data points / MIN_PTS justification~~ CUT -- a methods sentence
  plus supplementary table is sufficient.

### Fig S1: ancillary model justification

- s1a: coupon collector plot (formula matches simulations)
  - Done: `coupon_collector.ipynb`, SVG generated
- s1b: collision rate barchart
  - Done: `estimated_percent_conflict.ipynb`, SVG generated
- s1c: bounds plots (transfection + library models from extracted bounds)
  - Done: all 10 presets extracted, SVGs in `abstract_bounds/output/{cohen,shendure,seelig}/`

### Fig 2: statistical tests & re-analysis of existing datasets

- 2a: PRC and ROC curves of different stat tests (MWU, Wald OPG, Wald
  sandwich, permutation) for 3x datasets, activity hypotheses, "30% active"
  - Infrastructure exists (shendure + cohen calibration/power notebooks)
  - Needs re-running on phantom-zero orthos
- 2b: volcano plots from best test applied to real data
  - `hypothesis_testing_SHENDURE_empirical.ipynb` exists for shendure
  - Cohen + seelig volcano plots needed
  - BLOCKED on Wald fix for phantom / NB
- 2c: comparison of our calls to existing calls from original papers
  - Limit to shendure initially (most public results available)
  - Pull significance calls from Lalanne et al. 2024

### Fig 3: error control -- "the new stuff we discovered is real"

Discovery count + empirical FDR argument:
1. Show calibration on null sims (p-values uniform, FPR controlled at alpha)
2. Apply to real data, report D discoveries at chosen alpha
3. Given calibrated FPR of Z%, expect at most K = D * Z false positives

Pair the calibration plot (from simulation infrastructure) with the volcano
plot from 2b. This is the "we found X new things, at least Y are real" story.

### Fig 4: variant effect suitability across assay designs

Pairwise power heatmaps showing power as a function of baseline activity and
fold change, across datasets.

Existing:
- `cohen_pairwise_power_mwu.ipynb` + `replot_cohen_power.py`
- `shendure_pairwise_power_mwu.ipynb` + `replot_shendure_power.py`
- Output SVGs in respective `output/` directories

Needs re-running on phantom-zero orthos.

### Fig 5: design space exploration

Demonstrate the software's utility by simulating regimes not yet tested
experimentally. Minimal scope for v1: a single panel showing how power varies
with cell number (or another key design parameter).

Could incorporate the transfection reporter evaluation (see above) if that
analysis is completed in time. Otherwise, a straightforward power-vs-N curve
from `de_novo_simulation`.

### Supplementary: QC / convergence diagnostics

Model convergence diagnostics across fits. Infrastructure exists (`run_qc.py`,
`wrap_qc.sh`, seelig QC plots). Could strengthen methods section.

---

## Caveats

- **Wald SE computation requires GPU.** The TF2 GradientTape implementation
  (which replaced TF1 graph mode to fix a memory leak after ~400 models) has
  high per-observation dispatch overhead on CPU (~60s/model vs ~1.2s on GPU for
  cohen). CPU runs time out before completing even small datasets. The old TF1
  graph-mode code was CPU-viable (~1s/model) but leaked memory and segfaulted.
  Maintaining both code paths is not justified for v1. Users without GPU access
  should use MWU or bootstrap tests instead of Wald.

- **Seelig by_cell_type Wald is slow.** The 2 by_cell_type models have
  1344-column design matrices, making per-observation Jacobian computation
  ~3000s/model on H200 (~100 min total for both). Fine for a single empirical
  analysis, but expensive for power analysis (20 replicates x 2 models x
  3000s = 33h). MWU/bootstrap may be more practical for seelig cell-type
  power simulations. Seelig by_cre Wald (1344 models, 2 columns each) is fast.

---

## Dream features (post-publication)

Ideas that are out of scope for v1 but worth recording for a future version.

- **Empirical Bayes shrinkage.** Estimate the cross-CRE (or cross-cell-type)
  distribution of effect sizes from MLE results, then shrink noisy per-CRE
  estimates toward the population mean -- analogous to limma/DESeq2 dispersion
  shrinkage. Purely a post-processing step on existing ortho coefficients, no
  new fitting engine needed. Main benefit: stabilizes inference for low-count
  CREs without discarding them via MIN_PTS.

- **Full hierarchical Bayesian ZINB.** Partial pooling across CREs via
  hierarchical priors on regression coefficients and dispersion. Posterior
  contrasts replace Wald/MWU/bootstrap. Would unify the three testing
  frameworks and give principled uncertainty on the ZI mixing weight. Likely
  engine: NumPyro (JAX, GPU-native MCMC). Abandoned early in development due
  to scale: MCMC on thousands of coupled ZINB models is 10-100x slower than
  MLE, and consider_missing at Bayesian scale is probably intractable.
  Empirical Bayes (above) captures most of the benefit at a fraction of the
  cost.

- **Cell-type-specific variant effects.** Testing whether a variant's effect
  (CRE_A vs CRE_B) differs across cell types -- a CRE x cell-type interaction.
  The current stratified approach decouples these axes: by_cell_type gives
  variant effects within a cell type, by_cre gives cell-type effects within a
  CRE, but neither directly tests the interaction. For v1, MWU on the paired
  observations is sufficient (prelim results show Wald doesn't outperform MWU
  anyway). A future version could use a post-hoc contrast-of-contrasts approach:
  compute delta_k = beta_A - beta_B within each cell-type model (with proper
  within-model covariance from Wald), then test heterogeneity of delta_k across
  cell types via Cochran's Q. Cross-model independence is exact by construction
  (cell-type models are fitted independently), so no joint covariance estimation
  needed. Would pair well with empirical Bayes shrinkage for low-count CREs.
