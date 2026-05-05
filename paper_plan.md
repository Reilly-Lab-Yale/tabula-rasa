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
| Cohen     | obs (coarse)    | 54%                | 82.5%            | 2.87%              |
| Cohen     | obs (obsingle)  | 0%                 | 0%               | ~0%                |
| Shendure  | obs             | 6.7%               | 25%              | 12.3%              |
| Takeshi   | obs             | 0%                 | 0%               | ~0%                |

### Conclusions and canonical model choices

- **Seelig (MOIB)**: NB is the right model. MOIB improves QC and LRT over
  plain CM. Zero models significant, mu parameters identical to 4 decimal places.
  **Canonical: `SEELIG_BOUNDS` = MOIB NB** (`seelig_cm_moib_nb_phantom`)

- **Cohen (obs, obsingle)**: Under CRE-coarse reporter expansion, ZINB appeared
  warranted (54% Bonferroni sig). However, coarse expansion floods zeros
  (N barcodes per U6 detection), crushing mu estimates (median 0.001 vs true
  0.88) and distorting the ZI component. Under obsingle expansion (1 zero per
  U6 detection), mu estimates recover (median 0.88), QC r-values improve
  dramatically (by_cell_type: 0.59-0.72 vs negative under coarse), and LRT/AIC
  confirm NB is sufficient (0% significant, 0% AIC).
  **Canonical: `COHEN_BOUNDS` = obsingle NB** (`cohen_obsingle_nb_phantom`)

- **Shendure (obs)**: NB is sufficient. Only 7% survive Bonferroni,
  BIC agrees (3% prefer ZINB by_cre, 0% by_cell_type).
  **Canonical: `SHENDURE_BOUNDS` = obs NB** (`shendure_obs_nb_phantom`)

- **Takeshi (obs)**: NB is the right model. Zero significant under any test
  (LRT, AIC, BIC), both by_cre (150 models) and by_cell_type (3 models).
  ZINB lambda negative for by_cell_type -- NB fits strictly better.
  **Canonical: `TAKESHI_BOUNDS` = obs NB** (`takeshi_obs_nb_phantom`)

These canonical choices are set as default aliases in `scMPRAforge.core` and
will be used to parameterize all downstream simulations.

(generally, ZI has a very hard time pulling out structural zeroes because of the identifiability issue: nb is too flexible, and the means are too low to allow separation. we see this in how CM shendure does OK for ZINB because shendure has high averages. we can ignore this and just do nb for shendure obs becuase it has a fine reporter. we can probably do the same (simple nb) for cohen obsingle, which is better than cohen obs. for seelig, we don't have the luxury of a tfection reporter AND we have low counts causing the identifyability problem, so we introduce the MOI prior).

### Cohen obsingle expansion -- impact on power analysis

Switching from CRE-coarse to obsingle reporter expansion dramatically changed
Cohen simulation results (5x5 activity, 70% inactive):

| Test | Coarse ZINB AUROC | Obsingle NB AUROC | Delta |
|------|-------------------|-------------------|-------|
| MWU  | 0.592             | 0.880             | +0.29 |
| Wald | 0.739             | 0.751             | +0.01 |

Key observations:
- MWU jumps ~29 AUROC points -- coarse expansion was poisoning rank-based tests
  by flooding identical zeros across all barcodes of a CRE
- Wald is similar but noisier under obsingle (std=0.13 vs 0.07), with GT draw 0
  at 0.89 AUROC vs draws 1-4 at 0.58-0.65. Needs investigation.
- MWU now beats Wald under obsingle NB (reversal from coarse ZINB where Wald
  was better). This is consistent with the other datasets (shendure, seelig)
  where MWU also outperforms Wald for activity hypotheses.

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
- [x] Fix _classifier_summary empty-array crash on seelig CM metrics
  - roc_auc_score throws on empty input when all hypotheses are NA after MIN_PTS filtering
  - Need guard upstream of sklearn call; subset-based evaluation (only score hypotheses valid in all tests)
  - Wald+MWU test results save fine, only the summary/metrics phase crashes
- [x] MOI-based phantom zero downweighting (see section below)
- [x] Strip out consider_missing memory heuristics in core.py (no longer needed with phantom efficiency gains)
- [ ] Add bulk to seelig 5x5
- [ ] Wire MOIB into Wald (currently not implemented for CM phantom weights in
  `precompute_wald` -- seelig Wald results are not trustworthy until this lands;
  treat current seelig Wald numbers as placeholder)
- [ ] Create all figures
- [ ] Takeshi data analysis (4th dataset)
  - [x] Re-fitting all 4 orthos (obs/CM x NB/ZINB) with corrected negative controls
    (was 18, now 20 -- added 10:57784083:A:G:R:wC and 4:108065576:T:C:R:wC)
  - [ ] power analysis
- [ ] Code style corrections
- [ ] (v2) Port Wald SE computation from TF1 graph mode to TF2 eager + GradientTape for GPU acceleration
  - Current: TF1 tf.hessians + tf.map_fn per-obs scores, CPU-only (~1s/model cohen, ~0.04s shendure)
  - Target: tf.GradientTape.jacobian for vectorized per-obs scores on GPU, est. 10-50x speedup
  - Enables dense simulation power sweeps at scale

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

### Design (revised)

Simulate with reporter (shendure-like or cohen-like bounds), then fit the
SAME simulated data twice:

| Condition | Fit strategy | Model |
|-----------|-------------|-------|
| A: with reporter | obs-condition, reporter-informed zeros | NB |
| B: without reporter | strip reporter, consider_missing + MOIB | NB |

Identical ground truth, identical cells, identical simulated counts. The only
variable is whether the reporter information is used at fit time. Any power
difference is purely attributable to the reporter.

This is cleaner than the old design (which simulated two separate datasets
with different parameters). No confounds from different simulation settings.

Seelig is already the "no reporter" case -- serves as empirical validation
that MOIB produces usable results without a reporter.

### Prerequisites

- NB vs ZINB evaluation complete (done -- see table above)
- MOIB implementation and validation (done -- shendure convergence r=0.994)
- Existing `de_novo_simulation` framework

### Metrics

- P-value calibration (QQ under null)
- Power at alpha=0.05 vs effect size
- AUROC / AUPRC
- Precision-recall curves at median-AUPRC replicate
- Discovery count comparison (how many additional positives/negatives
  detected with vs without reporter)

### Open questions

- How to handle shendure's clonotype bottlenecking? `two_thirds_inactive`
  notebook handles this -- review and reuse.
- Number of simulation replicates: 5 is fast but noisy, 10-20 better.
- The MOIB 1.8x mu offset for shendure may or may not affect power --
  parameter accuracy != discrimination ability. Must test empirically.
- **Effect-size priors for contextualizing power curves**: The FC-at-80%-power
  summary shows what effect sizes are *detectable*, but doesn't say whether
  those effect sizes *exist* in real biology. Extracting empirical FC
  distributions from the same datasets is circular (observed FCs are already
  truncated by the power limitations being characterized). An external prior
  on CRE effect sizes -- e.g. from new experimental data with known
  ground-truth perturbations -- would let us say "X% of real effects are
  detectable with/without reporter." PIN: Mackenzie is generating wet-lab
  data that could provide this prior (PCR-based, pending).

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
    truth. Existing: `analyses/model_selection/bias/vsground_subsampling.ipynb`
    (synthetic data, 100/1000/10000 cells). MSE drops ~50x, bias near-zero at
    all sizes -- the interesting finding is how few cells suffice.
  - Needs: updated paths (old McCleary /gpfs/gibbs/ paths), re-run, saved SVGs.
  - Source data: `/nfs/roberts/project/pi_skr2/shared/tabula_data/simulated/fake_cres_v3_*.tsv`
- ~~1e: minimum data points / MIN_PTS justification~~ CUT -- a methods sentence
  plus supplementary table is sufficient.

### Fig S1: ancillary model justification

- s1a: coupon collector plot (formula matches simulations, validating assumption that barcode complexity overwhelms multi-transfection. worst for cohen not because of hi moi but because of low barcode complexity)
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

- **MOIB not wired into Wald yet.** `precompute_wald` reconstructs CM design
  matrices without applying MOI-based phantom downweighting, so Wald on
  CM-fitted orthos uses unconditional (deflated) mu and inflated phantom
  weights. This breaks seelig Wald specifically (CM is mandatory there).
  Treat current seelig Wald p-values, AUROCs, and power numbers as
  placeholders -- they are not real results and must be regenerated after
  MOIB is plumbed through `precompute_wald`. Shendure/cohen Wald (obs
  condition) is unaffected.

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

- **Per-cell-type MOI for MOIB.** Currently MOIB uses a single global MOI
  (from `describe_transfection`) to compute P(transfected) for all groups.
  Different cell types may have different transfection efficiencies, yielding
  different effective MOIs. A future version could estimate per-cell-type MOI
  and use cell-type-specific P(transfected) when computing phantom zero weights
  in `_cm_group_totals`. Would require parameterizing `moi_correction` as a
  dict keyed by cell_type rather than a single scalar.
