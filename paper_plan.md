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

### Conclusions

- **Seelig (CM)**: NB is the right model. Zero models significant, zero AIC
  wins, mu parameters identical to 4 decimal places.

- **Cohen (obs)**: ZINB is warranted. Majority significant under Bonferroni,
  mu stable (~3% shift).

- **Shendure (obs)**: NB is probably sufficient. Only 7% survive Bonferroni,
  BIC agrees (3% prefer ZINB), 12% mu shift.

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

- [ ] Consolidate md files & revise figure plan to match analysis goals (this task)
- [ ] Map analysis flow, make sure everything is organized
- [ ] Create new data directory
- [ ] Organize pre-processing
- [ ] Re-fit models with phantom zeroes (all datasets, obs + cm)
- [ ] Re-extract bounds (requires from_ortho bug fix first -- see below)
- [ ] Fix Wald for phantom and Wald for NB
- [x] Delete pre-phantom run_stats/run_summary files (obsolete after phantom refit)
- [ ] Strip out consider_missing memory heuristics in core.py (no longer needed with phantom efficiency gains)
- [ ] Review simulation notebooks (shendure + cohen calibration, power, pairwise)
- [ ] Create all figures
- [ ] Code style corrections

---

## Known bugs blocking progress

### Bounds.from_ortho() fails on seelig_ortho_20260320

```
TypeError: float() argument must be a string or a real number, not 'Series'
```

In `Bounds.from_ortho()`, `current.max()` returns a Series (one max per column)
when `current` is a DataFrame, but `float()` expects a scalar. The NB parameter
storage format changed; earlier orthos store 1-D structures, seelig stores 2-D
DataFrames.

**Impact**: Cannot extract bounds from seelig ortho (or any new ortho with the
updated format). Existing shendure/cohen presets unaffected.

**Workaround used**: Collision rate computed directly from raw TSV
(`seelig_collision_rate.py`).

**Fix**: In `Bounds.from_ortho()` min/max loop, replace `float(current.max())`
with something robust to both 1-D and 2-D, e.g.
`float(np.nanmax(current.values))`. But confirm parameter storage format first.

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
  - Plotting routines exist but BLOCKED on `from_ortho` bug fix
  - Need to re-extract bounds for all datasets post-phantom

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
