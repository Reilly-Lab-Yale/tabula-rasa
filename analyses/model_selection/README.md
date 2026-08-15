# Model selection

Which count distribution describes scMPRA UMIs. `lrt_nb_vs_zinb.py` and
`zi_zero_decomposition.py` cover the NB-vs-ZINB comparison across datasets;
the two below are cited directly by the manuscript.

## plot_nb_vs_zinb_shendure.py

Generates the manuscript's worked-example figure: per-fit
`dAIC = AIC_ZINB - AIC_NB` for Lalanne et al. in both stratification
directions, ten by-cell-type fits shown individually and 208 by-CRE fits as a
distribution. Positive favours NB.

Reads the saved orthos directly (`shendure_obs_nb` and its zinb
counterpart) with a tolerant unpickler, so it needs neither `scMPRAforge` nor
a refit -- only the recorded per-fit scalars.

    /opt/anaconda3/envs/data-analysis-generic/bin/python \
        analyses/model_selection/plot_nb_vs_zinb_shendure.py

Writes `output/nb_vs_zinb_shendure.{svg,png}`; the manuscript pulls the SVG in
via its sync manifest.

Result: NB preferred in 8 of 10 cell types and 76% of 208 CREs. The by-CRE
panel uses a symmetric log axis because a handful of CREs favour ZINB by more
than 1000 AIC units, which on a linear axis collapses the bulk into a spike.
Those outliers are the point of keeping ZINB support in the package, so they
are shown rather than clipped.

## overdispersion.py

Whether scMPRA counts require an overdispersed model at all -- the evidence
behind Methods "Overdispersion of scMPRA counts". Deliberately independent of
`scMPRAforge`: it reads a published count table and fits with `statsmodels`
only, so the conclusion cannot be an artifact of the manuscript's own fitting
choices.

Per cell type, a Poisson GLM of `umis_mpra_bc ~ C(cre_id)` gives the Pearson
dispersion `phi = X2/(n-p)`, which is 1 when counts are Poisson about their
fitted means. A negative binomial fit to the same design corroborates by AIC.

Outputs `overdispersion.tsv`.

Result: `phi` runs 6.1 to 28.4 across the ten Lalanne et al. cell types
(median 16.4), against a maximum of 1.06 in 2,000 draws simulated from the
fitted Poisson. NB is preferred in 10/10, with dispersion 2.5-4.8.

Run with an environment that has statsmodels:

    /opt/anaconda3/envs/data-analysis-generic/bin/python analysis/overdispersion.py

### Data

Reads the published count tables under
`/nfs/roberts/project/pi_skr2/shared/tabula_data/`. On a machine without that
mount, `/etc/synthetic.conf` can map `/nfs` to a local mirror (macOS creates
the link at boot); only the Lalanne et al. table is needed for the default run.

### Caveats

- **Lalanne et al. only, by necessity.** The three published tables do not
  share a zero convention: Lalanne et al. materialises zero rows (85% of the
  table), while Zhao et al. and Yin et al. list observed rows only, with a
  minimum count of 1. Fitting an untruncated model to a table that omits its
  zeros conditions on detection -- for Zhao et al. that inverts the result,
  giving variance/mean around 0.2 and apparent *under*dispersion, which is an
  artifact of the file format and not a property of the assay. Under a
  zero-truncated likelihood both are overdispersed as expected (spot-checked
  on random CRE panels, negative binomial preferred in 12/12), so the
  conclusion generalises; it simply cannot be shown with this analysis.
  `frac_zero` is recorded per row so the distinction is visible.
- **`phi` conditions on CRE; the marginal variance/mean does not.** The
  marginal ratio (median 317 here) also contains genuine expression
  differences between elements, which the Poisson mean structure already
  absorbs. It is reported as `fano_marginal` for context only, and always
  exceeds `phi`.
- **The null is simulated, not taken from chi-square.** About 85% of fitted
  means fall below 1, where the chi-square reference for `X2` is not
  calibrated -- simulation puts the null near 0.8 rather than 1.
- **Negative binomial fitting is fragile.** `sm.NegativeBinomial` returns a
  `nan` likelihood from a cold start with a few hundred dummy columns, so the
  fit starts from the Poisson solution with a method-of-moments dispersion and
  tries several optimizers. Failures are recorded as `NaN` and counted in the
  summary rather than dropped. Nothing in the headline claim depends on it.
- **`freq_weights` is silently ignored** by `statsmodels`' truncated count
  models: it is accepted through `**kwargs` and never applied, returning a
  plausible but wrong likelihood. Collapsing these tables to per-`(cre, count)`
  weights is therefore not available as a shortcut.
