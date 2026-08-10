"""Are scMPRA UMI counts overdispersed relative to Poisson?

This is the analysis that motivates using negative binomial models at all. It
is deliberately self-contained: it reads the published count tables directly
and fits with statsmodels, using no part of scMPRAforge. Nothing here depends
on the manuscript's fitting choices, so the conclusion cannot be an artifact
of them.

For each dataset and each cell type we fit two GLMs to the same design,

    umis_mpra_bc ~ C(cre_id)          log link

one Poisson and one negative binomial, and compare by AIC. The design
conditions on CRE identity, so any overdispersion reported here is dispersion
*within* a (cell type, CRE) group -- not the spread between CREs, which a
Poisson model would also have to explain but which is not the point at issue.

The headline statistic is the Pearson dispersion of the Poisson fit,
phi = X2/(n-p), which is 1 under Poisson and is the factor by which the counts
are more variable than it allows. Its null is simulated rather than taken from
chi-square, because most fitted means here are below 1 and the asymptotic
reference is not calibrated in that regime. The NB fit is reported alongside
as corroboration, but nothing depends on it converging.

Deliberately excluded: transfection-reporter handling, phantom or structural
zero modelling, and the per-dataset expansions used for the canonical fits.
Those matter for effect estimates but not for whether the count process is
overdispersed. Only the observed counts are used.

SCOPE -- shendure only. The published tables do not share a zero convention:

    shendure   materialises zeros (85% of rows)   -> analysed here
    cohen      observed rows only, min count 1    -> excluded
    seelig     observed rows only, min count 1    -> excluded

Fitting an untruncated model to a table that omits its zeros conditions on
detection. For cohen that inverts the result: counts concentrate on 1, giving
var/mean ~ 0.2 -- apparent *under*dispersion, an artifact of the table format
rather than a property of the assay. Under a zero-truncated likelihood cohen
is overdispersed as expected, so the conclusion generalises; it simply cannot
be shown with this analysis, which is why the scope is one dataset.
frac_zero is recorded per row so the distinction stays visible.

    python analyses/model_selection/overdispersion.py [--datasets shendure cohen]
"""
import argparse
import pathlib
import sys
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm

HERE = pathlib.Path(__file__).resolve().parent
DATA = pathlib.Path("/nfs/roberts/project/pi_skr2/shared/tabula_data")

# Shendure only. The other two published tables omit their zero rows (see the
# note above), so an untruncated fit to them measures the wrong thing. Their
# paths are kept for reference but are not analysed here.
TABLES = {
    "shendure": DATA / "shendure" / "shendure_processed.tsv",
}
TRUNCATED_TABLES = {
    "cohen": DATA / "cohen" / "retina_single_counting_u6.tsv",
    "seelig": DATA / "seelig" / "seelig_scmpra_umiwise.tsv.gz",
}
COLS = ["cre_id", "cell_type", "umis_mpra_bc"]

# A cell type needs enough CREs and observations for the comparison to mean
# anything; below this the NB dispersion is not usefully estimable.
MIN_ROWS = 200
MIN_CRES = 5


def dispersion(y, X, n_boot, rng):
    """Pearson dispersion of a Poisson fit, with a simulated null.

    phi = X2/(n-p) is 1 when the counts really are Poisson about the fitted
    means, and is the factor by which they are more variable when they are
    not. It needs only the Poisson fit, which converges reliably, so the
    headline claim does not depend on the fragile NB optimisation.

    The chi-square reference for X2 is not trusted here: most fitted means are
    well below 1, and under those conditions phi does not centre on 1 (it runs
    low). The null is therefore simulated -- draw Poisson counts at the fitted
    means, recompute phi -- which is both properly calibrated and a direct
    answer to "how overdispersed could this look by chance".
    """
    pois = sm.GLM(y, X, family=sm.families.Poisson()).fit()
    mu = np.asarray(pois.fittedvalues, dtype=float)
    n, p = len(y), int(pois.df_model) + 1
    phi = float(pois.pearson_chi2 / pois.df_resid)

    null = np.array([
        float(np.sum((rng.poisson(mu) - mu) ** 2 / np.maximum(mu, 1e-12))
              / pois.df_resid)
        for _ in range(n_boot)
    ]) if n_boot else np.array([])
    return pois, mu, phi, null


def fit_nb(y, X, pois):
    """Fit NB to the same design as an already-fitted Poisson.

    statsmodels' NB maximum likelihood diverges from a cold start with this
    many dummy columns, so the mean coefficients start at the Poisson solution
    and the dispersion at its method-of-moments estimate. Several optimizers
    are tried because convergence is design-dependent.
    """
    mu = np.asarray(pois.fittedvalues, dtype=float)
    alpha0 = float(((y - mu) ** 2 - mu).sum() / max((mu ** 2).sum(), 1e-12))
    start = np.append(np.asarray(pois.params, dtype=float), max(alpha0, 1e-3))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for method in ("lbfgs", "bfgs", "nm"):
            try:
                nb = sm.NegativeBinomial(y, X).fit(
                    start_params=start, method=method, maxiter=500, disp=0)
            except Exception:
                continue
            if nb.mle_retvals.get("converged") and np.isfinite(nb.llf):
                return nb, method, alpha0
    return None, None, alpha0


def analyse(name, path, rows, n_boot, rng):
    print(f"\n=== {name} ===", flush=True)
    df = pd.read_csv(path, sep="\t", usecols=COLS)
    df["umis_mpra_bc"] = pd.to_numeric(df["umis_mpra_bc"], errors="coerce").fillna(0)
    print(f"{len(df):,} rows, {df.cre_id.nunique():,} CREs, "
          f"{df.cell_type.nunique()} cell types", flush=True)

    for ct, sub in df.groupby("cell_type", observed=True):
        if len(sub) < MIN_ROWS or sub.cre_id.nunique() < MIN_CRES:
            print(f"  {str(ct)[:24]:26s} skipped (n={len(sub)}, "
                  f"cres={sub.cre_id.nunique()})", flush=True)
            continue
        y = sub["umis_mpra_bc"].to_numpy(float)
        X = sm.add_constant(
            pd.get_dummies(sub["cre_id"], drop_first=True, dtype=float),
            has_constant="add").to_numpy(float)

        pois, mu, phi, null = dispersion(y, X, n_boot, rng)
        nb, method, alpha0 = fit_nb(y, X, pois)
        # fano is the marginal variance/mean of the raw counts. It does NOT
        # condition on CRE, so it also contains genuine expression differences
        # between elements, which a Poisson GLM is entitled to explain through
        # its mean structure. phi is the honest effect size; fano is reported
        # only for context and always exceeds it.
        fano = float(y.var() / max(y.mean(), 1e-12))
        rec = dict(dataset=name, cell_type=str(ct), n=len(y),
                   n_cre=int(sub.cre_id.nunique()), mean=float(y.mean()),
                   var=float(y.var()), fano_marginal=fano, phi=phi,
                   phi_null_mean=float(null.mean()) if null.size else np.nan,
                   phi_null_max=float(null.max()) if null.size else np.nan,
                   phi_null_sd=float(null.std()) if null.size else np.nan,
                   frac_zero=float((y == 0).mean()),
                   frac_mu_below_1=float((mu < 1).mean()),
                   aic_pois=float(pois.aic), ll_pois=float(pois.llf),
                   alpha_mom=alpha0)
        if nb is not None:
            alpha = float(np.asarray(nb.params)[-1])
            rec.update(aic_nb=float(nb.aic), ll_nb=float(nb.llf),
                       alpha_nb=alpha, nb_method=method,
                       delta_aic=float(pois.aic - nb.aic))
            print(f"  {str(ct)[:24]:26s} n={len(y):>8,} phi={phi:>8.2f} "
                  f"(null max {null.max():.2f})  alpha={alpha:6.3f}  "
                  f"dAIC={pois.aic - nb.aic:>13,.0f}", flush=True)
        else:
            rec.update(aic_nb=np.nan, ll_nb=np.nan, alpha_nb=np.nan,
                       nb_method="did not converge", delta_aic=np.nan)
            print(f"  {str(ct)[:24]:26s} n={len(y):>8,} phi={phi:>8.2f} "
                  f"(null max {null.max():.2f})   NB did not converge", flush=True)
        rows.append(rec)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="*", default=list(TABLES))
    ap.add_argument("--n-boot", type=int, default=200,
                    dest="n_boot", help="Poisson null draws per cell type")
    ap.add_argument("--out", type=pathlib.Path, default=HERE / "overdispersion.tsv")
    args = ap.parse_args()

    rng = np.random.default_rng(0)
    rows = []
    for name in args.datasets:
        path = TABLES[name]
        if not path.exists():
            print(f"skip {name}: {path} not found", file=sys.stderr)
            continue
        analyse(name, path, rows, args.n_boot, rng)

    assert rows, "no cell types analysed"
    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"\nwrote {args.out.name} ({len(out)} rows)")

    ok = out.dropna(subset=["delta_aic"])
    print(f"\nNB preferred over Poisson in {int((ok.delta_aic > 0).sum())} "
          f"of {len(ok)} converged fits")
    print(f"Pearson dispersion phi: min {out.phi.min():.2f}, "
          f"median {out.phi.median():.2f}, max {out.phi.max():.2f}")
    if out.phi_null_max.notna().any():
        print(f"largest phi seen in {args.n_boot} Poisson simulations: "
              f"{out.phi_null_max.max():.2f}")
    print(f"marginal var/mean (does not condition on CRE): "
          f"median {out.fano_marginal.median():.1f}")
    if len(ok) < len(out):
        print(f"NB failed to converge for {len(out) - len(ok)} cell types "
              "(reported as NaN, not dropped)")


if __name__ == "__main__":
    main()
