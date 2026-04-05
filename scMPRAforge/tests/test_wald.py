"""
Tests for the Wald test pipeline: phantom-zero weights + NB-only support.

Tests are written first (TDD). They exercise the numerical core bottom-up:
  1. Log-likelihood (TF graph)
  2. Covariance / SE estimation
  3. Worker function (_build_wald_precomp_for_subset)

No Dask required for unit tests -- these are pure numpy/TF.
"""
import numpy as np
import pandas as pd
import pytest

tf = pytest.importorskip("tensorflow")
pytest.importorskip("tensorzinb.tensorzinb")

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from scMPRAforge.core import (
    _zinb_loglik_tf,
    _setup_params_from_fit,
    _estimate_cov_se,
    _slice_se,
)


# ---------------------------------------------------------------------------
# Fixtures: synthetic data generators
# ---------------------------------------------------------------------------

def _make_zinb_data(n=300, k_nb=3, k_zi=2, seed=42):
    """
    Generate synthetic ZINB data with known parameters.
    Returns exog, infl, endog (as numpy), plus true params vector.
    """
    rng = np.random.default_rng(seed)

    # Design matrices
    exog = np.column_stack([
        np.ones(n),
        rng.standard_normal((n, k_nb - 1)),
    ])
    infl = np.column_stack([
        np.ones(n),
        rng.standard_normal((n, k_zi - 1)),
    ])

    # True parameters
    beta_mu = np.array([1.0] + [0.3] * (k_nb - 1))
    beta_pi = np.array([-1.0] + [0.2] * (k_zi - 1))
    log_theta = np.array([1.0])

    theta = np.exp(log_theta[0])
    mu = np.exp(exog @ beta_mu)
    pi_prob = 1.0 / (1.0 + np.exp(-(infl @ beta_pi)))

    # Sample from ZINB
    from scipy.stats import nbinom
    p_nb = theta / (theta + mu)
    y_nb = nbinom.rvs(theta, p_nb, random_state=rng)
    zi_mask = rng.random(n) < pi_prob
    y = np.where(zi_mask, 0, y_nb)

    params = np.concatenate([beta_mu, beta_pi, log_theta])
    endog = y.reshape(-1, 1).astype(np.float64)

    return exog, infl, endog, params


def _make_nb_data(n=300, k_nb=3, seed=42):
    """
    Generate synthetic NB data (no zero-inflation).
    Returns exog, endog (as numpy), plus true params vector [x_mu, log_theta].
    """
    rng = np.random.default_rng(seed)

    exog = np.column_stack([
        np.ones(n),
        rng.standard_normal((n, k_nb - 1)),
    ])

    beta_mu = np.array([1.0] + [0.3] * (k_nb - 1))
    log_theta = np.array([1.0])

    theta = np.exp(log_theta[0])
    mu = np.exp(exog @ beta_mu)

    from scipy.stats import nbinom
    p_nb = theta / (theta + mu)
    y = nbinom.rvs(theta, p_nb, random_state=rng)

    params = np.concatenate([beta_mu, log_theta])
    endog = y.reshape(-1, 1).astype(np.float64)

    return exog, endog, params


def _compress_with_weights(exog, infl, endog, seed=123):
    """
    Simulate phantom-zero compression: duplicate each row 1-5 times,
    then compress back by summing weights for identical rows.
    Returns compressed (exog, infl, endog, weights) that should give
    equivalent SEs to the expanded version.
    """
    rng = np.random.default_rng(seed)
    n = exog.shape[0]
    repeats = rng.integers(1, 6, size=n)

    exog_exp = np.repeat(exog, repeats, axis=0)
    infl_exp = np.repeat(infl, repeats, axis=0)
    endog_exp = np.repeat(endog, repeats, axis=0)

    # Now compress: since each original row is repeated identically,
    # compressing just gives back originals with weights = repeats
    weights = repeats.astype(np.float64)

    return exog, infl, endog, weights, exog_exp, infl_exp, endog_exp


# ---------------------------------------------------------------------------
# Phase 0: Baseline tests (should pass on current code)
# ---------------------------------------------------------------------------

class TestZINBLoglikBaseline:
    """Baseline: verify _zinb_loglik_tf produces finite values."""

    def test_gradient_finite(self):
        exog, infl, endog, params = _make_zinb_data(n=50, k_nb=3, k_zi=2)

        params_var = tf.Variable(params, dtype=tf.float64)
        exog_t = tf.constant(exog, dtype=tf.float64)
        infl_t = tf.constant(infl, dtype=tf.float64)
        endog_t = tf.constant(endog, dtype=tf.float64)

        with tf.GradientTape() as tape:
            tape.watch(params_var)
            ll_sum, ll_vec = _zinb_loglik_tf(
                params_var, exog_t, infl_t, endog_t, return_per_obs=True
            )
        grad_val = tape.gradient(ll_sum, params_var).numpy()
        ll_val = ll_sum.numpy()
        ll_v = ll_vec.numpy()

        assert np.isfinite(ll_val), f"ll_sum not finite: {ll_val}"
        assert np.all(np.isfinite(ll_v)), "some per-obs ll not finite"
        assert np.all(np.isfinite(grad_val)), "some gradient elements not finite"


class TestEstimateCovSeBaseline:
    """Baseline: _estimate_cov_se on ZINB data, no weights, no nb_only."""

    def test_ses_finite_positive(self):
        exog, infl, endog, params = _make_zinb_data(n=300, k_nb=3, k_zi=2)

        se, cov, diag = _estimate_cov_se(
            params, exog, infl, endog,
            ridge_h=1e-8, ridge_j=1e-12, opg_only=False,
        )

        assert np.all(np.isfinite(se)), f"non-finite SEs: {se}"
        assert np.all(se >= 0), f"negative SEs: {se}"
        # SEs should be in a reasonable range for this synthetic data
        assert np.all(se < 100), f"implausibly large SEs: {se}"


# ---------------------------------------------------------------------------
# Phase 1: NB-only tests
# ---------------------------------------------------------------------------

class TestZINBLoglikNBOnly:
    """_zinb_loglik_tf with nb_only=True."""

    def test_nb_only_finite(self):
        exog, endog, params = _make_nb_data(n=50, k_nb=3)
        infl = np.ones((50, 1), dtype=np.float64)

        params_var = tf.Variable(params, dtype=tf.float64)
        with tf.GradientTape() as tape:
            tape.watch(params_var)
            ll_sum, ll_vec = _zinb_loglik_tf(
                params_var,
                tf.constant(exog, dtype=tf.float64),
                tf.constant(infl, dtype=tf.float64),
                tf.constant(endog, dtype=tf.float64),
                return_per_obs=True, nb_only=True,
            )
        grad_val = tape.gradient(ll_sum, params_var).numpy()
        ll_val = ll_sum.numpy()

        assert np.isfinite(ll_val), f"NB-only ll not finite: {ll_val}"
        assert np.all(np.isfinite(grad_val)), "NB-only gradient not finite"
        assert grad_val.shape[0] == params.size

    def test_nb_only_differs_from_zinb(self):
        """NB-only ll should differ from ZINB ll on same data (ZI terms absent)."""
        exog, infl, endog, params_zinb = _make_zinb_data(n=50, k_nb=3, k_zi=2)

        k_nb = 3
        params_nb = np.concatenate([params_zinb[:k_nb], params_zinb[-1:]])
        infl_dummy = np.ones((50, 1), dtype=np.float64)

        exog_t = tf.constant(exog, dtype=tf.float64)
        endog_t = tf.constant(endog, dtype=tf.float64)

        ll_z = _zinb_loglik_tf(
            tf.constant(params_zinb, dtype=tf.float64),
            exog_t,
            tf.constant(infl, dtype=tf.float64),
            endog_t, return_per_obs=False,
        ).numpy()

        ll_n = _zinb_loglik_tf(
            tf.constant(params_nb, dtype=tf.float64),
            exog_t,
            tf.constant(infl_dummy, dtype=tf.float64),
            endog_t, return_per_obs=False, nb_only=True,
        ).numpy()

        assert np.isfinite(ll_z) and np.isfinite(ll_n)
        assert ll_z != ll_n, "ZINB and NB-only should give different ll values"


class TestSetupParamsNBOnly:
    """_setup_params_from_fit with nb_only=True model dict."""

    def test_nb_only_no_xpi(self):
        model_dict = {
            "nb_only": True,
            "weights": {
                "x_mu": pd.Series([1.0, 0.5, -0.3], index=["Intercept", "A", "B"]),
                "theta": np.array([1.0]),
            },
        }
        params = _setup_params_from_fit(model_dict)
        # Should be [x_mu(3), log_theta(1)] = 4 elements
        assert params.shape == (4,), f"Expected 4 params, got {params.shape}"
        np.testing.assert_allclose(params[:3], [1.0, 0.5, -0.3])
        np.testing.assert_allclose(params[3], 1.0)

    def test_zinb_has_xpi(self):
        model_dict = {
            "nb_only": False,
            "weights": {
                "x_mu": pd.Series([1.0, 0.5], index=["Intercept", "A"]),
                "x_pi": pd.Series([-1.0, 0.2], index=["Intercept", "B"]),
                "theta": np.array([1.0]),
            },
        }
        params = _setup_params_from_fit(model_dict)
        # Should be [x_mu(2), x_pi(2), log_theta(1)] = 5 elements
        assert params.shape == (5,), f"Expected 5 params, got {params.shape}"


class TestSliceSeNBOnly:
    """_slice_se with nb_only flag."""

    def test_nb_only_no_zi_block(self):
        # 3 NB + 1 theta = 4 total SEs
        se = np.array([0.1, 0.2, 0.3, 0.05])
        exog_cols = ["Intercept", "A", "B"]
        infl_cols = []  # NB-only has no ZI cols

        se_mu, se_pi, se_theta = _slice_se(se, exog_cols, infl_cols, nb_only=True)

        np.testing.assert_array_equal(se_mu, [0.1, 0.2, 0.3])
        assert len(se_pi) == 0
        np.testing.assert_array_equal(se_theta, [0.05])

    def test_zinb_standard(self):
        # 3 NB + 2 ZI + 1 theta = 6 total SEs
        se = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.05])
        exog_cols = ["Intercept", "A", "B"]
        infl_cols = ["Intercept_zi", "C"]

        se_mu, se_pi, se_theta = _slice_se(se, exog_cols, infl_cols)

        np.testing.assert_array_equal(se_mu, [0.1, 0.2, 0.3])
        np.testing.assert_array_equal(se_pi, [0.4, 0.5])
        np.testing.assert_array_equal(se_theta, [0.05])


# ---------------------------------------------------------------------------
# Phase 2: Weighted covariance tests
# ---------------------------------------------------------------------------

class TestEstimateCovSeWeighted:
    """_estimate_cov_se with phantom weights."""

    def test_weighted_matches_expanded(self):
        """Weighted SEs on compressed data should match unweighted on expanded data."""
        exog, infl, endog, params = _make_zinb_data(n=100, k_nb=3, k_zi=2)

        (exog_c, infl_c, endog_c, weights,
         exog_e, infl_e, endog_e) = _compress_with_weights(exog, infl, endog)

        # Unweighted on expanded data
        se_exp, _, _ = _estimate_cov_se(
            params, exog_e, infl_e, endog_e,
            ridge_h=1e-8, ridge_j=1e-12, opg_only=False,
        )

        # Weighted on compressed data
        se_wt, _, _ = _estimate_cov_se(
            params, exog_c, infl_c, endog_c,
            weights_np=weights,
            ridge_h=1e-8, ridge_j=1e-12, opg_only=False,
        )

        assert np.all(np.isfinite(se_wt)), f"weighted SEs not finite: {se_wt}"
        assert np.all(np.isfinite(se_exp)), f"expanded SEs not finite: {se_exp}"
        # Should match within 5% relative tolerance
        np.testing.assert_allclose(se_wt, se_exp, rtol=0.05,
            err_msg="Weighted SEs diverge from expanded unweighted SEs")


class TestEstimateCovSeNBOnly:
    """_estimate_cov_se with nb_only=True."""

    def test_nb_only_ses_finite(self):
        exog, endog, params = _make_nb_data(n=300, k_nb=3)
        infl = np.ones((300, 1), dtype=np.float64)

        se, cov, diag = _estimate_cov_se(
            params, exog, infl, endog,
            nb_only=True,
            ridge_h=1e-8, ridge_j=1e-12, opg_only=False,
        )

        assert np.all(np.isfinite(se)), f"NB-only SEs not finite: {se}"
        assert np.all(se >= 0), f"negative SEs: {se}"
        # k_nb + 1 (theta) = 4 SEs total
        assert se.shape[0] == 4, f"Expected 4 SEs, got {se.shape[0]}"

    def test_nb_only_cross_check_statsmodels(self):
        """Cross-check NB-only SEs against statsmodels NB2."""
        import statsmodels.discrete.discrete_model as smd

        exog, endog, params = _make_nb_data(n=500, k_nb=3, seed=99)
        infl = np.ones((500, 1), dtype=np.float64)

        # Fit with statsmodels to get "ground truth" SEs
        sm_nb = smd.NegativeBinomial(endog.ravel(), exog).fit(disp=False)
        if not sm_nb.converged:
            pytest.skip("statsmodels NB did not converge")

        # Use statsmodels MLE params for our SE computation
        sm_params = np.concatenate([sm_nb.params[:-1], [sm_nb.params[-1]]])
        # Note: statsmodels parameterizes alpha, we use log(theta) = log(1/alpha)
        # Adjust: our log_theta ~ -log(alpha)
        sm_alpha = sm_nb.params[-1]
        our_params = np.concatenate([
            sm_nb.params[:-1],  # x_mu coefficients
            [np.log(1.0 / sm_alpha)],  # log_theta = log(1/alpha)
        ])

        se, cov, diag = _estimate_cov_se(
            our_params, exog, infl, endog,
            nb_only=True,
            ridge_h=1e-8, ridge_j=1e-12, opg_only=False,
        )

        # Compare NB coefficient SEs (not theta -- different parameterization)
        k_nb = exog.shape[1]
        our_se_mu = se[:k_nb]
        sm_se_mu = sm_nb.bse[:k_nb]

        # Should be in same ballpark (within 50% -- different estimators)
        ratio = our_se_mu / sm_se_mu
        assert np.all(ratio > 0.3) and np.all(ratio < 3.0), (
            f"SE ratio out of range: {ratio}\n"
            f"Our SEs: {our_se_mu}\nSM SEs: {sm_se_mu}"
        )


class TestEstimateCovSeNBOnlyWeighted:
    """Combined: NB-only + weighted."""

    def test_nb_only_weighted_finite(self):
        exog, endog, params = _make_nb_data(n=100, k_nb=3)
        infl = np.ones((100, 1), dtype=np.float64)
        weights = np.random.default_rng(77).integers(1, 6, size=100).astype(np.float64)

        se, cov, diag = _estimate_cov_se(
            params, exog, infl, endog,
            weights_np=weights,
            nb_only=True,
            ridge_h=1e-8, ridge_j=1e-12, opg_only=False,
        )

        assert np.all(np.isfinite(se)), f"NB+weighted SEs not finite: {se}"
        assert np.all(se >= 0), f"negative SEs: {se}"


# ---------------------------------------------------------------------------
# Simulation-based end-to-end tests
# ---------------------------------------------------------------------------

def _simulate_two_group_nb(n_per_group=500, beta_treatment=0.0, seed=42):
    """
    Simulate NB data from a two-group design (reference + treatment).
    Returns exog, endog, true_params, and group labels.

    Design: y ~ Intercept + treatment
      intercept (log-mu for reference) = 2.0
      treatment effect = beta_treatment (log-fold-change)
      theta = exp(1.0) ~ 2.72
    """
    rng = np.random.default_rng(seed)
    from scipy.stats import nbinom

    n = 2 * n_per_group
    # Design matrix: intercept + treatment indicator
    exog = np.zeros((n, 2), dtype=np.float64)
    exog[:, 0] = 1.0  # intercept
    exog[n_per_group:, 1] = 1.0  # treatment group

    beta_mu = np.array([2.0, beta_treatment])
    log_theta = 1.0
    theta = np.exp(log_theta)

    mu = np.exp(exog @ beta_mu)
    p_nb = theta / (theta + mu)
    y = nbinom.rvs(theta, p_nb, random_state=rng)

    params = np.concatenate([beta_mu, [log_theta]])
    endog = y.reshape(-1, 1).astype(np.float64)

    return exog, endog, params


class TestWaldSimulation:
    """End-to-end: simulate known effects, verify Wald detects them."""

    def _run_wald_z_test(self, exog, endog, params):
        """Helper: compute SEs and return z-score + p-value for treatment coef."""
        from scipy.stats import chi2
        infl = np.ones((exog.shape[0], 1), dtype=np.float64)

        se, cov, diag = _estimate_cov_se(
            params, exog, infl, endog,
            nb_only=True,
            ridge_h=1e-8, ridge_j=1e-12, opg_only=False,
        )

        se_mu, _, _ = _slice_se(se, ["Intercept", "treatment"], [], nb_only=True)
        beta_treatment = params[1]
        se_treatment = se_mu[1]
        z = beta_treatment / se_treatment
        p = float(chi2.sf(z * z, 1))
        return z, p, se_treatment

    def test_strong_effect_detected(self):
        """Large treatment effect (2x fold change) should be significant."""
        exog, endog, params = _simulate_two_group_nb(
            n_per_group=500, beta_treatment=np.log(2.0), seed=42
        )
        z, p, se = self._run_wald_z_test(exog, endog, params)

        assert np.isfinite(z), f"z not finite: {z}"
        assert p < 0.01, f"Strong effect not detected: p={p:.4f}, z={z:.2f}"

    def test_null_effect_not_detected(self):
        """Zero treatment effect should NOT be significant."""
        exog, endog, params = _simulate_two_group_nb(
            n_per_group=500, beta_treatment=0.0, seed=99
        )
        z, p, se = self._run_wald_z_test(exog, endog, params)

        assert np.isfinite(z), f"z not finite: {z}"
        assert p > 0.05, f"Null effect falsely detected: p={p:.4f}, z={z:.2f}"

    def test_weighted_detects_same_as_unweighted(self):
        """Weighted (compressed) data should give same significance call."""
        from scipy.stats import chi2

        exog, endog, params = _simulate_two_group_nb(
            n_per_group=300, beta_treatment=np.log(1.5), seed=77
        )
        infl = np.ones((exog.shape[0], 1), dtype=np.float64)

        # Unweighted
        se_uw, _, _ = _estimate_cov_se(
            params, exog, infl, endog,
            nb_only=True, ridge_h=1e-8, ridge_j=1e-12, opg_only=False,
        )
        z_uw = params[1] / se_uw[1]
        p_uw = float(chi2.sf(z_uw * z_uw, 1))

        # Compress with weights
        (exog_c, infl_c, endog_c, weights,
         exog_e, infl_e, endog_e) = _compress_with_weights(exog, infl, endog)

        se_wt, _, _ = _estimate_cov_se(
            params, exog_c, infl_c, endog_c,
            weights_np=weights, nb_only=True,
            ridge_h=1e-8, ridge_j=1e-12, opg_only=False,
        )
        z_wt = params[1] / se_wt[1]
        p_wt = float(chi2.sf(z_wt * z_wt, 1))

        # Both should agree on significance
        assert (p_uw < 0.05) == (p_wt < 0.05), (
            f"Significance mismatch: unweighted p={p_uw:.4f}, weighted p={p_wt:.4f}"
        )
