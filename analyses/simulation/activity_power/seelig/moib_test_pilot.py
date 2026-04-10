#!/usr/bin/env python3
"""
Pilot: MOIB at hypothesis-test time.

STATUS: NEGATIVE RESULT -- do not promote this approach to scMPRAforge core.

Summary of finding (2026-04-10):
    Augmenting observed nonzero counts with MOI-prior-derived synthetic zeros
    INFLATES Type I error 2-5x at alpha=0.05. The inflation grows with
    baseline mu (null FPR ~0.10 at mu=0.10 up to ~0.25 at mu=3.0 for both
    Welch t-test and MWU). Apparent power gains over the naive (drop-zeros)
    baseline are therefore at least partially spurious.

Mechanism:
    Augmenting with E[n_transfected] zeros pushes effective N from ~10s-100s
    of nonzero observations to ~thousands. Welch's SE collapses with large N,
    so any sampling-noise mean difference between groups becomes "significant."
    Specifically: realized n_nonzero fluctuates binomially between groups even
    when mu_a == mu_b; that random difference shifts the augmented-array means
    by O(mu * sigma_binom / N), tiny in absolute terms but large relative to
    the now-collapsed SE. Inflation grows linearly in mu, matching the
    observed null FPR pattern. MWU is no better than t-test, ruling out a
    pure variance-assumption fix.

Why naive fixes don't work:
  - Augmenting with NB(mu_hat, theta) draws instead of pure zeros: circular,
    you'd be testing whether your own estimate of mu differs between groups.
  - Increasing N_REPS or alpha correction: doesn't address the bias, just
    moves the threshold.
  - The principled fix is fit-time MOIB (NB ortho with phantom weights ->
    Wald), which is already implemented in scMPRAforge core. Test-time
    augmentation cannot recover the variance properties of a proper Wald
    test on the conditional mu.

Conclusion:
    MOIB-at-test-time is rejected as a kwarg flag for mwu/ttest/ks/pseudobulk.
    For Seelig power analysis, use the naive drop-zeros condition (the same
    `has_reporter=False` path used for Shendure/Cohen deflated comparisons),
    and accept the reduced power as the cost of having no transfection
    reporter. This script is retained as a record of the negative result.

----------------------------------------------------------------------------
ORIGINAL DESIGN (kept for context):

Question: For a Seelig-like no-reporter dataset, does augmenting the observed
counts with MOI-prior-derived synthetic zeros recover statistical power
relative to the naive (drop-zeros) baseline?

Setup (keep it simple -- direct NB draws, no scMPRAforge simulator):
  - Two CREs: "A" (reference, mu = baseline) and "B" (variant, mu = baseline * 2^log2_fc)
  - Each CRE has a fixed number of barcodes (n_bc_per_cre)
  - For each (cell, barcode) pair, independently draw Bernoulli(p_transfected)
    where p_transfected = 1 - (1 - 1/n_lib)^moi (per-barcode transfection prob)
  - For transfected pairs, draw UMI count from NB(mu, theta)
  - "Observed" = rows with count > 0 (this mimics what Seelig actually sees)

Two conditions on the same simulated observed data:
  1. NAIVE: t-test / MWU on the observed nonzero values directly
  2. MOIB: compute n_transfected = n_cells * n_barcodes_per_cre * p_transfected
     per CRE, subtract n_nonzero, append that many zeros to the array, then
     run the test.

For each (baseline_mu, log2_fc) grid cell we repeat n_reps times and compute
TPR at alpha = 0.05. Comparing naive vs MOIB TPR tells us whether MOIB
augmentation is net-helpful.

The Seelig parameters are pulled from SEELIG_BOUNDS so the pilot runs on a
realistic cell count, MOI, library size, and dispersion. The MOI passed to
the MOIB step is the SAME MOI the simulator used -- this is "oracle MOI" but
NOT "oracle zero count." The MOIB step still has to estimate the zero count
by inference, it just isn't penalized for wrong MOI. A misspecified-MOI
sensitivity sweep is a reasonable follow-up if MOIB passes this pilot.

Usage:
    python moib_test_pilot.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_repo_root = Path(__file__).resolve().parents[4]
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

import scMPRAforge as scm

OUTPUT_DIR = Path(__file__).resolve().parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Seelig parameters (realistic)
# ---------------------------------------------------------------------------
SEELIG = scm.SEELIG_BOUNDS
MOI = float(SEELIG.get_effective_moi())          # ~73.57
N_LIB = int(SEELIG.total_uniq_mpra_bc)           # 6646
THETA = float(SEELIG.theta)                       # 1.628
# Cell count: use K562 (4478) as the test cell type
N_CELLS = int(SEELIG.cells_per_cell_type["K562"])
# Typical barcodes per CRE for Seelig -- rough estimate
N_BC_PER_CRE = 3

P_TRANSFECTED = 1.0 - (1.0 - 1.0 / N_LIB) ** MOI

print(f"Seelig parameters:")
print(f"  MOI               = {MOI:.2f}")
print(f"  n_library_bc      = {N_LIB}")
print(f"  theta (NB disp)   = {THETA:.3f}")
print(f"  n_cells (K562)    = {N_CELLS}")
print(f"  n_bc per CRE      = {N_BC_PER_CRE}")
print(f"  p_transfected/bc  = {P_TRANSFECTED:.5f}")
print(f"  n_transfected/CRE = {N_CELLS * N_BC_PER_CRE * P_TRANSFECTED:.1f}")
print()

# ---------------------------------------------------------------------------
# Sim parameters
# ---------------------------------------------------------------------------
# Baseline mus log-spaced across Seelig's observed activity range, skipping
# the tails.
BASELINE_MUS = np.geomspace(0.1, 3.0, 5)
# Positive log2 FCs -- test is symmetric so we skip negatives.
# log2_fc = 0 is included as a null calibration check (Type I error at alpha).
LOG2_FCS = np.array([0.0, 0.1, 0.2, 0.3, 0.5, 1.0])
N_REPS = 200
ALPHA = 0.05
RNG = np.random.default_rng(20260410)


def draw_nb_counts(mu, theta, n):
    """Draw n counts from NB with mean mu and dispersion theta.

    Uses the numpy parameterization p = theta / (theta + mu), n_param = theta.
    """
    if n <= 0:
        return np.empty(0, dtype=np.int64)
    p = theta / (theta + mu)
    return RNG.negative_binomial(theta, p, size=int(n))


def simulate_cre(mu, n_cells=N_CELLS, n_bc=N_BC_PER_CRE):
    """Simulate one CRE: return the observed (nonzero) UMI counts.

    Each (cell, barcode) pair is independently transfected with
    probability P_TRANSFECTED. Transfected pairs draw a count from NB.
    Observed = nonzero counts.

    Returns
    -------
    observed : np.ndarray
        Nonzero UMI counts (length = n_nonzero)
    n_transfected : int
        Number of transfected (cell, barcode) pairs (used by MOIB at
        the test step to estimate n_zeros to add back)
    """
    n_pairs = n_cells * n_bc
    # Deterministic (expected) n_transfected -- MOIB uses the expectation,
    # not the realization, because at test time we only know MOI and library
    # sizes, not which specific pairs transfected.
    n_transfected_expected = int(round(n_pairs * P_TRANSFECTED))
    # Actual realization for simulating observations
    n_transfected_actual = RNG.binomial(n_pairs, P_TRANSFECTED)
    counts = draw_nb_counts(mu, THETA, n_transfected_actual)
    observed = counts[counts > 0]
    return observed, n_transfected_expected


def naive_ttest(obs_a, obs_b):
    """t-test on observed nonzero values only."""
    if len(obs_a) < 2 or len(obs_b) < 2:
        return np.nan
    return ttest_ind(obs_a, obs_b, equal_var=False).pvalue


def moib_ttest(obs_a, obs_b, n_t_a, n_t_b):
    """Drop-zeros + MOIB augmentation + t-test.

    Append max(0, n_transfected - n_nonzero) zeros to each group, then ttest.
    """
    n_zeros_a = max(0, n_t_a - len(obs_a))
    n_zeros_b = max(0, n_t_b - len(obs_b))
    aug_a = np.concatenate([obs_a, np.zeros(n_zeros_a)])
    aug_b = np.concatenate([obs_b, np.zeros(n_zeros_b)])
    if len(aug_a) < 2 or len(aug_b) < 2:
        return np.nan
    return ttest_ind(aug_a, aug_b, equal_var=False).pvalue


def naive_mwu(obs_a, obs_b):
    if len(obs_a) < 2 or len(obs_b) < 2:
        return np.nan
    return mannwhitneyu(obs_a, obs_b, alternative="two-sided").pvalue


def moib_mwu(obs_a, obs_b, n_t_a, n_t_b):
    n_zeros_a = max(0, n_t_a - len(obs_a))
    n_zeros_b = max(0, n_t_b - len(obs_b))
    aug_a = np.concatenate([obs_a, np.zeros(n_zeros_a)])
    aug_b = np.concatenate([obs_b, np.zeros(n_zeros_b)])
    if len(aug_a) < 2 or len(aug_b) < 2:
        return np.nan
    return mannwhitneyu(aug_a, aug_b, alternative="two-sided").pvalue


# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------

rows = []
for baseline_mu in BASELINE_MUS:
    for log2_fc in LOG2_FCS:
        mu_a = baseline_mu
        mu_b = baseline_mu * 2 ** log2_fc

        rejects = {"naive_t": 0, "moib_t": 0, "naive_mwu": 0, "moib_mwu": 0}
        n_nan = {"naive_t": 0, "moib_t": 0, "naive_mwu": 0, "moib_mwu": 0}

        for rep in range(N_REPS):
            obs_a, n_t_a = simulate_cre(mu_a)
            obs_b, n_t_b = simulate_cre(mu_b)

            for label, fn in [
                ("naive_t",   lambda: naive_ttest(obs_a, obs_b)),
                ("moib_t",    lambda: moib_ttest(obs_a, obs_b, n_t_a, n_t_b)),
                ("naive_mwu", lambda: naive_mwu(obs_a, obs_b)),
                ("moib_mwu",  lambda: moib_mwu(obs_a, obs_b, n_t_a, n_t_b)),
            ]:
                p = fn()
                if np.isnan(p):
                    n_nan[label] += 1
                elif p < ALPHA:
                    rejects[label] += 1

        for label in rejects:
            n_valid = N_REPS - n_nan[label]
            tpr = rejects[label] / n_valid if n_valid > 0 else np.nan
            rows.append({
                "baseline_mu": baseline_mu,
                "log2_fc": log2_fc,
                "test": label,
                "tpr": tpr,
                "n_valid": n_valid,
                "n_nan": n_nan[label],
            })

        # Also compute the "typical" n_nonzero and n_zeros for context
        obs_a_sample, n_t_a_sample = simulate_cre(mu_a)
        print(
            f"baseline={baseline_mu:.3f} log2fc={log2_fc:.2f} "
            f"mu_a={mu_a:.3f} mu_b={mu_b:.3f}  "
            f"n_obs_a~{len(obs_a_sample)}/{n_t_a_sample}  "
            f"naive_t={rejects['naive_t']/N_REPS:.2f}  "
            f"moib_t={rejects['moib_t']/N_REPS:.2f}",
            flush=True,
        )

res = pd.DataFrame(rows)
res.to_parquet(OUTPUT_DIR / "moib_pilot_results.parquet")
print(f"\nSaved: {OUTPUT_DIR / 'moib_pilot_results.parquet'}")

# ---------------------------------------------------------------------------
# Plot: TPR curves per baseline_mu, one subplot per baseline, naive vs MOIB
# ---------------------------------------------------------------------------

fig, axes = plt.subplots(
    1, len(BASELINE_MUS), figsize=(4 * len(BASELINE_MUS), 4), sharey=True
)
if len(BASELINE_MUS) == 1:
    axes = [axes]

colors = {"naive_t": "coral", "moib_t": "steelblue",
          "naive_mwu": "lightsalmon", "moib_mwu": "lightsteelblue"}
linestyles = {"naive_t": "-", "moib_t": "-",
              "naive_mwu": "--", "moib_mwu": "--"}

for ax, baseline_mu in zip(axes, BASELINE_MUS):
    sub = res[np.isclose(res["baseline_mu"], baseline_mu)]
    for label in ["naive_t", "moib_t", "naive_mwu", "moib_mwu"]:
        d = sub[sub["test"] == label].sort_values("log2_fc")
        ax.plot(
            d["log2_fc"], d["tpr"],
            marker="o", markersize=4,
            color=colors[label], linestyle=linestyles[label],
            label=label,
        )
    ax.axhline(0.8, color="grey", linestyle=":", lw=0.8)
    ax.axhline(ALPHA, color="black", linestyle=":", lw=0.5)
    ax.set_title(f"baseline mu = {baseline_mu:.2f}", fontsize=10)
    ax.set_xlabel("log2 FC")
    ax.set_ylim(0, 1.05)
    ax.grid(alpha=0.3)

axes[0].set_ylabel("Power (TPR) @ alpha=0.05")
axes[0].legend(fontsize=7, loc="upper left")
fig.suptitle(
    f"MOIB test-time pilot -- Seelig K562 (n_cells={N_CELLS}, "
    f"MOI={MOI:.1f}, n_bc/CRE={N_BC_PER_CRE}, N_REPS={N_REPS})",
    fontsize=11, y=1.02,
)
plt.tight_layout()
out = OUTPUT_DIR / "moib_pilot_power_curves.svg"
fig.savefig(out, format="svg", bbox_inches="tight")
print(f"Saved: {out}")
