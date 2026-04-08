#!/usr/bin/env python
"""
Likelihood Ratio Test: NB vs ZINB for consider_missing orthos.

For each (dataset, model_type) pair, computes per-model LRT statistics:
  Λ = 2 × (LL_ZINB - LL_NB)
  p-value from χ²(q), where q = df_ZINB - df_NB (number of ZI parameters)

This is the conservative test — the true null distribution is a 50:50 mixture
of χ²(0) and χ²(q), so χ²(q) p-values are upper bounds.  If significant even
under the conservative test, ZI is definitively needed.

All pairs use phantom-zero compressed orthos from tabula_data_new.

Available pairs:
  - shendure obs: ZINB shendure_obs_zinb_phantom  vs NB shendure_obs_nb_phantom
  - shendure cm:  ZINB shendure_cm_zinb_phantom   vs NB shendure_cm_nb_phantom
  - cohen obs:    ZINB cohen_obs_zinb_phantom_20260401  vs NB cohen_obs_nb_phantom_20260401
  - cohen cm:     ZINB cohen_cm_zinb_phantom_20260401   vs NB cohen_cm_nb_phantom_20260401
  - seelig cm:    ZINB seelig_cm_zinb_phantom     vs NB seelig_cm_nb_phantom
"""
import pickle
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

data_root = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new")

pairs = {
    "shendure_obs": {
        "zinb": data_root / "shendure" / "shendure_obs_zinb_phantom",
        "nb":   data_root / "shendure" / "shendure_obs_nb_phantom",
    },
    "shendure_cm": {
        "zinb": data_root / "shendure" / "shendure_cm_zinb_phantom",
        "nb":   data_root / "shendure" / "shendure_cm_nb_phantom",
    },
    "cohen_obs": {
        "zinb": data_root / "cohen" / "cohen_obs_zinb_phantom_20260401",
        "nb":   data_root / "cohen" / "cohen_obs_nb_phantom_20260401",
    },
    "cohen_cm": {
        "zinb": data_root / "cohen" / "cohen_cm_zinb_phantom_20260401",
        "nb":   data_root / "cohen" / "cohen_cm_nb_phantom_20260401",
    },
    "seelig_cm": {
        "zinb": data_root / "seelig" / "seelig_cm_zinb_phantom",
        "nb":   data_root / "seelig" / "seelig_cm_nb_phantom",
    },
    "takeshi_obs": {
        "zinb": data_root / "takeshi" / "takeshi_obs_zinb_phantom",
        "nb":   data_root / "takeshi" / "takeshi_obs_nb_phantom",
    },
    "takeshi_cm": {
        "zinb": data_root / "takeshi" / "takeshi_cm_zinb_phantom",
        "nb":   data_root / "takeshi" / "takeshi_cm_nb_phantom",
    },
    "cohen_obsingle": {
        "zinb": data_root / "cohen" / "cohen_obsingle_zinb_phantom",
        "nb":   data_root / "cohen" / "cohen_obsingle_nb_phantom",
    },
    "seelig_cm_moib": {
        "zinb": data_root / "seelig" / "seelig_cm_moib_zinb_phantom",
        "nb":   data_root / "seelig" / "seelig_cm_moib_phantom",
    },
}

def load_models(path, filename):
    with open(path / filename, "rb") as f:
        return pickle.load(f)

def lrt_table(zinb_models, nb_models, label):
    """Compute per-model LRT statistics."""
    keys = sorted(set(zinb_models.model.keys()) & set(nb_models.model.keys()))
    if set(zinb_models.model.keys()) != set(nb_models.model.keys()):
        missing_zinb = set(nb_models.model.keys()) - set(zinb_models.model.keys())
        missing_nb = set(zinb_models.model.keys()) - set(nb_models.model.keys())
        if missing_zinb:
            print(f"  WARNING: {len(missing_zinb)} models in NB but not ZINB")
        if missing_nb:
            print(f"  WARNING: {len(missing_nb)} models in ZINB but not NB")

    rows = []
    for k in keys:
        zm = zinb_models.model[k]
        nm = nb_models.model[k]

        ll_zinb = zm["llf_total"]
        ll_nb = nm["llf_total"]
        df_zinb = zm["df_model_total"]
        df_nb = nm["df_model_total"]
        q = df_zinb - df_nb  # number of ZI parameters
        n = zm["num_sample"]

        lambda_lr = 2 * (ll_zinb - ll_nb)
        p_conservative = stats.chi2.sf(lambda_lr, df=q) if lambda_lr > 0 else 1.0

        # AIC / BIC for comparison
        aic_zinb = zm.get("aic_total", 2 * df_zinb - 2 * ll_zinb)
        aic_nb = nm.get("aic_total", 2 * df_nb - 2 * ll_nb)
        bic_zinb = df_zinb * np.log(n) - 2 * ll_zinb
        bic_nb = df_nb * np.log(n) - 2 * ll_nb

        rows.append({
            "model": k,
            "n": n,
            "ll_zinb": ll_zinb,
            "ll_nb": ll_nb,
            "df_zinb": df_zinb,
            "df_nb": df_nb,
            "q": q,
            "lambda_lr": lambda_lr,
            "p_conservative": p_conservative,
            "aic_zinb": aic_zinb,
            "aic_nb": aic_nb,
            "delta_aic": aic_nb - aic_zinb,  # positive = ZINB preferred
            "bic_zinb": bic_zinb,
            "bic_nb": bic_nb,
            "delta_bic": bic_nb - bic_zinb,  # positive = ZINB preferred
        })

    return pd.DataFrame(rows)


def summarize(df, label):
    """Print summary statistics for an LRT table."""
    n_models = len(df)
    sig_005 = (df["p_conservative"] < 0.05).sum()
    sig_001 = (df["p_conservative"] < 0.01).sum()
    sig_bonf = (df["p_conservative"] < 0.05 / n_models).sum()
    zinb_aic_wins = (df["delta_aic"] > 0).sum()
    zinb_bic_wins = (df["delta_bic"] > 0).sum()

    print(f"\n{'='*70}")
    print(f"  {label}  (n={n_models} models)")
    print(f"{'='*70}")
    print(f"  LRT (conservative χ² test, H₀: NB is sufficient):")
    print(f"    Significant at α=0.05:          {sig_005:>5d} / {n_models}  ({100*sig_005/n_models:.1f}%)")
    print(f"    Significant at α=0.01:          {sig_001:>5d} / {n_models}  ({100*sig_001/n_models:.1f}%)")
    print(f"    Significant (Bonferroni 0.05):  {sig_bonf:>5d} / {n_models}  ({100*sig_bonf/n_models:.1f}%)")
    print(f"  Λ (test statistic):")
    print(f"    median={df['lambda_lr'].median():.2f}  mean={df['lambda_lr'].mean():.2f}  "
          f"max={df['lambda_lr'].max():.2f}  min={df['lambda_lr'].min():.2f}")
    print(f"  ΔAIC (positive = ZINB preferred):")
    print(f"    median={df['delta_aic'].median():.2f}  mean={df['delta_aic'].mean():.2f}")
    print(f"    ZINB wins: {zinb_aic_wins}/{n_models} ({100*zinb_aic_wins/n_models:.1f}%)")
    print(f"  ΔBIC (positive = ZINB preferred):")
    print(f"    median={df['delta_bic'].median():.2f}  mean={df['delta_bic'].mean():.2f}")
    print(f"    ZINB wins: {zinb_bic_wins}/{n_models} ({100*zinb_bic_wins/n_models:.1f}%)")


if __name__ == "__main__":
    all_results = {}

    for pair_name, paths in pairs.items():
        print(f"\n--- {pair_name} ---")
        print(f"  ZINB: {paths['zinb']}")
        print(f"  NB:   {paths['nb']}")

        for model_type in ["by_cre", "by_cell_type"]:
            label = f"{pair_name} / {model_type}"
            print(f"\n  Loading {model_type}...")
            zinb = load_models(paths["zinb"], f"{model_type}.pkl")
            nb = load_models(paths["nb"], f"{model_type}.pkl")

            df = lrt_table(zinb, nb, label)
            summarize(df, label)
            all_results[label] = df

    # Save combined results
    out_path = Path(__file__).parent / "lrt_nb_vs_zinb_results.tsv"
    combined = pd.concat(
        {k: v for k, v in all_results.items()},
        names=["comparison", "idx"],
    ).reset_index(level="comparison").reset_index(drop=True)
    combined.to_csv(out_path, sep="\t", index=False)
    print(f"\nResults saved to {out_path}")
