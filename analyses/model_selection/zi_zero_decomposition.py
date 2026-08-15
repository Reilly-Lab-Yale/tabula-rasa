#!/usr/bin/env python
"""
Supplemental figure: Zero-inflation decomposition across datasets.

Shows that ZINB doesn't always add zeros -- in Shendure, the ZI component
absorbs structural zeros, allowing the NB mean to rise, which *reduces*
total predicted zeros relative to a plain NB fit.

For each dataset x cell type, computes:
  - Observed %zero (from training data)
  - NB-predicted %zero:   P(X=0|NB) = (theta/(theta+mu))^theta, averaged over CREs
  - ZINB-predicted %zero:  P(X=0|ZINB) = pi + (1-pi) * (theta_z/(theta_z+mu_z))^theta_z
  - Excess = observed - NB predicted (positive = real ZI signal)
  - ZINB shift = ZINB predicted - NB predicted (negative = ZINB predicts fewer zeros)

Requires obs NB and obs ZINB orthos for each dataset.
"""
import pickle, sys, numpy as np, pandas as pd
from pathlib import Path
from scipy.special import expit

sys.path.insert(0, '/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa-takeshi')
import scMPRAforge.core as scm
from dask.distributed import Client, LocalCluster

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
data_root = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new")
OUT_DIR = Path(__file__).resolve().parent

datasets = {
    "Shendure (obs)": {
        "nb":   data_root / "shendure" / "shendure_obs_nb",
        "zinb": data_root / "shendure" / "shendure_obs_zinb",
        "dir":  "shendure",
    },
    "Shendure (CM)": {
        "nb":   data_root / "shendure" / "shendure_cm_nb_phantom",
        "zinb": data_root / "shendure" / "shendure_cm_zinb_phantom",
        "dir":  "shendure",
    },
    "Takeshi (obs)": {
        "nb":   data_root / "takeshi" / "takeshi_obs_nb",
        "zinb": data_root / "takeshi" / "takeshi_obs_zinb",
        "dir":  "takeshi",
    },
    "Takeshi (CM)": {
        "nb":   data_root / "takeshi" / "takeshi_cm_nb_phantom",
        "zinb": data_root / "takeshi" / "takeshi_cm_zinb_phantom",
        "dir":  "takeshi",
    },
}

# ---------------------------------------------------------------------------
cluster = LocalCluster(n_workers=1, threads_per_worker=1,
                       processes=False, memory_limit="30GB")
client = Client(cluster)

def get_theta(params, key):
    t = params.theta[key]
    t = t.result() if hasattr(t, "result") else t
    return float(np.array(t).flatten()[0])

def get_nb_params(params, key):
    nb = params.nb[key]
    nb_df = nb.result() if hasattr(nb, "result") else nb
    return nb_df["mu"].values, get_theta(params, key)

def get_zi_probs(params, key):
    """Extract per-rep ZI probabilities from ZINB model.
    params.zi already stores probability-scale values (post-expit),
    NOT logit-scale. Do not apply expit again.
    """
    zi = params.zi[key]
    zi_df = zi.result() if hasattr(zi, "result") else zi
    return zi_df.values.flatten()

def nb_p_zero(mu, theta):
    """P(X=0|NB) per CRE."""
    return (theta / (theta + mu)) ** theta

# ---------------------------------------------------------------------------
all_rows = []

for ds_label, paths in datasets.items():
    print(f"\n=== {ds_label} ===", flush=True)

    # Load orthos
    nb_ortho = scm.ortho.load(client, str(paths["nb"].parent), paths["nb"].name)
    zinb_ortho = scm.ortho.load(client, str(paths["zinb"].parent), paths["zinb"].name)

    # Get observed zero rates from training data
    td = nb_ortho.training_data
    data = td.get_data(include_missing=False)
    if hasattr(data["umis_mpra_bc"].dtype, "fill_value"):
        data["umis_mpra_bc"] = data["umis_mpra_bc"].astype("int64")
    pdf = data[["cell_type", "umis_mpra_bc"]].compute()

    nb_params = nb_ortho.by_cell_type_parameters
    zinb_params = zinb_ortho.by_cell_type_parameters

    for ct in sorted(pdf["cell_type"].unique()):
        sub = pdf[pdf["cell_type"] == ct]
        obs_pct_zero = 100 * (sub["umis_mpra_bc"] == 0).mean()

        # NB predicted zeros
        mu_nb, theta_nb = get_nb_params(nb_params, ct)
        p0_nb = nb_p_zero(mu_nb, theta_nb)
        nb_pct_zero = 100 * np.mean(p0_nb)

        # ZINB predicted zeros: pi + (1-pi) * P(X=0|NB_zinb)
        mu_zinb, theta_zinb = get_nb_params(zinb_params, ct)
        pi_vals = get_zi_probs(zinb_params, ct)
        p0_nb_zinb = nb_p_zero(mu_zinb, theta_zinb)
        # pi_vals may be scalar (1 per rep) or per-CRE
        # Average pi across reps if needed
        pi_mean = np.mean(pi_vals)
        p0_zinb = pi_mean + (1 - pi_mean) * p0_nb_zinb
        zinb_pct_zero = 100 * np.mean(p0_zinb)

        excess = obs_pct_zero - nb_pct_zero
        zinb_shift = zinb_pct_zero - nb_pct_zero
        mu_shift = np.mean(mu_zinb) - np.mean(mu_nb)

        row = {
            "dataset": ds_label,
            "cell_type": ct,
            "obs_%zero": obs_pct_zero,
            "NB_pred_%zero": nb_pct_zero,
            "ZINB_pred_%zero": zinb_pct_zero,
            "excess_pp": excess,
            "zinb_shift_pp": zinb_shift,
            "mean_pi": pi_mean,
            "mean_mu_NB": np.mean(mu_nb),
            "mean_mu_ZINB": np.mean(mu_zinb),
            "mu_shift": mu_shift,
            "theta_NB": theta_nb,
            "theta_ZINB": theta_zinb,
        }
        all_rows.append(row)
        print(f"  {ct:25s}  obs={obs_pct_zero:5.1f}%  NB={nb_pct_zero:5.1f}%  "
              f"ZINB={zinb_pct_zero:5.1f}%  excess={excess:+5.1f}pp  "
              f"zinb_shift={zinb_shift:+5.1f}pp  pi={pi_mean:.3f}  "
              f"mu_NB={np.mean(mu_nb):.4f}  mu_ZINB={np.mean(mu_zinb):.4f}",
              flush=True)

# ---------------------------------------------------------------------------
# By-CRE analysis: aggregate pi and predicted zeros across CREs
# ---------------------------------------------------------------------------
print("\n\n=== BY-CRE SUMMARY (where LRT matters) ===", flush=True)
print(f"{'dataset':20s} {'n_cre':>6s} {'med_pi':>8s} {'mean_pi':>8s} {'max_pi':>8s}  "
      f"{'med_mu_NB':>10s} {'med_mu_ZINB':>12s}", flush=True)
print("-" * 85)

cre_rows = []
for ds_label, paths in datasets.items():
    nb_ortho = scm.ortho.load(client, str(paths["nb"].parent), paths["nb"].name)
    zinb_ortho = scm.ortho.load(client, str(paths["zinb"].parent), paths["zinb"].name)

    nb_cre_params = nb_ortho.by_cre_parameters
    zinb_cre_params = zinb_ortho.by_cre_parameters

    all_pi = []
    all_mu_nb = []
    all_mu_zinb = []
    for key in zinb_cre_params.keys:
        # ZI probs
        zi = zinb_cre_params.zi[key]
        zi_df = zi.result() if hasattr(zi, "result") else zi
        if zi_df is not None:
            all_pi.extend(zi_df.values.flatten().tolist())

        # NB mu
        nb_nb = nb_cre_params.nb[key]
        nb_nb_df = nb_nb.result() if hasattr(nb_nb, "result") else nb_nb
        all_mu_nb.extend(nb_nb_df["mu"].values.tolist())

        # ZINB mu
        zinb_nb = zinb_cre_params.nb[key]
        zinb_nb_df = zinb_nb.result() if hasattr(zinb_nb, "result") else zinb_nb
        all_mu_zinb.extend(zinb_nb_df["mu"].values.tolist())

    pi_arr = np.array(all_pi)
    mu_nb_arr = np.array(all_mu_nb)
    mu_zinb_arr = np.array(all_mu_zinb)

    print(f"{ds_label:20s} {len(zinb_cre_params.keys):6d} {np.median(pi_arr):8.4f} "
          f"{np.mean(pi_arr):8.4f} {np.max(pi_arr):8.4f}  "
          f"{np.median(mu_nb_arr):10.4f} {np.median(mu_zinb_arr):12.4f}", flush=True)

    cre_rows.append({
        "dataset": ds_label,
        "n_cre": len(zinb_cre_params.keys),
        "median_pi": np.median(pi_arr),
        "mean_pi": np.mean(pi_arr),
        "max_pi": np.max(pi_arr),
        "pct_pi_gt_0.05": 100 * np.mean(pi_arr > 0.05),
        "pct_pi_gt_0.10": 100 * np.mean(pi_arr > 0.10),
        "median_mu_NB": np.median(mu_nb_arr),
        "median_mu_ZINB": np.median(mu_zinb_arr),
        "mean_mu_NB": np.mean(mu_nb_arr),
        "mean_mu_ZINB": np.mean(mu_zinb_arr),
    })

cre_df = pd.DataFrame(cre_rows)
cre_df.to_csv(OUT_DIR / "zi_zero_decomposition_by_cre.tsv", sep="\t", index=False)
print(f"\nSaved by-CRE table to {OUT_DIR / 'zi_zero_decomposition_by_cre.tsv'}")

client.close()
cluster.close()

result_df = pd.DataFrame(all_rows)
result_df.to_csv(OUT_DIR / "zi_zero_decomposition.tsv", sep="\t", index=False)
print(f"\nSaved table to {OUT_DIR / 'zi_zero_decomposition.tsv'}")

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
print("\nGenerating figure...", flush=True)

ds_names = list(datasets.keys())
n_ds = len(ds_names)

fig, axes = plt.subplots(1, n_ds, figsize=(5 * n_ds, 5), sharey=False)
if n_ds == 1:
    axes = [axes]

for ax, ds_label in zip(axes, ds_names):
    sub = result_df[result_df["dataset"] == ds_label].copy()
    sub = sub.sort_values("cell_type")
    ct_labels = sub["cell_type"].values
    x = np.arange(len(ct_labels))
    w = 0.25

    bars_obs  = ax.bar(x - w, sub["obs_%zero"],       w, label="Observed", color="gray", alpha=0.8)
    bars_nb   = ax.bar(x,     sub["NB_pred_%zero"],    w, label="NB predicted", color="steelblue", alpha=0.8)
    bars_zinb = ax.bar(x + w, sub["ZINB_pred_%zero"],  w, label="ZINB predicted", color="coral", alpha=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels(ct_labels, rotation=30, ha="right", fontsize=7)
    ax.set_ylabel("% zero observations")
    ax.set_title(ds_label)
    ax.legend(fontsize=7, loc="lower left")

    # Annotate excess and ZINB shift
    for i, (_, row) in enumerate(sub.iterrows()):
        # Small annotation above bars
        ymax = max(row["obs_%zero"], row["NB_pred_%zero"], row["ZINB_pred_%zero"])
        ax.text(i, ymax + 0.5, f'pi={row["mean_pi"]:.2f}', ha="center",
                fontsize=6, color="coral")

fig.suptitle("Zero-inflation decomposition: observed vs model-predicted zero rates\n"
             "ZINB can predict FEWER zeros than NB by allowing higher mu",
             fontsize=11)
plt.tight_layout()
fig.savefig(OUT_DIR / "zi_zero_decomposition_bars.svg", format="svg", bbox_inches="tight")
print(f"Saved bar chart to {OUT_DIR / 'zi_zero_decomposition_bars.svg'}")

# --- Panel 2: mu shift scatter ---
fig2, ax2 = plt.subplots(figsize=(6, 5))

for ds_label in ds_names:
    sub = result_df[result_df["dataset"] == ds_label]
    ax2.scatter(sub["mean_pi"], sub["zinb_shift_pp"],
                label=ds_label, s=80, alpha=0.8, zorder=3)

ax2.axhline(0, color="black", lw=0.8, ls="--", zorder=0)
ax2.set_xlabel("Mean ZI probability (pi)")
ax2.set_ylabel("ZINB shift in predicted %zero (pp)\n(ZINB - NB; negative = fewer zeros)")
ax2.set_title("Zero-inflation paradox: higher pi can mean fewer predicted zeros\n"
              "ZI absorbs structural zeros, freeing NB to fit higher mu")
ax2.legend()

# Annotate cell types
for _, row in result_df.iterrows():
    ct_short = row["cell_type"][:12]
    ax2.annotate(ct_short, (row["mean_pi"], row["zinb_shift_pp"]),
                 fontsize=5, textcoords="offset points", xytext=(5, 3))

plt.tight_layout()
fig2.savefig(OUT_DIR / "zi_zero_decomposition_scatter.svg", format="svg", bbox_inches="tight")
print(f"Saved scatter to {OUT_DIR / 'zi_zero_decomposition_scatter.svg'}")

plt.close("all")
print("\nDone.")
