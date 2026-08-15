#!/usr/bin/env python3
"""
Cross-dataset power regression: log10(FC@50%-power) ~ library/cell factors.

Naive first-pass: across all (dataset, cell_type) combinations from
shendure/cohen/seelig power-curve simulations, regress the per-CT
sensitivity (FC at which Welch's t-test power crosses 0.5) on four
factors:
    log10(n_cells)        per CT
    log10(n_cres)         per CT (matches what each sim actually used)
    log10(bcs_per_cre)    dataset-level (total_uniq_mpra_bc / n_cres_total)
    log10(moi_median)     per CT (median per-cell unique mpra_bc count;
                                  see compute_per_ct_moi.py)

No dataset intercept -- the predictors are meant to absorb dataset
differences. With 16 rows and 4 predictors plus intercept, residual df=11.
Tight but adequate for descriptive purposes.

Outputs (under output/):
    regression_data.parquet     -- the 16-row design matrix + response
    regression_summary.txt      -- statsmodels OLS summary + partial R^2
    regression_obs_vs_pred.svg  -- observed vs predicted, colored by dataset
    regression_partial.svg      -- 4-panel partial regression plots
"""

import sys
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_repo_root = Path(__file__).resolve().parents[3]
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

import scMPRAforge as scm

OUT_DIR = Path(__file__).resolve().parent / "output"
OUT_DIR.mkdir(exist_ok=True)
DATA = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new")
ANALYSES = Path("/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa-reporter-fig/analyses")

# ---------------------------------------------------------------------------
# Per-dataset config
#
# - power_df:  per-CT (reject_null, fc) for the canonical test condition that
#              each dataset's existing power sim uses (reporter for
#              shendure/cohen, deflated for seelig since seelig has no
#              transfection reporter).
# - ortho_dir: ortho whose by_cell_type_parameters matches what the power sim
#              actually drew n_cres from. Shendure uses obs (per-CT varies),
#              cohen hardcodes N_CRES=115 (essentially constant), seelig uses
#              the cm phantom (constant).
# - n_cres:    if not None, override the per-CT value with this constant.
#              Matches the simulation's actual N_CRES.
# - bounds:    the canonical Bounds object whose total_uniq_mpra_bc and
#              cells_per_cell_type drive the simulation.
# ---------------------------------------------------------------------------
DATASETS = {
    "shendure": dict(
        power_df=ANALYSES / "simulation/activity_power/shendure/output/power_df_reporter.parquet",
        ortho_dir=DATA / "shendure/shendure_obs_nb",
        n_cres=None,  # per-CT, varies
        bounds=scm.SHENDURE_BOUNDS,
    ),
    "cohen": dict(
        power_df=ANALYSES / "simulation/activity_power/cohen/output/cohen_power_df_reporter.parquet",
        ortho_dir=DATA / "cohen/cohen_obsingle_nb_phantom",
        n_cres=115,  # cohen power sim hardcodes this
        bounds=scm.COHEN_BOUNDS,
    ),
    "seelig": dict(
        power_df=ANALYSES / "simulation/activity_power/seelig/output/seelig_power_df_deflated.parquet",
        ortho_dir=DATA / "seelig/seelig_cm_nb_phantom",
        n_cres=None,
        bounds=scm.SEELIG_BOUNDS,
    ),
}

PER_CT_MOI = pd.read_parquet(OUT_DIR / "per_ct_moi.parquet")


# ---------------------------------------------------------------------------
# Response: log10(FC@50%) per (dataset, cell_type)
# ---------------------------------------------------------------------------

def fc_at_power(power_df_path, target_power=0.5, n_bins=50):
    """For each cell type, return the smallest binned FC > 1 at which the
    rejection rate reaches `target_power`. NaN if never reached."""
    df = pd.read_parquet(power_df_path)
    df = df[df["fc"] >= 1.0].copy()
    out = {}
    for ct, sub in df.groupby("cell_type"):
        sub = sub.sort_values("fc").copy()
        sub["bin"] = pd.cut(sub["fc"], bins=n_bins)
        binned = (
            sub.groupby("bin", observed=True)
            .agg(rej=("reject_null", "mean"), mid=("fc", "mean"))
            .reset_index()
            .dropna()
        )
        cross = binned[binned["rej"] >= target_power]
        out[ct] = float(cross["mid"].iloc[0]) if len(cross) > 0 else np.nan
    return out


# ---------------------------------------------------------------------------
# Predictors per (dataset, cell_type)
# ---------------------------------------------------------------------------

def per_ct_n_cres(ortho_dir, override=None):
    if override is not None:
        return ("__constant__", int(override))
    fp = ortho_dir / "by_cell_type_parameters.pkl"
    with open(fp, "rb") as f:
        params = pickle.load(f)
    return {ct: len(params.nb[ct]) for ct in params.nb}


# ---------------------------------------------------------------------------
# Build the design matrix
# ---------------------------------------------------------------------------

rows = []
for ds, cfg in DATASETS.items():
    bounds = cfg["bounds"]
    cells = bounds.cells_per_cell_type  # pd.Series indexed by cell_type
    total_bcs = bounds.total_uniq_mpra_bc

    n_cres = per_ct_n_cres(cfg["ortho_dir"], cfg["n_cres"])
    moi_ds = PER_CT_MOI[PER_CT_MOI["dataset"] == ds].set_index("cell_type")

    # bcs_per_cre uses the simulation's effective N_CRES (constant per ds for
    # cohen/seelig, the obs total CRE count for shendure)
    if cfg["n_cres"] is not None:
        n_cres_for_lib = cfg["n_cres"]
    else:
        # use the obs ortho's max-CT n_cres as the union (library size)
        n_cres_for_lib = max(n_cres.values())
    bcs_per_cre = total_bcs / n_cres_for_lib

    # Pluripotent / Rod / HepG2 are stored as "reference" in the bounds
    # cells series, matching how the orthos relabel it.
    for ct in cells.index:
        if ct not in moi_ds.index:
            print(f"  WARN: no MOI for {ds}/{ct}, skipping", file=sys.stderr)
            continue
        n = int(cells.loc[ct])
        nc = n_cres if isinstance(n_cres, dict) else n_cres[1]
        if isinstance(nc, dict):
            ct_n_cres = nc.get(ct)
        else:
            ct_n_cres = nc
        rows.append(dict(
            dataset=ds,
            cell_type=ct,
            n_cells=n,
            n_cres=int(ct_n_cres),
            bcs_per_cre=float(bcs_per_cre),
            moi_median=float(moi_ds.loc[ct, "moi_median"]),
            moi_mean=float(moi_ds.loc[ct, "moi_mean"]),
        ))

design = pd.DataFrame(rows)

# Response
response_per_ds = {ds: fc_at_power(cfg["power_df"]) for ds, cfg in DATASETS.items()}

design["fc50"] = [
    response_per_ds[r["dataset"]].get(r["cell_type"], np.nan)
    for _, r in design.iterrows()
]

# log10 transforms
for col in ["n_cells", "n_cres", "bcs_per_cre", "moi_median", "moi_mean", "fc50"]:
    design[f"log10_{col}"] = np.log10(design[col])

# Drop any rows with NaN response
n_dropped = design["fc50"].isna().sum()
if n_dropped:
    print(f"Dropping {n_dropped} rows with no FC@50% (never crossed)", file=sys.stderr)
design = design.dropna(subset=["fc50"]).reset_index(drop=True)

print("\n=== Design matrix ===")
print(design[[
    "dataset", "cell_type", "n_cells", "n_cres", "bcs_per_cre",
    "moi_median", "fc50",
]].to_string(index=False))

design.to_parquet(OUT_DIR / "regression_data.parquet")


# ---------------------------------------------------------------------------
# OLS fit
# ---------------------------------------------------------------------------

predictors = ["log10_n_cells", "log10_n_cres", "log10_bcs_per_cre", "log10_moi_median"]
X = sm.add_constant(design[predictors])
y = design["log10_fc50"]

model = sm.OLS(y, X).fit()

# Partial R^2 for each predictor: drop one at a time, refit, compare R^2
def partial_r2(model_full, X_full, y, drop_col):
    X_red = X_full.drop(columns=[drop_col])
    model_red = sm.OLS(y, X_red).fit()
    return (model_full.rsquared - model_red.rsquared) / (1.0 - model_red.rsquared)

partial = {p: partial_r2(model, X, y, p) for p in predictors}

summary_path = OUT_DIR / "regression_summary.txt"
with open(summary_path, "w") as f:
    f.write("Cross-dataset power regression\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"Response: log10(FC at which Welch t-test power crosses 0.5)\n")
    f.write(f"n = {len(design)} (dataset x cell_type)\n\n")
    f.write(str(model.summary()))
    f.write("\n\n")
    f.write("Partial R^2 (variance explained uniquely by each predictor):\n")
    for p, r2 in partial.items():
        f.write(f"  {p:25s} {r2:+.4f}\n")
    f.write(f"\nTotal R^2:    {model.rsquared:.4f}\n")
    f.write(f"Adj. R^2:     {model.rsquared_adj:.4f}\n")

print(f"\nSaved: {summary_path}")
print("\n--- OLS summary ---")
print(model.summary())
print("\nPartial R^2:")
for p, r2 in partial.items():
    print(f"  {p:25s} {r2:+.4f}")


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

ds_colors = {"shendure": "#377eb8", "cohen": "#e41a1c", "seelig": "#4daf4a"}

# Obs vs pred
fig, ax = plt.subplots(figsize=(5, 5))
y_pred = model.fittedvalues
for ds in design["dataset"].unique():
    mask = design["dataset"] == ds
    ax.scatter(y[mask], y_pred[mask], color=ds_colors[ds], label=ds, s=40, alpha=0.85)
lims = [min(y.min(), y_pred.min()) - 0.02, max(y.max(), y_pred.max()) + 0.02]
ax.plot(lims, lims, "k--", lw=0.8)
ax.set_xlim(lims); ax.set_ylim(lims)
ax.set_xlabel("observed log10(FC@50%)")
ax.set_ylabel("predicted log10(FC@50%)")
ax.set_title(
    f"Cross-dataset power regression\n"
    f"R^2 = {model.rsquared:.3f}, n = {len(design)}",
    fontsize=10,
)
ax.legend(fontsize=8, loc="upper left")
plt.tight_layout()
fig.savefig(OUT_DIR / "regression_obs_vs_pred.svg", format="svg", bbox_inches="tight")
print(f"Saved: {OUT_DIR / 'regression_obs_vs_pred.svg'}")

# Partial regression: for each predictor, plot residuals of y vs predictor
# after partialing out the other predictors
fig, axes = plt.subplots(2, 2, figsize=(9, 7))
for i, p in enumerate(predictors):
    ax = axes.flat[i]
    others = [q for q in predictors if q != p]
    X_oth = sm.add_constant(design[others])
    res_y = sm.OLS(y, X_oth).fit().resid
    res_x = sm.OLS(design[p], X_oth).fit().resid
    for ds in design["dataset"].unique():
        mask = (design["dataset"] == ds).values
        ax.scatter(res_x[mask], res_y[mask], color=ds_colors[ds], s=30, alpha=0.85, label=ds)
    # fit line
    slope = model.params[p]
    xs = np.array([res_x.min(), res_x.max()])
    ax.plot(xs, slope * xs, "k-", lw=0.8)
    ax.axhline(0, color="grey", lw=0.5, linestyle=":")
    ax.axvline(0, color="grey", lw=0.5, linestyle=":")
    ax.set_xlabel(f"{p} | others")
    ax.set_ylabel("log10(FC@50%) | others")
    ax.set_title(
        f"{p}\nbeta={slope:+.3f}, partial R^2={partial[p]:+.3f}",
        fontsize=9,
    )
    if i == 0:
        ax.legend(fontsize=7, loc="best")

plt.tight_layout()
fig.savefig(OUT_DIR / "regression_partial.svg", format="svg", bbox_inches="tight")
print(f"Saved: {OUT_DIR / 'regression_partial.svg'}")
