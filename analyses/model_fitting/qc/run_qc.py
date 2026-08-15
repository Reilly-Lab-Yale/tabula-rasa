#!/usr/bin/env python
"""
Generic ortho QC script.

Usage:
    ipython run_qc.py -- --dataset shendure --ortho shendure_obs_nb
    ipython run_qc.py -- --dataset cohen --ortho cohen_cm_zinb_phantom_20260401
    ipython run_qc.py -- --dataset seelig --ortho seelig_cm_nb_phantom

Detects NB-only models automatically (skips ZI plots).
Generates:  qc/{ortho_name}/plots/*.svg  and  qc/{ortho_name}/summary.txt
"""
import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from dask.distributed import Client, LocalCluster

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent))
import scMPRAforge.core as scm

# -- config per dataset -------------------------------------------------------
DATASET_CONFIG = {
    "shendure": {
        "data_root": "/nfs/roberts/project/pi_skr2/shared/tabula_data_new/shendure",
        "cell_type_labels": {
            "reference": "Pluripotent (reference)",
            "Cardiomyocytes": "Cardiomyocytes",
            "EpiblastPrimitiveStreak": "Epiblast/PrimStreak",
            "ExEndodermParietal": "ExEndoderm Parietal",
            "ExEndodermVisceral": "ExEndoderm Visceral",
            "Haematoendothelial": "Haematoendothelial",
            "Mesoderm": "Mesoderm",
            "NeuroectodermBrain": "Neuroectoderm Brain",
            "NeuroectodermRostral": "Neuroectoderm Rostral",
            "SurfaceEctoderm": "Surface Ectoderm",
        },
    },
    "cohen": {
        "data_root": "/nfs/roberts/project/pi_skr2/shared/tabula_data_new/cohen",
        "cell_type_labels": {
            "reference": "Retina (reference)",
            "Bipolar": "Bipolar",
            "Mueller Glia": "Mueller Glia",
            "Interneuron": "Interneuron",
        },
    },
    "seelig": {
        "data_root": "/nfs/roberts/project/pi_skr2/shared/tabula_data_new/seelig",
        "cell_type_labels": {"reference": "HepG2 (reference)", "K562": "K562"},
    },
    "takeshi": {
        "data_root": "/nfs/roberts/project/pi_skr2/shared/tabula_data_new/takeshi",
        "cell_type_labels": {"reference": "HepG2 (reference)", "K562": "K562", "SKNSH": "SK-N-SH"},
    },
}

# -- parse args ----------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--dataset", required=True, choices=DATASET_CONFIG.keys())
parser.add_argument("--ortho", required=True, help="Ortho directory name")
args = parser.parse_args()

cfg = DATASET_CONFIG[args.dataset]
DATA_ROOT = Path(cfg["data_root"])
ORTHO_NAME = args.ortho
CELL_TYPE_LABELS = cfg["cell_type_labels"]

WORK_DIR = Path(__file__).resolve().parent / ORTHO_NAME
PLOT_DIR = WORK_DIR / "plots"
PLOT_DIR.mkdir(parents=True, exist_ok=True)
SUMMARY_PATH = WORK_DIR / "summary.txt"

# -- cluster -------------------------------------------------------------------
cluster = LocalCluster(n_workers=4, threads_per_worker=1, memory_limit="50GB", processes=False)
client = Client(cluster)
print(client.dashboard_link, flush=True)

# -- load ----------------------------------------------------------------------
print(f"[+] Loading ortho: {DATA_ROOT / ORTHO_NAME}", flush=True)
primordial = scm.ortho.load(client, DATA_ROOT, ORTHO_NAME)

# -- detect NB-only ------------------------------------------------------------
first_cre_key = next(iter(primordial.by_cre.model))
first_model = primordial.by_cre.model[first_cre_key]
if hasattr(first_model, 'result'):
    first_model = first_model.result()
is_nb_only = first_model.get("nb_only", False)
print(f"[+] Model type: {'NB-only' if is_nb_only else 'ZINB'}", flush=True)

# -- detect fit_mode for axis labels -------------------------------------------
fit_mode = None
for design_dict in (primordial.by_cell_type_design, primordial.by_cre_design):
    if design_dict is not None:
        sample_key = next(iter(design_dict))
        sample = design_dict[sample_key]
        if hasattr(sample, 'result'):
            sample = sample.result()
        fit_mode = sample.get('fit_mode', None)
        if fit_mode is not None:
            break
is_moib = fit_mode == "cm_phantom_moib"
mean_label = "mean UMI (conditional)" if is_moib else "mean UMI (observed)"

# -- compute QC ----------------------------------------------------------------
print("[+] Computing model QC...", flush=True)
primordial.compute_model_qc()
by_cell = primordial.by_cell_qc
by_cre  = primordial.by_cre_qc

def label(ct):
    return CELL_TYPE_LABELS.get(ct, ct)

def save(fig, name):
    path = PLOT_DIR / name
    fig.savefig(path, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {path}", flush=True)

LOW_R_THRESH = 0.8

# -- summary.txt ---------------------------------------------------------------
lines = [
    "=" * 72,
    f"  {ORTHO_NAME} -- Model QC Summary",
    f"  Model type: {'NB-only' if is_nb_only else 'ZINB'}",
    "=" * 72, "",
]

lines += ["-- Convergence ---------------------------------------------------------", ""]
for lbl, qc in [("by_cell_type", by_cell), ("by_cre", by_cre)]:
    total = len(qc)
    success = sum(1 for v in qc.values() if v["success"])
    failed = [k for k, v in qc.items() if not v["success"]]
    lines.append(f"  {lbl}: {success}/{total} converged")
    if failed:
        lines.append(f"    Failed levels:")
        for f in failed:
            short = f[:60] + "..." if len(f) > 60 else f
            lines.append(f"      {short}")
lines.append("")

lines += ["-- Mu range (NB mean estimate) -----------------------------------------", ""]
for lbl, qc in [("by_cell_type", by_cell), ("by_cre", by_cre)]:
    all_mu = pd.concat(
        [v["dat"]["mu"] for v in qc.values() if v["success"] and "mu" in v["dat"].columns],
        ignore_index=True,
    )
    lines.append(f"  {lbl}: min={all_mu.min():.4f}  max={all_mu.max():.4f}  "
                 f"median={all_mu.median():.4f}")
    if all_mu.min() <= 0:
        lines.append(f"    WARNING: mu <= 0 detected")
    if all_mu.max() > 1000:
        lines.append(f"    WARNING: mu > 1000 detected")
lines.append("")

lines += ["-- Pearson r (mu vs mean UMI) ------------------------------------------", ""]
for lbl, qc in [("by_cell_type", by_cell), ("by_cre", by_cre)]:
    r_vals = {k: v["r_value"] for k, v in qc.items() if v["success"] and v["r_value"] is not None}
    arr = np.array(list(r_vals.values()))
    n_nan = int(np.isnan(arr).sum())
    arr_clean = arr[~np.isnan(arr)]
    low = {k: r for k, r in r_vals.items() if not np.isnan(r) and r < LOW_R_THRESH}
    # Summarize NaN anti-levels dropped before regression
    total_dropped = sum(v.get("n_dropped_nan", 0) for v in qc.values() if v["success"])
    total_anti = sum(len(v["dat"]) for v in qc.values() if v["success"])
    if len(arr_clean) > 0:
        lines.append(f"  {lbl}: mean r={np.nanmean(arr):.3f}  min={np.nanmin(arr):.3f}  "
                     f"max={np.nanmax(arr):.3f}  n(r<{LOW_R_THRESH})={len(low)}  "
                     f"n_nan_r={n_nan}")
    else:
        lines.append(f"  {lbl}: no valid r-values  n_nan_r={n_nan}")
    lines.append(f"    NaN anti-levels dropped before regression: "
                 f"{total_dropped}/{total_anti} across all models")
    if low:
        lines.append(f"    Low-r levels (r < {LOW_R_THRESH}):")
        for k, r in sorted(low.items(), key=lambda x: x[1]):
            short = k[:60] + "..." if len(k) > 60 else k
            lines.append(f"      r={r:.3f}  {short}")
lines.append("")

SUMMARY_PATH.write_text("\n".join(lines) + "\n")
print(f"[+] Summary written: {SUMMARY_PATH}", flush=True)

# -- plot 1: theta distributions -----------------------------------------------
print("[+] Plotting theta distributions...", flush=True)

ct_params  = primordial.by_cell_type_parameters
cre_params = primordial.by_cre_parameters

ct_thetas = {
    label(k): np.array(ct_params.theta[k].result()).flatten()
    for k in ct_params.theta
}

cre_thetas_all = np.concatenate([
    np.array(cre_params.theta[k].result()).flatten()
    for k in cre_params.theta
])

fig, axes = plt.subplots(1, 2, figsize=(10, 5))
fig.suptitle(f"{ORTHO_NAME} -- Theta (dispersion) distributions", fontsize=13)

axes[0].set_title("by_cell_type models")
axes[0].violinplot(list(ct_thetas.values()), showmedians=True)
axes[0].set_xticks(range(1, len(ct_thetas) + 1))
axes[0].set_xticklabels(list(ct_thetas.keys()), rotation=20, ha="right", fontsize=7)
axes[0].set_ylabel("theta")
axes[0].set_yscale("log")

n_cre = len(cre_params.theta)
axes[1].set_title(f"by_cre models (all {n_cre} CREs)")
axes[1].violinplot([cre_thetas_all], showmedians=True)
axes[1].set_xticks([1])
axes[1].set_xticklabels(["all CREs"])
axes[1].set_ylabel("theta")
axes[1].set_yscale("log")

fig.tight_layout()
save(fig, "theta_distributions.svg")

# -- plot 2: r-value distributions ---------------------------------------------
print("[+] Plotting r-value distributions...", flush=True)

fig, axes = plt.subplots(1, 2, figsize=(10, 5))
fig.suptitle(f"{ORTHO_NAME} -- Pearson r (mu vs mean UMI)", fontsize=13)

for ax, (lbl, qc) in zip(axes, [("by_cell_type", by_cell), ("by_cre", by_cre)]):
    r_vals = [v["r_value"] for v in qc.values()
              if v["success"] and v["r_value"] is not None and not np.isnan(v["r_value"])]
    n_nan = sum(1 for v in qc.values()
                if v["success"] and v["r_value"] is not None and np.isnan(v["r_value"]))
    if r_vals:
        ax.violinplot([r_vals], showmedians=True)
        ax.scatter(
            np.ones(len(r_vals)) + np.random.uniform(-0.05, 0.05, len(r_vals)),
            r_vals, alpha=0.3, s=6, color="steelblue"
        )
    ax.axhline(LOW_R_THRESH, color="red", linestyle="--", alpha=0.7, label=f"r={LOW_R_THRESH}")
    n_low = sum(1 for r in r_vals if r < LOW_R_THRESH)
    ax.set_title(f"{lbl}\n(n_low={n_low}, n_nan={n_nan})")
    ax.set_ylabel("r")
    ax.set_ylim(-0.1, 1.05)
    ax.set_xticks([])
    ax.legend(fontsize=8)

fig.tight_layout()
save(fig, "r_values.svg")

# -- plot 3: mu vs mean UMI scatter -- by_cell_type ----------------------------
print("[+] Plotting mu vs mean UMI (by_cell_type)...", flush=True)

n_ct = len(by_cell)
ncols = min(n_ct, 6)
nrows = int(np.ceil(n_ct / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows), squeeze=False)
fig.suptitle(f"{ORTHO_NAME} -- mu vs mean(UMI) by cell type", fontsize=13)

for i, (ct, qc_entry) in enumerate(by_cell.items()):
    ax = axes[i // ncols][i % ncols]
    dat = qc_entry["dat"].dropna(subset=["mu", "mean(umis_mpra_bc)"])
    n_drop = len(qc_entry["dat"]) - len(dat)
    ax.scatter(dat["mean(umis_mpra_bc)"], dat["mu"], alpha=0.4, s=8, color="steelblue")
    r = qc_entry["r_value"]
    slope = qc_entry["slope"]
    intercept = qc_entry["intercept"]
    if r is not None and not np.isnan(r):
        x = np.linspace(dat["mean(umis_mpra_bc)"].min(), dat["mean(umis_mpra_bc)"].max(), 100)
        ax.plot(x, slope * x + intercept, color="tomato", linewidth=1.5, label=f"r={r:.3f}")
        ax.legend(fontsize=9)
    subtitle = label(ct)
    if n_drop > 0:
        subtitle += f" ({n_drop} NaN dropped)"
    ax.set_xlabel(mean_label)
    ax.set_ylabel("mu (NB estimate)")
    ax.set_title(subtitle)

for i in range(len(by_cell), nrows * ncols):
    axes[i // ncols][i % ncols].set_visible(False)

fig.tight_layout()
save(fig, "mu_vs_mean_by_cell_type.svg")

# -- plot 4: mu vs mean UMI -- by_cre (worst 20 + random 20) ------------------
print("[+] Plotting mu vs mean UMI (by_cre, worst + sample)...", flush=True)

r_sorted = sorted(
    [(k, v) for k, v in by_cre.items()
     if v["success"] and v["r_value"] is not None and not np.isnan(v["r_value"])],
    key=lambda x: x[1]["r_value"]
)
n_nan_cre = sum(1 for v in by_cre.values()
                if v["success"] and v["r_value"] is not None and np.isnan(v["r_value"]))
worst_20 = r_sorted[:20]
rng = np.random.default_rng(42)
rest = r_sorted[20:]
sample_20 = [rest[i] for i in rng.choice(len(rest), min(20, len(rest)), replace=False)]
to_plot = worst_20 + sample_20

ncols_p = 8
nrows_p = max(1, int(np.ceil(len(to_plot) / ncols_p)))
fig, axes = plt.subplots(nrows_p, ncols_p, figsize=(ncols_p * 2.5, nrows_p * 2.5))
nan_note = f" [{n_nan_cre} CREs had NaN r, excluded]" if n_nan_cre else ""
fig.suptitle(f"{ORTHO_NAME} -- mu vs mean(UMI) by CRE\n"
             f"(worst 20 r-values + 20 random; red=worst){nan_note}", fontsize=11)
axes_flat = np.atleast_1d(axes).flatten()

for i, (cre_id, qc_entry) in enumerate(to_plot):
    ax = axes_flat[i]
    dat = qc_entry["dat"].dropna(subset=["mu", "mean(umis_mpra_bc)"])
    color = "tomato" if i < 20 else "steelblue"
    ax.scatter(dat["mean(umis_mpra_bc)"], dat["mu"], alpha=0.7, s=20, color=color)
    r = qc_entry["r_value"]
    ax.set_title(f"r={r:.2f}", fontsize=7)
    ax.tick_params(labelsize=6)
    ax.set_xlabel(mean_label, fontsize=6)
    ax.set_ylabel("mu", fontsize=6)

for ax in axes_flat[len(to_plot):]:
    ax.set_visible(False)

fig.tight_layout()
save(fig, "mu_vs_mean_by_cre_worst.svg")

# -- plot 5: ZI parameter by cell type (ZINB only) ----------------------------
if not is_nb_only:
    print("[+] Plotting ZI parameters...", flush=True)

    zi_data = {}
    for ct in ct_params.zi:
        df = ct_params.zi[ct].result()
        if df is None:
            continue
        zi_col = ct if ct in df.columns else df.columns[0]
        zi_data[label(ct)] = df[zi_col].values

    if zi_data:
        fig, ax = plt.subplots(figsize=(6, 5))
        fig.suptitle(f"{ORTHO_NAME} -- Zero-inflation (ZI) by cell type", fontsize=13)

        labels_zi = list(zi_data.keys())
        data_zi = list(zi_data.values())
        ax.violinplot(data_zi, showmedians=True)
        for i, (lbl_zi, vals) in enumerate(zip(labels_zi, data_zi)):
            ax.scatter(
                np.ones(len(vals)) * (i + 1) + np.random.uniform(-0.05, 0.05, len(vals)),
                vals, alpha=0.3, s=6, color="steelblue"
            )
        ax.set_xticks(range(1, len(labels_zi) + 1))
        ax.set_xticklabels(labels_zi, rotation=15, ha="right")
        ax.set_ylabel("ZI probability")
        ax.set_title("Per-CRE zero-inflation estimates")

        fig.tight_layout()
        save(fig, "zi_by_cell_type.svg")
    else:
        print("[+] No ZI data available.", flush=True)
else:
    print("[+] Skipping ZI plot (NB-only model).", flush=True)

# -- plot 6: r-value vs missing data (consider_missing datasets) ---------------
has_missing = False
for cre_id, qc_entry in by_cre.items():
    if qc_entry["success"] and qc_entry["dat"]["mean(umis_mpra_bc)"].isna().any():
        has_missing = True
        break

if has_missing:
    print("[+] Plotting r-value vs missing data...", flush=True)

    r_vals_cre = []
    n_missing = []
    for cre_id, qc_entry in by_cre.items():
        if not qc_entry["success"] or qc_entry["r_value"] is None:
            continue
        dat = qc_entry["dat"]
        n_nan = dat["mean(umis_mpra_bc)"].isna().sum()
        r_vals_cre.append(qc_entry["r_value"])
        n_missing.append(n_nan)

    fig, ax = plt.subplots(figsize=(6, 5))
    fig.suptitle(f"{ORTHO_NAME} -- r-value vs missing cell types (by_cre)", fontsize=12)
    ax.scatter(n_missing, r_vals_cre, alpha=0.3, s=8, color="steelblue")
    ax.axhline(LOW_R_THRESH, color="red", linestyle="--", alpha=0.7, label=f"r={LOW_R_THRESH}")
    ax.set_xlabel("n cell types with missing observations")
    ax.set_ylabel("Pearson r (mu vs mean UMI)")
    ax.legend(fontsize=9)
    fig.tight_layout()
    save(fig, "r_vs_missing.svg")
else:
    print("[+] Skipping r-vs-missing plot (no missing data detected).", flush=True)

# -- done ----------------------------------------------------------------------
print(f"[+] All plots saved to: {PLOT_DIR}", flush=True)
print(f"[+] Summary: {SUMMARY_PATH}", flush=True)
client.close()
cluster.close()
