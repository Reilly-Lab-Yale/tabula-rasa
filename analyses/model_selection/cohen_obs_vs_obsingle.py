"""
Cohen obs (CRE-coarse reporter expansion) vs obsingle (single-zero expansion)
comparison figure.

Demonstrates that CRE-coarse expansion floods zeros and crushes mu estimates,
while obsingle preserves the conditional mean and yields good QC.

Generates: cohen_obs_vs_obsingle.svg
"""
import sys
import re
import pickle
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
import scMPRAforge.core as scm

DATA_ROOT = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new/cohen")
QC_ROOT = Path(__file__).resolve().parent.parent / "model_fitting" / "qc"

# Load by_cell_type parameters for both conditions
def load_params(ortho_name):
    with open(DATA_ROOT / ortho_name / "by_cell_type_parameters.pkl", "rb") as f:
        return pickle.load(f)

obs_params = load_params("cohen_obs_nb_phantom_20260401")
obsingle_params = load_params("cohen_obsingle_nb_phantom")

# Extract mu values per (cell_type, cre_id) from each
def extract_mu_table(params):
    rows = []
    for level_key in params.keys:
        df = params.nb[level_key]  # DataFrame with 'mu' column, indexed by cre_id
        for cre_id, mu in df["mu"].items():
            rows.append({"level": level_key, "regressor": cre_id, "mu": float(mu)})
    return pd.DataFrame(rows)

obs_mu = extract_mu_table(obs_params)
obsingle_mu = extract_mu_table(obsingle_params)

# Cached so the manuscript panel can be drawn without the orthos: the coarse
# expansion is a fit nobody keeps a local copy of, and plotting should not
# require unpickling one. plot_cohen_expansion.py reads this.
MU_TSV = Path(__file__).resolve().parent / "cohen_expansion_mu.tsv"
_cached = pd.concat([obs_mu.assign(expansion="coarse"),
                     obsingle_mu.assign(expansion="fine")], ignore_index=True)
assert set(_cached.expansion) == {"coarse", "fine"}, "both expansions must be present"
assert (_cached.mu >= 0).all(), "fitted means cannot be negative"
_cached.to_csv(MU_TSV, sep="\t", index=False)
print(f"wrote {MU_TSV} ({len(_cached)} rows)")

# Parse r-values from QC summary files
def parse_r_values(summary_path, section="by_cell_type"):
    """Extract per-level r-values from a QC summary.txt file."""
    text = summary_path.read_text()
    # Find the section
    in_section = False
    r_values = {}
    for line in text.split("\n"):
        if f"{section}: mean r=" in line:
            in_section = True
            continue
        if in_section:
            m = re.match(r"\s+r=([\-0-9.nan]+)\s+(.*)", line)
            if m:
                r_str, name = m.group(1), m.group(2).strip()
                try:
                    r_values[name] = float(r_str)
                except ValueError:
                    r_values[name] = np.nan
            elif line.strip().startswith("by_") or line.strip().startswith("--"):
                break  # next section
    return r_values

obs_r = parse_r_values(
    QC_ROOT / "cohen_obs_nb_phantom_20260401" / "summary.txt", "by_cell_type")
obsingle_r = parse_r_values(
    QC_ROOT / "cohen_obsingle_nb_phantom" / "summary.txt", "by_cell_type")

# For by_cell_type, the "levels" are cell types and there are only 4.
# The summary only lists levels with r < 0.8. For obsingle, all 4 are < 0.8
# so they're all listed. For obs, 3/4 are listed (the 4th has r=NaN).
# Let's also parse the full r line for any levels NOT listed (r >= 0.8).
def parse_r_summary_line(summary_path, section="by_cell_type"):
    text = summary_path.read_text()
    for line in text.split("\n"):
        if f"{section}: mean r=" in line:
            return line.strip()
    return ""

obs_r_line = parse_r_summary_line(
    QC_ROOT / "cohen_obs_nb_phantom_20260401" / "summary.txt")
obsingle_r_line = parse_r_summary_line(
    QC_ROOT / "cohen_obsingle_nb_phantom" / "summary.txt")
print("Obs r-line:", obs_r_line)
print("Obsingle r-line:", obsingle_r_line)

# --- Figure: 2x2 layout ---
fig, axes = plt.subplots(2, 2, figsize=(11, 9), constrained_layout=True)

# Panel A: mu distribution comparison (by_cell_type)
ax = axes[0, 0]
obs_vals = obs_mu["mu"].values
obsingle_vals = obsingle_mu["mu"].values

bins = np.logspace(np.log10(max(1e-5, obs_vals[obs_vals > 0].min())),
                   np.log10(max(obs_vals.max(), obsingle_vals.max())),
                   50)
ax.hist(obs_vals, bins=bins, alpha=0.6, label="obs (CRE-coarse)", color="firebrick")
ax.hist(obsingle_vals, bins=bins, alpha=0.6, label="obsingle", color="steelblue")
ax.set_xscale("log")
ax.set_xlabel("NB mu (fitted mean)", fontsize=9)
ax.set_ylabel("Count", fontsize=9)
ax.set_title("A. Mu distribution (by cell type)", fontsize=10)
ax.legend(fontsize=8)
ax.axvline(np.median(obs_vals), color="firebrick", ls="--", lw=1, alpha=0.7)
ax.axvline(np.median(obsingle_vals), color="steelblue", ls="--", lw=1, alpha=0.7)
ax.text(np.median(obs_vals), ax.get_ylim()[1] * 0.9,
        f"  med={np.median(obs_vals):.4f}", fontsize=7, color="firebrick",
        va="top")
ax.text(np.median(obsingle_vals), ax.get_ylim()[1] * 0.8,
        f"  med={np.median(obsingle_vals):.3f}", fontsize=7, color="steelblue",
        va="top")

# Panel B: matched mu scatter (obs vs obsingle, same level+regressor)
ax = axes[0, 1]
merged = obs_mu.merge(obsingle_mu, on=["level", "regressor"],
                      suffixes=("_obs", "_obsingle"))
ax.scatter(merged["mu_obsingle"], merged["mu_obs"], s=8, alpha=0.5,
           color="gray", edgecolors="none")
lims = [1e-5, max(merged["mu_obs"].max(), merged["mu_obsingle"].max()) * 2]
ax.plot(lims, lims, "k--", lw=0.8, alpha=0.5)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Obsingle mu", fontsize=9)
ax.set_ylabel("Obs (CRE-coarse) mu", fontsize=9)
ax.set_title("B. Matched mu: obs vs obsingle", fontsize=10)
ax.set_xlim(lims)
ax.set_ylim(lims)
n_crushed = (merged["mu_obs"] < 0.01).sum()
n_total = len(merged)
ax.text(0.05, 0.95, f"{n_crushed}/{n_total} obs mu < 0.01",
        transform=ax.transAxes, fontsize=8, va="top", color="firebrick")

# Panel C: per-cell-type r-value comparison
ax = axes[1, 0]
# Build matched r-value table
all_cts = sorted(set(list(obs_r.keys()) + list(obsingle_r.keys())))
r_data = []
for ct in all_cts:
    r_data.append({
        "cell_type": ct,
        "r_obs": obs_r.get(ct, np.nan),
        "r_obsingle": obsingle_r.get(ct, np.nan),
    })
r_df = pd.DataFrame(r_data)

x = np.arange(len(r_df))
w = 0.35
bars_obs = ax.bar(x - w/2, r_df["r_obs"], w, label="obs (CRE-coarse)",
                  color="firebrick", alpha=0.7)
bars_osingle = ax.bar(x + w/2, r_df["r_obsingle"], w, label="obsingle",
                      color="steelblue", alpha=0.7)
ax.axhline(0, color="black", lw=0.5)
ax.axhline(0.8, color="gray", lw=0.8, ls=":", alpha=0.5)
ax.text(len(r_df) - 0.5, 0.82, "r = 0.8", fontsize=7, color="gray",
        ha="right")
ax.set_xticks(x)
ax.set_xticklabels(r_df["cell_type"], fontsize=8, rotation=30, ha="right")
ax.set_ylabel("Pearson r (mu vs mean UMI)", fontsize=9)
ax.set_title("C. Model-data correlation (by cell type)", fontsize=10)
ax.set_ylim(-0.15, 1.0)
ax.legend(fontsize=7, loc="upper left")
# Annotate bar values
for xi, row in enumerate(r_df.itertuples()):
    if not np.isnan(row.r_obs):
        ypos = max(row.r_obs, 0) + 0.02
        ax.text(xi - w/2, ypos, f"{row.r_obs:.2f}", ha="center", fontsize=6,
                color="firebrick")
    if not np.isnan(row.r_obsingle):
        ypos = max(row.r_obsingle, 0) + 0.02
        ax.text(xi + w/2, ypos, f"{row.r_obsingle:.2f}", ha="center",
                fontsize=6, color="steelblue")

# Panel D: per-cell-type median mu comparison
ax = axes[1, 1]
ct_summary = []
for ct in sorted(obs_mu["level"].unique()):
    obs_sub = obs_mu[obs_mu["level"] == ct]["mu"]
    osingle_sub = obsingle_mu[obsingle_mu["level"] == ct]["mu"]
    ct_summary.append({
        "cell_type": ct,
        "obs_median_mu": obs_sub.median(),
        "obsingle_median_mu": osingle_sub.median(),
    })
ct_df = pd.DataFrame(ct_summary)

x = np.arange(len(ct_df))
w = 0.35
ax.bar(x - w/2, ct_df["obs_median_mu"], w, label="obs (CRE-coarse)",
       color="firebrick", alpha=0.7)
ax.bar(x + w/2, ct_df["obsingle_median_mu"], w, label="obsingle",
       color="steelblue", alpha=0.7)
ax.set_yscale("log")
ax.set_xticks(x)
ax.set_xticklabels(ct_df["cell_type"], fontsize=8, rotation=30, ha="right")
ax.set_ylabel("Median mu", fontsize=9)
ax.set_title("D. Median mu by cell type", fontsize=10)
ax.legend(fontsize=7)

fig.suptitle("Cohen: CRE-coarse vs single-zero reporter expansion",
             fontsize=12, fontweight="bold")

out = Path(__file__).parent / "cohen_obs_vs_obsingle.svg"
fig.savefig(out, format="svg", bbox_inches="tight")
print(f"Saved {out}")
