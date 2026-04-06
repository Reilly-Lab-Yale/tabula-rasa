"""
MWU hypothesis testing for Takeshi emVar allelic comparisons.

For each locus x cell type, tests A-allele vs R-allele UMI counts
using Mann-Whitney U, then compares to bulk MPRA emVar calls.

No ortho fitting required -- operates directly on scMPRA_data.
"""

import sys
sys.path.insert(0, "/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa-takeshi")

import pandas as pd
import numpy as np
import scMPRAforge as scm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DATA_ROOT = "/nfs/roberts/project/pi_skr2/shared/tabula_data_new/takeshi"
SCMPRA_PATH = f"{DATA_ROOT}/takeshi_scmpra_umiwise.scmpra"
BULK_EMVAR = f"{DATA_ROOT}/bulk_emvar.xls"
BULK_ACTIVITY = f"{DATA_ROOT}/bulk_activity.xls"
OUT_DIR = "/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa-takeshi/notebooks/results/takeshi"

# ---------------------------------------------------------------------------
# Load scMPRA data
# ---------------------------------------------------------------------------
print("Loading scMPRA data...")
dat = scm.scMPRA_data.from_parquet(SCMPRA_PATH)
df = dat.data.compute() if hasattr(dat.data, 'compute') else dat.data

cell_types = sorted(df["cell_type"].unique())
print(f"  Cell types: {cell_types}")

# ---------------------------------------------------------------------------
# Build R/A locus pairs from CRE names
# ---------------------------------------------------------------------------
cre_ids = sorted(df["cre_id"].unique())
# Variant CREs follow pattern: chr:pos:ref:alt:R:wC or chr:pos:ref:alt:A:wC
variant_cres = [c for c in cre_ids if c.endswith(":wC")]

# Extract locus (strip :R:wC or :A:wC)
locus_map = {}  # locus -> {"R": cre_id, "A": cre_id}
for cre in variant_cres:
    parts = cre.rsplit(":", 2)  # [locus, R_or_A, wC]
    if len(parts) == 3 and parts[1] in ("R", "A"):
        locus = parts[0]
        allele = parts[1]
        locus_map.setdefault(locus, {})[allele] = cre

# Keep only loci with both R and A
paired_loci = {loc: alleles for loc, alleles in locus_map.items()
               if "R" in alleles and "A" in alleles}
print(f"  Paired R/A loci: {len(paired_loci)}")

orphan_loci = {loc: alleles for loc, alleles in locus_map.items()
               if not ("R" in alleles and "A" in alleles)}
if orphan_loci:
    print(f"  Orphan loci (missing one allele): {orphan_loci}")

# ---------------------------------------------------------------------------
# Build hypothesis set: A vs R within each cell type
# ---------------------------------------------------------------------------
hyp_rows = []
for locus, alleles in sorted(paired_loci.items()):
    for ct in cell_types:
        hyp_rows.append({
            "comparison_CRE": alleles["A"],
            "comparison_cell_type": ct,
            "reference_CRE": alleles["R"],
            "reference_cell_type": ct,
            "meta": locus,
        })

hyp_df = pd.DataFrame(hyp_rows)
hs = scm.HypothesisSet.from_dataframe(hyp_df)
print(f"  Hypotheses: {len(hyp_df)} ({len(paired_loci)} loci x {len(cell_types)} cell types)")

# ---------------------------------------------------------------------------
# Run MWU
# ---------------------------------------------------------------------------
print("Running MWU...")
runner = scm.HypothesisTester("mwu", alternative="two-sided")
mwu_res = runner.run(hs, dat)
mwu_df = mwu_res.to_dataframe()

# Log fold change
mwu_df["log2_fc"] = np.log2(mwu_df["fold_change"])

print(f"  Results: {len(mwu_df)} rows")
print(f"  Significant (bh_p < 0.05): {(mwu_df['bh_p'] < 0.05).sum()}")
print(f"  NaN p-values: {mwu_df['p_value'].isna().sum()}")

# Save MWU results
mwu_out = f"{OUT_DIR}/mwu_emvar_results.tsv"
mwu_df.to_csv(mwu_out, sep="\t", index=False)
print(f"  Saved to {mwu_out}")

# ---------------------------------------------------------------------------
# Load bulk emVar data
# ---------------------------------------------------------------------------
print("\nLoading bulk emVar data...")
bulk = pd.read_csv(BULK_EMVAR, sep="\t", index_col=0)
print(f"  Bulk loci: {len(bulk)}")

# Melt bulk to long form: one row per (locus, cell_type)
bulk_long_rows = []
for ct_raw, ct_clean in [("K562", "K562"), ("HepG2", "HepG2"), ("SKNSH", "SKNSH")]:
    for locus in bulk.index:
        skew = bulk.loc[locus, f"{ct_raw}_Log2Skew"]
        logp = bulk.loc[locus, f"{ct_raw}_logPadj_BF"]
        bulk_long_rows.append({
            "locus": locus,
            "cell_type": ct_clean,
            "bulk_log2skew": skew,
            "bulk_logPadj_BF": logp,
        })

bulk_long = pd.DataFrame(bulk_long_rows)
# bulk uses logPadj_BF = -log10(padj) with BF correction; 0 means not significant
# significant if logPadj_BF > -log10(0.05) ~ 1.301, but convention may vary
# Let's use logPadj_BF > 0 as "tested" and a threshold for significance
bulk_long["bulk_sig"] = bulk_long["bulk_logPadj_BF"] > 1.301  # p < 0.05

# ---------------------------------------------------------------------------
# Merge sc MWU with bulk
# ---------------------------------------------------------------------------
mwu_df["locus"] = mwu_df["meta"]
mwu_df["cell_type"] = mwu_df["comparison_cell_type"]
merged = mwu_df.merge(bulk_long, on=["locus", "cell_type"], how="inner")
print(f"  Merged rows: {len(merged)}")

merged["sc_sig"] = merged["bh_p"] < 0.05

# ---------------------------------------------------------------------------
# Concordance analysis
# ---------------------------------------------------------------------------
print("\n=== CONCORDANCE: sc MWU vs bulk emVar ===\n")

# Confusion matrix
for ct in cell_types:
    sub = merged[merged["cell_type"] == ct]
    ct_tab = pd.crosstab(sub["sc_sig"], sub["bulk_sig"],
                         rownames=["sc_MWU_sig"], colnames=["bulk_sig"])
    print(f"--- {ct} ---")
    print(ct_tab)
    agree = (sub["sc_sig"] == sub["bulk_sig"]).mean()
    print(f"Agreement: {agree:.3f}\n")

# Overall
print("--- Overall ---")
overall_tab = pd.crosstab(merged["sc_sig"], merged["bulk_sig"],
                          rownames=["sc_MWU_sig"], colnames=["bulk_sig"])
print(overall_tab)
agree_all = (merged["sc_sig"] == merged["bulk_sig"]).mean()
print(f"Agreement: {agree_all:.3f}\n")

# Direction concordance (among loci significant in both)
both_sig = merged[(merged["sc_sig"]) & (merged["bulk_sig"])]
if len(both_sig) > 0:
    both_sig = both_sig.copy()
    both_sig["sc_direction"] = np.sign(both_sig["log2_fc"])
    both_sig["bulk_direction"] = np.sign(both_sig["bulk_log2skew"])
    direction_agree = (both_sig["sc_direction"] == both_sig["bulk_direction"]).mean()
    print(f"Direction concordance (both significant, n={len(both_sig)}): {direction_agree:.3f}")

    for ct in cell_types:
        sub = both_sig[both_sig["cell_type"] == ct]
        if len(sub) > 0:
            da = (sub["sc_direction"] == sub["bulk_direction"]).mean()
            print(f"  {ct}: {da:.3f} (n={len(sub)})")
else:
    print("No loci significant in both sc and bulk.")

# ---------------------------------------------------------------------------
# Correlation of effect sizes
# ---------------------------------------------------------------------------
print("\n=== EFFECT SIZE CORRELATION ===\n")
from scipy.stats import pearsonr, spearmanr

for ct in cell_types:
    sub = merged[merged["cell_type"] == ct].dropna(subset=["log2_fc", "bulk_log2skew"])
    if len(sub) < 3:
        print(f"  {ct}: too few points")
        continue
    # filter infinite
    sub = sub[np.isfinite(sub["log2_fc"]) & np.isfinite(sub["bulk_log2skew"])]
    r_p, p_p = pearsonr(sub["log2_fc"], sub["bulk_log2skew"])
    r_s, p_s = spearmanr(sub["log2_fc"], sub["bulk_log2skew"])
    print(f"  {ct}: Pearson r={r_p:.3f} (p={p_p:.2e}), "
          f"Spearman rho={r_s:.3f} (p={p_s:.2e}), n={len(sub)}")

# All cell types pooled
sub_all = merged.dropna(subset=["log2_fc", "bulk_log2skew"])
sub_all = sub_all[np.isfinite(sub_all["log2_fc"]) & np.isfinite(sub_all["bulk_log2skew"])]
r_p, p_p = pearsonr(sub_all["log2_fc"], sub_all["bulk_log2skew"])
r_s, p_s = spearmanr(sub_all["log2_fc"], sub_all["bulk_log2skew"])
print(f"  Pooled: Pearson r={r_p:.3f} (p={p_p:.2e}), "
      f"Spearman rho={r_s:.3f} (p={p_s:.2e}), n={len(sub_all)}")

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
print("\nGenerating plots...")

# 1. Scatter: sc log2FC vs bulk log2Skew, per cell type
fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=False)
for i, ct in enumerate(cell_types):
    ax = axes[i]
    sub = merged[merged["cell_type"] == ct].dropna(subset=["log2_fc", "bulk_log2skew"])
    sub = sub[np.isfinite(sub["log2_fc"]) & np.isfinite(sub["bulk_log2skew"])]

    colors = []
    for _, row in sub.iterrows():
        if row["sc_sig"] and row["bulk_sig"]:
            colors.append("green")
        elif row["sc_sig"]:
            colors.append("red")
        elif row["bulk_sig"]:
            colors.append("blue")
        else:
            colors.append("gray")

    ax.scatter(sub["bulk_log2skew"], sub["log2_fc"], c=colors, alpha=0.6, s=30)
    ax.axhline(0, color="k", linestyle="--", linewidth=0.5)
    ax.axvline(0, color="k", linestyle="--", linewidth=0.5)
    ax.set_xlabel("Bulk log2(Skew)")
    ax.set_ylabel("sc MWU log2(FC)")
    ax.set_title(ct)

    # Correlation annotation
    if len(sub) >= 3:
        r_s, _ = spearmanr(sub["log2_fc"], sub["bulk_log2skew"])
        ax.text(0.05, 0.95, f"rho={r_s:.2f}", transform=ax.transAxes,
                va="top", fontsize=10)

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='green', label='Both sig', markersize=8),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='red', label='sc only', markersize=8),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', label='Bulk only', markersize=8),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', label='Neither', markersize=8),
]
fig.legend(handles=legend_elements, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.02))
fig.suptitle("emVar: sc MWU vs Bulk MPRA", fontsize=14)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/emvar_sc_vs_bulk_scatter.png", dpi=150, bbox_inches="tight")
print(f"  Saved scatter to {OUT_DIR}/emvar_sc_vs_bulk_scatter.png")

# 2. Volcano: sc MWU
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
for i, ct in enumerate(cell_types):
    ax = axes[i]
    sub = merged[merged["cell_type"] == ct].copy()
    sub = sub[np.isfinite(sub["log2_fc"])]
    neg_log10_p = -np.log10(sub["bh_p"].clip(lower=1e-300))

    colors = ["red" if sig else "gray" for sig in sub["sc_sig"]]
    ax.scatter(sub["log2_fc"], neg_log10_p, c=colors, alpha=0.5, s=20)
    ax.axhline(-np.log10(0.05), color="k", linestyle="--", linewidth=0.5)
    ax.set_xlabel("log2(FC) A vs R")
    ax.set_ylabel("-log10(BH p)")
    ax.set_title(ct)

fig.suptitle("sc MWU emVar Volcano", fontsize=14)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/emvar_mwu_volcano.png", dpi=150, bbox_inches="tight")
print(f"  Saved volcano to {OUT_DIR}/emvar_mwu_volcano.png")

plt.close("all")
print("\nDone.")
