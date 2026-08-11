"""
Activity volcanoes for the three empirical datasets.

Each CRE vs the flattened negative-control CRE ("reference"), within each cell
type, tested by Mann-Whitney U. No filtering: TSVs are loaded as-is, so each
dataset carries whatever zeros its assay design produced (shendure: oBC
barcode-level, cohen: U6 CRE-coarse, seelig: none).

Output layout, per dataset:
    output/<dataset>/<dataset>_activity_mwu.tsv        full results
    output/<dataset>/<dataset>_activity_volcano.svg    all cell types pooled
    output/<dataset>/<dataset>_activity_volcano_by_cell_type.svg   same, hued by cell type
    output/<dataset>/cell_type_specific/*.svg          one panel per cell type
"""
import sys, os
sys.path.insert(0, "/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scMPRAforge as scm  # noqa: F401  (sets editable-text rcParams on import)

DATA = "/nfs/roberts/project/pi_skr2/shared/tabula_data_new"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
ALPHA = 0.05
# Shendure's p-values run past double precision (mannwhitneyu returns exactly 0),
# which piles points on one line and crushes the rest of the axis. Panels whose
# -log10(BH p) exceeds this are clipped here, off-scale points drawn as triangles.
Y_CAP = 50.0

# Categorical hues, fixed order, for the by-cell-type overall panel. Validated
# all-pairs (every series can neighbour every other in a scatter): worst CVD
# dE 13.0, worst normal-vision dE 16.3. Only five hues clear that bar, so past
# five the marker shape carries identity alongside the hue.
PALETTE = ["#2a78d6", "#eda100", "#e87ba4", "#008300", "#4a3aa7"]
MARKERS = ["o", "s", "D"]

DATASETS = {
    "shendure": (f"{DATA}/shendure/shendure_processed.tsv", ["minP", "noP"]),
    "cohen": (f"{DATA}/cohen/retina_single_counting_u6.tsv", ["wt_1", "wt_2"]),
    "seelig": (f"{DATA}/seelig/seelig_scmpra_umiwise.tsv.gz", [
        "AACGCCCTCCACGGATGGGCCGGCCAATAAGAAGCGTTAGCGGACTCATGCGTTACGCGCCTCCGAGTTATGGGGGGGGAGGCGCGTATCTCGTGGAGAAGAAGCGATGTAACGCTTGGGCGATAAGCTTATAAGGAAGATATTT",
        "CCCTCGGAGTTAATAAGATACGCGGATCGATATCGGCTTGAAGAAGCGTATCTTATCTTCAGATGGGGATGTCGCGCATCCACCCAGTGGGCACCGCCGCTATAGAAGGGTGATAACGCTTCTCAGCCTTCAGGCTCTGGGTCTT",
    ]),
}


def volcano(sub, title, path, cap):
    """cap: -log10(BH p) ceiling, or None. Points above are drawn as triangles at cap."""
    fig, ax = plt.subplots(figsize=(4, 4))
    y = sub["neg_log10_bh_p"]
    off = y > cap if cap else pd.Series(False, index=sub.index)
    sig = sub["bh_p"] < ALPHA
    for mask, color in [(~sig, "lightgray"), (sig, "crimson")]:
        ax.scatter(sub.loc[mask & ~off, "log2_fc"], y[mask & ~off],
                   c=color, s=12, alpha=0.8, linewidths=0, rasterized=True)
        ax.scatter(sub.loc[mask & off, "log2_fc"], np.full(int((mask & off).sum()), cap),
                   c=color, s=22, marker="^", alpha=0.8, linewidths=0, clip_on=False)
    ax.axhline(-np.log10(ALPHA), color="k", ls="--", lw=0.6)
    ax.axvline(0, color="k", ls="--", lw=0.6)
    ax.set_xlabel("log2(FC) vs negative control")
    ax.set_ylabel("-log10(BH p)")
    if off.any():
        ax.set_ylim(top=cap)
        ax.set_ylabel(f"-log10(BH p)   (^ = BH p < 1e-{cap:g}, n={int(off.sum())})")
    ax.set_title(f"{title}\n{int(sig.sum())}/{len(sub)} sig at BH<{ALPHA}", fontsize=9)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    print(f"  saved {path}", flush=True)


def volcano_by_cell_type(sub, title, path, cap):
    """As volcano(), but significant points take a per-cell-type hue/marker."""
    fig, ax = plt.subplots(figsize=(5.8, 4))
    y = sub["neg_log10_bh_p"]
    off = y > cap if cap else pd.Series(False, index=sub.index)
    yc = y.clip(upper=cap) if cap else y
    sig = sub["bh_p"] < ALPHA

    ax.scatter(sub.loc[~sig, "log2_fc"], yc[~sig], c="lightgray", s=12,
               alpha=0.7, linewidths=0, rasterized=True, label="n.s.")
    cts = sorted(sub["comparison_cell_type"].unique())
    assert len(cts) <= len(PALETTE) * len(MARKERS), f"{len(cts)} cell types exceeds palette"
    for i, ct in enumerate(cts):
        m = sig & (sub["comparison_cell_type"] == ct)
        ax.scatter(sub.loc[m, "log2_fc"], yc[m], c=PALETTE[i % len(PALETTE)],
                   marker=MARKERS[i // len(PALETTE)], s=14, alpha=0.8,
                   linewidths=0, rasterized=True, label=str(ct))

    ax.axhline(-np.log10(ALPHA), color="k", ls="--", lw=0.6)
    ax.axvline(0, color="k", ls="--", lw=0.6)
    ax.set_xlabel("log2(FC) vs negative control")
    ax.set_ylabel("-log10(BH p)")
    if off.any():
        ax.set_ylim(top=cap)
        ax.set_ylabel(f"-log10(BH p)   (top row: BH p < 1e-{cap:g}, n={int(off.sum())})")
    ax.set_title(f"{title}\n{int(sig.sum())}/{len(sub)} sig at BH<{ALPHA}", fontsize=9)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=7,
              frameon=False, markerscale=1.6, handletextpad=0.2)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    print(f"  saved {path}", flush=True)


def main():
    for name, (tsv, neg_controls) in DATASETS.items():
        print(f"=== {name} ===", flush=True)
        ds_dir = f"{OUT}/{name}"
        ct_dir = f"{ds_dir}/cell_type_specific"
        os.makedirs(ct_dir, exist_ok=True)

        dat = scm.scMPRA_data.from_tsv(tsv)
        dat.set_negative_controls(neg_controls)

        hs = scm.make_all_by_celltype_hypotheses(counts=dat, reference_cre="reference")
        print(f"  {len(hs)} hypotheses", flush=True)

        df = scm.HypothesisTester("mwu").run(hs, dat).to_dataframe()
        assert len(df) == len(hs), f"MWU changed row count: {len(hs)} -> {len(df)}"
        assert (df["comparison_CRE"] != "reference").all(), "reference CRE tested against itself"

        df["log2_fc"] = np.log2(df["fold_change"])
        df["neg_log10_bh_p"] = -np.log10(df["bh_p"].clip(lower=1e-300))
        n_na = int(df["p_value"].isna().sum())
        print(f"  NaN p (empty group on one side): {n_na}", flush=True)
        print(f"  significant at BH<{ALPHA}: {int((df['bh_p'] < ALPHA).sum())}", flush=True)

        df.to_csv(f"{ds_dir}/{name}_activity_mwu.tsv", sep="\t", index=False)

        plottable = df[np.isfinite(df["log2_fc"]) & df["bh_p"].notna()]
        assert len(plottable) > 0, f"{name}: nothing plottable"

        # One ceiling per dataset so all its panels share a y scale.
        cap = Y_CAP if plottable["neg_log10_bh_p"].max() > Y_CAP else None
        print(f"  y cap: {cap}", flush=True)

        volcano(plottable, f"{name}: all cell types",
                f"{ds_dir}/{name}_activity_volcano.svg", cap)
        volcano_by_cell_type(plottable, f"{name}: all cell types",
                             f"{ds_dir}/{name}_activity_volcano_by_cell_type.svg", cap)
        for ct, sub in plottable.groupby("comparison_cell_type", observed=True):
            safe = str(ct).replace(" ", "_").replace("/", "_")
            volcano(sub, f"{name}: {ct}",
                    f"{ct_dir}/{name}_{safe}_activity_volcano.svg", cap)

    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
