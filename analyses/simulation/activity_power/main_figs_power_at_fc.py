#!/usr/bin/env python3
"""
Main-text figures: reporter effect (Fig A) and dataset comparison (Fig B).

Both use power at FC=1.5 (upregulation), computed from existing parquets.
No re-simulation. FC=1.5 is in the sim range for all three datasets (Cohen
caps at 1.6) and gives a clean 3-tier visual under the "bigger=better"
convention: Cohen near-saturated, Shendure mid-range, Seelig floored.

Fig A: Reporter dumbbell for cohen + shendure (the two datasets with
       reporters). Per cell type, -reporter dot and +reporter dot connected
       by a line; line length is the reporter's power contribution.
Fig B: Best-available power@FC=1.5 bars across all three datasets.
       cohen+rep, shendure+rep, seelig-rep. Horizontal bars, bigger=better.

Full butterflies (FC@80% power, both directions) live in per-dataset
subdirs and are supplemental.

Usage:
    python main_figs_power_at_fc.py
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

TARGET_FC = 1.5
FC_TOL = 0.05  # +/- 5% window around target

SCRIPT_DIR = Path(__file__).resolve().parent
OUT_DIR = SCRIPT_DIR / "output"
OUT_DIR.mkdir(exist_ok=True)

DATASETS = {
    "cohen": {
        "reporter": SCRIPT_DIR / "cohen/output/cohen_power_df_reporter.parquet",
        "deflated": SCRIPT_DIR / "cohen/output/cohen_power_df_deflated.parquet",
        "ref_display": "Rod",
        "title": "Cohen (episomal, CRE-reporter)",
    },
    "shendure": {
        "reporter": SCRIPT_DIR / "shendure/output/power_df_reporter.parquet",
        "deflated": SCRIPT_DIR / "shendure/output/power_df_deflated.parquet",
        "ref_display": "Pluripotent",
        "title": "Shendure (piggyBac, oBC reporter)",
    },
    "seelig": {
        "deflated": SCRIPT_DIR / "seelig/output/seelig_power_df_deflated.parquet",
        "ref_display": "HepG2",
        "title": "Seelig (episomal, no reporter)",
    },
}

COLOR_REP = "steelblue"
COLOR_DEFL = "coral"


def power_at_fc(df, target_fc=TARGET_FC, tol=FC_TOL):
    """Mean reject_null over a window of fc around target, per cell type."""
    lo = target_fc * (1 - tol)
    hi = target_fc * (1 + tol)
    window = df[(df["fc"] >= lo) & (df["fc"] <= hi)]
    return (
        window.groupby("cell_type")["reject_null"]
        .mean()
        .rename("power")
        .reset_index()
    )


def rename_reference(df, display):
    df = df.copy()
    df.loc[df["cell_type"] == "reference", "cell_type"] = display
    return df


def load_powers():
    out = {}
    for name, d in DATASETS.items():
        out[name] = {}
        if "reporter" in d:
            rep = pd.read_parquet(d["reporter"])
            out[name]["reporter"] = rename_reference(
                power_at_fc(rep), d["ref_display"]
            )
        defl = pd.read_parquet(d["deflated"])
        out[name]["deflated"] = rename_reference(
            power_at_fc(defl), d["ref_display"]
        )
    return out


def fig_a_reporter_dumbbell(powers):
    """Fig A: +rep vs -rep at FC=1.5, cohen + shendure panels."""
    dsets = ["cohen", "shendure"]
    n_cts = [len(powers[d]["reporter"]) for d in dsets]
    fig, axes = plt.subplots(
        2, 1, figsize=(6, 5),
        gridspec_kw={"height_ratios": n_cts},
        sharex=True,
    )

    for ax, name in zip(axes, dsets):
        rep = powers[name]["reporter"].set_index("cell_type")["power"]
        defl = powers[name]["deflated"].set_index("cell_type")["power"]
        # Sort by +reporter power ascending
        cts = rep.sort_values().index.tolist()
        y = np.arange(len(cts))
        rep_vals = rep.loc[cts].values
        defl_vals = defl.loc[cts].values

        for yi, dv, rv in zip(y, defl_vals, rep_vals):
            ax.plot([dv, rv], [yi, yi], color="gray", lw=1.5, zorder=1)
        ax.scatter(
            defl_vals, y, color=COLOR_DEFL, s=55,
            label="-reporter", zorder=3, edgecolor="white", lw=0.6,
        )
        ax.scatter(
            rep_vals, y, color=COLOR_REP, s=55,
            label="+reporter", zorder=3, edgecolor="white", lw=0.6,
        )
        ax.set_yticks(y)
        ax.set_yticklabels(cts, fontsize=9)
        ax.set_title(DATASETS[name]["title"], fontsize=10, loc="left")
        ax.axvline(0.8, color="black", ls="--", lw=0.6, alpha=0.4)
        ax.grid(axis="x", alpha=0.3)
        ax.set_xlim(-0.02, 1.02)

    axes[0].legend(
        loc="lower left", fontsize=8, frameon=True,
        bbox_to_anchor=(0.0, 1.15), ncol=2,
    )
    axes[-1].set_xlabel(f"Power at FC = {TARGET_FC}", fontsize=10)
    plt.tight_layout()
    out = OUT_DIR / "fig_a_reporter_dumbbell.svg"
    fig.savefig(out, format="svg", bbox_inches="tight")
    print(f"Saved: {out}")
    plt.close(fig)


def fig_b_dataset_bars(powers):
    """Fig B: best-available power@FC=1.5 across all three datasets."""
    best_mode = {
        "cohen": "reporter",
        "shendure": "reporter",
        "seelig": "deflated",
    }
    mode_label = {"reporter": "+reporter", "deflated": "-reporter"}
    mode_color = {"reporter": COLOR_REP, "deflated": COLOR_DEFL}

    dsets = ["cohen", "shendure", "seelig"]
    n_cts = [len(powers[d][best_mode[d]]) for d in dsets]
    heights = [max(n, 1) for n in n_cts]

    fig, axes = plt.subplots(
        3, 1, figsize=(6, 5.5),
        gridspec_kw={"height_ratios": heights},
        sharex=True,
    )

    for ax, name in zip(axes, dsets):
        mode = best_mode[name]
        df = powers[name][mode].sort_values("power")
        cts = df["cell_type"].tolist()
        vals = df["power"].values
        y = np.arange(len(cts))
        ax.barh(
            y, vals, color=mode_color[mode],
            edgecolor="white", lw=0.6, height=0.7,
        )
        ax.set_yticks(y)
        ax.set_yticklabels(cts, fontsize=9)
        ax.set_xlim(0, 1.02)
        ax.axvline(0.8, color="black", ls="--", lw=0.6, alpha=0.4)
        ax.set_title(
            f"{DATASETS[name]['title']}  ({mode_label[mode]})",
            fontsize=9, loc="left",
        )
        ax.grid(axis="x", alpha=0.3)
        ax.invert_yaxis()  # strongest at top

    axes[-1].set_xlabel(f"Power at FC = {TARGET_FC}", fontsize=10)
    plt.tight_layout()
    out = OUT_DIR / "fig_b_dataset_bars.svg"
    fig.savefig(out, format="svg", bbox_inches="tight")
    print(f"Saved: {out}")
    plt.close(fig)


def main():
    powers = load_powers()

    print("=== Power at FC = {:.2f} (+/- {:.0%}) ===".format(TARGET_FC, FC_TOL))
    for name, modes in powers.items():
        for mode, df in modes.items():
            print(f"\n{name} ({mode}):")
            print(df.to_string(index=False))
    print()

    fig_a_reporter_dumbbell(powers)
    fig_b_dataset_bars(powers)


if __name__ == "__main__":
    main()
