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
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
from pathlib import Path

TARGET_FC = 1.5
FC_TOL = 0.05  # +/- 5% window around target

# Per-dataset FC triplets for multi-FC variants. Each spans the dataset's
# informative regime. Cohen caps at 1.6; seelig needs higher FCs to show any
# power; shendure is in the middle.
FC_TRIPLET = {
    "cohen":    [1.1, 1.2, 1.3],
    "shendure": [1.2, 1.5, 2.0],
    "seelig":   [1.5, 2.0, 3.0],
}

SCRIPT_DIR = Path(__file__).resolve().parent
OUT_DIR = SCRIPT_DIR / "output"
OUT_DIR.mkdir(exist_ok=True)

# The reporter/deflated pair is MWU. Measuring the reporter's power benefit
# with the t-test confounds it: the t-test over-rejects in exactly the deflated
# condition, so part of the apparent benefit is that miscalibration rather than
# the reporter. The 2026-08-13 run added the MWU arms for this reason -- see the
# ARMS comment in shendure/shendure_power_ttest_all_cell_types.py.
# Yin et al. is still t-test: no MWU arm has been run on its power sims.
DATASETS = {
    "cohen": {
        "reporter": SCRIPT_DIR / "output/cohen_power_df_2026-08-13_mwu.parquet",
        "deflated": SCRIPT_DIR / "output/cohen_power_df_2026-08-13_mwu_deflated.parquet",
        "ref_display": "Rod",
        "title": "Cohen (episomal, CRE-reporter)",
    },
    "shendure": {
        "reporter": SCRIPT_DIR / "shendure/output/power_df_mwu_2026-08-13.parquet",
        "deflated": SCRIPT_DIR / "shendure/output/power_df_mwu_deflated_2026-08-13.parquet",
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


def load_multi_fc_powers():
    """Return {dataset: {mode: {fc_val: power_df}}} using FC_TRIPLET."""
    out = {}
    for name, d in DATASETS.items():
        out[name] = {}
        triplet = FC_TRIPLET[name]
        for parquet_key in ("reporter", "deflated"):
            if parquet_key not in d:
                continue
            raw = pd.read_parquet(d[parquet_key])
            per_fc = {}
            for fc in triplet:
                per_fc[fc] = rename_reference(
                    power_at_fc(raw, target_fc=fc, tol=FC_TOL),
                    d["ref_display"],
                )
            out[name][parquet_key] = per_fc
    return out


def _shades(base_color, n=3):
    """Build n shades from light to dark for a base matplotlib color."""
    import matplotlib.colors as mcolors
    rgb = np.array(mcolors.to_rgb(base_color))
    # fractions: larger -> closer to base (darker).
    fracs = np.linspace(0.30, 0.90, n)
    return [tuple(1.0 - f * (1.0 - rgb)) for f in fracs]


def fig_b_telescope(mpows):
    """Telescope (nested) bars: 3 overlaid bars per CT, one per FC."""
    best_mode = {
        "cohen": "reporter",
        "shendure": "reporter",
        "seelig": "deflated",
    }
    base_color = {"reporter": COLOR_REP, "deflated": COLOR_DEFL}

    dsets = ["cohen", "shendure", "seelig"]
    n_cts = [len(mpows[d][best_mode[d]][FC_TRIPLET[d][0]]) for d in dsets]
    heights = [max(n, 1) for n in n_cts]

    fig, axes = plt.subplots(
        3, 1, figsize=(6.5, 6),
        gridspec_kw={"height_ratios": heights},
        sharex=True,
    )

    for ax, name in zip(axes, dsets):
        mode = best_mode[name]
        triplet = FC_TRIPLET[name]  # ascending
        light, mid, dark = _shades(base_color[mode], n=3)

        # Sort CTs by the highest-FC power (strongest at top).
        sort_df = mpows[name][mode][triplet[-1]].sort_values("power")
        cts = sort_df["cell_type"].tolist()
        y = np.arange(len(cts))

        fc_lo, fc_mi, fc_hi = triplet
        hi = mpows[name][mode][fc_hi].set_index("cell_type").loc[cts, "power"].values
        mi = mpows[name][mode][fc_mi].set_index("cell_type").loc[cts, "power"].values
        lo = mpows[name][mode][fc_lo].set_index("cell_type").loc[cts, "power"].values

        # Draw longest first so shorter bars paint on top.
        ax.barh(y, hi, color=light, edgecolor="white", lw=0.4,
                height=0.75, label=f"FC={fc_hi}")
        ax.barh(y, mi, color=mid, edgecolor="white", lw=0.4,
                height=0.75, label=f"FC={fc_mi}")
        ax.barh(y, lo, color=dark, edgecolor="white", lw=0.4,
                height=0.75, label=f"FC={fc_lo}")

        ax.set_yticks(y)
        ax.set_yticklabels(cts, fontsize=9)
        ax.set_xlim(0, 1.02)
        ax.axvline(0.8, color="black", ls="--", lw=0.6, alpha=0.4)
        ax.set_title(DATASETS[name]["title"], fontsize=9, loc="left")
        ax.grid(axis="x", alpha=0.3)
        ax.invert_yaxis()
        ax.legend(
            loc="lower right", fontsize=7, frameon=True, ncol=3,
            handlelength=1.2, handletextpad=0.4, columnspacing=0.8,
        )

    axes[-1].set_xlabel("Power", fontsize=10)
    plt.tight_layout()
    out = OUT_DIR / "fig_b_telescope.svg"
    fig.savefig(out, format="svg", bbox_inches="tight")
    print(f"Saved: {out}")
    plt.close(fig)


def fig_b_slope(mpows):
    """Slope plot: one panel per dataset, x=FC triplet, y=power, lines per CT."""
    best_mode = {
        "cohen": "reporter",
        "shendure": "reporter",
        "seelig": "deflated",
    }
    base_color = {"reporter": COLOR_REP, "deflated": COLOR_DEFL}
    ref_display = {d: DATASETS[d]["ref_display"] for d in DATASETS}

    dsets = ["cohen", "shendure", "seelig"]
    fig, axes = plt.subplots(1, 3, figsize=(9, 3.6), sharey=True)

    for ax, name in zip(axes, dsets):
        mode = best_mode[name]
        triplet = FC_TRIPLET[name]
        color = base_color[mode]
        ref_ct = ref_display[name]

        cts = mpows[name][mode][triplet[0]]["cell_type"].tolist()
        vals_by_ct = {
            ct: [
                float(mpows[name][mode][fc].set_index("cell_type").loc[ct, "power"])
                for fc in triplet
            ]
            for ct in cts
        }

        x = np.arange(len(triplet))
        for ct, v in vals_by_ct.items():
            is_ref = (ct == ref_ct)
            ax.plot(
                x, v,
                marker="o",
                color=("black" if is_ref else color),
                alpha=(1.0 if is_ref else 0.55),
                lw=(1.8 if is_ref else 1.0),
                markersize=(6 if is_ref else 4),
                zorder=(3 if is_ref else 2),
                label=(ref_ct if is_ref else None),
            )

        ax.set_xticks(x)
        ax.set_xticklabels([f"{fc}" for fc in triplet], fontsize=9)
        ax.set_xlabel("Fold change", fontsize=9)
        ax.set_xlim(-0.25, len(triplet) - 0.75)
        ax.set_ylim(0, 1.05)
        ax.axhline(0.8, color="black", ls="--", lw=0.6, alpha=0.4)
        ax.set_title(DATASETS[name]["title"], fontsize=9)
        ax.grid(axis="y", alpha=0.3)
        ax.legend(
            loc="lower right", fontsize=7, frameon=True,
            title="reference CT", title_fontsize=7,
        )

    axes[0].set_ylabel("Power", fontsize=10)
    plt.tight_layout()
    out = OUT_DIR / "fig_b_slope.svg"
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

    mpows = load_multi_fc_powers()
    fig_b_telescope(mpows)
    fig_b_slope(mpows)


if __name__ == "__main__":
    main()
