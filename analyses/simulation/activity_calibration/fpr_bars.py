"""FPR@0.05 summary bar charts per dataset.

Two two-panel views, each chosen via --axis:
- axis=reporter: panels=tests (t-test, MWU). Within each panel, +reporter vs
  -reporter (deflated) bars per cell type. Shows that t-test miscalibration
  is concentrated in the -reporter condition.
- axis=test:     panels=conditions (+reporter, -reporter). Within each
  panel, t-test vs MWU bars per cell type.

Reads the cached null-pvals parquets ({ds}_null_pvals_{[mwu_]}{cond}.parquet).
Seelig has no transfection reporter; its axis=reporter figure is skipped,
and its axis=test figure has only the -reporter panel.

Usage:
    python fpr_bars.py --dataset shendure --axis reporter
    python fpr_bars.py --dataset shendure --axis test
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import scMPRAforge  # noqa: F401  (sets editable-text rcParams on import)

HERE = Path(__file__).resolve().parent

REFERENCE_DISPLAY = {
    "shendure": "Pluripotent",
    "cohen": "Rod",
    "seelig": "HepG2",
    "takeshi": "HepG2",
}

REPORTER_COLOR = {"reporter": "steelblue", "deflated": "coral"}
REPORTER_LABEL = {"reporter": "+reporter", "deflated": "-reporter (deflated)"}
TEST_COLOR = {"ttest": "#5e3c99", "mwu": "#e66101"}
TEST_LABEL = {"ttest": "t-test", "mwu": "MWU"}


def _pvals_path(dataset: str, test: str, condition: str) -> Path:
    test_tag = "" if test == "ttest" else f"{test}_"
    return HERE / dataset / "output" / f"{dataset}_null_pvals_{test_tag}{condition}.parquet"


def _fpr_per_ct(path: Path) -> dict[str, float]:
    df = pd.read_parquet(path)
    return {
        ct: float(np.mean(g["p_value"] < 0.05))
        for ct, g in df.groupby("cell_type", observed=True)
    }


def _ordered_cts(fpr_dicts: list[dict[str, float]]) -> list[str]:
    seen: list[str] = []
    sset: set[str] = set()
    for fd in fpr_dicts:
        for ct in fd:
            if ct not in sset:
                sset.add(ct)
                seen.append(ct)
    return ["reference"] + [c for c in seen if c != "reference" and c in sset]


def _draw_panel(ax, cts, ref_display, groups, title):
    """groups: list of (color, label, fpr_dict)"""
    n_ct = len(cts)
    n_group = len(groups)
    width = 0.8 / n_group
    x = np.arange(n_ct)
    for i, (color, label, fd) in enumerate(groups):
        h = [fd.get(ct, np.nan) for ct in cts]
        ax.bar(x + (i - (n_group - 1) / 2) * width, h, width=width,
               color=color, label=label, edgecolor="black", linewidth=0.3)
    ax.axhline(0.05, color="red", linestyle="--", lw=0.8, zorder=0)
    ax.set_xticks(x)
    ax.set_xticklabels([ref_display if c == "reference" else c for c in cts],
                       rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("FPR @ p<0.05")
    ax.set_title(title, fontsize=10)
    ax.set_ylim(0, max(0.08, ax.get_ylim()[1]))
    ax.legend(fontsize=8, loc="upper right", frameon=False)
    ax.grid(axis="y", alpha=0.3, lw=0.4)


def _figure_axis_reporter(dataset, ref_display):
    """Panels=tests; within each, +reporter vs -reporter (deflated) bars."""
    panels = []  # (panel_title, groups, fpr_dicts)
    for test in ("ttest", "mwu"):
        groups = []
        fpr_dicts = []
        for cond in ("reporter", "deflated"):
            p = _pvals_path(dataset, test, cond)
            if not p.exists():
                continue
            fd = _fpr_per_ct(p)
            fpr_dicts.append(fd)
            groups.append((REPORTER_COLOR[cond], REPORTER_LABEL[cond], fd))
        if len(groups) >= 2:
            panels.append((TEST_LABEL[test], groups, fpr_dicts))
    if not panels:
        return None
    all_dicts = [fd for _, _, fds in panels for fd in fds]
    cts = _ordered_cts(all_dicts)
    fig, axes = plt.subplots(len(panels), 1,
                             figsize=(max(5, len(cts) * 0.55), 3.5 * len(panels)),
                             sharex=True, squeeze=False)
    axes = axes[:, 0]
    for ax, (panel_title, groups, _) in zip(axes, panels):
        _draw_panel(ax, cts, ref_display, groups, panel_title)
    fig.suptitle(f"{dataset} -- null FPR by cell type: +reporter vs -reporter",
                 fontsize=11, y=1.005)
    plt.tight_layout()
    return fig


def _figure_axis_test(dataset, ref_display):
    """Panels=conditions; within each, t-test vs MWU bars."""
    panels = []
    for cond in ("reporter", "deflated"):
        groups = []
        fpr_dicts = []
        for test in ("ttest", "mwu"):
            p = _pvals_path(dataset, test, cond)
            if not p.exists():
                continue
            fd = _fpr_per_ct(p)
            fpr_dicts.append(fd)
            groups.append((TEST_COLOR[test], TEST_LABEL[test], fd))
        if len(groups) >= 2:
            panels.append((REPORTER_LABEL[cond], groups, fpr_dicts))
    if not panels:
        return None
    all_dicts = [fd for _, _, fds in panels for fd in fds]
    cts = _ordered_cts(all_dicts)
    fig, axes = plt.subplots(len(panels), 1,
                             figsize=(max(5, len(cts) * 0.55), 3.5 * len(panels)),
                             sharex=True, squeeze=False)
    axes = axes[:, 0]
    for ax, (panel_title, groups, _) in zip(axes, panels):
        _draw_panel(ax, cts, ref_display, groups, panel_title)
    fig.suptitle(f"{dataset} -- null FPR by cell type: t-test vs MWU",
                 fontsize=11, y=1.005)
    plt.tight_layout()
    return fig


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True, choices=list(REFERENCE_DISPLAY))
    parser.add_argument("--axis", required=True, choices=["reporter", "test"])
    args = parser.parse_args()

    ref_display = REFERENCE_DISPLAY[args.dataset]
    if args.axis == "reporter":
        fig = _figure_axis_reporter(args.dataset, ref_display)
    else:
        fig = _figure_axis_test(args.dataset, ref_display)

    if fig is None:
        sys.exit(f"No qualifying panels for dataset={args.dataset} axis={args.axis} "
                 f"(seelig has no reporter condition).")

    out = HERE / args.dataset / "output" / f"{args.dataset}_calibration_fpr_summary_{args.axis}.svg"
    fig.savefig(out, format="svg", bbox_inches="tight")
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
