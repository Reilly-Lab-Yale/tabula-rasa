#!/usr/bin/env python3
"""Zhao et al. (code name cohen) pairwise MWU power heatmaps, one SVG per cell
type into output/panels_mwu/.

The reference-cell-type panel is one third of manuscript Fig 3B; all four are
supplementary Fig S6. No simulation -- a replot from the cached aggregate.

Drawing lives in ../panel_style.py, shared with the other two datasets.

    python replot_panels_mwu.py
"""

import sys
from pathlib import Path

import pandas as pd

_repo_root = Path(__file__).resolve().parents[4]
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import scMPRAforge as scm
import panel_style

DATASET = "cohen"
OUTPUT_DIR = Path(__file__).resolve().parent / "output"
PANEL_DIR = OUTPUT_DIR / "panels_mwu"


def main():
    power_df = pd.read_parquet(OUTPUT_DIR / "cohen_pairwise_power_df.parquet")
    panel_style.render_dataset(
        power_df, DATASET,
        scm.COHEN_BOUNDS.cells_per_cell_type,
        PANEL_DIR,
    )


if __name__ == "__main__":
    main()
