#!/usr/bin/env python3
"""
Synthetic factorial sweep over the design-space encompassing the empirical
datasets shendure, cohen, and takeshi.

Latin Hypercube sample over 7 axes, simulate at each sample, run MWU on the
activity hypothesis set, then plot smoothed marginals (plus pairwise heatmaps
and empirical overlays). Empirical anchors: shendure-Pluripotent and cohen-Rod
(DEFAULT_ANCHORS, the methods-paper figures), plus takeshi-HepG2 (the
supplemental with_takeshi marginal variants and the takeshi attribution pairs).

The 1D sweep (analyses/simulation/design_space/shendure_1d) anchored at
shendure-Pluripotent and varied one axis at a time. This script removes the
anchor: each LHS sample is a fresh 7-tuple of axis values, with all other
axes at the LHS-drawn value rather than a fixed reference. The empirical
datasets become predictions of the meta-surface rather than starting points.

Axes and bounds vary by mode (see the *_AXIS_BOUNDS dicts). The canonical
sweep is `union`: a single 5000-point LHS over UNION_AXIS_BOUNDS, log-uniform
on every axis. The box is drawn wide enough to contain the empirical anchors
in EMPIRICAL, so each dataset is interior to the sweep rather than an
extrapolation from it. Read bounds and anchor values off those two dicts.

Three of the SVGs this writes are placed in the manuscript, and are drawn at
the width the page gives them, since the manuscript sync applies no font
scaling:

    pairwise_heatmaps_union.svg          Fig 5A, at \textwidth (6.90 in)
    marginals_auc_union_published.svg    Fig 5B, at 0.95\linewidth (6.56 in)
    pairwise_heatmaps_all_union.svg      Fig S10, at \linewidth (6.90 in)

Fig 5A shows that the axes act close to independently, which is what licenses
reading each of them on its own in Fig 5B; Fig S10 is the exhaustive version
of 5A. Same type scale and greys as activity_power/replot_overlay_ct.py and
variant_power/panel_style.py. The `manuscript` phase renders just these three.

Two axis names need glossing:
    lib_alpha_nb    NB2 alpha for the library; theta = 1/alpha
    minP            the reference element's fitted mean, not a promoter
                    strength -- see the comment on AXIS_DISPLAY

The earlier full/topup/topup3 modes (AXIS_BOUNDS / TOPUP*_AXIS_BOUNDS) were
exploratory sub-boxes used to calibrate this union box; see methods_sweeps.md.

Single test cell type. NB only. Reporter regime fixed at "obs" (best case).
MWU is the canonical test.

Usage:
    python synthetic_factorial.py <phase> [mode] [slice n_slices]
        phase:    samples | simulate | plot | replot | manuscript | all
        mode:     smoke | pilot | full | topup | topup3 | union
                  (default: full; canonical analysis sweep: union)
        slice:    0..n_slices-1   (default: 0; only meaningful for simulate)
        n_slices: total slice count (default: from MODES[mode]["n_slices"])

    samples  -- write the LHS samples_<mode>.parquet table; do nothing else.
    simulate -- run the simulation phase for one slice. Each slice is an
                independent sbatch job with its own dask scheduler+workers.
                Slice i processes rows where (idx % n_slices) == i.
    plot     -- aggregate all slices' outputs and produce SVGs.
    replot   -- re-render SVGs from the cached samples_power<mode>.parquet
                without walking scratch or starting a dask client. Fails if
                the cache is missing.
    manuscript -- as replot, but only the three SVGs the manuscript places:
                pairwise_heatmaps (Fig 5A), marginals_auc..._published
                (Fig 5B) and pairwise_heatmaps_all (Fig S10).
    all      -- samples then simulate (slice 0 only) then plot. Useful only
                when n_slices=1 (smoke). For pilot/full use launch.sh.

Why slices?  A single dask scheduler chokes when 8+ threads concurrently
submit graph work (observed 2026-05-05 pilot deadlocked at ~3-5 samples).
Multiple independent driver jobs, each with its own scheduler over its own
worker pool, sidestep all nested-dask and thread-contention pitfalls. Each
slice runs serially internally; parallelism is achieved across slices.

Sims are written under SIM_ROOT (envvar SYNTHETIC_FACTORIAL_SIM_ROOT to
override the default), so the directory can be moved (e.g. to cold storage)
without breaking the plot phase. The samples.parquet table is the source
of truth for which LHS draws exist; its sample_id values index into
SIM_ROOT/<sample_id>/.
"""

import sys
import os
import json
import gc
import pickle
import copy
import shutil
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import qmc

# File lives at <repo>/analyses/simulation/design_space/synthetic_factorial/<this>.py.
# parents[4] = repo root (containing the scMPRAforge/ package).
_repo_root = Path(__file__).resolve().parents[4]
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

import scMPRAforge as scm
from dask.distributed import Client, LocalCluster
from dask_jobqueue import SLURMCluster

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

OUT = Path(__file__).resolve().parent / "output"
OUT.mkdir(exist_ok=True)

SIM_DATE = "2026-05-05"
DEFAULT_SIM_ROOT = (
    Path("/nfs/roberts/scratch/pi_skr2/mcn26/synthetic_factorial_sims")
    / f"{SIM_DATE}_synthetic_factorial"
)
# Override default with envvar SYNTHETIC_FACTORIAL_SIM_ROOT for cold-storage
# relocation. The sim format is location-portable (de_novo_simulation re-reads
# from `location` at load time), so as long as the directory tree is intact
# (sample_id subdirs each containing 0/, 1/, ... sim dirs) the plot phase will
# pick it up wherever it lives.
# SIM_ROOT and SAMPLES_PATH are mode-suffixed (set in _apply_mode_paths) so
# smoke / pilot / full runs do not collide.
SIM_ROOT: Path = DEFAULT_SIM_ROOT
SAMPLES_PATH: Path = OUT / "samples.parquet"


def _apply_mode_paths(mode: str):
    """Mode-suffixed sim and samples paths so different scales coexist."""
    global SIM_ROOT, SAMPLES_PATH
    base = Path(os.environ.get("SYNTHETIC_FACTORIAL_SIM_ROOT", DEFAULT_SIM_ROOT))
    SIM_ROOT = base if mode == "full" else base.parent / f"{base.name}_{mode}"
    SIM_ROOT.mkdir(parents=True, exist_ok=True)
    SAMPLES_PATH = OUT / (f"samples_{mode}.parquet" if mode != "full" else "samples.parquet")
    print(f"SIM_ROOT={SIM_ROOT}", flush=True)
    print(f"SAMPLES_PATH={SAMPLES_PATH}", flush=True)

# ---------------------------------------------------------------------------
# Axis definitions
# ---------------------------------------------------------------------------

CELL_TYPE = "reference"

# (low, high, log_scale)
AXIS_BOUNDS = {
    # n_cells upper gives roughly 3x headroom over the largest empirical cell
    # type without blowing past the MWU sweet spot -- an order of magnitude
    # higher caused a multi-hour MWU stall in pilot.
    "n_cells":           (5.0e2, 5.0e4, True),
    "n_cres":            (5.0e1, 2.0e3, True),
    "bcs_per_cre":       (3.0,   5.0e2, True),
    "moi":               (0.5,   3.0e1, True),
    "lib_alpha_nb":      (0.02,  2.0,   True),
    "minP":              (0.02,  2.0,   True),
    "activity_max_mult": (2.0,   8.0,   False),
}

# Top-up sweep bounds. The original LHS box left the cohen-Rod anchor outside
# on bcs_per_cre and moi, by more than an order of magnitude on the former.
# Top-up draws 100 samples in the previously unexplored high-bcs/high-moi
# corner. Other axes use the same ranges as the main sweep so the top-up
# surfaces still combine cleanly with the original samples.
TOPUP_AXIS_BOUNDS = {
    "n_cells":           (5.0e2, 5.0e4, True),
    "n_cres":            (5.0e1, 2.0e3, True),
    # Start where the main sweep ended, extend past cohen.
    "bcs_per_cre":       (5.0e2, 5.0e4, True),
    "moi":               (3.0e1, 2.0e2, True),
    "lib_alpha_nb":      (0.02,  2.0,   True),
    "minP":              (0.02,  2.0,   True),
    "activity_max_mult": (2.0,   8.0,   False),
}

# Second top-up bounds. Extends moi, minP, and activity_max_mult so the
# takeshi-HepG2 anchor falls inside, which incidentally closes outstanding
# gaps on activity_max_mult for the other anchors too. Other axes use the
# FULL combined range from prior sweeps so topup3 samples span the 7-D
# space broadly.
TOPUP3_AXIS_BOUNDS = {
    "n_cells":           (5.0e2, 5.0e4, True),
    "n_cres":            (5.0e1, 2.0e3, True),
    "bcs_per_cre":       (3.0,   5.0e4, True),
    "moi":               (2.0e2, 3.5e2, True),
    "lib_alpha_nb":      (0.02,  2.0,   True),
    "minP":              (3.0e-3, 2.0e-2, True),
    "activity_max_mult": (1.0,   1.2e2, True),    # log, unlike the full-mode box
}

# The box enclosing the full + topup + topup3 ranges, swept in one draw.
# "Union" here means the union of the axis *ranges*, not a union of sample
# sets: a single LHS over this box is what the methods paper uses. Pooling
# separate sweeps instead induces cross-axis correlation, because a top-up
# that extends two axes at once makes them dependent in the pooled sample
# -- a high-moi draw could only have come from the top-up, and so carries
# its high bcs_per_cre as well.
# activity_max_mult is log-scaled here (vs linear in full) because the
# union range spans two orders of magnitude.
UNION_AXIS_BOUNDS = {
    "n_cells":           (5.0e2, 5.0e4, True),
    "n_cres":            (5.0e1, 2.0e3, True),
    "bcs_per_cre":       (3.0,   5.0e4, True),
    "moi":               (0.5,   3.5e2, True),
    "lib_alpha_nb":      (0.02,  2.0,   True),
    "minP":              (3.0e-3, 2.0,  True),
    "activity_max_mult": (1.0,   1.2e2, True),    # log over full union
}

AXIS_NAMES = list(AXIS_BOUNDS.keys())

# Sample / rep / sim / worker budget per mode. Default is "full"; the wrapper
# selects the mode via CLI arg.
# n_workers is per-slice (each slice has its own dask cluster).
# Total worker pool = n_slices * n_workers.
MODES = {
    "smoke": dict(n_lhs=4,    n_library_reps=2, n_sims=2, n_workers=4,  worker_mem="16G",  n_slices=1),
    "pilot": dict(n_lhs=30,   n_library_reps=3, n_sims=5, n_workers=20, worker_mem="48G",  n_slices=3),
    "full":  dict(n_lhs=1000, n_library_reps=5, n_sims=5, n_workers=50, worker_mem="48G",  n_slices=10),
    "topup": dict(n_lhs=100,  n_library_reps=5, n_sims=5, n_workers=50, worker_mem="128G", n_slices=10),
    "topup3":dict(n_lhs=120,  n_library_reps=5, n_sims=5, n_workers=50, worker_mem="128G", n_slices=10),
    # Brute-force LHS over the full union of axis ranges. n_lhs=5000 gives
    # ~5x the density of the original full sweep across a ~10x larger box.
    # worker_mem=128G to handle the heavy-tail corner without OOM (see
    # TODO.md 2026-05-10 incident).
    # n_workers=5 (not 50): n_sims=5 is the real intra-rep parallelism,
    # more workers per driver are wasted AND when many drivers run in
    # parallel as a job array they overwhelm YCRC's per-hour sbatch rate
    # limit (200/hr/user) with worker submissions. With 100 concurrent
    # array tasks x 5 workers each = 500 worker sbatches total over the
    # run, vs 100 x 50 = 5000 which hit the limit hard on 2026-05-11.
    "union": dict(n_lhs=5000, n_library_reps=5, n_sims=5, n_workers=5,  worker_mem="128G", n_slices=50),
}

# Mutable globals -- set in main() based on chosen mode. Defaults match "full"
# so the script is importable without arg parsing.
N_LHS = MODES["full"]["n_lhs"]
N_LIBRARY_REPS = MODES["full"]["n_library_reps"]
N_SIMS = MODES["full"]["n_sims"]
N_WORKERS = MODES["full"]["n_workers"]
WORKER_MEM = MODES["full"]["worker_mem"]
N_SLICES = MODES["full"]["n_slices"]
CURRENT_MODE = "full"


def _apply_mode(mode: str):
    global N_LHS, N_LIBRARY_REPS, N_SIMS, N_WORKERS, WORKER_MEM, N_SLICES, CURRENT_MODE
    if mode not in MODES:
        raise ValueError(f"Unknown mode {mode!r}; pick one of {list(MODES)}")
    cfg = MODES[mode]
    N_LHS = cfg["n_lhs"]
    N_LIBRARY_REPS = cfg["n_library_reps"]
    N_SIMS = cfg["n_sims"]
    N_WORKERS = cfg["n_workers"]
    WORKER_MEM = cfg["worker_mem"]
    N_SLICES = cfg["n_slices"]
    CURRENT_MODE = mode
    print(f"Mode: {mode!r}: {cfg}", flush=True)


def active_bounds(mode: str = None) -> dict:
    """The axis bounds the given mode draws from.

    Plotting must ask for these rather than reading AXIS_BOUNDS directly.
    The two disagree on activity_max_mult: the union sweep draws it
    log-uniformly over a much wider range than the linear one AXIS_BOUNDS
    declares. Binning a union draw on the linear declaration crowds most of
    the samples into the first bin and leaves the top bins nearly empty,
    which reads as a flat response where there is none.
    """
    mode = CURRENT_MODE if mode is None else mode
    return {
        "topup": TOPUP_AXIS_BOUNDS,
        "topup3": TOPUP3_AXIS_BOUNDS,
        "union": UNION_AXIS_BOUNDS,
    }.get(mode, AXIS_BOUNDS)

# Activity floor: scale with minP rather than absolute. Empirically
# min_mpra_umi/minP ~ 0.02-0.08 across datasets; pick a single ratio.
MIN_ACT_FRAC_OF_MINP = 0.05

SEED = 20260505

# Empirical anchor values for plot overlay, and the source of truth for every
# per-dataset number in this module -- do not restate them in comments.
#
# activity_max_mult is the data-derived ratio p95(mu) / minP, i.e. how broad
# the assayed CRE activity dynamic range is, robust to single-CRE outliers.
# The spread across datasets is real and large: Lalanne et al.'s library
# carries strong promoters producing a heavy tail, while Zhao et al.'s is
# variants of a single element, all clustered near baseline. That spread is
# why the union box draws this axis log-uniformly; the superseded full-mode
# box was linear and too narrow to contain either anchor.
EMPIRICAL = {
    "shendure-Pluripotent": dict(
        n_cells=8201, n_cres=207, bcs_per_cre=136.06, moi=18.01,
        lib_alpha_nb=0.2241, minP=0.0414, activity_max_mult=99.21,
    ),
    "cohen-Rod": dict(
        n_cells=18633, n_cres=116, bcs_per_cre=17244.5, moi=149.15,
        lib_alpha_nb=1.3873, minP=0.9363, activity_max_mult=1.11,
    ),
    # Yin et al. synthesised its barcodes rather than adding them by PCR, so
    # its library is near-deterministic and the NB fit returns alpha at the
    # optimiser floor. That value sits below UNION_AXIS_BOUNDS, so this anchor
    # is the one that clamps on lib_alpha_nb; attribution against it should
    # treat that axis as not applicable rather than clamped.
    "seelig-HepG2": dict(
        n_cells=5858, n_cres=1344, bcs_per_cre=4.945, moi=73.57,
        lib_alpha_nb=1e-06, minP=0.2092, activity_max_mult=15.04,
    ),
    # Unpublished as of 2026-05-07; the with_takeshi figure variant is
    # supplemental, and the methods-paper version uses no_takeshi
    # (DEFAULT_ANCHORS below).
    "takeshi-HepG2": dict(
        n_cells=1690, n_cres=149, bcs_per_cre=394.07, moi=265.91,
        lib_alpha_nb=1.7174, minP=0.0116, activity_max_mult=37.22,
    ),
}

# Methods-paper figures use only shendure + cohen overlays. with_takeshi
# variants (supplemental) overlay all three. Anchor names index into
# EMPIRICAL above.
# Display names for plots. The internal keys are the column names in the
# samples parquets and are left alone so cached sweeps stay readable.
#
# "minP" is renamed on sight because the name misleads. It is the fitted mean
# of the reference element, i.e. the construct's output with no CRE, and the
# simulation draws every element's activity as a multiple of it
# (activity_max_mult * minP). It is therefore the level the whole assay
# operates at, not a promoter one chooses independently of the elements. That
# is not an extra assumption of the sweep: the GLM adds the intercept and the
# CRE coefficient in log space, so basal level and fold change are
# multiplicative on the count scale throughout this work.
#
# "moi" is the transfection NB mean, i.e. the number of library constructs
# delivered per cell (see make_synthetic_bounds -> set_effective_moi), so the
# gloss names the count rather than repeating the acronym.
AXIS_DISPLAY = {
    "minP": "basal expression",
    "lib_alpha_nb": "library skew (NB alpha)",
    "bcs_per_cre": "barcodes per element",
    "n_cres": "elements",
    "n_cells": "cells",
    "moi": "MOI (constructs per cell)",
    "activity_max_mult": "dynamic range (max / basal)",
}

# Short forms for the heatmap grids, where an axis label is rotated into a
# panel barely an inch tall and the long form runs off both ends of it. The
# two that are still too long once the gloss is dropped are wrapped rather
# than abbreviated further.
AXIS_DISPLAY_SHORT = {
    "lib_alpha_nb": "library skew",
    "bcs_per_cre": "barcodes\nper element",
    "moi": "MOI",
    "minP": "basal\nexpression",
    "activity_max_mult": "dynamic range",
}


def axis_label(axis: str) -> str:
    return AXIS_DISPLAY.get(axis, axis)


def axis_label_short(axis: str) -> str:
    return AXIS_DISPLAY_SHORT.get(axis, axis_label(axis))


DEFAULT_ANCHORS = ["shendure-Pluripotent", "cohen-Rod"]
ALL_ANCHORS = list(EMPIRICAL.keys())
# The three published designs, which is what the manuscript overlays.
PUBLISHED_ANCHORS = ["shendure-Pluripotent", "cohen-Rod", "seelig-HepG2"]

# The house per-dataset colours, the same three Okabe-Ito hues and the same
# assignment as model_selection/cross_family_agreement.py, so a dataset keeps
# its colour across the manuscript. The triple passes every categorical check
# on white with all pairs in play (worst CVD dE 11.0, normal-vision 18.7, all
# above 3:1). Takeshi is supplemental only; adding it puts the worst pair in
# the 6-8 CVD band, which the single-letter tag beside each rule covers.
EMPIRICAL_COLORS = {
    "shendure-Pluripotent": "#0072b2",
    "cohen-Rod":            "#d55e00",
    "seelig-HepG2":         "#009e73",
    "takeshi-HepG2":        "#cc79a7",
}

# Author-citation form plus the one-letter tag that marks the anchor's rule
# inside a panel. Written here rather than left to the manuscript's
# dataset_names SVG transform, which only rewrites text that survived as
# <text> nodes. The tag carries the identity that colour alone would not.
ANCHOR_DISPLAY = {
    "shendure-Pluripotent": ("L", "Lalanne et al."),
    "cohen-Rod":            ("Z", "Zhao et al."),
    "seelig-HepG2":         ("Y", "Yin et al."),
    "takeshi-HepG2":        ("T", "Takeshi (unpublished)"),
}

# ---------------------------------------------------------------------------
# Figure style
#
# Both manuscript panels are drawn at the width the page gives them -- Fig 5A
# at \textwidth of a 6.90 in text block, Fig 5B at 0.95\linewidth of it -- and
# the manuscript sync applies no font scaling, so a point size set here is the
# point size a reader gets. Same type scale and same ink/muted/hairline greys
# as activity_power/replot_overlay_ct.py and variant_power/panel_style.py.
# ---------------------------------------------------------------------------

INK, MUTED, HAIR = "#1a1a1a", "#666666", "#c9c9c9"
# Scatter is background texture behind the smoother, not a series.
DOT = "#c2c2c2"

# Power is a magnitude on a fixed 0-1 scale, so the ramp is sequential and
# single-hue. YlOrRd is what the Fig 3B power heatmaps already use for the
# same quantity.
POWER_CMAP = "YlOrRd"

plt.rcParams.update({
    "svg.fonttype": "none",       # keep text editable downstream
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.size": 7,
    "axes.labelsize": 7.5,
    "axes.titlesize": 7.5,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5,
    "axes.unicode_minus": False,  # the repo is plain ASCII
})


def _fmt_decade(k: int) -> str:
    """A power of ten as a plain decimal.

    Plain rather than mathtext: svg.fonttype="none" leaves mathtext as text
    nodes in a maths font the manuscript's font transform does not know about,
    and the repo is plain ASCII besides.
    """
    return f"{10.0 ** k:g}" if k >= 0 else f"{10.0 ** k:.{-k}f}"


def _decade_ticks(lo: float, hi: float, max_ticks: int = 4) -> list:
    """Decades spanning [lo, hi] in log10 units, thinned to at most max_ticks.

    Thinned by stride so the ladder stays evenly spaced; the ends are what a
    reader uses to place the axis, so the first decade in range is always kept.
    """
    decades = list(range(int(np.ceil(lo)), int(np.floor(hi)) + 1))
    if len(decades) > max_ticks:
        stride = int(np.ceil(len(decades) / max_ticks))
        decades = decades[::stride]
    return decades


def _assert_labels_fit(fig, what):
    """No axis label may be longer than the panel edge it names.

    Nothing is cropped when one is -- the figure is saved at a fixed size, so
    an over-long label simply runs into its neighbour, which is exactly the
    collision that is hard to notice in a grid of twenty-one panels.
    """
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    for ax in fig.axes:
        box = ax.get_window_extent(r)
        for label, along in ((ax.yaxis.label, box.height),
                             (ax.xaxis.label, box.width)):
            if not label.get_text():
                continue
            ext = label.get_window_extent(r)
            run = ext.height if label is ax.yaxis.label else ext.width
            assert run <= along, (
                f"{what}: label {label.get_text()!r} runs {run / fig.dpi:.2f} "
                f"in along a {along / fig.dpi:.2f} in panel edge")


def _assert_text_fits(fig, texts, limit_in, what):
    """No header line may run past the width it was written to fit."""
    fig.canvas.draw()
    for t in texts:
        right = t.get_window_extent(fig.canvas.get_renderer()).x1 / fig.dpi
        assert right <= limit_in, (
            f"{what}: {t.get_text()[:40]!r} ends at {right:.2f} in, past the "
            f"{limit_in:.2f} in it has")


# ---------------------------------------------------------------------------
# LHS sampling
# ---------------------------------------------------------------------------

def draw_lhs(n: int, seed: int = SEED, bounds: dict = None,
             sample_id_prefix: str = "s") -> pd.DataFrame:
    """LHS draw over the supplied axis bounds.

    bounds: dict[axis] -> (low, high, log_scale). Defaults to AXIS_BOUNDS.
            For the bounds-expansion top-up sweep, pass TOPUP_AXIS_BOUNDS.
    sample_id_prefix: one letter per sweep ("s" full, "t" topup, "u"
            topup3, "v" union), so a parquet can be traced to the draw that
            produced it. The manuscript uses "v" exclusively.
    """
    if bounds is None:
        bounds = AXIS_BOUNDS
    sampler = qmc.LatinHypercube(d=len(AXIS_NAMES), seed=seed)
    u = sampler.random(n=n)  # (n, d) in [0, 1)
    rows = {}
    for i, ax in enumerate(AXIS_NAMES):
        lo, hi, logsc = bounds[ax]
        if logsc:
            vals = np.exp(np.log(lo) + u[:, i] * (np.log(hi) - np.log(lo)))
        else:
            vals = lo + u[:, i] * (hi - lo)
        rows[ax] = vals
    df = pd.DataFrame(rows)
    # Round integer-valued axes
    df["n_cells"] = df["n_cells"].round().astype(int)
    df["n_cres"] = df["n_cres"].round().astype(int)
    df.insert(0, "sample_id", [f"{sample_id_prefix}{i:04d}" for i in range(n)])
    return df


# ---------------------------------------------------------------------------
# Bounds construction
# ---------------------------------------------------------------------------

def make_synthetic_bounds(sample: pd.Series):
    """Mutate a deep copy of SHENDURE_BOUNDS to match this sample's axes."""
    b = scm.SHENDURE_BOUNDS.copy()

    # Single test cell type
    if CELL_TYPE not in b.cells_per_cell_type.index:
        # SHENDURE_BOUNDS uses "reference" already, but be defensive
        ct = b.cells_per_cell_type.index[0]
        b.cells_per_cell_type = b.cells_per_cell_type.loc[[ct]].copy()
        b.cells_per_cell_type.index = [CELL_TYPE]
    else:
        b.cells_per_cell_type = b.cells_per_cell_type.loc[[CELL_TYPE]].copy()
    b.cells_per_cell_type.loc[CELL_TYPE] = int(sample["n_cells"])

    # Library NB
    b.library_model.mu_nb = float(sample["bcs_per_cre"])
    b.library_model.alpha_nb = float(sample["lib_alpha_nb"])
    b.library_model.update_alt_nb_param()

    # MOI / transfection NB. set_effective_moi rescales mu while preserving
    # dispersion shape (see core.py).
    b.set_effective_moi(float(sample["moi"]))

    # Reference activity
    b.reference_activity = float(sample["minP"])
    b.min_mpra_umi = MIN_ACT_FRAC_OF_MINP * float(sample["minP"])

    return b


# ---------------------------------------------------------------------------
# Phase: simulate
# ---------------------------------------------------------------------------

def _sample_dir(sample_id: str) -> Path:
    return SIM_ROOT / sample_id


def _count_done(pt_dir: Path, test_type: str = "mwu") -> int:
    n = 0
    for d in sorted(pt_dir.iterdir()) if pt_dir.exists() else []:
        if not d.is_dir():
            continue
        tt_dir = d / "tests" / "hs_act" / test_type
        if tt_dir.exists() and any(tt_dir.iterdir()):
            n += 1
    return n


def _run_one_sample(row: pd.Series, client) -> None:
    """Run all reps for one LHS sample on the given client.

    Each driver job runs this serially in a for-loop -- no thread pool. Across
    slices, multiple drivers run their own copies in parallel, each with its
    own dask cluster. This avoids nested-dask deadlock and shared-scheduler
    contention seen in the 2026-05-05 threaded approach.
    """
    sid = row["sample_id"]
    pt_dir = _sample_dir(sid)
    # Skip if already cached + pruned: cached_results.parquet is the
    # ground-truth marker that this sample's reps were all completed and
    # the raw sim_<hash> dirs were collapsed. _count_done can't see that
    # state because it looks for tests/hs_act/mwu which gets pruned.
    if (pt_dir / "cached_results.parquet").exists():
        print(f"SKIP   {sid} (cached)", flush=True)
        return
    pt_dir.mkdir(parents=True, exist_ok=True)
    n_cres = int(row["n_cres"])
    minP = float(row["minP"])
    max_act = float(row["activity_max_mult"]) * minP
    min_act = MIN_ACT_FRAC_OF_MINP * minP

    n_done = _count_done(pt_dir)
    head = (
        f"{sid}: ncells={row['n_cells']} ncres={n_cres} "
        f"bcs={row['bcs_per_cre']:.1f} moi={row['moi']:.2f} "
        f"alpha={row['lib_alpha_nb']:.3f} minP={minP:.3f} "
        f"maxmult={row['activity_max_mult']:.1f} (done={n_done})"
    )
    if n_done >= N_LIBRARY_REPS:
        print(f"SKIP   {head}", flush=True)
        return
    print(f"START  {head}", flush=True)

    bound = make_synthetic_bounds(row)
    n_remaining = N_LIBRARY_REPS - n_done
    hs = None
    for i in range(n_remaining):
        _, sim = scm.one_library_replicate(
            root=pt_dir,
            n_sims=N_SIMS,
            client=client,
            flatten_overtransfection=True,
            bound=bound,
            n_cres=n_cres,
            min=min_act,
            max=max_act,
            minP=minP,
            cell_type=CELL_TYPE,
        )
        sim.save()
        if hs is None:
            example = scm.scMPRA_data.from_parquet(sim.scmpradatp / "0.scmpra")
            hs = scm.make_all_by_celltype_hypotheses(
                counts=example, reference_cre="reference"
            )
        sim.add_hypothesis_set("hs_act", hs)
        sim.mwu("hs_act")
        sim.save()
        del sim
        gc.collect()

    # Cache + prune inline. Scratch is tight (pi_skr2 group quota is 89%
    # full as of 2026-05-10), so each completed sample collapses its
    # ~10-15 GB of raw sim dirs to a ~10 KB cached_results.parquet
    # before the next sample starts. Without inline pruning peak disk
    # would be n_samples x ~10 GB, which exceeds available scratch.
    if _cache_sample(pt_dir, sid, client):
        try:
            freed = _prune_sample(pt_dir)
            print(f"DONE   {sid}  pruned {freed/1e9:.1f} GB", flush=True)
        except Exception as e:
            print(f"DONE   {sid}  prune skipped: {type(e).__name__}: {e}", flush=True)
    else:
        print(f"DONE   {sid}  cache failed; not pruning", flush=True)


def _slice_samples(samples: pd.DataFrame, slice_idx: int, n_slices: int) -> pd.DataFrame:
    """Return rows whose 0-indexed position satisfies (i % n_slices) == slice_idx.

    Modulo round-robin (rather than contiguous chunking) gives each slice a
    representative spread of axis values, so per-slice runtime is roughly
    balanced even when the LHS rows happen to be sorted in a way that
    correlates with sample size.
    """
    if n_slices <= 1:
        return samples
    if slice_idx < 0 or slice_idx >= n_slices:
        raise ValueError(f"slice_idx={slice_idx} must be in [0,{n_slices})")
    mask = (np.arange(len(samples)) % n_slices) == slice_idx
    return samples.loc[mask].reset_index(drop=True)


def phase_simulate(samples: pd.DataFrame, client, slice_idx: int = 0, n_slices: int = 1):
    """Serial driver over a slice of samples.

    Each invocation runs ONE slice; cross-slice parallelism comes from
    submitting n_slices independent sbatch driver jobs (see launch.sh).
    """
    sliced = _slice_samples(samples, slice_idx, n_slices)
    print(f"=== synthetic factorial slice {slice_idx}/{n_slices}: "
          f"{len(sliced)}/{len(samples)} samples on {N_WORKERS} workers ===",
          flush=True)
    completed = 0
    failed = 0
    n_total = len(sliced)
    for _, row in sliced.iterrows():
        try:
            _run_one_sample(row, client)
            completed += 1
            print(f"[slice {slice_idx}: {completed+failed}/{n_total}] ok", flush=True)
        except Exception as e:
            failed += 1
            print(f"[slice {slice_idx}: {completed+failed}/{n_total}] FAILED "
                  f"{row['sample_id']}: {type(e).__name__}: {e}", flush=True)
            traceback.print_exc()
    print(f"slice {slice_idx}: {completed} completed, {failed} failed", flush=True)


# ---------------------------------------------------------------------------
# Phase: plot
# ---------------------------------------------------------------------------

def _power_metrics(cat: pd.DataFrame) -> dict:
    """Compute multiple power summaries from a per-sample (fc, reject_null)
    dataframe. We compute several so that downstream plots aren't trapped
    by P@FC=2 saturation -- ~40% of LHS draws sit at the 0 or 1 ceiling
    when summarized by P@FC=2 alone (2026-05-06 audit).

    Metrics:
      power_at_fc1p5 : window-averaged power at FC=1.5 (+/- 0.10).
                      Less ceiling-saturated than P@FC=2.
      power_at_fc2  : window-averaged power at FC=2.0 (+/- 0.15). Legacy.
      power_at_fc3  : window-averaged power at FC=3.0 (+/- 0.30). Sparse
                      because LHS draws fewer activities at high FC.
      power_auc_1to3: mean of power across 20 bins in FC ∈ [1.0, 3.0].
                      Continuous, never saturates.
      fc_at_p50    : smallest FC where smoothed power crosses 0.5.
                      "Minimum detectable effect"; smaller is better.
                      NaN when power never reaches 0.5.
    """
    cat = cat.dropna(subset=["fc"])
    cat = cat[cat["fc"] > 0]
    if len(cat) == 0:
        return dict.fromkeys(
            ["power_at_fc1p5", "n_obs_fc1p5",
             "power_at_fc2", "n_obs_fc2",
             "power_at_fc3", "n_obs_fc3",
             "power_auc_1to3", "fc_at_p50"], np.nan)

    out = {}
    for cut, half, key in [(1.5, 0.10, "1p5"), (2.0, 0.15, "2"), (3.0, 0.30, "3")]:
        win = cat[(cat["fc"] >= cut - half) & (cat["fc"] <= cut + half)]
        out[f"power_at_fc{key}"] = float(win["reject_null"].mean()) if len(win) > 0 else np.nan
        out[f"n_obs_fc{key}"] = int(len(win))

    # AUC: mean of binwise power across FC in [1, 3]
    bin_edges = np.linspace(1.0, 3.0, 21)
    bin_powers = []
    for lo, hi in zip(bin_edges[:-1], bin_edges[1:]):
        m = cat[(cat["fc"] >= lo) & (cat["fc"] < hi)]
        bin_powers.append(m["reject_null"].mean() if len(m) > 0 else np.nan)
    bp = np.array(bin_powers, dtype=float)
    out["power_auc_1to3"] = float(np.nanmean(bp)) if not np.all(np.isnan(bp)) else np.nan

    # FC at 50% power: rolling smooth across fc-sorted observations
    sd = cat.sort_values("fc")
    if len(sd) >= 20:
        wsize = max(20, len(sd) // 30)
        roll = sd["reject_null"].astype(float).rolling(wsize, min_periods=10, center=True).mean()
        cross = sd[(roll >= 0.5)]
        out["fc_at_p50"] = float(cross["fc"].iloc[0]) if len(cross) > 0 else np.nan
    else:
        out["fc_at_p50"] = np.nan
    return out


def _cache_sample(pt_dir: Path, sample_id: str, client,
                   test_type: str = "mwu",
                   hypothesis_set: str = "hs_act") -> bool:
    """Walk the sim_<hash> dirs of one sample, merge ground truth with
    test results, and write a single <pt_dir>/cached_results.parquet.

    Columns: sample_id, sim_id, rep_idx, reject_null, fc,
             comparison_truth, reference_truth.

    No-op if the cache already exists. Returns True on success (cache
    present afterward), False if no valid sim dirs were found.
    """
    cache_path = pt_dir / "cached_results.parquet"
    if cache_path.exists():
        return True
    if not pt_dir.exists():
        return False
    rows = []
    for d in sorted(pt_dir.iterdir()):
        if not d.is_dir():
            continue
        res_dir = d / "tests" / hypothesis_set / test_type
        if not (res_dir.exists() and any(res_dir.iterdir())):
            continue
        try:
            sim = scm.de_novo_simulation(location=pt_dir, name=d.name, client=client)
            n_sims = sim.get_state_field("n_sims")
        except Exception:
            continue
        for i in range(n_sims):
            try:
                m = sim._merge_in_ground_truth(
                    hypothesis_set_name=hypothesis_set, test_type=test_type, index=i)
            except Exception:
                continue
            m["sample_id"] = sample_id
            m["sim_id"] = d.name
            m["rep_idx"] = i
            m["comparison_truth"] = m["comparison_truth"].astype(float)
            m["reference_truth"] = m["reference_truth"].astype(float)
            m["fc"] = m["comparison_truth"] / m["reference_truth"]
            rows.append(m[["sample_id", "sim_id", "rep_idx",
                            "reject_null", "fc",
                            "comparison_truth", "reference_truth"]].copy())
    if not rows:
        return False
    pd.concat(rows, ignore_index=True).to_parquet(cache_path)
    return True


def _prune_sample(pt_dir: Path) -> int:
    """Delete everything in pt_dir except cached_results.parquet.

    Refuses to prune if no cache exists -- the cache is the only thing
    that lets _aggregate run later. Returns approximate bytes freed.
    """
    cache_path = pt_dir / "cached_results.parquet"
    if not cache_path.exists():
        raise RuntimeError(f"refusing to prune {pt_dir}: no cached_results.parquet")
    freed = 0
    for item in pt_dir.iterdir():
        if item.name == "cached_results.parquet":
            continue
        if item.is_dir():
            for root, _, files in os.walk(item):
                for f in files:
                    try:
                        freed += os.path.getsize(os.path.join(root, f))
                    except OSError:
                        pass
            shutil.rmtree(item, ignore_errors=True)
        else:
            try:
                freed += item.stat().st_size
            except OSError:
                pass
            try:
                item.unlink()
            except OSError:
                pass
    return freed


def _aggregate(samples: pd.DataFrame, client, test_type: str = "mwu") -> pd.DataFrame:
    """Walk scratch, per sample compute several power summaries, return one
    row per sample.

    Prefers <pt_dir>/cached_results.parquet (cheap read) when present;
    falls back to walking the sim_<hash> dirs for samples that haven't
    been cached + pruned yet.
    """
    rows = []
    skipped = 0
    for _, row in samples.iterrows():
        sid = row["sample_id"]
        pt_dir = _sample_dir(sid)
        if not pt_dir.exists():
            continue
        cache_path = pt_dir / "cached_results.parquet"
        if cache_path.exists():
            try:
                cat = pd.read_parquet(cache_path)
            except Exception as e:
                skipped += 1
                print(f"  skip {sid} cache: {type(e).__name__}: {e}", flush=True)
                continue
            if "reject_null" not in cat.columns or "fc" not in cat.columns:
                skipped += 1
                continue
            cat = cat[["reject_null", "fc"]].dropna(subset=["fc"])
        else:
            sim_dirs = [d for d in sorted(pt_dir.iterdir()) if d.is_dir()]
            per_sample = []
            for d in sim_dirs:
                res_dir = d / "tests" / "hs_act" / test_type
                if not (res_dir.exists() and any(res_dir.iterdir())):
                    continue
                required = [
                    d / "ground_truth.tsv.gz",
                    d / "futures.tsv.gz",
                    d / "state.parquet",
                ] + [res_dir / f"{i}_results.tsv" for i in range(N_SIMS)]
                if not all(p.exists() for p in required):
                    skipped += 1
                    continue
                try:
                    sim = scm.de_novo_simulation(location=pt_dir, name=d.name, client=client)
                    reps = sim.get_state_field("n_sims")
                    for i in range(reps):
                        mergy = sim._merge_in_ground_truth(
                            hypothesis_set_name="hs_act", test_type=test_type, index=i
                        )
                        mergy["comparison_truth"] = mergy["comparison_truth"].astype(float)
                        mergy["reference_truth"] = mergy["reference_truth"].astype(float)
                        mergy["fc"] = mergy["comparison_truth"] / mergy["reference_truth"]
                        per_sample.append(mergy[["reject_null", "fc"]].copy())
                except Exception as e:
                    skipped += 1
                    print(f"  skip {sid}/{d.name}: {type(e).__name__}: {e}", flush=True)
                    continue
            if not per_sample:
                continue
            cat = pd.concat(per_sample, ignore_index=True)
        out = row.to_dict()
        out.update(_power_metrics(cat))
        rows.append(out)
    if skipped:
        print(f"  aggregate: skipped {skipped} unloadable items", flush=True)
    return pd.DataFrame(rows)


# Memoized results of _loess_band for the marginal panels, keyed by
# (metric, axis). Cleared at the top of each phase_plot so a second
# mode in the same process cannot read the first mode's curves.
_LOESS_CACHE: dict = {}


def _loess_band(x, y, n_grid=100, frac=0.4, n_boot=200, seed=42):
    """LOESS-like smoother via locally-weighted linear regression with
    bootstrap 95% bands. Returns (xg, yhat, lo, hi)."""
    from statsmodels.nonparametric.smoothers_lowess import lowess
    xg = np.linspace(np.min(x), np.max(x), n_grid)
    yhat = lowess(y, x, frac=frac, xvals=xg, return_sorted=False)
    rng = np.random.default_rng(seed)
    boots = np.empty((n_boot, n_grid))
    for b in range(n_boot):
        idx = rng.integers(0, len(x), size=len(x))
        try:
            boots[b] = lowess(y[idx], x[idx], frac=frac, xvals=xg, return_sorted=False)
        except Exception:
            boots[b] = np.nan
    lo = np.nanpercentile(boots, 2.5, axis=0)
    hi = np.nanpercentile(boots, 97.5, axis=0)
    return xg, yhat, lo, hi


# --- Fig 5B geometry, in inches at the size the page gives the panel --------
# Three columns of seven axes: at 6.56 in a fourth column would leave each
# panel too narrow for a decade ladder underneath it, and two columns would
# make the figure taller than the page has left once Fig 5A is above it.
MARG_FIG_W = 6.56            # 0.95\linewidth of a 6.90 in text block
MARG_NCOL = 3
# The left margin is narrow because the y axis is labelled once for the whole
# grid: it holds a three-rung ladder and nothing else. It is also what puts
# this panel's first column roughly under Fig 5A's, once the 0.95\linewidth
# inclusion has centred the panel in the text block.
MARG_LEFT_IN, MARG_RIGHT_IN = 0.36, 0.06
MARG_COLGAP_IN = 0.28
MARG_PLOT_H_IN = 1.15
MARG_HEADROOM_IN = 0.20      # anchor tags sit above each panel's top spine
MARG_XFURN_IN = 0.42         # tick ladder and axis name under each row
MARG_TOP_IN = 0.12           # headroom only; the panel carries no title
MARG_BOTTOM_IN = 0.06
MARG_PLOT_W_IN = (MARG_FIG_W - MARG_LEFT_IN - MARG_RIGHT_IN
                  - (MARG_NCOL - 1) * MARG_COLGAP_IN) / MARG_NCOL
MARG_ROW_IN = MARG_HEADROOM_IN + MARG_PLOT_H_IN + MARG_XFURN_IN
MARG_NROW = int(np.ceil(len(AXIS_NAMES) / MARG_NCOL))
MARG_FIG_H = MARG_TOP_IN + MARG_NROW * MARG_ROW_IN + MARG_BOTTOM_IN

# Two anchor tags closer than this in axes fraction go on separate lines.
ANCHOR_TAG_SEP = 0.06


def _marg_rect(i: int) -> list:
    """Figure-fraction rect of the i-th marginal panel, row-major."""
    r, c = divmod(i, MARG_NCOL)
    x = MARG_LEFT_IN + c * (MARG_PLOT_W_IN + MARG_COLGAP_IN)
    top = MARG_TOP_IN + r * MARG_ROW_IN + MARG_HEADROOM_IN
    return [x / MARG_FIG_W, (MARG_FIG_H - top - MARG_PLOT_H_IN) / MARG_FIG_H,
            MARG_PLOT_W_IN / MARG_FIG_W, MARG_PLOT_H_IN / MARG_FIG_H]


def _style_marg_axes(ax, show_yticks: bool):
    """Hairline y grid, no box, ticks without tick marks."""
    ax.grid(axis="y", color=HAIR, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(HAIR)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(length=0, colors=INK, pad=2)
    if not show_yticks:
        ax.tick_params(labelleft=False)


def _draw_anchor_rules(ax, axis: str, anchors: list, data_lo: float,
                       data_hi: float):
    """Dashed rule and one-letter tag where each published design sits.

    Call after the x limits are final: the tags are laid out in axes fraction
    so that two designs at nearby values end up on separate lines rather than
    on top of one another, and that fraction depends on the limits.

    An anchor outside the swept range is marked at the edge rather than
    plotted at its value: extending the axis to reach it would imply the
    smoother extends there too, and would squash the region that carries the
    data. Yin et al. on lib_alpha_nb is the case in hand -- see the note on
    that axis in EMPIRICAL.
    """
    placed = []
    for name in anchors:
        coord = EMPIRICAL.get(name)
        if coord is None or coord.get(axis) is None:
            continue
        xv = float(coord[axis])
        color = EMPIRICAL_COLORS.get(name, INK)
        tag = ANCHOR_DISPLAY[name][0]
        if xv < data_lo or xv > data_hi:
            off_low = xv < data_lo
            edge = data_lo if off_low else data_hi
            ax.plot([edge], [1.015], marker="<" if off_low else ">",
                    color=color, markersize=3, clip_on=False,
                    transform=ax.get_xaxis_transform())
            ax.text(edge, 1.015, f"  {tag} off scale" if off_low
                    else f"{tag} off scale  ",
                    fontsize=6, color=color, va="baseline",
                    ha="left" if off_low else "right",
                    transform=ax.get_xaxis_transform())
            continue
        frac = ax.transAxes.inverted().transform(
            ax.transData.transform((xv, 0)))[0]
        row = 0
        while any(abs(frac - f) < ANCHOR_TAG_SEP and r == row
                  for f, r in placed):
            row += 1
        ax.axvline(xv, color=color, linestyle=(0, (3, 2)), lw=0.8, zorder=3)
        ax.text(xv, 1.015 + 0.085 * row, tag, fontsize=6, fontweight="bold",
                ha="center", va="baseline", color=color,
                transform=ax.get_xaxis_transform())
        placed.append((frac, row))
    return placed


def _draw_anchor_key(fig, rect, anchors: list):
    """Name the dashed rules, in the grid slot the seven axes leave empty."""
    x = rect[0] + 0.01
    y = rect[1] + rect[3] - 0.02
    dy = 0.115 / MARG_FIG_H
    fig.text(x, y, "published designs", fontsize=6.5, color=INK,
             ha="left", va="top")
    for k, name in enumerate(anchors):
        color = EMPIRICAL_COLORS.get(name, INK)
        tag, full = ANCHOR_DISPLAY[name]
        yk = y - (k + 1.35) * dy
        fig.add_artist(plt.Line2D([x, x + 0.10 / MARG_FIG_W], [yk] * 2,
                                  color=color, lw=0.8, ls=(0, (3, 2))))
        fig.text(x + 0.13 / MARG_FIG_W, yk, f"{tag}   {full}", fontsize=6,
                 color=MUTED, ha="left", va="center")


def _plot_marginals_for_metric(df: pd.DataFrame, metric: str, ylim: "tuple[float,float]",
                                ylabel: str, title: str, out_path: Path,
                                subtitle: str = "",
                                captioned: bool = True,
                                hline: "float | None" = None,
                                invert_y: bool = False,
                                force_linear_x: bool = False,
                                anchors: "list[str] | None" = None):
    """Render the 7-axis marginals figure for a chosen power metric column.

    One panel per design axis: the sampled designs as a light scatter, the
    LOESS smoother and its bootstrap band on top, and a dashed rule where each
    published design sits. The smoother is what the figure is about, so it is
    the only thing drawn at full contrast; the scatter is the residual spread
    the other six axes leave behind, and reads as texture behind it.

    For axes that have log_scale=True in the mode's bounds, the smoother runs
    in log10-space (otherwise it is dominated by the dense end of the
    distribution) but the axis is drawn on a log scale, so the ladder shows
    raw values rather than their logarithms.

    force_linear_x=True overrides the per-axis log scale and sets all panels
    to linear x. The smoother *still runs in log space* for those axes, and
    the curve is drawn back on the linear scale. Useful when reviewers want a
    "true scale" view alongside the log view.

    The y scale is shared across panels and labelled once, on the left.
    Repeating it seven times spends a seventh of the figure restating that
    every panel plots the same quantity.
    """
    from matplotlib import patheffects as pe
    from matplotlib.ticker import (FixedFormatter, FixedLocator, MaxNLocator,
                                   NullLocator)

    if anchors is None:
        anchors = DEFAULT_ANCHORS
    y_all = df[metric].values.astype(float)
    if ylim is None:
        # A shared range even when the caller has no fixed one: per-panel
        # autoscaling would draw the same curve at seven different heights.
        finite = y_all[np.isfinite(y_all)]
        assert len(finite), f"{metric}: no finite values in {len(y_all)} samples"
        ylim = (float(np.min(finite)), float(np.nanpercentile(finite, 99)))

    fig = plt.figure(figsize=(MARG_FIG_W, MARG_FIG_H))
    for ax_i, axis in enumerate(AXIS_NAMES):
        ax = fig.add_axes(_marg_rect(ax_i))
        x = df[axis].values.astype(float)
        y = y_all
        log_scale = active_bounds()[axis][2]
        display_log = log_scale and not force_linear_x
        x_smooth = np.log10(x) if log_scale else x

        ax.scatter(x, y, s=1.2, alpha=0.32, color=DOT, linewidths=0,
                   zorder=1, rasterized=True)

        # Limits before the anchors, which need them to lay their tags out,
        # and computed from the data rather than read back off the axes: at
        # this point the axis is still linear, so its left edge can be ~0 and
        # a later set_xscale("log") would open a multi-decade empty region.
        x_finite = x[np.isfinite(x)]
        assert len(x_finite), f"{axis}: no finite values in {len(x)} samples"
        data_lo, data_hi = float(x_finite.min()), float(x_finite.max())
        if display_log:
            ax.set_xscale("log")
            pad = 10 ** (0.03 * (np.log10(data_hi) - np.log10(data_lo)))
            ax.set_xlim(data_lo / pad, data_hi * pad)
            decades = _decade_ticks(np.log10(data_lo), np.log10(data_hi))
            ax.xaxis.set_major_locator(FixedLocator([10.0 ** k for k in decades]))
            ax.xaxis.set_major_formatter(
                FixedFormatter([_fmt_decade(k) for k in decades]))
            ax.xaxis.set_minor_locator(NullLocator())
        else:
            pad = 0.03 * (data_hi - data_lo)
            ax.set_xlim(data_lo - pad, data_hi + pad)
        ax.set_ylim(*ylim)
        if ylim == (0, 1):
            # Three rungs on a 1.15 in panel. A six-rung 0.0-1.0 ladder is
            # more axis furniture than the curve it is there to read.
            ax.set_yticks([0.0, 0.5, 1.0], ["0", "0.5", "1"])
        else:
            ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
        if invert_y:
            ax.invert_yaxis()

        valid = np.isfinite(y) & np.isfinite(x_smooth)
        if valid.sum() >= 20:
            try:
                # The smoother and its 200 bootstraps depend only on the
                # (metric, axis) pair: anchor overlays are drawn on top of it
                # and force_linear_x only changes the display scale, since
                # x_smooth keys off log_scale rather than display_log. Without
                # this cache the same curve is recomputed once per variant --
                # 210 fits where 35 suffice, about an hour of the run.
                # _loess_band is seeded, so caching cannot change the output.
                ck = (metric, axis)
                if ck not in _LOESS_CACHE:
                    _LOESS_CACHE[ck] = _loess_band(x_smooth[valid], y[valid])
                xg, yhat, lo, hi = _LOESS_CACHE[ck]
                xg_disp = (10.0 ** xg) if log_scale else xg
                ax.fill_between(xg_disp, lo, hi, color=INK, alpha=0.16, lw=0,
                                zorder=4)
                # Cased in white so the curve stays legible where it crosses
                # the densest part of the scatter.
                ax.plot(xg_disp, yhat, color=INK, lw=1.3, zorder=5,
                        solid_capstyle="round",
                        path_effects=[pe.Stroke(linewidth=2.8,
                                                foreground="white"),
                                      pe.Normal()])
            except Exception as e:
                print(f"  smoother failed for {axis}/{metric}: {e}", flush=True)

        _draw_anchor_rules(ax, axis, anchors, data_lo, data_hi)
        if hline is not None:
            ax.axhline(hline, color=MUTED, lw=0.6, ls=(0, (3, 2)), zorder=2)
        ax.set_xlabel(axis_label(axis), color=INK, labelpad=2)
        _style_marg_axes(ax, show_yticks=(ax_i % MARG_NCOL == 0))

    _draw_anchor_key(fig, _marg_rect(len(AXIS_NAMES)), anchors)

    # captioned=False for the manuscript panel: its figure legend says what
    # this is, and a title repeating the legend is wasted space. The
    # exploratory variants are read on their own and keep theirs.
    if captioned:
        fig.text(MARG_LEFT_IN / MARG_FIG_W, 1 - 0.115 / MARG_FIG_H,
                 title + ("  [linear x]" if force_linear_x else ""),
                 ha="left", va="center", fontsize=7.5, color=INK)
        if subtitle:
            fig.text(MARG_LEFT_IN / MARG_FIG_W, 1 - 0.255 / MARG_FIG_H,
                     subtitle, ha="left", va="center", fontsize=6.5,
                     color=MUTED)
    fig.text(0.075 / MARG_FIG_W, 0.5, ylabel, rotation=90, ha="left",
             va="center", fontsize=7.5, color=INK)
    _assert_labels_fit(fig, f"marginals/{metric}")
    # No bbox_inches="tight": the panel is placed at a fixed width on the page
    # and the margins above are what put its plot areas where they are.
    fig.savefig(out_path, format="svg", dpi=600)
    plt.close(fig)
    print(f"Saved: {out_path} ({MARG_FIG_W:.2f} x {MARG_FIG_H:.2f} in)",
          flush=True)


# One marginals SVG per metric. The legacy filename "marginals.svg" remains
# the P@FC=2 plot for backwards compat with existing references.
# Titles name the quantity and the reading; the y label names the metric, once.
METRIC_SPECS = [
    # (metric, ylim, ylabel, hline, invert_y, fname_suffix, title)
    ("power_at_fc1p5", (0, 1), "power at 1.5x fold change", 0.8, False,
     "_p15", "Power at 1.5x against each design axis"),
    ("power_at_fc2",   (0, 1), "power at 2x fold change",   0.8, False,
     "", "Power at 2x against each design axis"),
    ("power_at_fc3",   (0, 1), "power at 3x fold change",   0.8, False,
     "_p3", "Power at 3x against each design axis"),
    ("power_auc_1to3", (0, 1), "mean power, fold change 1-3x", None, False,
     "_auc", "Mean power against each design axis"),
    ("fc_at_p50",      None,   "fold change at 50% power", None, True,
     "_fc50", "Minimum detectable fold change against each design axis"),
]

# published: the three assays the manuscript reports on, which is the variant
# its figures use. The takeshi variants predate Yin et al. being available as
# an anchor and are kept for the unpublished-data comparison.
ANCHOR_VARIANTS = [
    ("_published",    PUBLISHED_ANCHORS),
    ("_no_takeshi",   DEFAULT_ANCHORS),
    ("_with_takeshi", ALL_ANCHORS),
]

# The manuscript's Fig 5B, and the two Fig 5A / Fig S10 heatmap grids.
MANUSCRIPT_METRIC = "power_auc_1to3"
MANUSCRIPT_ANCHORS = "_published"


def _marginals_subtitle(df: pd.DataFrame) -> str:
    return (f"{len(df):,} designs sampled at random over all seven axes; "
            f"LOESS fit with 95% bootstrap band")


def _plot_manuscript_marginals(df: pd.DataFrame, suffix: str):
    """Fig 5B only: the AUC metric with the three published designs overlaid."""
    metric, ylim, ylabel, hline, invert, fsuf, title = next(
        s for s in METRIC_SPECS if s[0] == MANUSCRIPT_METRIC)
    anchors = dict(ANCHOR_VARIANTS)[MANUSCRIPT_ANCHORS]
    _plot_marginals_for_metric(
        df, metric=metric, ylim=ylim, ylabel=ylabel, title=title,
        subtitle=_marginals_subtitle(df), captioned=False,
        hline=hline, invert_y=invert, anchors=anchors,
        out_path=OUT / f"marginals{fsuf}{suffix}{MANUSCRIPT_ANCHORS}.svg")


def phase_manuscript(df: pd.DataFrame, mode: str):
    """Just the three SVGs the manuscript places, from the cached sweep.

    phase_plot renders thirty marginal variants on top of these; this is the
    entry point for iterating on the figures themselves.
    """
    _LOESS_CACHE.clear()
    suffix = "" if mode == "full" else f"_{mode}"
    _plot_manuscript_marginals(df, suffix)
    _pairwise_heatmaps(df, suffix=suffix, title_suffix="")
    _pairwise_heatmaps_all(df, suffix=suffix)


def phase_plot(samples_with_power: pd.DataFrame, mode: str):
    _LOESS_CACHE.clear()
    df = samples_with_power.copy()
    suffix = "" if mode == "full" else f"_{mode}"
    df.to_parquet(OUT / f"samples_power{suffix}.parquet")
    print(f"plotting {len(df)} samples; metrics: power_at_fc1p5/2/3, "
          f"power_auc_1to3, fc_at_p50", flush=True)

    for metric, ylim, ylabel, hline, invert, fsuf, title in METRIC_SPECS:
        for asuf, alist in ANCHOR_VARIANTS:
            for linx_suf, linx in [("", False), ("_linearx", True)]:
                out_path = OUT / f"marginals{fsuf}{suffix}{asuf}{linx_suf}.svg"
                _plot_marginals_for_metric(
                    df, metric=metric, ylim=ylim, ylabel=ylabel, title=title,
                    subtitle=_marginals_subtitle(df),
                    out_path=out_path, hline=hline, invert_y=invert,
                    force_linear_x=linx, anchors=alist)

    # Pairwise heatmaps for top axes by smoother range (= biggest LOESS swing).
    # One LHS draw in, one set of heatmaps out. Fusing separate sweeps is not
    # an option here: sample sets welded together are not a Latin hypercube of
    # the combined box, and the earlier full+topup+topup3 patchwork induced
    # cross-axis correlation because its top-ups extended two axes at once.
    # The "union" modes below draw once over the enclosing box instead.
    _pairwise_heatmaps(df, suffix=suffix, title_suffix="")
    _pairwise_heatmaps_all(df, suffix=suffix)


def _heat_grid(df, a, b, nb=8):
    """Mean power in an nb x nb grid over axes (a, b), and the bin edges.

    Bins on whichever scale the mode drew the axis on; binning a log-drawn
    axis linearly piles most of the samples into the first bin and reads as a
    flat response. Edges come back in that same (log10, where log-drawn)
    space, which is what the panel is drawn in.
    """
    xa = df[a].values.astype(float)
    xb = df[b].values.astype(float)
    if active_bounds()[a][2]:
        xa = np.log10(xa)
    if active_bounds()[b][2]:
        xb = np.log10(xb)
    # Mask non-finite power before the weighted histogram: np.histogram2d sums
    # weights, so one NaN weight poisons its whole bin to NaN (imshow then
    # renders it transparent). power_auc_1to3 is NaN-free, but mask
    # defensively -- the earlier power_at_fc2 weighting was ~15% NaN and
    # blanked the entire figure.
    pw = df["power_auc_1to3"].values.astype(float)
    m = np.isfinite(xa) & np.isfinite(xb) & np.isfinite(pw)
    bins_a = np.linspace(xa[m].min(), xa[m].max(), nb + 1)
    bins_b = np.linspace(xb[m].min(), xb[m].max(), nb + 1)
    H_sum, _, _ = np.histogram2d(xa[m], xb[m], bins=[bins_a, bins_b], weights=pw[m])
    H_n, _, _ = np.histogram2d(xa[m], xb[m], bins=[bins_a, bins_b])
    assert H_n.sum() == m.sum(), (
        f"{a} x {b}: {int(H_n.sum())} of {int(m.sum())} samples landed in a "
        f"bin; the rest fell outside the edges")
    with np.errstate(invalid="ignore"):
        H = np.where(H_n > 0, H_sum / H_n, np.nan)
    filled = np.isfinite(H)
    assert filled.any(), f"{a} x {b}: every one of {H.size} bins is empty"
    assert (H[filled] >= 0).all() and (H[filled] <= 1).all(), (
        f"{a} x {b}: mean power outside [0, 1], range "
        f"{np.nanmin(H)} to {np.nanmax(H)}")
    return H, bins_a, bins_b, int(H.size - filled.sum())


def _heat_panel(ax, df, a, b, nb=8, xlabel=True, ylabel=True, labelsize=None):
    """One (a, b) power surface. Returns the image handle for a colourbar.

    Drawn in log10 units where the axis was drawn there, but the ladder is
    labelled with raw values: "log10 cells" makes a reader do the arithmetic
    the axis could have done for them.
    """
    from matplotlib.ticker import FixedFormatter, FixedLocator

    H, bins_a, bins_b, n_empty = _heat_grid(df, a, b, nb=nb)
    im = ax.imshow(H.T, origin="lower", aspect="auto", cmap=POWER_CMAP,
                   vmin=0, vmax=1,
                   extent=[bins_a[0], bins_a[-1], bins_b[0], bins_b[-1]])
    for edges, mpl_axis, is_log in (
            (bins_a, ax.xaxis, active_bounds()[a][2]),
            (bins_b, ax.yaxis, active_bounds()[b][2])):
        if is_log:
            # Three decades at most: these panels are barely an inch wide, and
            # a four-rung ladder runs its labels into one another.
            decades = _decade_ticks(edges[0], edges[-1], max_ticks=3)
            mpl_axis.set_major_locator(FixedLocator(decades))
            mpl_axis.set_major_formatter(
                FixedFormatter([_fmt_decade(k) for k in decades]))
    ax.tick_params(length=0, colors=INK, pad=2)
    # A label on a tick hard against either end of the axis would hang into
    # the neighbouring panel, so those two are aligned inwards. Only those
    # two: nudging a tick that sits well inside the axis pulls it towards its
    # neighbour instead, which is the collision it was meant to prevent.
    x_lo, x_hi = ax.get_xlim()
    for tick, label in zip(ax.get_xticks(), ax.get_xticklabels()):
        frac = (tick - x_lo) / (x_hi - x_lo)
        if frac < 0.1:
            label.set_ha("left")
        elif frac > 0.9:
            label.set_ha("right")
    for spine in ax.spines.values():
        spine.set_color(HAIR)
        spine.set_linewidth(0.6)
    if xlabel:
        ax.set_xlabel(axis_label_short(a), color=INK, labelpad=2,
                      fontsize=labelsize)
    if ylabel:
        ax.set_ylabel(axis_label_short(b), color=INK, labelpad=2,
                      fontsize=labelsize)
    return im, n_empty


STRIP_W_IN, STRIP_H_IN, STRIP_RIGHT_IN = 1.30, 0.05, 0.06
# Two lines: on one it is wider than the strip it names, and widening the
# strip to match would eat into the room the title has.
STRIP_LABEL = "mean power,\nfold change 1-3x"


def _power_strip(fig, im, fig_w, fig_h, top_in=0.26):
    """The 0-1 power scale as a slim horizontal strip in the top margin.

    One strip for the whole figure, at its top right. Every panel is on the
    same fixed 0-1 scale, so a bar beside each of them is the same key drawn
    once per panel, and in the margin it costs no panel width at all.

    Returns the left edge of the strip in inches, which is where the header
    text on the other side of the margin has to stop.
    """
    x0_in = fig_w - STRIP_RIGHT_IN - STRIP_W_IN
    cax = fig.add_axes([x0_in / fig_w, (fig_h - top_in - STRIP_H_IN) / fig_h,
                        STRIP_W_IN / fig_w, STRIP_H_IN / fig_h])
    cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
    cbar.set_ticks([0.0, 0.5, 1.0])
    cbar.set_ticklabels(["0", "0.5", "1"])
    cbar.outline.set_visible(False)
    cax.tick_params(length=0, labelsize=6, colors=MUTED, pad=1.5)
    label = fig.text(x0_in / fig_w, 1 - (top_in - 0.03) / fig_h, STRIP_LABEL,
                     ha="left", va="bottom", fontsize=6.5, color=MUTED,
                     linespacing=1.35)
    _assert_text_fits(fig, [label], fig_w, "power strip label")
    return x0_in


# --- Fig S10 geometry, in inches at the size the page gives the figure ------
# Twenty-one pairs laid out as the strict lower triangle of a 7-axis matrix.
# A flat grid of 21 would repeat every axis label five or six times; in the
# triangle each axis is named once, on the outer edge, which is what buys the
# panels enough width to be worth looking at.
ALL_FIG_W = 6.90             # \linewidth of the text block
ALL_LEFT_IN, ALL_RIGHT_IN = 0.56, 0.06
# Top margin holds the power strip and its two-line label, nothing else.
ALL_TOP_IN, ALL_BOTTOM_IN = 0.38, 0.50
# Panels butt up against one another without this and the triangle reads as
# one continuous field rather than as a grid of separate surfaces.
ALL_GAP_IN = 0.05
# Rows are shorter than columns are wide. Six square rows plus this figure's
# own caption is more than a page holds, and nothing about the cells needs to
# be square -- their two axes are different quantities.
ALL_ROW_FRAC = 0.84
ALL_LABEL_PT = 6.5           # a 7.5 pt "basal expression" is taller than a row


def _pairwise_heatmaps_all(df, suffix: str):
    """Every axis pair, for the supplement.

    The main figure shows the three slices among the highest-swing axes; this
    is the exhaustive version, so a reader can check that the additive
    decomposition is not hiding structure in a pair nobody chose to look at.

    Row i against column j for j < i, so the pair at a given cell is read off
    the axis names on the bottom row and the left column. Panels share their
    row's and column's axis, and only the outer edge carries furniture.
    """
    n = len(AXIS_NAMES)
    pitch = (ALL_FIG_W - ALL_LEFT_IN - ALL_RIGHT_IN) / (n - 1)
    row_pitch = ALL_ROW_FRAC * pitch
    cell_w, cell_h = pitch - ALL_GAP_IN, row_pitch - ALL_GAP_IN
    fig_h = ALL_TOP_IN + (n - 1) * row_pitch + ALL_BOTTOM_IN
    fig = plt.figure(figsize=(ALL_FIG_W, fig_h))

    im, n_pairs, empty = None, 0, 0
    for i in range(1, n):          # row axis, drawn up the y
        for j in range(i):         # column axis, drawn along the x
            rect = [(ALL_LEFT_IN + j * pitch) / ALL_FIG_W,
                    (fig_h - ALL_TOP_IN - i * row_pitch + ALL_GAP_IN) / fig_h,
                    cell_w / ALL_FIG_W, cell_h / fig_h]
            ax = fig.add_axes(rect)
            im, n_empty = _heat_panel(ax, df, AXIS_NAMES[j], AXIS_NAMES[i],
                                      xlabel=(i == n - 1), ylabel=(j == 0),
                                      labelsize=ALL_LABEL_PT)
            empty += n_empty
            n_pairs += 1
            if i != n - 1:
                ax.tick_params(labelbottom=False)
            if j != 0:
                ax.tick_params(labelleft=False)
    assert n_pairs == n * (n - 1) // 2, (
        f"drew {n_pairs} panels for {n} axes, expected {n * (n - 1) // 2}")

    # See the note in _pairwise_heatmaps: caption's job, not the figure's.
    _power_strip(fig, im, ALL_FIG_W, fig_h)
    _assert_labels_fit(fig, "Fig S10")
    out = OUT / f"pairwise_heatmaps_all{suffix}.svg"
    fig.savefig(out, format="svg")
    plt.close(fig)
    print(f"Saved: {out} ({ALL_FIG_W:.2f} x {fig_h:.2f} in, {n_pairs} pairs, "
          f"{empty} empty bins)", flush=True)


# --- Fig 5A geometry, in inches at the size the page gives the panel --------
# Three pairs in one row. Each panel carries its own pair of axes, so each
# needs its own left gutter; the right-hand edge of the last one is the right
# edge of the text block.
PAIR_FIG_W = 6.90            # \textwidth of the text block
PAIR_GUTTER_IN, PAIR_RIGHT_IN = 0.47, 0.06
# Enough for the power strip, the only furniture in the top margin.
PAIR_TOP_IN, PAIR_BOTTOM_IN = 0.38, 0.44


def _pairwise_heatmaps(df, suffix: str, title_suffix: str = ""):
    ranges = {}
    for axis in AXIS_NAMES:
        x = df[axis].values.astype(float)
        y = df["power_auc_1to3"].values.astype(float)
        log_scale = active_bounds()[axis][2]
        x_plot = np.log10(x) if log_scale else x
        try:
            _, yhat, _, _ = _loess_band(x_plot, y, n_boot=20)
            ranges[axis] = float(np.nanmax(yhat) - np.nanmin(yhat))
        except Exception:
            ranges[axis] = 0.0
    top3 = sorted(ranges.items(), key=lambda kv: -kv[1])[:3]
    print(f"Top axes by marginal swing ({suffix or 'full'}): "
          + ", ".join(f"{a} {r:.3f}" for a, r in top3), flush=True)
    # Every pair among the three highest-swing axes, and nothing else. The
    # selection rule is stated in full by that sentence, which is what makes
    # the panel answerable to a reader; a hand-picked fourth slice would not
    # be. Fig S9 draws all 21 pairs for anyone who wants the rest.
    pairs = [(top3[0][0], top3[1][0]), (top3[0][0], top3[2][0]),
             (top3[1][0], top3[2][0])]

    n = len(pairs)
    cell = (PAIR_FIG_W - n * PAIR_GUTTER_IN - PAIR_RIGHT_IN) / n
    fig_h = PAIR_TOP_IN + cell + PAIR_BOTTOM_IN
    fig = plt.figure(figsize=(PAIR_FIG_W, fig_h))

    im, empty = None, 0
    for k, (a, b) in enumerate(pairs):
        rect = [(PAIR_GUTTER_IN + k * (cell + PAIR_GUTTER_IN)) / PAIR_FIG_W,
                PAIR_BOTTOM_IN / fig_h, cell / PAIR_FIG_W, cell / fig_h]
        ax = fig.add_axes(rect)
        im, n_empty = _heat_panel(ax, df, a, b)
        empty += n_empty

    # No title or explanatory subtitle: this is a manuscript panel and the
    # caption carries what it shows. The top margin stays as it is because the
    # power strip lives in it.
    _power_strip(fig, im, PAIR_FIG_W, fig_h)
    _assert_labels_fit(fig, "Fig 5A")
    out_pair = OUT / f"pairwise_heatmaps{suffix}.svg"
    fig.savefig(out_pair, format="svg")
    plt.close(fig)
    print(f"Saved: {out_pair} ({PAIR_FIG_W:.2f} x {fig_h:.2f} in, "
          f"{n} pairs, {empty} empty bins)", flush=True)


# ---------------------------------------------------------------------------
# Cluster setup
# ---------------------------------------------------------------------------

def _make_slurm_client():
    cluster = SLURMCluster(
        cores=1,
        memory=WORKER_MEM,
        processes=1,
        job_script_prologue=[f"export PYTHONPATH={_repo_root}:$PYTHONPATH"],
        job_extra_directives=[
            "-p priority",
            "--account=prio_skr2",
            "--job-name=synfac_worker",
            "--time=8:00:00",
            "--output=worker_%j.out",
        ],
    )
    # Use adapt(minimum, maximum) instead of scale(jobs=N) so worker sbatches
    # are issued gradually as work appears, rather than burst-submitting all
    # N at once. Burst submission of 10 slices x 50 workers triggered the
    # YCRC per-hour sbatch rate limit on 2026-05-05, leaving 6 of 10 slices
    # with zero workers ("No valid workers found" failures). Adaptive scaling
    # spreads submissions and avoids the cliff.
    cluster.adapt(minimum=1, maximum=N_WORKERS, interval="30s")
    client = Client(cluster, timeout="600s", heartbeat_interval="20s")
    print(f"Dask dashboard: {client.dashboard_link}", flush=True)
    return cluster, client


def _make_local_client():
    import psutil
    mem_mb = os.environ.get("SLURM_MEM_PER_NODE")
    if mem_mb:
        mem_limit = f"{int(int(mem_mb) * 0.9)}MB"
    else:
        mem_limit = f"{int(psutil.virtual_memory().total * 0.9 / 1e9)}GB"
    cluster = LocalCluster(
        n_workers=1, threads_per_worker=2, processes=False, memory_limit=mem_limit
    )
    return cluster, Client(cluster)


# ---------------------------------------------------------------------------
# Sample-table I/O
# ---------------------------------------------------------------------------


def get_samples() -> pd.DataFrame:
    """Load LHS samples table; create on first call so simulate/plot agree.

    For mode=="topup", draws from TOPUP_AXIS_BOUNDS with a different seed and
    prefixes sample_ids with "t" so they cannot collide with the original
    sweep's "s"-prefixed ids in any combined plot or aggregate.
    """
    if SAMPLES_PATH.exists():
        return pd.read_parquet(SAMPLES_PATH)
    if CURRENT_MODE == "topup":
        df = draw_lhs(N_LHS, seed=20260506, bounds=TOPUP_AXIS_BOUNDS,
                      sample_id_prefix="t")
    elif CURRENT_MODE == "topup3":
        # 'u' prefix keeps sample_ids disjoint from 's' (full) and 't' (topup).
        df = draw_lhs(N_LHS, seed=20260507, bounds=TOPUP3_AXIS_BOUNDS,
                      sample_id_prefix="u")
    elif CURRENT_MODE == "union":
        # 'v' prefix keeps sample_ids disjoint from prior sweeps.
        df = draw_lhs(N_LHS, seed=20260511, bounds=UNION_AXIS_BOUNDS,
                      sample_id_prefix="v")
    else:
        df = draw_lhs(N_LHS, seed=SEED)
    df.to_parquet(SAMPLES_PATH)
    print(f"Wrote {SAMPLES_PATH}: {len(df)} samples", flush=True)
    return df


def write_anchor(mode: str):
    meta = dict(
        mode=mode,
        n_lhs=N_LHS,
        n_library_reps=N_LIBRARY_REPS,
        n_sims=N_SIMS,
        n_workers=N_WORKERS,
        worker_mem=WORKER_MEM,
        seed=SEED,
        axis_bounds={k: list(v) for k, v in AXIS_BOUNDS.items()},
        cell_type=CELL_TYPE,
        sim_root=str(SIM_ROOT),
        empirical=EMPIRICAL,
        sim_date=SIM_DATE,
    )
    out_path = OUT / (f"anchor_{mode}.json" if mode != "full" else "anchor.json")
    with open(out_path, "w") as f:
        json.dump(meta, f, indent=2, default=str)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    phase = sys.argv[1] if len(sys.argv) > 1 else "all"
    mode = sys.argv[2] if len(sys.argv) > 2 else "full"
    _apply_mode(mode)
    _apply_mode_paths(mode)
    write_anchor(mode)

    # Parse optional slice/n_slices args (positional, after mode)
    slice_idx = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    n_slices_arg = int(sys.argv[4]) if len(sys.argv) > 4 else N_SLICES

    if phase in ("replot", "manuscript"):
        suffix = "" if mode == "full" else f"_{mode}"
        cache_path = OUT / f"samples_power{suffix}.parquet"
        if not cache_path.exists():
            print(f"{phase}: cache not found at {cache_path}; "
                  f"run `plot` first to aggregate from scratch.", flush=True)
            sys.exit(1)
        print(f"\n{'='*60}\nPhase: {phase} ({mode}) from {cache_path.name}"
              f"\n{'='*60}", flush=True)
        df_pow = pd.read_parquet(cache_path)
        if phase == "manuscript":
            phase_manuscript(df_pow, mode)
        else:
            phase_plot(df_pow, mode)
        print("\nDone.", flush=True)
        sys.exit(0)

    if phase == "samples":
        # Materialize the samples table and exit. Used by launch.sh to ensure
        # the table exists before parallel sbatch slice jobs race on it.
        get_samples()
        print("samples table ready", flush=True)
        sys.exit(0)

    if phase == "cache_prune":
        # Backfill cache + prune for samples in this mode's table. For each
        # sample with valid sim dirs but no cached_results.parquet, build
        # the cache then rm -rf the raw sim_<hash> contents. Existing
        # caches and missing samples are skipped silently.
        # Optional slice args (positional 3,4) partition the samples table
        # the same way phase_simulate does, so cache_prune can be sbatched
        # in parallel across slices. Default is one slice = all samples,
        # since cache_prune is fast and usually not worth slicing.
        cp_slice_idx = int(sys.argv[3]) if len(sys.argv) > 3 else 0
        cp_n_slices = int(sys.argv[4]) if len(sys.argv) > 4 else 1
        samples = get_samples()
        sliced = _slice_samples(samples, cp_slice_idx, cp_n_slices)
        print(f"\n{'='*60}\nPhase: cache_prune ({mode}, slice "
              f"{cp_slice_idx}/{cp_n_slices}: {len(sliced)}/{len(samples)})"
              f"\n{'='*60}", flush=True)
        cluster, client = _make_local_client()
        n_cached = n_pruned = n_skipped = 0
        total_freed = 0
        try:
            for _, row in sliced.iterrows():
                sid = row["sample_id"]
                pt_dir = _sample_dir(sid)
                if not pt_dir.exists():
                    n_skipped += 1
                    continue
                cache_path = pt_dir / "cached_results.parquet"
                already_cached = cache_path.exists()
                ok = _cache_sample(pt_dir, sid, client)
                if not ok:
                    n_skipped += 1
                    print(f"  {sid}: no valid sim dirs; skipping", flush=True)
                    continue
                if not already_cached:
                    n_cached += 1
                try:
                    freed = _prune_sample(pt_dir)
                    total_freed += freed
                    n_pruned += 1
                    print(f"  {sid}: {'cached + ' if not already_cached else ''}"
                          f"pruned {freed/1e9:.1f} GB", flush=True)
                except Exception as e:
                    print(f"  {sid}: prune failed: {type(e).__name__}: {e}",
                          flush=True)
        finally:
            client.close(); cluster.close()
        print(f"\nSummary: cached={n_cached} pruned={n_pruned} "
              f"skipped={n_skipped} freed={total_freed/1e12:.2f} TB",
              flush=True)
        sys.exit(0)

    samples = get_samples()

    if phase in ("simulate", "all"):
        print(f"\n{'='*60}\nPhase: simulate ({mode}, slice {slice_idx}/{n_slices_arg})\n{'='*60}",
              flush=True)
        cluster, client = _make_slurm_client()
        try:
            phase_simulate(samples, client, slice_idx=slice_idx, n_slices=n_slices_arg)
        finally:
            client.close(); cluster.close()

    if phase in ("plot", "all"):
        print(f"\n{'='*60}\nPhase: plot ({mode})\n{'='*60}", flush=True)
        cluster, client = _make_local_client()
        try:
            df_pow = _aggregate(samples, client, test_type="mwu")
            phase_plot(df_pow, mode)
        finally:
            client.close(); cluster.close()

    print("\nDone.", flush=True)
