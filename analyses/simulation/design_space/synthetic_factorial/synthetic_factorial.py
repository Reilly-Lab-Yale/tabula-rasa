#!/usr/bin/env python3
"""
Synthetic factorial sweep over the design-space encompassing shendure + cohen.

Latin Hypercube sample over 7 axes, simulate at each sample, run MWU on the
activity hypothesis set, then plot smoothed marginals (plus pairwise heatmaps
and the shendure/cohen empirical overlay).

The 1D sweep (analyses/simulation/design_space/shendure_1d) anchored at
shendure-Pluripotent and varied one axis at a time. This script removes the
anchor: each LHS sample is a fresh 7-tuple of axis values, with all other
axes at the LHS-drawn value rather than a fixed reference. The empirical
datasets become predictions of the meta-surface rather than starting points.

Axes (log-uniform unless noted):
    n_cells              (500 - 50000)   -- shendure max ~8.2k, cohen max ~18.6k
    n_cres               (50  - 2000)    -- shendure ~207, cohen ~1140
    bcs_per_cre          (3   - 500)
    moi                  (0.5 - 30)
    lib_alpha_nb         (0.02 - 2.0)   -- NB2 alpha for library; theta = 1/alpha
    minP                 (0.02 - 2.0)   -- reference_activity
    activity_max_mult    (linear, 2 - 8)

Single test cell type. NB only. Reporter regime fixed at "obs" (best case).
MWU is the canonical test.

Usage:
    python synthetic_factorial.py <phase> [mode] [slice n_slices]
        phase:    samples | simulate | plot | all
        mode:     smoke | pilot | full   (default: full)
        slice:    0..n_slices-1   (default: 0; only meaningful for simulate)
        n_slices: total slice count (default: from MODES[mode]["n_slices"])

    samples  -- write the LHS samples_<mode>.parquet table; do nothing else.
    simulate -- run the simulation phase for one slice. Each slice is an
                independent sbatch job with its own dask scheduler+workers.
                Slice i processes rows where (idx % n_slices) == i.
    plot     -- aggregate all slices' outputs and produce SVGs.
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
    # n_cells upper at 50k: cohen-Rod has 18.6k cells, shendure-Pluripotent
    # 8.2k. 50k gives ~3x headroom over the largest empirical CT without
    # blowing past the MWU sweet spot (200k caused a >1h MWU stall in pilot).
    "n_cells":           (5.0e2, 5.0e4, True),
    # n_cres upper at 2000: cohen 1140, shendure 207.
    "n_cres":            (5.0e1, 2.0e3, True),
    "bcs_per_cre":       (3.0,   5.0e2, True),
    "moi":               (0.5,   3.0e1, True),
    "lib_alpha_nb":      (0.02,  2.0,   True),
    "minP":              (0.02,  2.0,   True),
    "activity_max_mult": (2.0,   8.0,   False),
}

# Top-up sweep bounds. The original LHS box left cohen-Rod outside on
# bcs_per_cre (cohen=17244 vs LHS max 500, 34x past) and moi (cohen=149 vs
# LHS max 30, 5x past). Top-up draws 100 samples in the previously
# unexplored high-bcs/high-moi corner. Other axes use the same ranges as
# the main sweep so the top-up surfaces still combine cleanly with the
# original samples.
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
AXIS_NAMES = list(AXIS_BOUNDS.keys())

# Sample / rep / sim / worker budget per mode. Default is "full"; the wrapper
# selects the mode via CLI arg.
# n_workers is per-slice (each slice has its own dask cluster).
# Total worker pool = n_slices * n_workers.
MODES = {
    "smoke": dict(n_lhs=4,    n_library_reps=2, n_sims=2, n_workers=4,  worker_mem="16G", n_slices=1),
    "pilot": dict(n_lhs=30,   n_library_reps=3, n_sims=5, n_workers=20, worker_mem="48G", n_slices=3),
    "full":  dict(n_lhs=1000, n_library_reps=5, n_sims=5, n_workers=50, worker_mem="48G", n_slices=10),
    "topup": dict(n_lhs=100,  n_library_reps=5, n_sims=5, n_workers=50, worker_mem="48G", n_slices=10),
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

# Activity floor: scale with minP rather than absolute. Empirically
# min_mpra_umi/minP ~ 0.02-0.08 across datasets; pick a single ratio.
MIN_ACT_FRAC_OF_MINP = 0.05

SEED = 20260505

# Empirical anchors for plot overlay.
# Empirical anchor values for plot overlay. activity_max_mult is the
# data-derived ratio p95(mu) / minP -- i.e. how broad the assayed CRE
# activity dynamic range is, robust to single-CRE outliers. Shendure's
# library has strong promoters (eef1aP at 272x baseline) producing a
# heavy tail; cohen's library is variants of a single element, all
# clustered near baseline. Both empirical values fall outside the LHS
# range [2, 8] -- shendure above (99x), cohen below (1.11x).
EMPIRICAL = {
    "shendure-Pluripotent": dict(
        n_cells=8201, n_cres=207, bcs_per_cre=136.06, moi=18.01,
        lib_alpha_nb=0.2241, minP=0.0414, activity_max_mult=99.21,
    ),
    "cohen-Rod": dict(
        n_cells=18633, n_cres=116, bcs_per_cre=17244.5, moi=149.15,
        lib_alpha_nb=1.3873, minP=0.9363, activity_max_mult=1.11,
    ),
}


# ---------------------------------------------------------------------------
# LHS sampling
# ---------------------------------------------------------------------------

def draw_lhs(n: int, seed: int = SEED, bounds: dict = None,
             sample_id_prefix: str = "s") -> pd.DataFrame:
    """LHS draw over the supplied axis bounds.

    bounds: dict[axis] -> (low, high, log_scale). Defaults to AXIS_BOUNDS.
            For the bounds-expansion top-up sweep, pass TOPUP_AXIS_BOUNDS.
    sample_id_prefix: "s" for the original sweep, "t" for the top-up. Keeps
            sample_ids disjoint across sweeps so a combined plot can index by
            sample_id without collisions.
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
    print(f"DONE   {sid}", flush=True)


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


def _aggregate(samples: pd.DataFrame, client, test_type: str = "mwu") -> pd.DataFrame:
    """Walk scratch, per sample compute several power summaries, return one
    row per sample."""
    rows = []
    skipped = 0
    for _, row in samples.iterrows():
        sid = row["sample_id"]
        pt_dir = _sample_dir(sid)
        if not pt_dir.exists():
            continue
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
        print(f"  aggregate: skipped {skipped} unloadable sim dirs", flush=True)
    return pd.DataFrame(rows)


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


def _plot_marginals_for_metric(df: pd.DataFrame, metric: str, ylim: "tuple[float,float]",
                                ylabel: str, title: str, out_path: Path,
                                hline: "float | None" = None,
                                invert_y: bool = False):
    """Render the 7-axis marginals figure for a chosen power metric column."""
    fig, axes = plt.subplots(3, 3, figsize=(13, 11))
    axes = axes.ravel()
    empirical_colors = {"shendure-Pluripotent": "darkorange",
                        "cohen-Rod": "purple"}
    for ax_i, axis in enumerate(AXIS_NAMES):
        ax = axes[ax_i]
        x = df[axis].values.astype(float)
        y = df[metric].values.astype(float)
        log_scale = AXIS_BOUNDS[axis][2]
        x_plot = np.log10(x) if log_scale else x
        # Drop NaN y values from the LOESS input but keep them in scatter
        # (matplotlib already handles NaN).
        ax.scatter(x_plot, y, s=8, alpha=0.4, color="steelblue")
        valid = np.isfinite(y) & np.isfinite(x_plot)
        if valid.sum() >= 20:
            try:
                xg, yhat, lo, hi = _loess_band(x_plot[valid], y[valid])
                ax.plot(xg, yhat, color="firebrick", lw=2)
                ax.fill_between(xg, lo, hi, color="firebrick", alpha=0.2)
            except Exception as e:
                print(f"  smoother failed for {axis}/{metric}: {e}", flush=True)
        # Empirical overlays.
        cur_lo, cur_hi = ax.get_xlim()
        for name, coord in EMPIRICAL.items():
            v = coord.get(axis)
            if v is None:
                continue
            xv = np.log10(v) if log_scale else v
            cur_lo = min(cur_lo, xv)
            cur_hi = max(cur_hi, xv)
            color = empirical_colors.get(name, "black")
            ax.axvline(xv, color=color, linestyle="--", lw=1.0, alpha=0.8)
            ax.text(xv, 1.02, name.split("-")[0][:4],
                    fontsize=7, ha="center", color=color,
                    transform=ax.get_xaxis_transform())
        span = cur_hi - cur_lo
        if span > 0:
            ax.set_xlim(cur_lo - 0.03 * span, cur_hi + 0.03 * span)
        ax.set_xlabel(("log10 " if log_scale else "") + axis)
        ax.set_ylabel(ylabel)
        if ylim is not None:
            ax.set_ylim(*ylim)
        if invert_y:
            ax.invert_yaxis()
        if hline is not None:
            ax.axhline(hline, color="grey", lw=0.5, ls="--")
    for j in range(len(AXIS_NAMES), len(axes)):
        axes[j].axis("off")
    fig.suptitle(title)
    plt.tight_layout()
    fig.savefig(out_path, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}", flush=True)


def phase_plot(samples_with_power: pd.DataFrame, mode: str):
    df = samples_with_power.copy()
    suffix = "" if mode == "full" else f"_{mode}"
    df.to_parquet(OUT / f"samples_power{suffix}.parquet")
    print(f"plotting {len(df)} samples; metrics: power_at_fc1p5/2/3, "
          f"power_auc_1to3, fc_at_p50", flush=True)

    # One marginals SVG per metric. The legacy filename "marginals.svg"
    # remains the P@FC=2 plot for backwards compat with existing references.
    metric_specs = [
        # (metric, ylim, ylabel, hline, invert_y, fname_suffix, title_suffix)
        ("power_at_fc1p5", (0, 1), "power @ FC=1.5", 0.8, False,
         "_p15", "power@FC=1.5 (less ceiling-saturated than P@FC=2)"),
        ("power_at_fc2",   (0, 1), "power @ FC=2",   0.8, False,
         "", "power@FC=2 (legacy summary; ~40% of samples saturate at edges)"),
        ("power_at_fc3",   (0, 1), "power @ FC=3",   0.8, False,
         "_p3", "power@FC=3 (sparse: few CREs sampled at high FC)"),
        ("power_auc_1to3", (0, 1), "AUC of power(FC) over [1,3]", None, False,
         "_auc", "power-curve AUC over FC in [1,3] (continuous, no ceiling)"),
        ("fc_at_p50",      None,   "FC at 50% power (lower=better)", None, True,
         "_fc50", "FC at 50% power (minimum detectable effect; NaN if never reached)"),
    ]
    for metric, ylim, ylabel, hline, invert, fsuf, tsuf in metric_specs:
        out_path = OUT / f"marginals{fsuf}{suffix}.svg"
        _plot_marginals_for_metric(df, metric=metric, ylim=ylim, ylabel=ylabel,
                                   title=f"Synthetic factorial: marginal {tsuf} (LOESS w/ 95% boot band)",
                                   out_path=out_path, hline=hline, invert_y=invert)

    # Combined plot: when running in topup mode, fuse with the original full
    # sweep's samples_power.parquet so the marginals span the full
    # bcs_per_cre and moi range covered by the union. Both DataFrames must
    # have the same metric columns -- they will because both were aggregated
    # by the same _power_metrics function.
    if mode == "topup":
        full_pp = OUT / "samples_power.parquet"
        if full_pp.exists():
            full_df = pd.read_parquet(full_pp)
            # sanity: sample_ids should be disjoint by prefix (s vs t)
            shared = set(full_df["sample_id"]) & set(df["sample_id"])
            if shared:
                print(f"WARNING: {len(shared)} overlapping sample_ids; "
                      "skipping combined plot to avoid double-counting", flush=True)
            else:
                combined = pd.concat([full_df, df], ignore_index=True)
                combined.to_parquet(OUT / "samples_power_combined.parquet")
                print(f"\nCombined plot: {len(full_df)} (full) + {len(df)} (topup) "
                      f"= {len(combined)} samples", flush=True)
                for metric, ylim, ylabel, hline, invert, fsuf, tsuf in metric_specs:
                    out_path = OUT / f"marginals{fsuf}_combined.svg"
                    _plot_marginals_for_metric(
                        combined, metric=metric, ylim=ylim, ylabel=ylabel,
                        title=f"Combined factorial (full + topup): {tsuf} (LOESS w/ 95% boot band)",
                        out_path=out_path, hline=hline, invert_y=invert)
        else:
            print(f"NOTE: {full_pp} not found, skipping combined plot", flush=True)

    # Pairwise heatmaps for top axes by smoother range (= biggest LOESS swing)
    ranges = {}
    for axis in AXIS_NAMES:
        x = df[axis].values.astype(float)
        y = df["power_at_fc2"].values.astype(float)
        log_scale = AXIS_BOUNDS[axis][2]
        x_plot = np.log10(x) if log_scale else x
        try:
            _, yhat, _, _ = _loess_band(x_plot, y, n_boot=20)
            ranges[axis] = float(np.nanmax(yhat) - np.nanmin(yhat))
        except Exception:
            ranges[axis] = 0.0
    top3 = sorted(ranges.items(), key=lambda kv: -kv[1])[:3]
    print(f"Top axes by marginal swing: {top3}", flush=True)
    pairs = [(top3[0][0], top3[1][0]), (top3[0][0], top3[2][0]), (top3[1][0], top3[2][0])]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    for ax, (a, b) in zip(axes, pairs):
        xa = df[a].values.astype(float)
        xb = df[b].values.astype(float)
        if AXIS_BOUNDS[a][2]:
            xa = np.log10(xa); a_lab = f"log10 {a}"
        else:
            a_lab = a
        if AXIS_BOUNDS[b][2]:
            xb = np.log10(xb); b_lab = f"log10 {b}"
        else:
            b_lab = b
        # 8x8 binning, mean power per bin
        nb = 8
        bins_a = np.linspace(xa.min(), xa.max(), nb + 1)
        bins_b = np.linspace(xb.min(), xb.max(), nb + 1)
        H_sum, _, _ = np.histogram2d(xa, xb, bins=[bins_a, bins_b],
                                      weights=df["power_at_fc2"].values)
        H_n, _, _ = np.histogram2d(xa, xb, bins=[bins_a, bins_b])
        with np.errstate(invalid="ignore"):
            H = np.where(H_n > 0, H_sum / H_n, np.nan)
        im = ax.imshow(H.T, origin="lower", aspect="auto", cmap="viridis",
                       vmin=0, vmax=1,
                       extent=[bins_a[0], bins_a[-1], bins_b[0], bins_b[-1]])
        ax.set_xlabel(a_lab)
        ax.set_ylabel(b_lab)
        plt.colorbar(im, ax=ax, label="power @ FC=2")
    fig.suptitle("Pairwise heatmaps (top-3 axes by marginal swing)")
    plt.tight_layout()
    out_pair = OUT / f"pairwise_heatmaps{suffix}.svg"
    fig.savefig(out_pair, format="svg", bbox_inches="tight")
    print(f"Saved: {out_pair}", flush=True)


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

    if phase == "samples":
        # Materialize the samples table and exit. Used by launch.sh to ensure
        # the table exists before parallel sbatch slice jobs race on it.
        get_samples()
        print("samples table ready", flush=True)
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
