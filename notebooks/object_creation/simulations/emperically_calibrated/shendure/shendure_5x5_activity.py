"""
Shendure activity hypothesis power analysis -- 5 GT draws x 5 sim reps.

Scriptified from two_thirds_inactive.ipynb with the following changes:
  1. Paths updated to tabula_data_new (phantom orthos)
  2. flatten_overtransfection=True (paper_plan requirement)
  3. nb_only=True (canonical model choice for shendure)
  4. phantom_compress=False (simulated data has no reporter table)
  5. Single Wald with cov_method="auto" (auto fallback chain)
  6. 5 independent GT draws to reduce bias from a single draw
  7. Script format for sbatch submission
  8. Stale McCleary symlink fixed

Each GT draw:
  - Loads the canonical phantom ortho to extract real mu values
  - Shuffles mu values into a synthetic CRE x cell-type grid
  - Sets 70% to minP (inactive), 30% retain real empirical values
  - Generates 5 fresh libraries and 5 simulated replicates
  - Fits NB models, precomputes Wald SEs, runs Wald + MWU tests

The shendure dataset has extremely heterogeneous transfection: different CRE
sets appear in different cell types due to differential clonotype contribution.
We cannot use real (cre_id, cell_type) assignments directly because the
simulator allows any CRE in any cell type. Instead, we pool all empirical mu
values, shuffle them independently per cell type to fill the Cartesian product,
then subsample to max_tfection CREs (the number of unique CREs observed in the
most-covered cell type). This preserves the empirical effect size distribution
while avoiding missing-combination artifacts.

See commit 58a65def6964cea9247c15937dd60714489f1750 and 2026-10-23/2026-01-29
notes for further discussion on the clonotype bottleneck workaround.
"""

import sys
import os
import scMPRAforge as scm
from dask.distributed import Client, LocalCluster
from pathlib import Path
import pandas as pd
import numpy as np

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

N_GT_DRAWS = 5
N_SIMS_PER_DRAW = 5
DATA_ROOT = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new")
SIM_ROOT = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new/simulated")

# Canonical phantom ortho for shendure (obs NB)
ORTHO_DIR = DATA_ROOT / "shendure"
ORTHO_NAME = "shendure_obs_nb_phantom"


def sim_name(gt_idx):
    return f"shendure_5x5_activity_gt{gt_idx}"


# ---------------------------------------------------------------------------
# Phase: create -- GT construction + gamut for each draw
# ---------------------------------------------------------------------------

def phase_create(client):
    # Load the canonical phantom ortho to extract real mu values.
    # We use by_cell_type models for more robust per-CRE estimates.
    primordial = scm.ortho.load(client, ORTHO_DIR, ORTHO_NAME)
    primordial.compute_model_qc()

    vals = []
    for key in primordial.by_cell_qc.keys():
        working = primordial.by_cell_qc[key]["dat"].reset_index().drop(
            columns=["mean(umis_mpra_bc)"]
        )
        working["cell_type"] = key
        vals.append(working)
    gt_cell_type = pd.concat(vals)

    # Cast away from sparse dtype
    gt_cell_type["mu"] = gt_cell_type["mu"].astype(float)
    gt_cell_type["cre_id"] = gt_cell_type["cre_id"].astype(str)
    gt_cell_type["cell_type"] = gt_cell_type["cell_type"].astype(str)

    # Count how many CREs each cell type sees (for subsampling later).
    # max_tfection = the most CREs any single cell type has.
    combo_counts = (
        gt_cell_type[gt_cell_type["cre_id"] != "reference"]
        .drop(columns=["mu"])
        .groupby("cell_type")
        .nunique()
    )
    max_tfection = int(max(combo_counts["cre_id"]))

    # Strip identifying info -- just keep mu values as a pool
    gt_cell_type = gt_cell_type[gt_cell_type["cre_id"] != "reference"]
    real_means = gt_cell_type.drop(columns=["cre_id", "cell_type"])
    num_cell_types = len(gt_cell_type["cell_type"].unique())

    minP = scm.SHENDURE_BOUNDS.reference_activity

    # Build modified bounds: start from NB bounds (correct model parameters),
    # but NB bounds have zi=None which breaks the simulation's rep_id loop.
    # Pull in rep_ids from the ZINB bounds to preserve the 6-replicate
    # experiment structure, and construct a zero-valued zi Series so the
    # simulation generates pure NB data.
    bound_template = scm.SHENDURE_BOUNDS.copy()
    zinb_zi = scm.SHENDURE_OBS_ZINB_BOUNDS.zi
    bound_template.zi = zinb_zi * 0.0  # same index, all zeros
    bound_template.rep_ids = list(zinb_zi.index)
    s = bound_template.cells_per_cell_type
    s.index = [f"ct_{i}" for i in range(1, num_cell_types)] + ["reference"]
    bound_template.cells_per_cell_type = s

    for gt_idx in range(N_GT_DRAWS):
        print(f"--- GT draw {gt_idx} ---", flush=True)
        np.random.seed(gt_idx * 1000)  # reproducible but distinct per draw

        # Shuffle mu values independently per cell type to fill the Cartesian
        # product. Each cell type gets a random permutation of all observed mus,
        # assigned to synthetic CRE names.
        synth_cre_names = [f"synthcre_{i}" for i in range(len(real_means))]

        parts = []
        for i in range(num_cell_types):
            working = real_means.sample(frac=1).reset_index(drop=True)
            working["cre_id"] = synth_cre_names
            working["cell_type"] = f"ct_{i}"
            parts.append(working)

        cartesian = pd.concat(parts).sample(frac=1).reset_index(drop=True)

        # Set 70% of entries to minP (inactive). After the shuffle above this
        # is effectively random. Note: applied to the full Cartesian frame
        # before CRE subsampling, so the realized active fraction per CRE may
        # deviate slightly from 30%.
        num_inactive = int(len(cartesian) * 0.7)
        cartesian.loc[cartesian.index[:num_inactive], "mu"] = minP

        # Subsample to max_tfection CREs (the number observed in the
        # most-covered cell type in the real data).
        chosen_cre_id = np.random.choice(
            cartesian["cre_id"].unique(), max_tfection, replace=False
        )
        sampled = cartesian[cartesian["cre_id"].isin(chosen_cre_id)]

        # Add universal reference CRE (minP in every cell type)
        ref_df = pd.DataFrame(
            {"cell_type": [f"ct_{i}" for i in range(num_cell_types)]}
        )
        ref_df["mu"] = minP
        ref_df["cre_id"] = "reference"

        # Stack and rename ct_0 -> "reference" (one cell type serves as ref)
        final_gt = pd.concat([ref_df, sampled])
        final_gt["cell_type"] = final_gt["cell_type"].replace(
            {"ct_0": "reference"}
        )

        # Simulate 5 independent libraries from shendure library model
        libraries = [
            scm.simulate_library(
                CREs=final_gt["cre_id"],
                library_model=scm.SHENDURE_BOUNDS.library_model,
                reference_pooling=scm.SHENDURE_BOUNDS.n_negative_controls,
            )
            for _ in range(N_SIMS_PER_DRAW)
        ]

        bound = bound_template.copy()

        sim = scm.de_novo_simulation(
            location=SIM_ROOT,
            name=sim_name(gt_idx),
            client=client,
            libraries=libraries,
            library_mapping="corresponding",
            flatten_overtransfection=True,
            n_sims=N_SIMS_PER_DRAW,
            experiment_bounds=bound,
            ground_truth=final_gt,
        )

        sim.gamut()
        sim.save()
        print(f"GT draw {gt_idx}: gamut done, saved.", flush=True)


# ---------------------------------------------------------------------------
# Phase: fit -- NB model fitting for each sim
# ---------------------------------------------------------------------------

def phase_fit(client):
    for gt_idx in range(N_GT_DRAWS):
        print(f"--- Fitting GT draw {gt_idx} ---", flush=True)
        sim = scm.de_novo_simulation(
            location=SIM_ROOT, name=sim_name(gt_idx), client=client
        )
        # phantom_compress=False: simulated data has no coarse reporter table
        # (phantom compression is for empirical data with transfection reporters)
        # direction="by_cell_type": activity hypotheses only need by_cell_type
        # models (CRE as treatment within cell type). Skip by_cre.
        # serial_orthos=True required: with LocalCluster, fit_ortho_helper
        # submits sub-tasks (standard_fit) that deadlock if the single worker
        # is blocked on the outer task.
        sim.fit_orthos(
            direction="by_cell_type",
            serial_orthos=True,
            nb_only=True,
            phantom_compress=False,
            gpu=True,
        )
        sim.save()
        print(f"GT draw {gt_idx}: fit done.", flush=True)


# ---------------------------------------------------------------------------
# Phase: wald_precomp -- Wald SE precomputation (sandwich with auto fallback)
# ---------------------------------------------------------------------------

def phase_wald_precomp(client):
    # Bypass sim.precompute_wald() to avoid nested Dask task deadlock.
    # Load each ortho directly and call ortho.precompute_wald(client).
    # Then inject dummy futures so sim.wald() sees the "auto" column.
    for gt_idx in range(N_GT_DRAWS):
        print(f"--- Wald precomp GT draw {gt_idx} ---", flush=True)
        sim = scm.de_novo_simulation(
            location=SIM_ROOT, name=sim_name(gt_idx), client=client
        )
        autod = sim.fullp / "auto"
        autod.mkdir(exist_ok=True)
        n_sims = sim.get_state_field("n_sims")
        done_futures = []
        for idx in range(n_sims):
            print(f"  rep {idx}...", flush=True)
            orth = scm.ortho.load(client, path=sim.orthod, name=str(idx))
            orth.precompute_wald(client, cov_method="sandwich")
            orth.save(path=autod, name=str(idx))
            done_futures.append(client.submit(lambda: True))
            print(f"  rep {idx} done.", flush=True)
        # Inject completed futures so sim.wald("hs_all_ct", cov_method="auto")
        # finds the "auto" column in sim.futures.
        sim.futures["auto"] = done_futures
        sim.save()
        print(f"GT draw {gt_idx}: wald precomp done.", flush=True)


# ---------------------------------------------------------------------------
# Phase: test -- hypothesis testing (Wald + MWU)
# ---------------------------------------------------------------------------

def phase_test(client):
    for gt_idx in range(N_GT_DRAWS):
        print(f"--- Testing GT draw {gt_idx} ---", flush=True)
        sim = scm.de_novo_simulation(
            location=SIM_ROOT, name=sim_name(gt_idx), client=client
        )

        # Add activity hypotheses (CRE vs reference, per cell type)
        example_data = scm.scMPRA_data.from_parquet(
            sim.scmpradatp / "0.scmpra"
        )
        hs_all_ct = scm.make_all_by_celltype_hypotheses(
            counts=example_data,
            reference_cre="reference",
        )
        sim.add_hypothesis_set("hs_all_ct", hs_all_ct)
        sim.save()

        sim.wald("hs_all_ct", cov_method="auto")
        sim.save()
        print(f"GT draw {gt_idx}: wald done.", flush=True)

        sim.mwu("hs_all_ct")
        sim.save()
        print(f"GT draw {gt_idx}: mwu done.", flush=True)


# ---------------------------------------------------------------------------
# Phase: metrics -- compute and print AUROC/AUPRC summaries
# ---------------------------------------------------------------------------

def phase_metrics(client):
    all_summaries = []
    for gt_idx in range(N_GT_DRAWS):
        sim = scm.de_novo_simulation(
            location=SIM_ROOT, name=sim_name(gt_idx), client=client
        )
        summary = sim._all_classifier_summary("hs_all_ct")
        summary["gt_draw"] = gt_idx
        all_summaries.append(summary)
        print(f"GT draw {gt_idx}:", flush=True)
        print(summary.to_string(), flush=True)

    combined = pd.concat(all_summaries, ignore_index=True)
    print("\n=== Aggregated across all GT draws ===", flush=True)
    print(
        combined.groupby("test")[["auroc", "auprc"]]
        .agg(["mean", "std", "min", "max"])
        .to_string(),
        flush=True,
    )

    # Save combined results
    outdir = Path(__file__).parent / "output"
    outdir.mkdir(exist_ok=True)
    combined.to_csv(outdir / "shendure_5x5_activity_summary.tsv", sep="\t", index=False)
    print(f"\nSaved to {outdir / 'shendure_5x5_activity_summary.tsv'}", flush=True)


# ---------------------------------------------------------------------------
# Main dispatch
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    phase = sys.argv[1] if len(sys.argv) > 1 else "all"

    # Detect available memory from SLURM allocation (or system total).
    import psutil
    _slurm_mem_mb = os.environ.get("SLURM_MEM_PER_NODE")
    if _slurm_mem_mb:
        _mem_limit = f"{int(int(_slurm_mem_mb) * 0.9)}MB"
    else:
        _mem_limit = f"{int(psutil.virtual_memory().total * 0.9 / 1e9)}GB"

    # Single worker: fit uses serial_orthos=True, wald_precomp is bypassed
    # (direct Python loop). No nested-task deadlock risk.
    # processes=False: avoids TF multi-thread/process conflicts.
    cluster = LocalCluster(
        n_workers=1,
        threads_per_worker=2,
        processes=False,
        memory_limit=_mem_limit,
        resources={"GPU": 1},
    )
    client = Client(cluster)
    print(client.dashboard_link, flush=True)

    phases = {
        "create": phase_create,
        "fit": phase_fit,
        "wald_precomp": phase_wald_precomp,
        "test": phase_test,
        "metrics": phase_metrics,
    }

    if phase == "all":
        for name, fn in phases.items():
            print(f"\n{'='*60}", flush=True)
            print(f"PHASE: {name}", flush=True)
            print(f"{'='*60}", flush=True)
            fn(client)
    elif phase in phases:
        phases[phase](client)
    else:
        print(f"Unknown phase '{phase}'. Valid: {list(phases.keys()) + ['all']}")
        sys.exit(1)

    client.close()
    cluster.close()
    print("DONE", flush=True)
