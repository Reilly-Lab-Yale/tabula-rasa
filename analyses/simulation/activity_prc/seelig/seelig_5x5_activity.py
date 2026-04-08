"""
Seelig activity hypothesis power analysis -- 5 GT draws x 5 sim reps.

Seelig-specific notes:
  - Canonical model: NB (CM condition, no zero inflation detected)
  - Only 2 cell types: K562 + reference
  - 1 biological replicate (R1) -- zi splice from ZINB bounds needed
  - Fixed-count library: exactly 5 barcodes per CRE (synthesized, not PCR'd)
  - Seelig has no transfection reporter, so the empirical data uses
    consider_missing. For simulation, this doesn't matter -- we generate
    synthetic data directly.
  - With only 2 cell types, the by_cell_type models have ~1344 CREs as
    regressors. This makes Wald precompute slow (~3000s/model on GPU).
    For 5x5=25 orthos this is feasible but not fast.

See shendure_5x5_activity.py for detailed comments on the general approach.
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

# Canonical phantom ortho for seelig (MOIB NB)
ORTHO_DIR = DATA_ROOT / "seelig"
ORTHO_NAME = "seelig_cm_moib_phantom"


def sim_name(gt_idx):
    return f"seelig_5x5_activity_gt{gt_idx}"


# ---------------------------------------------------------------------------
# Phase: create -- GT construction + gamut for each draw
# ---------------------------------------------------------------------------

def phase_create(client):
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

    gt_cell_type["mu"] = gt_cell_type["mu"].astype(float)
    gt_cell_type["cre_id"] = gt_cell_type["cre_id"].astype(str)
    gt_cell_type["cell_type"] = gt_cell_type["cell_type"].astype(str)

    combo_counts = (
        gt_cell_type[gt_cell_type["cre_id"] != "reference"]
        .drop(columns=["mu"])
        .groupby("cell_type")
        .nunique()
    )
    max_tfection = int(max(combo_counts["cre_id"]))

    gt_cell_type = gt_cell_type[gt_cell_type["cre_id"] != "reference"]
    real_means = gt_cell_type.drop(columns=["cre_id", "cell_type"])
    num_cell_types = len(gt_cell_type["cell_type"].unique())

    minP = scm.SEELIG_BOUNDS.reference_activity

    # NB bounds -- splice in zero-valued zi from ZINB bounds for rep structure
    # Seelig has only 1 biological replicate (R1)
    bound_template = scm.SEELIG_BOUNDS.copy()
    zinb_zi = scm.SEELIG_CM_ZINB_BOUNDS.zi
    bound_template.zi = zinb_zi * 0.0  # same index (R1), zero ZI
    bound_template.rep_ids = list(zinb_zi.index)
    s = bound_template.cells_per_cell_type
    # Seelig only has 2 cell types: K562 + reference
    s.index = [f"ct_{i}" for i in range(1, num_cell_types)] + ["reference"]
    bound_template.cells_per_cell_type = s

    for gt_idx in range(N_GT_DRAWS):
        print(f"--- GT draw {gt_idx} ---", flush=True)
        np.random.seed(gt_idx * 3000)  # distinct from shendure/cohen seeds

        synth_cre_names = [f"synthcre_{i}" for i in range(len(real_means))]

        parts = []
        for i in range(num_cell_types):
            working = real_means.sample(frac=1).reset_index(drop=True)
            working["cre_id"] = synth_cre_names
            working["cell_type"] = f"ct_{i}"
            parts.append(working)

        cartesian = pd.concat(parts).sample(frac=1).reset_index(drop=True)

        # 70% inactive
        num_inactive = int(len(cartesian) * 0.7)
        cartesian.loc[cartesian.index[:num_inactive], "mu"] = minP

        chosen_cre_id = np.random.choice(
            cartesian["cre_id"].unique(), max_tfection, replace=False
        )
        sampled = cartesian[cartesian["cre_id"].isin(chosen_cre_id)]

        ref_df = pd.DataFrame(
            {"cell_type": [f"ct_{i}" for i in range(num_cell_types)]}
        )
        ref_df["mu"] = minP
        ref_df["cre_id"] = "reference"

        final_gt = pd.concat([ref_df, sampled])
        final_gt["cell_type"] = final_gt["cell_type"].replace(
            {"ct_0": "reference"}
        )

        # Seelig has a fixed-count library (5 barcodes per CRE)
        libraries = [
            scm.simulate_library(
                CREs=final_gt["cre_id"],
                library_model=scm.SEELIG_BOUNDS.library_model,
                reference_pooling=scm.SEELIG_BOUNDS.n_negative_controls,
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

        # consider_missing=True: seelig has no transfection reporter, so the
        # canonical analysis uses CM condition. This sets consider_missing_enabled
        # on each simulated scMPRA_data, which persists through save/load and is
        # picked up by the fitting code.
        sim.gamut(consider_missing=True)
        sim.save()
        print(f"GT draw {gt_idx}: gamut done, saved.", flush=True)


# ---------------------------------------------------------------------------
# Phase: fit -- NB model fitting
# ---------------------------------------------------------------------------

def phase_fit(client):
    for gt_idx in range(N_GT_DRAWS):
        print(f"--- Fitting GT draw {gt_idx} ---", flush=True)
        sim = scm.de_novo_simulation(
            location=SIM_ROOT, name=sim_name(gt_idx), client=client
        )
        # Only by_cell_type for activity hypotheses
        # serial_orthos=True: single worker avoids deadlock and TF GPU conflicts.
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
# Phase: wald_precomp
# ---------------------------------------------------------------------------

def phase_wald_precomp(client):
    # Bypass sim.precompute_wald() to avoid nested Dask task deadlock.
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
# Phase: metrics
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

    outdir = Path(__file__).parent / "output"
    outdir.mkdir(exist_ok=True)
    combined.to_csv(outdir / "seelig_5x5_activity_summary.tsv", sep="\t", index=False)
    print(f"\nSaved to {outdir / 'seelig_5x5_activity_summary.tsv'}", flush=True)


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
