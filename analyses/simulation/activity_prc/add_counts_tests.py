"""
Add count-based tests (KS, pseudobulk) to existing 5x5 activity simulations.

Pseudobulk requires >= 2 biological replicates, so it is skipped for Seelig
(1 replicate). KS runs on all datasets.

Usage:
    python add_counts_tests.py [dataset]
    dataset: shendure, cohen, seelig, or all (default: all)
"""
import sys
import os
sys.path.insert(0, "/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa")
import scMPRAforge as scm
import pandas as pd
from pathlib import Path
from dask.distributed import Client, LocalCluster

SIM_ROOT = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new/simulated")
N_GT_DRAWS = 5
HYPOTHESIS_SET = "hs_all_ct"

DATASETS = {
    "shendure": {
        "prefix": "shendure_5x5_activity",
        "has_reporter": True,
        "has_replicates": True,
    },
    "cohen": {
        "prefix": "cohen_5x5_activity",
        "has_reporter": True,
        "has_replicates": True,
    },
    "seelig": {
        "prefix": "seelig_5x5_activity",
        "has_reporter": False,
        "has_replicates": False,
    },
}


def run_tests(client, dataset_name, cfg):
    prefix = cfg["prefix"]
    has_reporter = cfg["has_reporter"]
    has_replicates = cfg["has_replicates"]
    print(f"\n{'='*60}", flush=True)
    print(f"Dataset: {dataset_name} (has_reporter={has_reporter}, "
          f"has_replicates={has_replicates})", flush=True)
    print(f"{'='*60}", flush=True)

    # Which tests to run
    tests_to_run = [("ks", "ks")]
    if has_replicates:
        tests_to_run.append(("pseudobulk", "pseudobulk"))

    all_summaries = []
    for gt_idx in range(N_GT_DRAWS):
        name = f"{prefix}_gt{gt_idx}"
        print(f"  GT draw {gt_idx}: loading {name}...", flush=True)
        try:
            sim = scm.de_novo_simulation(
                location=SIM_ROOT, name=name, client=client
            )
        except Exception as e:
            print(f"  GT draw {gt_idx}: SKIP ({e})", flush=True)
            continue

        for test_label, test_method in tests_to_run:
            test_dir = sim.testd / HYPOTHESIS_SET / test_label
            if test_dir.is_dir():
                print(f"  GT draw {gt_idx}: {test_label} already exists, "
                      f"skipping.", flush=True)
            else:
                getattr(sim, test_method)(
                    HYPOTHESIS_SET, has_reporter=has_reporter
                )
                sim.save()
                print(f"  GT draw {gt_idx}: {test_label} done.", flush=True)

        # Collect metrics
        summary = sim._all_classifier_summary(HYPOTHESIS_SET)
        summary["gt_draw"] = gt_idx
        all_summaries.append(summary)

    if all_summaries:
        combined = pd.concat(all_summaries, ignore_index=True)
        print(f"\n=== {dataset_name} aggregated ===", flush=True)
        print(
            combined.groupby("test")[["auroc", "auprc"]]
            .agg(["mean", "std"])
            .to_string(),
            flush=True,
        )

        outdir = Path(__file__).parent / dataset_name / "output"
        outdir.mkdir(parents=True, exist_ok=True)
        outpath = outdir / f"{dataset_name}_5x5_activity_summary.tsv"
        combined.to_csv(outpath, sep="\t", index=False)
        print(f"Saved to {outpath}", flush=True)


if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else "all"

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
    )
    client = Client(cluster)
    print(f"Dask dashboard: {client.dashboard_link}", flush=True)

    if target == "all":
        for ds_name, cfg in DATASETS.items():
            run_tests(client, ds_name, cfg)
    elif target in DATASETS:
        run_tests(client, target, DATASETS[target])
    else:
        print(f"Unknown dataset '{target}'. "
              f"Choose from: {list(DATASETS.keys())} or 'all'")
        sys.exit(1)

    client.close()
    cluster.close()
