"""
Benchmark: Wald precomputation timing on real phantom orthos.

Measures per-model wall time for _build_wald_precomp_for_subset
across different dataset sizes and model types.

Results written to bench_wald_timing_results.txt in the same directory.

Usage (CPU):
    sbatch --mem=64G -c 4 -t 2:00:00 -p priority -A prio_skr2 \
        -J wald_bench_cpu \
        --wrap="module load miniconda; conda activate tz; \
                python3 scMPRAforge/tests/bench_wald_timing.py --device cpu"

Usage (GPU):
    sbatch --mem=64G -c 4 -t 2:00:00 -p priority_gpu -A prio_skr2 \
        --gpus=1 -J wald_bench_gpu \
        --wrap="module load miniconda; conda activate tz; \
                python3 scMPRAforge/tests/bench_wald_timing.py --device gpu"
"""
import sys, time, argparse
sys.path.insert(0, '/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa')

import numpy as np
import os

parser = argparse.ArgumentParser()
parser.add_argument("--device", choices=["cpu", "gpu"], default="cpu")
args = parser.parse_args()

# Force TF device placement
if args.device == "cpu":
    os.environ["CUDA_VISIBLE_DEVICES"] = ""
# For GPU, let TF find the device naturally

import tensorflow as tf
print(f"TF devices: {[d.name for d in tf.config.list_physical_devices()]}", flush=True)
print(f"Benchmark device: {args.device}", flush=True)

from dask.distributed import Client, LocalCluster
import scMPRAforge.core as scm
from pathlib import Path

DATA_NEW = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new")
OUT_FILE = Path(__file__).parent / f"bench_wald_timing_results_{args.device}.txt"

# Orthos to benchmark, covering different sizes and model types
BENCHMARKS = [
    # (dataset_dir, ortho_name, label)
    ("cohen", "cohen_obs_nb_phantom_20260401", "cohen_obs_nb"),
    ("cohen", "cohen_obs_zinb_phantom_20260401", "cohen_obs_zinb"),
    ("shendure", "shendure_obs_nb_phantom", "shendure_obs_nb"),
    ("shendure", "shendure_obs_zinb_phantom", "shendure_obs_zinb"),
    ("seelig", "seelig_cm_nb_phantom", "seelig_cm_nb"),
    ("seelig", "seelig_cm_zinb_phantom", "seelig_cm_zinb"),
]

cluster = LocalCluster(n_workers=1, threads_per_worker=2, processes=False, memory_limit="50GB")
client = Client(cluster)
print(f"Dashboard: {client.dashboard_link}", flush=True)

results = []

for dataset, ortho_name, label in BENCHMARKS:
    print(f"\n{'='*60}", flush=True)
    print(f"Benchmarking: {label}", flush=True)
    print(f"{'='*60}", flush=True)

    data_dir = str(DATA_NEW / dataset)
    ortho = scm.ortho.load(client, data_dir, ortho_name)

    n_by_ct = len(ortho.by_cell_type.model) if ortho.by_cell_type else 0
    n_by_cre = len(ortho.by_cre.model) if ortho.by_cre else 0
    n_total = n_by_ct + n_by_cre
    is_nb = "nb" in label

    print(f"  Models: {n_by_ct} by_ct + {n_by_cre} by_cre = {n_total} total", flush=True)
    print(f"  Type: {'NB' if is_nb else 'ZINB'}", flush=True)

    # Time precompute_wald
    t0 = time.time()
    ortho.precompute_wald(client, cov_method="sandwich")
    elapsed = time.time() - t0

    # Flatten and check results
    wp = ortho.wald_precomp.flattened_copy()

    n_ok = 0
    n_fail = 0
    for side in [wp.by_cell_type, wp.by_cre]:
        if side is None:
            continue
        for key, entry in side.items():
            if np.all(np.isfinite(entry.se_x_mu)) and np.all(entry.se_x_mu > 0):
                n_ok += 1
            else:
                n_fail += 1

    per_model = elapsed / max(n_total, 1)

    result = {
        "label": label,
        "device": args.device,
        "n_by_ct": n_by_ct,
        "n_by_cre": n_by_cre,
        "n_total": n_total,
        "elapsed_s": elapsed,
        "per_model_s": per_model,
        "n_ok": n_ok,
        "n_fail": n_fail,
    }
    results.append(result)

    print(f"  Elapsed: {elapsed:.1f}s ({elapsed/60:.1f} min)", flush=True)
    print(f"  Per model: {per_model:.2f}s", flush=True)
    print(f"  Success: {n_ok}/{n_total}, Failed: {n_fail}/{n_total}", flush=True)

    # Free memory
    del ortho, wp

client.close()
cluster.close()

# Write results
print(f"\n{'='*60}", flush=True)
print(f"SUMMARY ({args.device.upper()})", flush=True)
print(f"{'='*60}", flush=True)

lines = [
    f"Wald precompute_wald timing benchmark -- device={args.device}",
    f"Date: {time.strftime('%Y-%m-%d %H:%M:%S')}",
    f"TF devices: {[d.name for d in tf.config.list_physical_devices()]}",
    "",
    f"{'Label':30s} {'Models':>8s} {'Total(s)':>10s} {'Per-model(s)':>14s} {'OK':>5s} {'Fail':>5s}",
    "-" * 80,
]

for r in results:
    line = (f"{r['label']:30s} {r['n_total']:>8d} {r['elapsed_s']:>10.1f} "
            f"{r['per_model_s']:>14.2f} {r['n_ok']:>5d} {r['n_fail']:>5d}")
    lines.append(line)
    print(line, flush=True)

# Extrapolation for power analysis
lines.append("")
lines.append("Extrapolation for power analysis:")
lines.append("  Assuming 20 simulation replicates, each needing Wald precomp:")

for r in results:
    n_sims = 20
    est_time = r["per_model_s"] * r["n_total"] * n_sims
    lines.append(
        f"  {r['label']:30s}: {n_sims} sims x {r['n_total']} models x "
        f"{r['per_model_s']:.2f}s = {est_time:.0f}s ({est_time/3600:.1f}h)"
    )

OUT_FILE.write_text("\n".join(lines) + "\n")
print(f"\nResults written to: {OUT_FILE}", flush=True)
print("Done.", flush=True)
