#!/usr/bin/env python
"""
Shendure NB-only by-cell-type fit — GPU variant for head-to-head vs CPU run.

Races against fit_by_cell_type.py (CPU, 16×96GB) to benchmark GPU speedup.
Output saved to shendure_cm_nb_20260329_by_cell_type_gpu (separate path).

4 CPU workers handle setup (consider_missing expansion + design matrices),
then GPU workers are added just-in-time via pre_fit_hook for TensorZINB fitting.
"""
import signal
import subprocess
import tempfile
import os
import threading
import time
import scMPRAforge as scm
from dask_jobqueue import SLURMCluster
from dask.distributed import Client, performance_report
from pathlib import Path

signal.signal(signal.SIGTERM, lambda s, f: (_ for _ in ()).throw(KeyboardInterrupt()))

data_root = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data")
path = data_root / "shendure"
name = "shendure_cm_nb_20260329_by_cell_type_gpu"
work_dir = str(Path(__file__).resolve().parent)
report_path = Path(__file__).with_name(f"{name}_dask_performance_report.html")
JOB = "shend_cm_nb_ct_gpu"

def notify(msg, title=None, tags=None, priority=None):
    env = os.environ.copy()
    env["NTFY_TITLE"] = title or JOB
    if tags:
        env["NTFY_TAGS"] = tags
    if priority:
        env["NTFY_PRIORITY"] = priority
    subprocess.run(["/home/mcn26/.local/bin/notify-job", msg], env=env, check=False)

# ── CPU cluster for setup ──────────────────────────────────────────────────────
cluster = SLURMCluster(
    cores=1, memory="96G", processes=1,
    job_extra_directives=[
        "-p priority", "-A prio_skr2",
        "--job-name=shend_cm_nb_ct_gpu_cpu_w", "--time=2-23:40:00",
        f"--output={work_dir}/worker_cpu_gpu_%j.out"
    ]
)
cluster.scale(jobs=4)
client = Client(cluster, timeout=f"{5*60}s", heartbeat_interval="20s")
print(client.dashboard_link, flush=True)

scheduler_addr = cluster.scheduler_address
print(f"[+] Scheduler: {scheduler_addr}", flush=True)
print("[+] Waiting for 4 CPU workers...", flush=True)
client.wait_for_workers(n_workers=4, timeout=600)
print("[+] CPU workers connected.", flush=True)

gpu_script = f"""#!/bin/bash
#SBATCH -p priority_gpu
#SBATCH -A prio_skr2
#SBATCH --job-name=shend_cm_nb_ct_gpu_w
#SBATCH -t 1-00:00:00
#SBATCH -c 1
#SBATCH --mem=64G
#SBATCH --gres=gpu:h200:1
#SBATCH --output={work_dir}/worker_gpu_%j.out

CUDA_ROOT=/apps/software/2024a/software/CUDA/12.6.0
export LD_LIBRARY_PATH=${{CUDA_ROOT}}/targets/x86_64-linux/lib:${{CUDA_ROOT}}/extras/CUPTI/lib64:/home/mcn26/.conda/envs/tz/lib:$LD_LIBRARY_PATH
export PYTHONPATH=/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa:$PYTHONPATH
/home/mcn26/.conda/envs/tz/bin/dask-worker {scheduler_addr} --resources GPU=1 --nthreads 1 --memory-limit 60GiB
"""

gpu_job_ids = []
_success = False

def spin_up_gpu():
    """Called by standard_fit after setup completes, before _tensorzinb_fit."""
    print("[+] Setup complete. Submitting GPU workers...", flush=True)
    for _ in range(2):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sh', delete=False, dir='/tmp') as f:
            f.write(gpu_script)
            tmp = f.name
        result = subprocess.run(['sbatch', tmp], capture_output=True, text=True)
        os.unlink(tmp)
        job_id = result.stdout.strip().split()[-1]
        gpu_job_ids.append(job_id)
        print(f"[+] Submitted GPU worker: {job_id}", flush=True)

    notify(
        f"GPU workers submitted: {', '.join(gpu_job_ids)}. Waiting to connect...",
        title=f"{JOB} — GPU submitted", tags="satellite"
    )

    print("[+] Waiting for GPU workers to connect...", flush=True)
    client.wait_for_workers(n_workers=6, timeout=1200)
    print("[+] All 6 workers connected (4 CPU + 2 GPU). Fitting started.", flush=True)

    notify(
        f"All 6 workers connected (4 CPU + 2 GPU). Fitting started.\nGPU jobs: {', '.join(gpu_job_ids)}",
        title=f"{JOB} — fitting started", tags="zap"
    )

try:
    with performance_report(filename=str(report_path)):
        if (path / name).is_dir():
            print("[+] Partial result found. Loading...", flush=True)
            primordial = scm.ortho.load(client, path, name)
        else:
            print("[+] Creating...", flush=True)
            shendure = scm.scMPRA_data.from_tsv(str(path / "shendure_processed.tsv"))
            shendure.set_negative_controls(["minP", "noP"])
            shendure.set_reference_cell("Pluripotent")
            shendure.ortho_filter()
            shendure.set_consider_missing(True)

            primordial = scm.ortho()
            primordial.fit_by_cell_type_models(
                client=client, dat=shendure,
                fit_resources={"GPU": 1},
                pre_fit_hook=spin_up_gpu,
                nb_only=True)
            primordial.extract_params(client)
            primordial.save(path, name, client=client)
            print("[+] Done.", flush=True)
            _success = True
            notify(
                f"GPU by_cell_type fit complete. Saved to {path / name}",
                title=f"{JOB} — SUCCESS", tags="white_check_mark"
            )
finally:
    if not _success:
        notify(
            f"GPU by_cell_type FAILED or cancelled. Check {work_dir}/slurm-*.out",
            title=f"{JOB} — FAILED", tags="x,warning", priority="high"
        )
    for jid in gpu_job_ids:
        subprocess.run(['scancel', jid], capture_output=True)
    print("! Done, shutting down", flush=True)
    time.sleep(10)
    client.close()
    cluster.close()
    time.sleep(10)
