#!/usr/bin/env python
"""
Cohen ZINB by-cell-type fit (consider_missing=True, GPU).
New condition for NB vs ZINB comparison — tests the "ignore reporter" scenario.
Run in parallel with fit_by_cre.py, then combine with merge.py.

GPU workers are submitted just-in-time via pre_fit_hook.
"""
import signal
import subprocess
import tempfile
import os
import threading
import time
import scMPRAforge as scm
from dask_jobqueue import SLURMCluster
from dask.distributed import Client
from pathlib import Path

signal.signal(signal.SIGTERM, lambda s, f: (_ for _ in ()).throw(KeyboardInterrupt()))

data_root = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data")
path = data_root / "cohen"
name = "cohen_cm_zinb_20260329_by_cell_type"
work_dir = str(Path(__file__).resolve().parent)
JOB = "cohen_cm_zinb"

def notify(msg, title=None, tags=None, priority=None):
    env = os.environ.copy()
    env["NTFY_TITLE"] = title or JOB
    if tags:
        env["NTFY_TAGS"] = tags
    if priority:
        env["NTFY_PRIORITY"] = priority
    subprocess.run(["/home/mcn26/.local/bin/notify-job", msg], env=env)

def gpu_util_check(job_ids, delay=600):
    """Fire ~10 min after GPU workers connect; ssh to worker nodes for nvidia-smi."""
    time.sleep(delay)
    lines = []
    for jid in job_ids:
        r = subprocess.run(
            ['squeue', '-j', jid, '-h', '-o', '%N'],
            capture_output=True, text=True
        )
        node = r.stdout.strip()
        if not node or node in ('', 'None assigned', '(null)'):
            lines.append(f"job {jid}: node not yet assigned or job finished")
            continue
        r2 = subprocess.run(
            ['ssh', '-o', 'StrictHostKeyChecking=no', '-o', 'ConnectTimeout=10',
             node, 'nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total --format=csv,noheader'],
            capture_output=True, text=True
        )
        if r2.returncode == 0:
            lines.append(f"{node}: {r2.stdout.strip()}")
        else:
            lines.append(f"{node}: ssh failed ({r2.stderr.strip()[:60]})")
    msg = "GPU utilization check (~10 min post-start):\n" + "\n".join(lines) if lines else "GPU util: no data"
    notify(msg, title=f"{JOB} — GPU util check", tags="chart_with_upwards_trend")

# ── cluster ────────────────────────────────────────────────────────────────────
cluster = SLURMCluster(
    cores=1, memory="128G", processes=1,
    job_extra_directives=[
        "-p priority", "-A prio_skr2",
        "--job-name=cohen_cm_zi_ct_w", "--time=12:00:00",
        f"--output={work_dir}/worker_cpu_%j.out"
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
#SBATCH --job-name=cohen_cm_zi_ct_gpu
#SBATCH -t 1-00:00:00
#SBATCH -c 1
#SBATCH --mem=64G
#SBATCH --gres=gpu:h200:1
#SBATCH --output={work_dir}/worker_gpu_%j.out

CUDA_ROOT=/apps/software/2024a/software/CUDA/12.6.0
export LD_LIBRARY_PATH=${{CUDA_ROOT}}/targets/x86_64-linux/lib:${{CUDA_ROOT}}/extras/CUPTI/lib64:/home/mcn26/.conda/envs/tz/lib:$LD_LIBRARY_PATH
export PYTHONPATH=/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa-cohen-regen:$PYTHONPATH
/home/mcn26/.conda/envs/tz/bin/dask-worker {scheduler_addr} --resources GPU=1 --nthreads 1 --memory-limit 60GiB
"""

gpu_job_ids = []
_success = False

def spin_up_gpu():
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
        f"GPU workers submitted (jobs {', '.join(gpu_job_ids)}). Waiting to connect...",
        title=f"{JOB} — GPU submitted", tags="rocket"
    )

    print("[+] Waiting for GPU workers to connect...", flush=True)
    client.wait_for_workers(n_workers=6, timeout=1200)
    print("[+] All 6 workers connected. Submitting fits.", flush=True)

    notify(
        f"All 6 workers connected (4 CPU + 2 GPU). Fitting started.\nGPU jobs: {', '.join(gpu_job_ids)}",
        title=f"{JOB} — fitting started", tags="zap"
    )

    threading.Thread(
        target=gpu_util_check, args=(gpu_job_ids,), daemon=True
    ).start()

try:
    if (path / name).is_dir():
        print("[+] Partial result found. Loading...", flush=True)
        primordial = scm.ortho.load(client, path, name)
    else:
        print("[+] Creating...", flush=True)
        # Load pre-CRE-coarse-expansion data (filtered read-wise → UMI-wise, ~3.2M rows).
        # Do NOT load retina_single_counting_u6.scmpra here — that object contains the
        # full CRE-coarse reporter-informed zeros (1.49B rows), which OOMs workers on
        # even trivial operations like enumerating cell_type levels. CM fits must start
        # from compact pre-expansion data and let consider_missing handle zero inflation.
        cohen = scm.scMPRA_data.from_tsv(str(path / "unjoined/read_wise_mpra_retina_filtered.tsv"))
        cohen.read_wise_to_umi_wise()
        cohen.set_negative_controls(["wt_1", "wt_2"])
        cohen.set_reference_cell("Rod")
        cohen.ortho_filter()
        cohen.set_consider_missing(True)
        # The memory cap estimator assumes dense string-per-row representation and
        # doesn't account for sparse encoding of umis_mpra_bc (99.8% zeros).
        # Bypass it — actual memory usage will be far lower than the estimate.
        cohen.consider_missing_max_memory_gb = None

        primordial = scm.ortho()
        primordial.fit_by_cell_type_models(
            client=client, dat=cohen,
            fit_resources={"GPU": 1},
            pre_fit_hook=spin_up_gpu)
        primordial.extract_params(client)
        primordial.save(path, name, client=client)
        print("[+] Done.", flush=True)
        _success = True
        notify(
            f"by_cell_type fit complete. Saved to {path / name}",
            title=f"{JOB} — SUCCESS", tags="white_check_mark"
        )
finally:
    if not _success:
        notify(
            f"by_cell_type fit FAILED or was cancelled. Check {work_dir}/slurm-*.out",
            title=f"{JOB} — FAILED", tags="x,warning", priority="high"
        )
    for jid in gpu_job_ids:
        subprocess.run(['scancel', jid], capture_output=True)
    client.close()
    cluster.close()
