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
import scMPRAforge as scm
from dask_jobqueue import SLURMCluster
from dask.distributed import Client
from pathlib import Path

signal.signal(signal.SIGTERM, lambda s, f: (_ for _ in ()).throw(KeyboardInterrupt()))

data_root = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data")
path = data_root / "cohen"
name = "cohen_cm_zinb_20260329_by_cell_type"
work_dir = str(Path(__file__).resolve().parent)

cluster = SLURMCluster(
    cores=1, memory="64G", processes=1,
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
/home/mcn26/.conda/envs/tz/bin/dask-worker {scheduler_addr} --resources GPU=1 --nthreads 1 --memory-limit 60GiB
"""

gpu_job_ids = []

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
    print("[+] Waiting for GPU workers to connect...", flush=True)
    client.wait_for_workers(n_workers=6, timeout=1200)
    print("[+] All 6 workers connected. Submitting fits.", flush=True)

try:
    if (path / name).is_dir():
        print("[+] Partial result found. Loading...", flush=True)
        primordial = scm.ortho.load(client, path, name)
    else:
        print("[+] Creating...", flush=True)
        cohen = scm.scMPRA_data.from_parquet(str(path / "retina_single_counting_u6.scmpra"))
        cohen.set_negative_controls(["wt_1", "wt_2"])
        cohen.set_reference_cell("Rod")
        cohen.ortho_filter()
        cohen.set_consider_missing(True)

        primordial = scm.ortho()
        primordial.fit_by_cell_type_models(
            client=client, dat=cohen,
            fit_resources={"GPU": 1},
            pre_fit_hook=spin_up_gpu)
        primordial.extract_params(client)
        primordial.save(path, name, client=client)
        print("[+] Done.", flush=True)
finally:
    for jid in gpu_job_ids:
        subprocess.run(['scancel', jid], capture_output=True)
    client.close()
    cluster.close()
