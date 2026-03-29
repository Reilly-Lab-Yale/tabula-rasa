#!/usr/bin/env python
"""
Fits only the by-cell-type half of seelig_ortho.
Saves to seelig_ortho_20260320_by_cell_type/.
Run in parallel with fit_seelig_ortho_by_cre.py,
then combine with merge_seelig_ortho.py.

GPU workers are submitted just-in-time via pre_fit_hook, after all
setup (consider_missing expansion, design matrices, MoM init) has
completed on CPU workers.  This avoids idle-GPU preemption on
priority_gpu, which cancels jobs not using the GPU after ~1 hour.
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
path = data_root / "seelig"
name = "seelig_ortho_20260320_by_cell_type"
work_dir = "/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/notebooks/object_creation/orthos"

# CPU workers for setup (design matrices, MoM init, etc.)
cluster = SLURMCluster(
    cores=1, memory="64G", processes=1,
    job_extra_directives=[
        "-p priority", "-A prio_skr2",
        "--job-name=seelig_ct_worker", "--time=12:00:00",
        f"--output={work_dir}/worker_cpu_%j.out"
    ]
)
cluster.scale(jobs=4)
client = Client(cluster, timeout=f"{5*60}s", heartbeat_interval="20s")
print(client.dashboard_link, flush=True)

scheduler_addr = cluster.scheduler_address
print(f"[+] Scheduler: {scheduler_addr}", flush=True)

# Wait for CPU workers only
print(f"[+] Waiting for 4 CPU workers...", flush=True)
client.wait_for_workers(n_workers=4, timeout=600)
print(f"[+] CPU workers connected.", flush=True)

gpu_script = f"""#!/bin/bash
#SBATCH -p priority_gpu
#SBATCH -A prio_skr2
#SBATCH --job-name=seelig_ct_gpu_worker
#SBATCH -t 1-00:00:00
#SBATCH -c 1
#SBATCH --mem=64G
#SBATCH --gres=gpu:h200:1
#SBATCH --output={work_dir}/worker_gpu_%j.out

# Set CUDA and conda lib paths explicitly — avoids mixing module/conda systems.
# Equivalent to: module load CUDA/12.6.0 + conda activate tz lib path.
CUDA_ROOT=/apps/software/2024a/software/CUDA/12.6.0
export LD_LIBRARY_PATH=${{CUDA_ROOT}}/targets/x86_64-linux/lib:${{CUDA_ROOT}}/extras/CUPTI/lib64:/home/mcn26/.conda/envs/tz/lib:$LD_LIBRARY_PATH
export PYTHONPATH={work_dir}:$PYTHONPATH
/home/mcn26/.conda/envs/tz/bin/dask-worker {scheduler_addr} --resources GPU=1 --nthreads 1 --memory-limit 60GiB
"""

gpu_job_ids = []

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
    print(f"[+] Waiting for GPU workers to connect...", flush=True)
    client.wait_for_workers(n_workers=6, timeout=1200)
    print(f"[+] All 6 workers connected. Submitting fits.", flush=True)

try:
    if (path / name).is_dir():
        print("[+] Partial result found. Loading...", flush=True)
        primordial = scm.ortho.load(client, path, name)
    else:
        print("[+] Creating...", flush=True)

        seelig = scm.scMPRA_data.from_tsv(str(path / "seelig_scmpra_umiwise.tsv.gz"))
        seelig.set_negative_controls([
            "AACGCCCTCCACGGATGGGCCGGCCAATAAGAAGCGTTAGCGGACTCATGCGTTACGCGCCTCCGAGTTATGGGGGGGGAGGCGCGTATCTCGTGGAGAAGAAGCGATGTAACGCTTGGGCGATAAGCTTATAAGGAAGATATTT",
            "CCCTCGGAGTTAATAAGATACGCGGATCGATATCGGCTTGAAGAAGCGTATCTTATCTTCAGATGGGGATGTCGCGCATCCACCCAGTGGGCACCGCCGCTATAGAAGGGTGATAACGCTTCTCAGCCTTCAGGCTCTGGGTCTT"
        ])
        seelig.set_reference_cell("HepG2")
        seelig.ortho_filter()
        seelig.set_consider_missing(enabled=True)

        primordial = scm.ortho()
        primordial.fit_by_cell_type_models(
            client=client, dat=seelig,
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
