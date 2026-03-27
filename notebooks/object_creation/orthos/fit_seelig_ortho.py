#!/usr/bin/env python
import signal
import subprocess
import tempfile
import os
import scMPRAforge as scm
from dask_jobqueue import SLURMCluster
from dask.distributed import Client
from pathlib import Path

# Graceful shutdown on scancel so client/cluster close cleanly
signal.signal(signal.SIGTERM, lambda s, f: (_ for _ in ()).throw(KeyboardInterrupt()))

data_root = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data")
path = data_root / "seelig"
name = "seelig_ortho_20260320"
work_dir = "/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/notebooks/object_creation/orthos"

# Single-threaded CPU workers for by-CRE models (~17 min, 1345 models)
cluster = SLURMCluster(
    cores=1,
    memory="64G",
    processes=1,
    job_extra_directives=[
        "-p priority",
        "-A prio_skr2",
        "--job-name=seelig_ortho_worker",
        "--time=12:00:00",
        f"--output={work_dir}/worker_cpu_%j.out"
    ]
)
cluster.scale(jobs=8)
client = Client(cluster, timeout=f"{5*60}s", heartbeat_interval="20s")
print(client.dashboard_link, flush=True)

# Submit 2 GPU workers (one per cell-type model) pointing at the same scheduler.
# Uses Dask resource token GPU=1 so criss_cross(gpu=True) routes cell-type fits here.
scheduler_addr = cluster.scheduler_address
print(f"[+] Scheduler: {scheduler_addr}", flush=True)

gpu_script = f"""#!/bin/bash
#SBATCH -p priority_gpu
#SBATCH -A prio_skr2
#SBATCH --job-name=seelig_gpu_worker
#SBATCH -t 1-00:00:00
#SBATCH -c 1
#SBATCH --mem=64G
#SBATCH --gres=gpu:rtx_5000_ada:1
#SBATCH --output={work_dir}/worker_gpu_%j.out

# Use full path to avoid conda activation issues on GPU nodes.
# PYTHONPATH so the nanny-forked worker process can import scMPRAforge.
export PYTHONPATH={work_dir}:$PYTHONPATH
/home/mcn26/.conda/envs/tz/bin/dask-worker {scheduler_addr} --resources GPU=1 --nthreads 1 --memory-limit 60GiB
"""

gpu_job_ids = []
for _ in range(2):
    with tempfile.NamedTemporaryFile(mode='w', suffix='.sh', delete=False, dir='/tmp') as f:
        f.write(gpu_script)
        tmp = f.name
    result = subprocess.run(['sbatch', tmp], capture_output=True, text=True)
    os.unlink(tmp)
    job_id = result.stdout.strip().split()[-1]
    gpu_job_ids.append(job_id)
    print(f"[+] Submitted GPU worker: {job_id}", flush=True)

print(f"[+] Waiting for 10 workers (8 CPU + 2 GPU)...", flush=True)
client.wait_for_workers(n_workers=10, timeout=600)
print(f"[+] All workers connected.", flush=True)

try:
    if (path / name).is_dir():
        print("[+] Model found. Loading...", flush=True)
        primordial = scm.ortho.load(client, path, name)
    else:
        print("[+] Model not found. Creating...", flush=True)

        seelig = scm.scMPRA_data.from_tsv(str(path / "seelig_scmpra_umiwise.tsv.gz"))
        seelig.set_negative_controls([
            "AACGCCCTCCACGGATGGGCCGGCCAATAAGAAGCGTTAGCGGACTCATGCGTTACGCGCCTCCGAGTTATGGGGGGGGAGGCGCGTATCTCGTGGAGAAGAAGCGATGTAACGCTTGGGCGATAAGCTTATAAGGAAGATATTT",
            "CCCTCGGAGTTAATAAGATACGCGGATCGATATCGGCTTGAAGAAGCGTATCTTATCTTCAGATGGGGATGTCGCGCATCCACCCAGTGGGCACCGCCGCTATAGAAGGGTGATAACGCTTCTCAGCCTTCAGGCTCTGGGTCTT"
        ])
        seelig.set_reference_cell("HepG2")
        seelig.ortho_filter()
        # No transfection reporter — treat all unobserved (cell, CRE) combos as true zeroes
        seelig.set_consider_missing(enabled=True)

        primordial = scm.ortho()
        # gpu=True: cell-type model fits routed to GPU workers via resources={"GPU": 1}
        primordial.criss_cross(client=client, dat=seelig, gpu=True)
        primordial.extract_params(client)
        primordial.save(path, name)
        print("[+] Done.", flush=True)
finally:
    # Cancel GPU workers explicitly before closing cluster
    for jid in gpu_job_ids:
        subprocess.run(['scancel', jid], capture_output=True)
    client.close()
    cluster.close()
