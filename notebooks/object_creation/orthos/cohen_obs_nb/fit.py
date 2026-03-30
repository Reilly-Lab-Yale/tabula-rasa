#!/usr/bin/env python
"""
Cohen NB-only ortho (no consider_missing).
Counterpart to cohen_obs_zinb_20260204 for NB vs ZINB comparison.
"""
import signal
import scMPRAforge as scm
from dask_jobqueue import SLURMCluster
from dask.distributed import Client
from pathlib import Path

signal.signal(signal.SIGTERM, lambda s, f: (_ for _ in ()).throw(KeyboardInterrupt()))

data_root = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data")
path = data_root / "cohen"
name = "cohen_obs_nb_20260329"
work_dir = str(Path(__file__).resolve().parent)

cluster = SLURMCluster(
    cores=1, memory="32G", processes=1,
    job_extra_directives=[
        "-p priority", "-A prio_skr2",
        "--job-name=cohen_obs_nb_w", "--time=12:00:00",
        f"--output={work_dir}/worker_%j.out"
    ]
)
cluster.scale(jobs=4)
client = Client(cluster, timeout=f"{5*60}s", heartbeat_interval="20s")
print(client.dashboard_link, flush=True)

print("[+] Waiting for 4 workers...", flush=True)
client.wait_for_workers(n_workers=4, timeout=600)
print("[+] All workers connected.", flush=True)

try:
    if (path / name).is_dir():
        print("[+] Partial result found. Loading...", flush=True)
        primordial = scm.ortho.load(client, path, name)
    else:
        print("[+] Creating...", flush=True)
        cohen = scm.scMPRA_data.from_tsv(str(path / "retina_single_counting_u6.tsv"))
        cohen.set_negative_controls(["wt_1", "wt_2"])
        cohen.set_reference_cell("Rod")
        cohen.ortho_filter()

        primordial = scm.ortho()
        primordial.criss_cross(client=client, dat=cohen, nb_only=True)
        primordial.extract_params(client)
        primordial.save(path, name, client=client)
        print("[+] Done.", flush=True)
finally:
    client.close()
    cluster.close()
