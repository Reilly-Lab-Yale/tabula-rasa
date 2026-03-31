#!/usr/bin/env python
"""
Cohen ZINB by-CRE fit (consider_missing=True, CPU).
New condition for NB vs ZINB comparison — tests the "ignore reporter" scenario.
Run in parallel with fit_by_cell_type.py, then combine with merge.py.
"""
import signal
import scMPRAforge as scm
from dask_jobqueue import SLURMCluster
from dask.distributed import Client
from pathlib import Path

signal.signal(signal.SIGTERM, lambda s, f: (_ for _ in ()).throw(KeyboardInterrupt()))

data_root = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data")
path = data_root / "cohen"
name = "cohen_cm_zinb_20260329_by_cre"
work_dir = str(Path(__file__).resolve().parent)

cluster = SLURMCluster(
    cores=1, memory="64G", processes=1,
    job_extra_directives=[
        "-p priority", "-A prio_skr2",
        "--job-name=cohen_cm_zi_cre_w", "--time=12:00:00",
        f"--output={work_dir}/worker_cpu_%j.out"
    ]
)
cluster.scale(jobs=8)
client = Client(cluster, timeout=f"{5*60}s", heartbeat_interval="20s")
print(client.dashboard_link, flush=True)

print("[+] Waiting for 8 workers...", flush=True)
client.wait_for_workers(n_workers=8, timeout=600)
print("[+] All workers connected.", flush=True)

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
        # The memory cap estimator assumes dense string-per-row representation and
        # doesn't account for sparse encoding of umis_mpra_bc (99.8% zeros).
        # Bypass it — actual memory usage will be far lower than the estimate.
        cohen.consider_missing_max_memory_gb = None

        primordial = scm.ortho()
        primordial.fit_by_cre_models(client=client, dat=cohen)
        primordial.extract_params(client)
        primordial.save(path, name, client=client)
        print("[+] Done.", flush=True)
finally:
    client.close()
    cluster.close()
