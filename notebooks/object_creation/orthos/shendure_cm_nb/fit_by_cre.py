#!/usr/bin/env python
"""
Shendure NB-only by-CRE fit (consider_missing=True, CPU).
Counterpart to shendure_cm_zinb_20260320 for NB vs ZINB comparison.
Run in parallel with fit_by_cell_type.py, then combine with merge.py.

Matches the ZINB CM config (Attempt 11, job 2561074) for parity:
16 single-threaded workers × 96GB, CPU-only.
"""
import signal
import time
import scMPRAforge as scm
from dask_jobqueue import SLURMCluster
from dask.distributed import Client, performance_report
from pathlib import Path

signal.signal(signal.SIGTERM, lambda s, f: (_ for _ in ()).throw(KeyboardInterrupt()))

data_root = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data")
path = data_root / "shendure"
name = "shendure_cm_nb_20260329_by_cre"
work_dir = str(Path(__file__).resolve().parent)
report_path = Path(__file__).with_name(f"{name}_dask_performance_report.html")

# ── cluster ────────────────────────────────────────────────────────────────────
cluster = SLURMCluster(
    cores=1, memory="96G", processes=1,
    job_extra_directives=[
        "-p priority", "-A prio_skr2",
        "--job-name=shend_cm_nb_cre_w", "--time=2-23:40:00",
        f"--output={work_dir}/worker_cpu_%j.out"
    ]
)
cluster.scale(jobs=16)
client = Client(cluster, timeout=f"{5*60}s", heartbeat_interval="20s")
print(client.dashboard_link, flush=True)

print("[+] Waiting for 16 workers...", flush=True)
client.wait_for_workers(n_workers=16, timeout=600)
print("[+] All workers connected.", flush=True)

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
            primordial.fit_by_cre_models(client=client, dat=shendure, nb_only=True)
            primordial.extract_params(client)
            primordial.save(path, name, client=client)
            print("[+] Done.", flush=True)
finally:
    print("! Done, shutting down", flush=True)
    time.sleep(10)
    client.close()
    cluster.close()
    time.sleep(10)
