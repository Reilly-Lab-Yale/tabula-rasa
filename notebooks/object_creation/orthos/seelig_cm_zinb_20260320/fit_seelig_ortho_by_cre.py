#!/usr/bin/env python
"""
Fits only the by-CRE half of seelig_ortho (CPU, no GPU needed).
Saves to seelig_ortho_20260320_by_cre/.
Run in parallel with fit_seelig_ortho_by_cell_type.py,
then combine with merge_seelig_ortho.py.
"""
import signal
import scMPRAforge as scm
from dask_jobqueue import SLURMCluster
from dask.distributed import Client
from pathlib import Path

signal.signal(signal.SIGTERM, lambda s, f: (_ for _ in ()).throw(KeyboardInterrupt()))

data_root = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data")
path = data_root / "seelig"
name = "seelig_ortho_20260320_by_cre"
work_dir = "/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/notebooks/object_creation/orthos"

cluster = SLURMCluster(
    cores=1, memory="64G", processes=1,
    job_extra_directives=[
        "-p priority", "-A prio_skr2",
        "--job-name=seelig_cre_worker", "--time=12:00:00",
        f"--output={work_dir}/worker_cpu_%j.out"
    ]
)
cluster.scale(jobs=8)
client = Client(cluster, timeout=f"{5*60}s", heartbeat_interval="20s")
print(client.dashboard_link, flush=True)

print(f"[+] Waiting for 8 workers...", flush=True)
client.wait_for_workers(n_workers=8, timeout=600)
print(f"[+] All workers connected.", flush=True)

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
        primordial.fit_by_cre_models(client=client, dat=seelig)
        primordial.extract_params(client)
        primordial.save(path, name, client=client)
        print("[+] Done.", flush=True)
finally:
    client.close()
    cluster.close()
