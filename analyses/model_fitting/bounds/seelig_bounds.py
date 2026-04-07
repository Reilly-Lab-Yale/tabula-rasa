#!/usr/bin/env python
"""Extract Bounds object from seelig_ortho_20260320 and save preset."""
import sys
sys.path.insert(0, "/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa")

import scMPRAforge as scm
from dask.distributed import Client, LocalCluster
from pathlib import Path

data_root = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data")

cluster = LocalCluster(n_workers=4, memory_limit="6GB", processes=False)
client = Client(cluster)
print(client.dashboard_link, flush=True)

try:
    primordial = scm.ortho.load(client, path=data_root / "seelig", name="seelig_ortho_20260320")
    seelig_bounds = scm.Bounds.from_ortho(primordial)
    out = Path("/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/scMPRAforge/presets/seelig_bounds.tgz")
    seelig_bounds.to_tgz(str(out))
    print(f"[+] Saved to {out}", flush=True)
finally:
    client.close()
    cluster.close()
