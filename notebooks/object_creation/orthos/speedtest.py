import time
from datetime import timedelta

_start_time = time.time()

def timeprint(*args, **kwargs):
    """Print with timestamp since program start."""
    elapsed = time.time() - _start_time
    prefix = f"[+{timedelta(seconds=elapsed)}] "
    print(prefix, *args, **kwargs)


import scMPRAforge as scm


timeprint("imports done")


data_root="/home/mcn26/project_pi_skr2/shared/tabula_data"
path=f"{data_root}/cohen"
name="ortho_primordial"

cohen=scm.scMPRA_data.from_parquet(f"{path}/retina_single_counting_u6.scmpra")

cohen.set_negative_controls(["wt_1","wt_2"])
cohen.set_reference_cell("Rod")
cohen.ortho_filter()



"""
data=data.data
    levels=data[split].unique()

    mats_futures = {
        t: client.submit(
            _smart_matrix,
            data=data[data[split]==t],
            split=split
        )
        for t in levels
    }
"""
data=cohen.data[cohen.data["cell_type"]=="reference"]
mats=scm.smart_matrix(data=data,split="cell_type")

timeprint("matricies created")

print(mats)

import faulthandler, threading, time, os
import sys
from pathlib import Path

if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <dirname> <sparse/predense>")
        sys.exit(1)

if sys.argv[2]=="sparse":
        pass
elif sys.argv[2]=="predense":
        mats["nb_regressors"]=mats["nb_regressors"].sparse.to_dense()
        timeprint("finished matrix densification")


dirname = sys.argv[1]
path = Path(dirname)
path.mkdir(parents=True, exist_ok=False)

def start_rotating_stacks(interval=60, prefix=f"{dirname}/tracebacks"):
    faulthandler.enable()
    def _worker():
        while True:
            ts = time.strftime("%Y%m%d-%H%M%S")
            path = f"{prefix}_{ts}.log"
            with open(path, "w") as f:
                print(f"STACK DUMP @ {ts} (pid={os.getpid()})", file=f)
                import threading
                for t in threading.enumerate():
                    nid = getattr(t, "native_id", None)
                    print(f"  thread name={t.name!r} ident={t.ident} native_id={nid}", file=f)
                print("-"*80, file=f)
                faulthandler.dump_traceback(file=f, all_threads=True)
            time.sleep(interval)
    threading.Thread(target=_worker, name="faulthandler-rotator", daemon=True).start()

start_rotating_stacks()



model=scm.tensorzinb_fit(mats,"reference")

timeprint("model fitting done")

print(model)




