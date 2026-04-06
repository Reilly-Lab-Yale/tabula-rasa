#!/usr/bin/env python
"""
Takeshi CM ZINB fit -- full consider_missing, phantom compressed.
"""
import sys, time, os, subprocess
sys.path.insert(0, '/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa-takeshi')

from pathlib import Path

SCRIPT_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new/takeshi")
OUT_ROOT = DATA_DIR
OUT_NAME = "takeshi_cm_zinb_phantom"
SLURM_JOB_ID = os.environ.get("SLURM_JOB_ID", "unknown")

NEGATIVE_CONTROLS = [
    "1:113047064:NA:NA", "1:50938700:NA:NA", "1:65979308:NA:NA",
    "12:117944158:NA:NA", "12:69973147:NA:NA", "14:34259631:NA:NA",
    "18:73114109:NA:NA", "2:150206278:NA:NA", "2:153327926:NA:NA",
    "2:4341678:NA:NA", "2:79600326:NA:NA", "3:384644:NA:NA",
    "3:76838770:NA:NA", "5:136396708:NA:NA", "6:139910638:NA:NA",
    "7:24696737:NA:NA", "7:96368776:NA:NA", "9:121708667:NA:NA",
]


def main():
    from dask.distributed import Client, LocalCluster, performance_report
    import scMPRAforge.core as scm

    print("Setting up Dask local cluster...", flush=True)
    cluster = LocalCluster(n_workers=1, threads_per_worker=2,
                           processes=False, memory_limit="200GB")
    client = Client(cluster)
    print(f"Dashboard: {client.dashboard_link}", flush=True)

    print("Loading takeshi data...", flush=True)
    takeshi = scm.scMPRA_data.from_parquet(str(DATA_DIR / "takeshi_scmpra_umiwise.scmpra"))
    takeshi.set_negative_controls(NEGATIVE_CONTROLS)
    takeshi.set_reference_cell("HepG2")
    takeshi.ortho_filter()
    takeshi.set_consider_missing(True)
    print("Data loaded and filtered.", flush=True)

    ortho_obj = scm.ortho()

    report_path = str(SCRIPT_DIR / "dask_performance_report.html")
    t0 = time.time()
    with performance_report(filename=report_path):
        print("\nFitting by_cre (CM phantom)...", flush=True)
        t_cre = time.time()
        ortho_obj.fit_by_cre_models(client=client, dat=takeshi)
        print(f"by_cre done in {time.time()-t_cre:.1f}s", flush=True)

        print("\nFitting by_cell_type (CM phantom)...", flush=True)
        t_ct = time.time()
        ortho_obj.fit_by_cell_type_models(client=client, dat=takeshi)
        print(f"by_cell_type done in {time.time()-t_ct:.1f}s", flush=True)

    total_time = time.time() - t0
    print(f"\nTotal fit time: {total_time:.1f}s ({total_time/60:.1f}min)", flush=True)

    print("Extracting parameters...", flush=True)
    ortho_obj.extract_params(client)

    print(f"Saving to {OUT_ROOT / OUT_NAME}...", flush=True)
    ortho_obj.save(str(OUT_ROOT), OUT_NAME, client=client)
    print("Saved.", flush=True)

    client.close()
    cluster.close()

    print(f"\nTotal wall time: {total_time:.1f}s ({total_time/60:.1f}min)", flush=True)
    if SLURM_JOB_ID != "unknown":
        sacct_out = subprocess.run(
            ["sacct", "-j", SLURM_JOB_ID,
             "--format=JobID,State,Elapsed,MaxRSS,AveRSS,AllocCPUS,ReqMem",
             "--parsable2"],
            capture_output=True, text=True
        ).stdout
        print(f"\nsacct:\n{sacct_out}", flush=True)

        node = os.environ.get("SLURMD_NODENAME", "")
        if node:
            info = subprocess.run(
                ["scontrol", "show", "node", node],
                capture_output=True, text=True
            ).stdout
            for field in ["NodeName", "CPUTot", "RealMemory", "AvailableFeatures"]:
                for tok in info.split():
                    if tok.startswith(f"{field}="):
                        print(f"  {tok}", flush=True)

    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
