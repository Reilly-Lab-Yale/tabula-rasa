"""
Cohen obs NB fit -- reporter-informed zeros via coarse reporter table, NB-only.
Counterpart to cohen_obs_zinb_phantom for NB vs ZINB comparison.
"""
import sys, time, os, subprocess
sys.path.insert(0, '/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa-cohen-regen')

from pathlib import Path

SCRIPT_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data/cohen")
OUT_DIR = DATA_DIR / "cohen_obs_nb_phantom_20260401"
SLURM_JOB_ID = os.environ.get("SLURM_JOB_ID", "unknown")


def main():
    from dask.distributed import Client, LocalCluster, performance_report
    import scMPRAforge.core as scm

    print("Setting up Dask local cluster...", flush=True)
    cluster = LocalCluster(n_workers=1, threads_per_worker=2,
                           processes=False, memory_limit="60GB")
    client = Client(cluster)
    print(f"Dashboard: {client.dashboard_link}", flush=True)

    print("Loading Cohen data (TSV, nonzero only)...", flush=True)
    cohen = scm.scMPRA_data.from_tsv(str(DATA_DIR / "retina_single_counting_u6.tsv"))
    cohen.set_negative_controls(["wt_1", "wt_2"])
    cohen.set_reference_cell("Rod")
    cohen.ortho_filter()
    cohen.set_coarse_reporter(str(DATA_DIR / "unjoined" / "u6.tsv"))
    print("Data loaded and filtered.", flush=True)

    ortho_obj = scm.ortho()

    report_path = str(SCRIPT_DIR / "dask_performance_report.html")
    t0 = time.time()
    with performance_report(filename=report_path):
        print("\nFitting by_cre (obs phantom NB)...", flush=True)
        t_cre = time.time()
        ortho_obj.fit_by_cre_models(client=client, dat=cohen, phantom_compress=True, nb_only=True)
        print(f"by_cre done in {time.time()-t_cre:.1f}s", flush=True)

        print("\nFitting by_cell_type (obs phantom NB)...", flush=True)
        t_ct = time.time()
        ortho_obj.fit_by_cell_type_models(client=client, dat=cohen, phantom_compress=True, nb_only=True)
        print(f"by_cell_type done in {time.time()-t_ct:.1f}s", flush=True)

    total_time = time.time() - t0
    print(f"\nTotal fit time: {total_time:.1f}s ({total_time/60:.1f}min)", flush=True)

    print("Extracting parameters...", flush=True)
    ortho_obj.extract_params(client)

    print(f"Saving to {OUT_DIR}...", flush=True)
    ortho_obj.save(str(DATA_DIR), "cohen_obs_nb_phantom_20260401", client=client)
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
