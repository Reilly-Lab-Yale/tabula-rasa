import scMPRAforge as scm
import dask
from dask_jobqueue import SLURMCluster
from dask.distributed import Client, LocalCluster, performance_report
from pathlib import Path
import time
import os

def main():
    with dask.config.set({"distributed.worker.resources.FIT":1,"CELL_DESIGN":1}):
        cluster=LocalCluster()
        client=Client(cluster)
        
        job_id = os.environ.get('SLURM_JOB_ID', 'local')
        report_filename = f"dask-report-job-{job_id}.html"
        
        with performance_report(filename=report_filename):
            data_root=Path("/nfs/roberts/project/pi_skr2/shared/tabula_data")
            print(client.dashboard_link,flush=True)
            sim=scm.de_novo_simulation(location=data_root,
                                       name="twothird_pow_sim_2026-02-05",
                                       client=client)
            print("LOADED",flush=True)
            sim.fit_orthos(serial_orthos=False,direction="by_cell_type")
            print("SUBMITTED",flush=True)
            sim.save()
            print("DONE",flush=True)
            time.sleep(60)
            client.close()
            cluster.close()
            time.sleep(60)

if __name__ == "__main__":
    main()
