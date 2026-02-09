import scMPRAforge as scm
import dask
from dask_jobqueue import SLURMCluster
from dask.distributed import Client, LocalCluster
from pathlib import Path
import time

def main():
    cluster=LocalCluster()
    client=Client(cluster)
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
