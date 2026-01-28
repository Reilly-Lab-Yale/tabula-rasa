#!/usr/bin/env python
# coding: utf-8

# In[10]:


#imports

#from tensorzinb.tensorzinb import TensorZINB
import scMPRAforge as scm
from scMPRAforge.core import _smart_matrix, _mom_from_training_data, _matricies_to_order, _tensorzinb_fit


# In[11]:


import pandas as pd
import numpy as np
import time
import pickle
from formulaic import Formula
import seaborn as sns
import matplotlib.pyplot as plt


# In[ ]:


#create dask cluster

from dask_jobqueue import SLURMCluster
from dask.distributed import Client

cluster=SLURMCluster(
    cores=8,#cores per slurm job
    memory="512G",#memory per slurm job
    processes=1,#dask workers per slurm job
    job_extra_directives=["-p week", 
        f"--job-name=cell_type_1_worker_job",
        f"--time=7-00:00:00",
        f"--output=worker_cell_type_1_%j.out"]
)

cluster.scale(jobs=1)

client = Client(cluster,
        timeout=f"{5*60}s",   # Client <-> scheduler timeout 
        heartbeat_interval="20s"  # Worker heartbeat interval
    )

index = 1

#from dask.distributed import Client, LocalCluster
#cluster=LocalCluster(memory_limit='8GB')
#client = Client(cluster)


# In[ ]:


print(client.dashboard_link, flush=True)


# # Describe with Ortho
# 

# In[ ]:


data_root="/nfs/roberts/project/pi_skr2/shared/tabula_data"

def load_and_preprocess_data(data_root):
    """Load, preprocess the data and return processed data object."""
    seelig = scm.scMPRA_data.from_tsv(f"{data_root}/seelig/seelig_counts_grouped.txt")
    seelig.set_negative_controls(["AACGCCCTCCACGGATGGGCCGGCCAATAAGAAGCGTTAGCGGACTCATGCGTTACGCGCCTCCGAGTTATGGGGGGGGAGGCGCGTATCTCGTGGAGAAGAAGCGATGTAACGCTTGGGCGATAAGCTTATAAGGAAGATATTT",
    "CCCTCGGAGTTAATAAGATACGCGGATCGATATCGGCTTGAAGAAGCGTATCTTATCTTCAGATGGGGATGTCGCGCATCCACCCAGTGGGCACCGCCGCTATAGAAGGGTGATAACGCTTCTCAGCCTTCAGGCTCTGGGTCTT"])
    seelig.set_reference_cell("HEPG2")
    seelig.ortho_filter()
    return seelig



seelig = load_and_preprocess_data(data_root)
print('data loaded', flush=True)


# In[ ]:


def single_model_fit(data, split, index, client):
    data = data.data
    levels=data[split].unique()
    t = levels[index]

    #smart matrix
    t_future = client.submit(
            _smart_matrix,
            data=data[data[split]==t],
            split=split
        )
    print('smart matrix done', flush=True)
    # mom for training data, matrices to order
    if split=="cell_type":
        init_method="pass"
        init_vals=client.submit(_mom_from_training_data, 
                data=data,
                split="cell_type",
                subset=t,
                indicies=client.submit(_matricies_to_order, matricies=t_future)
                )
        print('MoM done', flush=True)

    else:
            init_method="nb"
            init_vals= None 

    #tensorzinb fit
    print('submited fitting', flush=True)
    tzinb_futures = client.submit(
                _tensorzinb_fit,
                t_future,
                t,
                init_method=init_method,
                init_vals=init_vals
            )

    print('fitting done', flush=True)
    return(scm.experiment_model(model={t:tzinb_futures},
                                split=split),
            {t:t_future})



# In[ ]:


model, design = single_model_fit(seelig, 'cell_type', index, client)


# In[ ]:


model.save(f"{data_root}/seelig/ortho_test_seelig/7_daycell_type_1_model.pkl")


# In[ ]:


print("finished!", flush=True)


# In[ ]:


# client.close()
# cluster.close()


# In[ ]:




