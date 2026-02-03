#!/usr/bin/env python
# coding: utf-8

# In[2]:


#imports

#from tensorzinb.tensorzinb import TensorZINB
import scMPRAforge as scm
from scMPRAforge.core import _smart_matrix, _mom_from_training_data, _matricies_to_order, _tensorzinb_fit


# In[3]:


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
    cores=2,#cores per slurm job
    memory="512G",#memory per slurm job
    processes=1,#dask workers per slurm job
    job_extra_directives=["-p week", 
        f"--job-name=seelig_fit_worker",
        f"--time=4-00:00:00",
        f"--output=worker_%j.out"]
)

cluster.scale(jobs=2)

client = Client(cluster,
        timeout=f"{5*60}s",   # Client <-> scheduler timeout 
        heartbeat_interval="20s"  # Worker heartbeat interval
    )

#from dask.distributed import Client, LocalCluster
#cluster=LocalCluster(memory_limit='8GB')
#client = Client(cluster)


# In[7]:


print(client.dashboard_link, flush=True)


# # Describe with Ortho
# 

# In[ ]:


data_root="/nfs/roberts/project/pi_skr2/shared/tabula_data"
path=f"{data_root}/seelig"
name="ortho_seelig_v1"

def load_and_preprocess_data(data_root):
    """Load, preprocess the data and return processed data object."""
    seelig = scm.scMPRA_data.from_tsv(f"{data_root}/seelig/seelig_counts_grouped.txt")
    seelig.set_negative_controls(["AACGCCCTCCACGGATGGGCCGGCCAATAAGAAGCGTTAGCGGACTCATGCGTTACGCGCCTCCGAGTTATGGGGGGGGAGGCGCGTATCTCGTGGAGAAGAAGCGATGTAACGCTTGGGCGATAAGCTTATAAGGAAGATATTT",
    "CCCTCGGAGTTAATAAGATACGCGGATCGATATCGGCTTGAAGAAGCGTATCTTATCTTCAGATGGGGATGTCGCGCATCCACCCAGTGGGCACCGCCGCTATAGAAGGGTGATAACGCTTCTCAGCCTTCAGGCTCTGGGTCTT"])
    seelig.set_reference_cell("HEPG2")
    seelig.ortho_filter()
    return seelig

import os
if os.path.isdir(path+"/"+name):
    print("[+] Model found. Loading...")
    primordial=scm.ortho.load(client,path,name)
    shendure=primordial.training_data
else:
    print("[+] Model not found. Creating...")

    #load data
    seelig = load_and_preprocess_data(data_root)
    print('data loaded', flush=True)

    primordial=scm.ortho()
    primordial.criss_cross(client=client,
                       dat=seelig)
    primordial.extract_params(client)
    primordial.save(path,name)


# In[ ]:


# def single_model_fit(data, split, design_only=False):
#     data = data.data
#     levels=data[split].unique()
#     t = levels[0]

#     #smart matrix
#     t_future = client.submit(
#             _smart_matrix,
#             data=data[data[split]==t],
#             split=split
#         )
#     if design_only:
#         return {t:t_future}
#     print('smart matrix done', flush=True)
#     # mom for training data, matrices to order
#     if split=="cell_type":
#         init_method="pass"
#         init_vals=client.submit(_mom_from_training_data, 
#                 data=data,
#                 split="cell_type",
#                 subset=t,
#                 indicies=client.submit(_matricies_to_order, matricies=t_future)
#                 )
#         print('MoM done', flush=True)

#     else:
#             init_method="nb"
#             init_vals= None 

#     #tensorzinb fit
#     print('submited fitting', flush=True)
#     tzinb_futures = client.submit(
#                 _tensorzinb_fit,
#                 t_future,
#                 t,
#                 init_method=init_method,
#                 init_vals=init_vals
#             )

#     print('fitting done', flush=True)
#     return(scm.experiment_model(model={t:tzinb_futures},
#                                 split=split),
#             {t:t_future})



# In[ ]:


# client.close()
# cluster.close()


# In[ ]:




