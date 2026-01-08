#!/usr/bin/env python
# coding: utf-8

# In[1]:


#imports

#from tensorzinb.tensorzinb import TensorZINB
import scMPRAforge as scm


# In[2]:


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
    cores=4,#cores per slurm job
    memory="512G",#memory per slurm job
    processes=1,#dask workers per slurm job
    job_extra_directives=["-p day", 
        f"--job-name=simclust_worker",
        f"--time=1-0",
        f"--output=worker_%j.out"]
)

cluster.scale(jobs=1)

client = Client(cluster,
        timeout=f"{5*60}s",   # Client <-> scheduler timeout 
        heartbeat_interval="20s"  # Worker heartbeat interval
    )

#from dask.distributed import Client, LocalCluster
#cluster=LocalCluster(memory_limit='8GB')
#client = Client(cluster)


# In[4]:


client.dashboard_link


# # Describe with Ortho
# 

# In[5]:


data_root="/home/sxl6/project_pi_skr2/sxl6/tabula_data/seelig"
path= data_root
name="ortho_test_seelig"

seelig=scm.scMPRA_data.from_tsv(f"{data_root}/susanna_seelig_counts_grouped.txt")


# In[6]:


seelig.set_negative_controls(["AACGCCCTCCACGGATGGGCCGGCCAATAAGAAGCGTTAGCGGACTCATGCGTTACGCGCCTCCGAGTTATGGGGGGGGAGGCGCGTATCTCGTGGAGAAGAAGCGATGTAACGCTTGGGCGATAAGCTTATAAGGAAGATATTT",
    "CCCTCGGAGTTAATAAGATACGCGGATCGATATCGGCTTGAAGAAGCGTATCTTATCTTCAGATGGGGATGTCGCGCATCCACCCAGTGGGCACCGCCGCTATAGAAGGGTGATAACGCTTCTCAGCCTTCAGGCTCTGGGTCTT"])
seelig.set_reference_cell("HEPG2")
seelig.ortho_filter()


# In[7]:


primordial=scm.ortho()
primordial.criss_cross(client=client,
                       dat=seelig)
primordial.extract_params(client)
primordial.save(path,name)


# In[ ]:


print("finished!")


# In[ ]:


client.close()
cluster.close()

