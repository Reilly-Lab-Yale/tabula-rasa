from typing import Dict
from dataclasses import dataclass, replace
import numpy as np
import pandas as pd
import json
import tarfile
import tempfile
import os

from scipy.stats import nbinom

import matplotlib.pyplot as plt
import seaborn as sns


@dataclass()
class Bounds:
    """
    Describes the boundaries of an experiment.
    Specifically the UMI portion.

    Bounds objects just average zero inflation parameters, as there
    is no reason to complicate simulations with multiple values.
    (If multiple ZI are expected for some reason, just adjust geometrically
    or run multiple simulation batches.)

    by_cell_type_theta and by_cell_type_zi are probably 
    more accurate for downstream simulation than their 
    by_cre counterparts, assuming that there are more 
    cres than cell types, which should be the case for most
    datasets.
    """

    def plot_transfection(self):
        self.update_transfection_params()
        #scale the plot to get 99% of the probability density
        k_max = int(np.ceil(nbinom.ppf(0.999, self.r, self.p)))
        k = np.arange(0, k_max + 1)
        pmf_nb = nbinom.pmf(k, self.r, self.p)
        
        fig, ax = plt.subplots()
        ax.plot(k, pmf_nb, marker='o', linestyle='-')#, label=f'NB fit (μ={mu_nb:.2f}, α={alpha_nb:.3f})'
        ax.set_xlabel('Unique MPRA barcodes per cell')
        ax.set_ylabel('Density')
        ax.legend()
        plt.tight_layout()
        plt.show()

    metadata:dict=None

    min_mpra_umi:float=None
    max_mpra_umi:float=None
    by_cre_theta:float=None
    by_cell_type_theta:float=None
    by_cre_zi:float=None
    by_cell_type_zi:float=None

    num_cres:int=None
    cells_per_cell_type:dict=None

    transfection_nb_mu:float=None
    transfection_nb_alpha:float=None

    def get_effective_moi(self):
        return self.transfection_nb_mu

    def set_effective_moi(self,moi):
        self.transfection_nb_mu=moi
        self.update_transfection_params()
    
    def update_transfection_params(self):
        """
        Adds an alternate parametrization of the NB model of
        transfection to the object. 
        """
        self.r = 1.0 / self.transfection_nb_alpha             # "size"
        self.p = self.r / (self.r + self.transfection_nb_mu)  # success prob

    def copy(self, **kwargs):
        return replace(self, **kwargs)

    def to_tgz(self, out_file):
        with tempfile.TemporaryDirectory() as tmpdir:
            # dump members as parquet files
            for name, value in self.__dict__.items():
                path = os.path.join(tmpdir, f"{name}.parquet")
                if isinstance(value, pd.DataFrame):
                    value.to_parquet(path, engine="pyarrow", index=True)
                elif isinstance(value, pd.Series):
                    value.to_frame(name).to_parquet(path, engine="pyarrow", index=True)
                else:
                    pd.DataFrame({name: [value]}).to_parquet(path, engine="pyarrow", index=False)

            # pack into a tgz
            with tarfile.open(out_file, "w:gz") as tar:
                tar.add(tmpdir, arcname="")

    @classmethod
    def from_tgz(cls, in_file):
        ret=Bounds()
        with tempfile.TemporaryDirectory() as tmpdir:
            # extract archive
            with tarfile.open(in_file, "r:gz") as tar:
                tar.extractall(tmpdir)

            # load parquet members back
            for fname in os.listdir(tmpdir):
                if fname.endswith(".parquet"):
                    name = fname[:-8]
                    path = os.path.join(tmpdir, fname)
                    df = pd.read_parquet(path, engine="pyarrow")
                    if df.shape == (1, 1):
                        val = df.iloc[0, 0]
                    elif df.shape[1] == 1:
                        val = df.iloc[:, 0]
                    else:
                        val = df
                    setattr(ret, name, val)
        ret.update_transfection_params()
        return ret
    
    @classmethod
    def from_ortho(cls,inp):
        """
        Takes an ortho object and abstracts out its bounds.

        This is an aggregation function and will hang if ortho is not done fitting yet. 
        
        This function requires that the ortho still have its training data 
        so we can extract things like "number of cells per cell-type" and 
        "number of MPRA barcodes per cell".

        Note that this function is totally replicate-agnostic. It averages estimated zero 
        inflation across replicates. 
        """
        
        ret=cls()

        #get the min & max nb parameters across all models
        mins=[]
        maxes=[]
        
        #for each model class
        for var in ["by_cell_type_parameters","by_cre_parameters"]:
            ## nb min & max ##
            for key in getattr(inp,var).nb:
                current=getattr(inp,var).nb[key].result()
                maxes.append(float(current.max()))
                mins.append(float(current.min()))
            
            thetas=[]
            for key in getattr(inp,var).theta:
                current=getattr(inp,var).theta[key].result()
                thetas.append(current)
            
            ## theta means ##
            #we could munge the strings & use setattr but i think this is more readable
            if var=="by_cell_type_parameters":
                ret.by_cell_type_theta=np.mean(thetas)
            elif var=="by_cre_parameters":
                ret.by_cre_theta=np.mean(thetas)
            
            
            ## zero inflation ##
            zis=[]
            for key in getattr(inp,var).zi:
                current=getattr(inp,var).zi[key].result()
                current=current.rename({'zi':key},axis=1)
                zis.append(current)
            zis=pd.concat(zis,axis=1)
            if var=="by_cell_type_parameters":
                ret.by_cell_type_zi=zis.mean(axis=1).mean()
            elif var=="by_cre_parameters":
                ret.by_cre_zi=zis.mean(axis=1).mean()

        
        #min & max nb and add to the return object
        ret.min_mpra_umi=min(mins)
        ret.max_mpra_umi=max(maxes)

        ret.theta=np.mean([ret.by_cell_type_theta,ret.by_cre_theta])

        ## Parameters from training data ##
        # transfection, as proxied by number of MPRA barcodes detected per cell
        tfection_params=inp.training_data.describe_transfection()
        ret.transfection_nb_mu=tfection_params["nb_mu"]
        ret.transfection_nb_alpha=tfection_params["nb_alpha"]

        #cells per cell type
        ret.cells_per_cell_type=inp.training_data.data.groupby("cell_type")["cell_bc"].nunique()

        #number of CREs
        ret.num_cres=inp.training_data.data["cre_id"].nunique()
        
        ret.update_transfection_params()
        return ret

#class wet_bounds
#   read_depth

from pathlib import Path
working_dir = Path(__file__).resolve().parent

SHENDURE_BOUNDS=Bounds.from_tgz(working_dir/"presets/shendure_bounds.tgz")