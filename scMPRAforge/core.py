#        1         2         3         4         5         6         7         8
#2345678901234567890123456789012345678901234567890123456789012345678901234567890

#all the main functions.

#external imports
import seaborn as sns
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np

import logging
import time
import pickle
from pathlib import Path
import copy

import pyarrow as pa
import pyarrow.parquet as pq
import json

from scipy.stats import linregress
import scipy

import statsmodels.api as sm
import statsmodels.discrete.discrete_model as smd

import patsy
from tensorzinb.tensorzinb import TensorZINB
from formulaic import Formula

from dask.distributed import Client
from dask.distributed import Future

import dask.dataframe as dd
import dask.array as da

import os

from enum import Enum
from typing import List

from dataclasses import dataclass, replace
import json
import tarfile
import tempfile

#internal imports
from .utils import unimplemented
from .utils import bcs_to_lut
from .utils import undo_one_hot_encoding
from .utils import dict_wrap, dict_unwrap
logger = logging.getLogger("scMPRAforge")

MIN_PTS=3
PARTITION_SIZE_MB=50

#functions
@unimplemented
def always_unfinished():
    """tests unimplemented decorator."""
    pass


def helloworld():
    print("hello world!")


def table_type(column_names):
    """
    Arguments:
        column_names <pd.Index>
    Returns:
        <str>

    Returns the putative table type, one of:
    - mpra_umiwise
    - mpra_readwise
    - hypotheses
    - results
    - malformed

    Is kind with respect to extra columns & optional columns. 
    
    (We could extend to type-checking as well, but that seems a tad draconian / unpythonic.)
    """
    #Performs subset checking so that extra columns are allowed.
    #Matching multiple definitions however is NOT allowed. 
    
    #If columns match no definitions the table is malformed: so this is the default. 
    ret='malformed'
    
    #check mandatory columns of read-wise MPRA
    if {'cell_bc', 'rep_id', 'cre_id', 'cell_type', 'mpra_bc', 'mpra_umi', 'reads_mpra'}<=set(column_names):
        if ret=='malformed':
            ret='mpra_readwise'
        else:
            return 'malformed'
    
    #check mandatory columns of umi-wise full MPRA
    if {'cell_bc', 'rep_id', 'cre_id', 'cell_type', 'umis_mpra_bc'}<=set(column_names):
        if ret=='malformed':
            ret='mpra_umiwise'
        else:
            return 'malformed'
    
    #check mandatory columns of umi-wise flattened MPRA
    #unimplemented

    return ret



def simple_spread(cell_types:List[str],min:int,max:int,fineness:int=10):
    """
    Create a ground truth dataframe tiling all cell-types.
    with synthetic CREs at a variety of strengths.

    Returns a tuple of (,None) 
    (none will be replaced by hypothesis object)

    Useful for simulation and power calculations.
    (see readme for ground truth dataframe specification)

    min, max are the min & max MPRA UMI / cell values.

    
    
    Also returns a hypothesis object (UNIMPLEMENTED)
    This has the cartesian product of CREs.
    In general, same CRE name 
    """
  
    
    #to generally tile the space
    even_tiling=np.linspace(start=min, stop=max, num=fineness)
    
    even_names=["reference"]+[f"CRE_even_{i}" for i in range(1,len(even_tiling))]
    even_df=pd.DataFrame({"cre_id":even_names,"true_mean":even_tiling})
    #duplicate df for each cell-type
    even_df=pd.concat([even_df]*len(cell_types))
    even_df["cell_type"]=np.repeat(cell_types,len(even_tiling))
    
    #to create a large number of low-expression CREs
    fractional_tiling=1/(2**np.linspace(start=1,stop=fineness,num=fineness))*(max-min)+min
    fractional_tiling_names=[f"CRE_fractional_{i}" for i in range(0,len(fractional_tiling))]
    fractional_tiling_df=pd.DataFrame({"cre_id":fractional_tiling_names,"true_mean":fractional_tiling})
    #duplicate df for each cell-type
    fractional_tiling_df=pd.concat([fractional_tiling_df]*len(cell_types))
    fractional_tiling_df["cell_type"]=np.repeat(cell_types,len(fractional_tiling))
    

    #Then we create some "cell-type-specific CREs". 
    #these are simply "high in one cell-type and low in every other".
    #where "high" is a value from even tiling+1. 
    #and "low" =ln(high)

    #cell-type specific
    ct_specific=even_df
    ct_specific=ct_specific.rename({"true_mean":"high","cell_type":"expressed_cell_type"},axis=1)
    ct_specific["high"]=ct_specific["high"]+1
    ct_specific["low"]=np.log(ct_specific["high"])

    #ok, so we've created a high and low value for each CRE. 
    #"expressed_cell_type" will be highest cell-type
    #so let's rename accordingly
    ct_specific["cre_id"]=\
        ct_specific["cre_id"].astype(str)+\
        "_high_in_"+\
        ct_specific["expressed_cell_type"].astype(str)
    
    #let's repeat the DF across all cell-types
    #every CRE in every cell-type, after all...
    ct_specific=pd.concat([ct_specific]*len(cell_types))
    ct_specific["cell_type"]=np.repeat(cell_types,int(len(ct_specific)/len(cell_types)))
    #now let's assign values...
    ct_specific=ct_specific.assign(
        true_mean=np.where(ct_specific["cell_type"]==ct_specific["expressed_cell_type"],
                           ct_specific["high"],
                           ct_specific["low"])
    )
    #drop extranious columns
    ct_specific=ct_specific[["cre_id","cell_type","true_mean"]]

    final_ground_truth=pd.concat([even_df,fractional_tiling_df,ct_specific])
    
    hypothesis_set=None
    return (final_ground_truth,hypothesis_set)

@unimplemented
def load_hypothesis_set(filepath):
    """
    Arguments
        filepath <str>
    Returns
        <pd.DataFrame>

    Loads a hypothesis or hypothesis+results set from disc.
    """
    #load table...
    #assert table_type(table.columns)=="hypothesis" or table_type(table.columns)=="results"
    #return table
    pass

class simple_count:
    """
    This class stores information pertaining to a simple negative binomial model.
    It is low performance and NOT used for primary modeling of RNA-sequencing data.
    Instead, it us used for small, discrete modeling tasks whch need flexibility but not 
    performance. 
    """
    def __init__(self,data=None):
        if not data is None:
            self.from_data(data)

    def from_data(self,data):
        """
        Initializes the object from count (not frequency) data.
        Computes poisson and nb using statsmodels.
        """
        if type(data) != list:
            data=data.tolist()

        self.data=data
        
        # intercept-only design
        X = np.ones((len(data), 1))

        #poisson
        pois = smd.Poisson(data, X).fit(disp=False)
        if not pois.converged:
            logger.warning("Poisson model of transfection failed to converge.")
        self.mu_pois = np.exp(pois.params[0])      # intercept-only mean

        # NB2
        # Var = mu + alpha*mu^2, alpha is estimated
        nb = smd.NegativeBinomial(data, X).fit(disp=False)
        self.mu_nb = np.exp(nb.params[0])          # intercept-only mean
        self.alpha_nb = nb.params[-1]              # dispersion (NB2 alpha)

        if not nb.converged:
            logger.error("Negative binomial model of transfection failed to converge.")
            raise(RuntimeError)
        
        self.update_alt_nb_param()
        self.original_fit=True

    def update_alt_nb_param(self):
        self.r = 1.0 / self.alpha_nb
        self.p = self.r / (self.r + self.mu_nb)
    
    def draw_nb(self,size):
        """
        returns a 1d numpy vector of draws from the nb
        model of the object. 
        """
        return scipy.stats.nbinom.rvs(self.r,self.p,size=size)
    
    def plot(self, max_bins: int = 25, binwidth: int | None = None):
        """
        Plot a histogram of the data with *integer-aligned coarse bins* and overlay
        NB and Poisson fits aggregated to the same bins.

        Parameters
        ----------
        max_bins : int
            Target maximum number of bars shown (used only if `binwidth` is None).
            The method will choose an integer bin width w >= 1 so that the span of
            the data is displayed in ~max_bins bars.
        binwidth : int | None
            If provided, forces that integer bin width (w). If None, a width is
            computed from `max_bins`.
        """
        # Prepare data and model PMFs over integer support
        data = np.asarray(self.data, dtype=int)
        dmin = int(data.min())
        dmax = int(data.max())

        k = np.arange(dmin, dmax + 1, dtype=int)
        pmf_nb = scipy.stats.nbinom.pmf(k, self.r, self.p)
        pmf_pois = scipy.stats.poisson.pmf(k, self.mu_pois)

        # Choose integer-aligned coarse bins
        span = dmax - dmin + 1
        if binwidth is None:
            w = max(1, int(np.ceil(span / max_bins)))
        else:
            w = max(1, int(binwidth))

        # Half-integer edges so each integer falls cleanly into one bin
        edges = np.arange(dmin - 0.5, dmax + 0.5 + w, w)
        nbins = len(edges) - 1

        # Aggregate model PMFs to the same coarse bins 
        def _bin_sums_for_pmf(k_arr, pmf_arr, dmin_val, width, n_bins):
            xs = []
            ys = []
            for i in range(n_bins):
                a = dmin_val + i * width
                b = a + width - 1
                mask = (k_arr >= a) & (k_arr <= b)
                ys.append(pmf_arr[mask].sum())
                xs.append((a + b) / 2.0)  # plot at bin center
            return np.asarray(xs), np.asarray(ys)

        x_nb, y_nb = _bin_sums_for_pmf(k, pmf_nb, dmin, w, nbins)
        x_pois, y_pois = _bin_sums_for_pmf(k, pmf_pois, dmin, w, nbins)

        # Plot
        fig, ax = plt.subplots()
        sns.histplot(
            data,
            bins=edges,
            stat="probability",   # bar height is probability mass in bin
            # Use discrete binning only for unit-width bins; otherwise seaborn forces 1 bar per integer.
            discrete=(w == 1),
            alpha=0.3,
            shrink=0.9,
            ax=ax,
        )
        ax.plot(x_nb, y_nb, marker='o', linestyle='-', label=f'NB fit (μ={self.mu_nb:.2f}, α={self.alpha_nb:.3f})')
        ax.plot(x_pois, y_pois, linestyle='--', label=f'Poisson fit (μ={self.mu_pois:.2f})')

        ax.set_xlabel('Count')
        ax.set_ylabel('Probability')
        ax.legend()
        plt.tight_layout()
        plt.show()

    def adjust_nb(self,new_nb):
        self.original_fit=False
        self.nb=new_nb
        self.update_alt_nb_param()
    
    def to_dataframe(self) -> pd.DataFrame:
        """
        Represent the object as a single-row DataFrame.
        All non-hidden attributes are included.
        """
        d = {
            k: v
            for k, v in self.__dict__.items()
            if not k.startswith("_")
        }
        return pd.DataFrame([d])
    
    @classmethod
    def from_dataframe(cls, df: pd.DataFrame) -> "simple_count":
        """
        Recreate an object from a single-row DataFrame.
        """
        if len(df) != 1:
            raise ValueError("Expected DataFrame with exactly one row.")
        row = df.iloc[0].to_dict()

        obj = cls.__new__(cls)  # bypass __init__
        for k, v in row.items():
            setattr(obj, k, v)

        return obj

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

    metadata:dict=None

    preferred:str=None
    
    min_mpra_umi:float=None
    max_mpra_umi:float=None
    by_cre_theta:float=None
    by_cell_type_theta:float=None
    zi:float=None
    by_cre_zi:float=None
    by_cell_type_zi:float=None

    cells_per_cell_type:dict=None

    transfection_model:simple_count=None
    library_model:simple_count=None

    def get_effective_moi(self):
        self.transfection_model.update_alt_nb_param()
        return self.transfection_model.mu_nb

    def set_effective_moi(self,moi):
        self.transfection_model.mu_nb=moi
        self.transfection_model.update_alt_nb_param()
    
    def plot_transfection(self,*args,**kwargs):
        self.transfection_model.plot(*args,**kwargs)

    def copy(self, **kwargs):
        return replace(copy.deepcopy(self), **kwargs)

    def to_tgz(self, out_file):
        with tempfile.TemporaryDirectory() as tmpdir:
            # dump members as parquet files
            for name, value in self.__dict__.items():
                path = os.path.join(tmpdir, f"{name}.parquet")
                if isinstance(value, pd.DataFrame):
                    value.to_parquet(path, engine="pyarrow", index=True)
                if isinstance(value,simple_count):
                    value.to_dataframe().to_parquet(path, engine="pyarrow", index=True)
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
        #re-initalize 
        ret.transfection_model=simple_count.from_dataframe(ret.transfection_model)
        ret.library_model=simple_count.from_dataframe(ret.library_model)

        ret.transfection_model.update_alt_nb_param()
        ret.library_model.update_alt_nb_param()

        return ret
    
    @classmethod
    def from_ortho(cls,inp,preferred="by_cell_type"):
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
        ret.transfection_model=inp.training_data.describe_transfection()

        # library model
        ret.library_model=inp.training_data.describe_library()
        

        #cells per cell type
        ret.cells_per_cell_type=inp.training_data.data.groupby("cell_type")["cell_bc"].nunique()
        
        ret.preferred=preferred
        
        if preferred=="by_cell_type":
            ret.zi=ret.by_cell_type_zi
        elif preferred=="by_cre":
            ret.zi=ret.by_cre_zi
        else:
            assert False, "Unrecognized direction."

        return ret

from pathlib import Path
working_dir = Path(__file__).resolve().parent
SHENDURE_BOUNDS=Bounds.from_tgz(working_dir/"presets/shendure_bounds.tgz")

class scMPRA_data:
    """
    Wrapper around a pandas dataframe of MPRA data. 
    The primary purpose of the object is to record what operations have been performed on the data
    (Pandas does not support metadata)

    Could possibly replace with an anndata object.
    Alternatively. also allow pass-through of pandas operations & record them... 
    Alternatively, just implement a couple common operations (subsetting & friends) manually
    """
    def __init__(self):
        self.data=None
        self.table_type=None
        self.source=None
        self.metadata={}
        self.operations=[]
        self.negative_controls=[]
        self.reference_cell_type=None
    
    def flag_synthetic(self):
        self.metadata["synthetic"]=True
    
    def flag_emperical(self):
        self.metadata["synthetic"]=False
    
    def describe_transfection(self):
        """
        Returns a simple_count object describing the number of transfections per cell, 
        as proxied by number of unique MPRA barcodes per cell.

        Returns a simple_count object
        Drawing from one of the distributions (or a similar distribution shifted to a different MOI) 
        can be used to simulate transfection.
        """
        assert self.table_type=="mpra_umiwise"

        unique_mpra_barcodes_per_cell=self.data.groupby("cell_bc")["mpra_bc"].nunique()
        return simple_count(data=unique_mpra_barcodes_per_cell)
        
    def describe_library(self):
        """
        Returns a simple_count object describing the number of unique MPRA
        barcodes for each 
        """
        y=self.data.groupby("cre_id")["mpra_bc"].nunique()
        return simple_count(data=y)

    
    def set_negative_controls(self,negative_controls:list[str]):
        """
        Takes a list of CRE names that we consider to be negative controls and give them all the name "negative_control", lumping all their data together.
        """

        if self.negative_controls==[]:
            #User has set no negative controls before now. Back up cre_id information before we mutate it.
            self.data["cre_id_original"]=self.data["cre_id"]
        
        #flatten all labels of negative controls
        for control in negative_controls:
            self.data["cre_id"]=self.data["cre_id"].replace(control, "reference")

        #record which names we have flattened.
        #(Can be deduced from difference between cre_id and 
        #cre_id_original, but this is more convienient).
        self.negative_controls=self.negative_controls+negative_controls
    
    
    def set_reference_cell(self,reference_cell_type):
        assert self.reference_cell_type is None, "Already set reference cell type."
        
        self.reference_cell_type=reference_cell_type

        self.data["cell_type"]=self.data["cell_type"].replace(reference_cell_type, "reference")
    
    def copy(self, exclude=()):
        """Return a deepcopy of the object, optionally excluding fields."""
        cls = self.__class__
        result = cls.__new__(cls)
        memo = {}

        for k, v in self.__dict__.items():
            if k in exclude:
                setattr(result, k, None)  # or preserve original: self.__dict__[k]
            else:
                setattr(result, k, copy.deepcopy(v, memo))

        return result
    
    def total_umi(self):
        #the same cell barcode in two different replicates is NOT the same cell. 
        umis_per_cell=self.data.groupby(["cell_bc","rep_id"],as_index=False)["umis_mpra_bc"].sum()
        mask=umis_per_cell["umis_mpra_bc"]<1

        total_cells=len(umis_per_cell[["cell_bc","rep_id"]].value_counts())
        uniq_dropped=umis_per_cell[mask][["cell_bc","rep_id"]].value_counts()
        num_cells_to_drop=len(uniq_dropped)

        logger.info(f"Dropping {num_cells_to_drop} cells with no MPRA UMIs, leaving {total_cells-num_cells_to_drop}.")

        umis_per_cell=umis_per_cell[~mask]
        umis_per_cell["ln_cell_umis_mpra"]=np.log(umis_per_cell["umis_mpra_bc"])

        self.data=self.data.merge(umis_per_cell[["cell_bc","rep_id","ln_cell_umis_mpra"]],on=["cell_bc","rep_id"],how="right")

        self.operations.append(('total_umi',uniq_dropped))
    
    @classmethod
    def from_tsv(cls, filepath):
        """
        Returns a <scMPRA_data> object with data loaded from `filepath`.
        """
        tab=pd.read_csv(filepath,sep="\t")
        tabtype=table_type(tab.columns)
        
        assert tabtype=="mpra_readwise" or tabtype=="mpra_umiwise", "Malformed table."
        
        ret=cls()
        ret.data=tab
        ret.table_type=tabtype
        ret.source=filepath

        return ret
    
    
    @classmethod
    def from_parquet(cls,path):
        """
        Returns a <scMPRA_data> object with data loaded from `path`.
        Takes full path, /path/to/data.scmpra.
        """
        #create return object
        ret=cls()
        
        pa_data_table=pq.read_table(path)
        data=pa_data_table.to_pandas(types_mapper=pd.ArrowDtype)

        ret.data=data
        
        #extract parquet metadata (bytes->bytes)
        pa_metadata = pa_data_table.schema.metadata or {}
        #extract & decode the item with members
        meta_dict = json.loads(pa_metadata.get(b"scMPRA_data.members", b"{}").decode("utf-8"))

        # Restore all saved metadata members
        for k, v in meta_dict.items():
            setattr(ret, k, v)

        return ret

    def to_parquet(self, path:str):
        """
        Saves to a parquet file using gzip compression.
        Takes full path, /path/to/data.scmpra
        WILL clobber existing files with the same path.
        """
        #create a parquet table from the scMPRA data
        pa_data_table=pa.Table.from_pandas(self.data,preserve_index=True)
        #extract parquet metadata created in above, defaulting to empty dict
        pa_metadata=dict(pa_data_table.schema.metadata or {})
        #get all members of the object other than the actual data
        nondata={key:val for key, val in self.__dict__.items() if key != "data"}
        #add the class members to parquet metadata
        pa_metadata[b"scMPRA_data.members"]=json.dumps(nondata, default=str).encode("utf-8")

        #dump
        pa_data_table=pa_data_table.replace_schema_metadata(pa_metadata)
        pq.write_table(pa_data_table,path,compression="gzip")

    
    def graph_chimeric(self, *args, **kwargs):
        """
        TODO: test again now that its moved to scMPRA data obj
        
        Arguments
            self
            *args
            **kwargs

        Takes `scmpra_data`, a pandas dataframe of read-wise MPRA data (see docs) 
        and plots a histogram of frequency of reads per UMI using seaborn.histplot. 

        All other arguments are passed to the histplot call to allow graph 
        customization. Particular useful are `bins`, `binrange`, and `log_scale`
        """
        assert table_type(self.data.columns) == "mpra_readwise"
        
        sns.histplot(self.data['reads'], *args, **kwargs)

        plt.xlabel('Reads')
        plt.ylabel('Frequency')
        plt.title('Histogram of Reads')
        plt.show()
    
    def read_wise_to_umi_wise(self,keep_reads=False):
        """
        Converts read-wise to UMI-wise (see readme for spec).

        TODO: test again now that its moved to scMPRA data obj
        """

        assert self.table_type == "mpra_readwise", "Wrong table type."
        
        grouping_columns = [col for col in self.data.columns if col not in ['umi', 'reads']]


        aggregations = {
            'umis': ('umi', 'nunique')  # Count unique UMIs
        }

        # Conditionally include 'reads' sum
        if keep_reads:
            aggregations['reads'] = ('reads', 'sum')

        self.data = self.data.groupby(grouping_columns).agg(**aggregations).reset_index()
        self.table_type="mpra_umiwise"
        self.operations.append("read_wise_to_umi_wise")
    
    def cut_chimeric_reads(self,threshold):
        """
        Arguments
            self
            threshold : <int>
        
        subsets to those UMIs which lie ABOVE the number-of-reads threshold, 
        removing chimeric reads. 
        """
        assert table_type(self.data.columns) == "mpra_readwise"
        assert threshold >=0, "threshold must be greater than zero."
        
        #Trim
        ret=self.data[self.data["reads"]>threshold]

        original_umi_count=len(self.data["umi"].unique())
        cut_umi_count=len(ret["umi"].unique())

        logger.info(f"Original={original_umi_count} UMIs, Cut={cut_umi_count} UMIs, Lost={original_umi_count-cut_umi_count} UMIs.")

        self.data=ret

        self.operations.append(f"cut_chimeric_reads, threshold={threshold}")

    
    def ortho_filter(self):
        """
        Removes combinations of cre_id, cell_type which have less than MIN_PTS non-zero observations. 
        This is much stricter than filter_low_umi_count
        """
        tabtype = table_type(self.data.columns)
        assert tabtype == "mpra_umiwise", "Malformed table."

        # Count non-zero values per (cell_type, cre_id) group
        nonzero_counts = (
            self.data[self.data['umis_mpra_bc'] > 0]
            .groupby(['cell_type', 'cre_id'])
            .size()
            .reset_index(name='nonzero_count')
        )

        valid_combos = nonzero_counts.query('nonzero_count >= @MIN_PTS')[['cell_type', 'cre_id']]
        all_combos = nonzero_counts[['cell_type', 'cre_id']]

        # Compute dropped combos
        dropped_combos = pd.merge(all_combos, valid_combos, on=['cell_type', 'cre_id'], how='outer', indicator=True)
        dropped_combos = dropped_combos[dropped_combos['_merge'] == 'left_only'][['cell_type', 'cre_id']]

        # Keep only rows matching valid (cell_type, cre_id) combos
        self.data = self.data.merge(valid_combos, on=['cell_type', 'cre_id'], how='inner')

        # Print stats
        n_total = len(all_combos)
        n_dropped = len(dropped_combos)
        logger.info(f"Dropped {n_dropped} of {n_total} (cell_type, cre_id) combos with fewer than {MIN_PTS} nonzero entries.")


        # Record that we performed this operation
        self.operations.append((f"filter_low_umi_count, threshold={MIN_PTS}",dropped_combos))


    @unimplemented
    def round_down_zeroes():
        """
        TODO: implement
        ON HOLD: not clear how much we really care about making comparisons to zeroed CREs..
        If, later on, we find that we want to make these comparisons, come back and implement this.

        This is a follow-up to "ortho_filter".

        When a CRE+Cell-type combination has too few cells with non-zero observations
        to be modeled, this could be for one of two reasons.
        1. The CRE (in that cell-type) is simply expressed at too low a level to be detected here.
        2. Something technical went wrong. e.g. CRE was transfected into that cell-type at a low efficiency. 

        In either case, we cannot model that CRE-cell-type combination. However, recognizing the difference
        can be important. e.g. in case 2, we can't make any statements that that CRE + cell-type is 
        sig different from another CRE+cell-type. In case 1, we can! We just can't use the wald test.

        Things that differentiate case 1 from case 2 : 
        - If the few measurement(s) which are actually present in the data are a high number of UMIs
            that suggests 2 over 1. If they are a low number, that suggests 1 over 2. 
        - If the data are post transfection reporter filtering, then a large number of zeroes indicates
        case 1

        This function will use that information to call each filtered CRE as either "not detected" (1) or "dropout" (2)
        The default is "dropout", and the presence of sufficient evidence can move a filtered combination
        to "not detected".
        
        """
        pass

#        1         2         3         4         5         6         7         8
#2345678901234567890123456789012345678901234567890123456789012345678901234567890
    

@unimplemented
def flatten_barcode_errors(df, barcode_column,*args,**kwargs):
    """
    Need to re-work to work with scMPRA data object
    Arguments
        df <pandas.DataFrame>
        barcode_column <str>
    Returns
        <pandas.DataFrame>

    Uses umitools to flatten different barcodes which are likely only different
    due to sequencing errors. Passes *args,**kwargs upstream to bcs_to_lut. 
    """
    ret=df.copy()
    
    lut=bcs_to_lut(ret[barcode_column].value_counts().to_dict(),*args,**kwargs)

    ret[barcode_column]=ret[barcode_column].map(lut)

    return ret


@unimplemented
def apply_deseq():
    """
    R quarantine zone. 
    """
    pass


@unimplemented
def hypothesis_tester(scmpra_models, hypotheses, flavor="wald"):
    """
    Arguments
        scmpra_models : <pd.DataFrame>
        hypotheses : <pd.Dataframe>
        flavor : <str>
    Returns
        <pd.DataFrame>

    Takes .. and a set of hypotheses and tests them all. Returns
    a results dataframe. Flavor selects the test type, and can be one of
    - wald : wald test
    - wilcox : wilcoxin-rank-sum
    - permute : permutation test
    - deseq : uses deseq2
    
    """
    #calls apply_deseq
    assert table_type(scmpra_data.columns) == "mpra_umiwise"
    assert table_type(hypotheses.columns) == "hypotheses"
    pass


"""
ortho
- contains multiple experimental models for the same dataset
- useful sets are 
    - criss-cross
    - criss-cross subsetted to hypothesis test...

fit()
- takes a client & model params
- internally defines 
- returns an experiment_model

experiment_model
- mostly just an elaborate struct : minimal code, mostly results...
- no more SM, only tzinb
- will contain 1x model
- will NOT keep data
- lazy eval
- parameter extraction pulls triples from model
- self-description
	- model decisions : each is a / string. Strings compared against supported
        - nb form
        - zi form
        - split by
    - str dataset name (incl. source)
    - metadata
		- additional unstructured metadata in a dictionary

extract_triples()

triples

"""

#suggested formulas. 
#SUGGESTED_NB=['reads_mpra_bc ~ C(cell_type)*C(cre_id)',
#    'reads_mpra_bc ~ C(cell_type)',
#    'reads_mpra_bc ~ umis_transfection_bc:C(cell_type) + umis_transfection_bc:C(cre_id) + umis_transfection_bc:C(cell_type):C(cre_id) -1']
#SUGGESTED_ZI=['C(replicate)']
#SUGGESTED_BREAKBY=['']



def abort_on_failure(future,client):
    """
    Call this function with a completed future that is strictly necessary for the task at hand.
    If it didn't work, we'll crash. 
    """
    status=None
    if type(future)==type([3]):
        #if it's a list
        if any(a_future.status=="error" for a_future in future):
            status="error"
    else:
        #it's not a list, presumably just one future
        status=future.status
    if future.status=="error":
        sprint("[!] Computation failed. Aborting.")
        sprint("[!] Check slave node logs for details.")
        client.shutdown()
        assert 1==2

#        1         2         3         4         5         6         7         8
#2345678901234567890123456789012345678901234567890123456789012345678901234567890

def _smart_matrix(data,split):
    """
    Takes data and split column & produces design matricies
    according to standard ortho schema.
    """
    anti=anti_split(split)

    ## Decide model-type ##
    #0 levels will have already been filtered
    num_levels=len(data[anti].unique())
    if num_levels ==0:
        assert False, "Data were probably not filtered correctly."
    if num_levels>1:
        #more than one level between which we can perform hypothesis testing...
        model_type="contrastable"
    else:
        #only one level. Useless for hypothesis testing, but we will keep it around
        #in case we want to use it for simulation
        model_type="simulation_only"
    
    #now, pick the reference level
    if "reference" in data[anti].values:
        #if there is already a user-specificed reference level, let's just use that. 
        reference="reference"
    else:
        #Otherwise, use the most frequent level as the reference
        level_frequency=data[anti].value_counts()
        level_frequency=level_frequency.sort_values(ascending=False)
        reference=level_frequency.index[0]

    zi_formula="C(rep_id)-1"
    nb_formula=f"umis_mpra_bc ~ C({anti}, contr.treatment(base='{reference}'))"
    y, X=Formula(nb_formula).get_model_matrix(data,output='pandas')
    Z=Formula(zi_formula).get_model_matrix(data,output='pandas')
    return {
        'nb_regressors':X,
        'regressand':y,
        'zi_regressors':Z,
        'model_type':model_type,
        'nb_formula':nb_formula,
        'zi_formula':zi_formula,
        'reference':reference,
        'split':split
    }


def _tensorzinb_fit(matricies,name):
    """
    Takes matricies & produces a single tensorzinb model
    """
    zinbo = TensorZINB(endog=matricies["regressand"]["umis_mpra_bc"].to_numpy().squeeze(),
                    exog=matricies["nb_regressors"].to_numpy(),
                    exog_infl=matricies["zi_regressors"].to_numpy())
    

    result = zinbo.fit(init_method="nb")

    del zinbo
    
    return result

def standard_fit(client,data,split):
    """
    Takes an scMPRA object and produces a set of models along one axis,
    specified by split.
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

    tzinb_futures = {
        t: client.submit(
                _tensorzinb_fit,
                mats_futures[t],
                t
            )
        for t in levels
    }
    
    return(experiment_model(model=tzinb_futures,
                                split=split),
            mats_futures)

class experiment_model:
    """
    Contrains one full model of a dataset.
    .model is a dict
    - keys are whatever the levels of "split" are.
    - values are futures of tensorzinb return dictionaries.
    """

    def __init__(self,
            model,
            split:str):
        self.split=split
        self.model=model
    
    @staticmethod
    def _label_tensorzinb_regressors(model,dm):
        """
        Takes one tensorzinb model (dict) and correspinding set of design matricies
        & labels the dict with the regressor and regressand names.
        Meant to be submitted to a dask cluster.
        """

        nb_regressor_names=list(dm["nb_regressors"].columns)
        zi_regressor_names=list(dm["zi_regressors"].columns)
        model["weights"]["x_mu"] = pd.Series(model["weights"]["x_mu"].flatten(),
                                            index=nb_regressor_names)
        model["weights"]["x_pi"] = pd.Series(model["weights"]["x_pi"].flatten(),
                                            index=zi_regressor_names)
        return model
    
    def label_regressors(self,client,design_matricies):
        """
        Takes design matricies used to generate the model &
        modifies self in-place to have 
        """
        for key in self.model:
            self.model[key]=client.submit(experiment_model._label_tensorzinb_regressors,
                model=self.model[key],
                dm=design_matricies[key])
    

    def _unflatten_futures(self,client):
        """
        For internal use only
        wraps all the models in futures
        """
        for key in self.model:
            self.model[key]=client.submit(lambda x: x, self.model[key])
        
    
    def flattened_copy(self):
        """
        Makes a copy where the members are not futures but just objects
        """
        #initalize a copy
        ret=experiment_model(
            model={key:None for key in self.model},
            split=self.split
        )
        #gather model results
        for key in ret.model:
            ret.model[key]=self.model[key].result()
        
        return ret
    
    def save(self,path):
        """
        Saves experimentmodel to a filepath.
        Will hang if computation is not done yet
        """
        with open(path,"wb") as f:
            pickle.dump(self.flattened_copy(),f)

    @staticmethod
    def load(client,path):
        """
        Loads an experimentmodel from a path
        & returns it. Requires a client to wrap the individual models
        in futures for use on a dask cluster. 
        """
        with open(path,"rb") as f:
            ret=pickle.load(f)
            if not ret is None:
                ret._unflatten_futures(client)
        return ret


class ortho:
    """
    Stores multiple models of the same data
    one set of by_cre models, and one set of by cell type models
    Not to be used with multiple datasets. 
    """
    def __init__(self):
        self.by_cre=None
        self.by_cre_parameters=None
        
        self.by_cell_type=None
        self.by_cell_type_parameters=None

        self.by_cre_design=None
        self.by_cell_type_design=None
        
        self.training_data=None

    def save(self,path,name,strip_training_data=False):
        """
        Simple pickle save.

        Will block & wait for results if not done computing

        creates directory 'name' in 'path'
        """
        #There are much nicer ways to structure this, but that level of effort
        #should be saved for non-pickle save/load
        full_path=Path(path)/name
        full_path.mkdir(parents=True)

        ## Function
        
        def simple_write(obj,filename):
            with open(full_path/filename,"wb") as f:
                pickle.dump(obj,f)
        

        ## Experiment models
        
        if self.by_cre is None:
            simple_write(self.by_cre,"by_cre.pkl")
        else:
            self.by_cre.save(full_path/"by_cre.pkl")
        
        if self.by_cell_type is None:
            simple_write(self.by_cell_type,"by_cell_type.pkl")
        else:
            self.by_cell_type.save(full_path/"by_cell_type.pkl")

        
        ## Parameters
        if self.by_cre_parameters is None:
            simple_write(self.by_cre_parameters,"by_cre_parameters.pkl")
        else:
            self.by_cre_parameters.save(full_path/"by_cre_parameters.pkl")

        if self.by_cell_type_parameters is None:
            simple_write(self.by_cell_type_parameters,"by_cell_type_parameters.pkl")
        else:
            self.by_cell_type_parameters.save(full_path/"by_cell_type_parameters.pkl")
        
        ## Design matricies
        if self.by_cre_design is None:
            simple_write(self.by_cre_design,"by_cre_design.pkl")
        else:
            simple_write(dict_unwrap(self.by_cre_design),"by_cre_design.pkl")
        
        if self.by_cell_type_design is None:
            simple_write(self.by_cell_type_design,"by_cell_type_design.pkl")
        else:
            simple_write(dict_unwrap(self.by_cell_type_design),"by_cell_type_design.pkl")

        ## training data
        if not strip_training_data:
            simple_write(self.training_data,"training_data.pkl")
        else:
            simple_write(None,"training_data.pkl")

    @classmethod
    def load(cls,client,path,name):
        """
        loads from a filepath, wrapping in futures on the provided cluster where appropriate
        """
        #There are much nicer ways to structure this, but that level of effort
        #should be saved for non-pickle save/load
        
        full_path=Path(path)/name

        ret_ortho=cls()

        # function
        def simple_load(filename):
            with open(full_path/filename,"rb") as f:
                ret=pickle.load(f)
            return ret
        
        
        ## Experiment models
        ret_ortho.by_cre=experiment_model.load(client,full_path/"by_cre.pkl")
        ret_ortho.by_cell_type=experiment_model.load(client,full_path/"by_cell_type.pkl")

        ## Parameters
        ret_ortho.by_cre_parameters=parameters.load(client,full_path/"by_cre_parameters.pkl")
        ret_ortho.by_cell_type_parameters=parameters.load(client,full_path/"by_cell_type_parameters.pkl")

        ## Design matricies
        ret_ortho.by_cell_type_design=dict_wrap(client,simple_load("by_cell_type_design.pkl"))
        ret_ortho.by_cre_design=dict_wrap(client,simple_load("by_cre_design.pkl"))

        ## Training data
        ret_ortho.training_data=simple_load("training_data.pkl")
   
        
        return ret_ortho

    @unimplemented
    def clean(self,kill_list="auto"):
        """
        Deletes intermediate values to save space. 

        `kill_list` is any or all of "training_data", "design_matricies", "models", "parameters"
        
        alternatively, "auto" is equivalent to ["training_data", "design_matricies"]
        """
        for target in kill_list:
            pass
    
    def criss_cross(self,client,dat):
        """
        Makes by_cre and by_cell_type models.

        Note: a little computationally intensive...
        retain_metadata will keep some information 'dat' in self.training_data
        The actual MPRA data will be stripped to save space, but metadata will be retained
        """
        self.by_cre, self.by_cre_design=standard_fit(client,
                                                     dat,
                                                     split="cre_id")
        
        self.by_cell_type, self.by_cell_type_design=standard_fit(client,
                                                        dat,
                                                        split="cell_type")
        
        self.training_data=dat.copy()
        self.annotate_models(client)

    def annotate_models(self,client):
        """
        Adds regressor names to each model
        """
        self.by_cre.label_regressors(client,self.by_cre_design)
        self.by_cell_type.label_regressors(client,self.by_cell_type_design)
    
    def extract_params(self,client):
        """Extracts parameters for all models in the object"""

        if not self.by_cre is None:
            self.by_cre_parameters=extract_parameters(
                client,
                self.by_cre,
                self.by_cre_design,
                "cre_id")
        
        if not self.by_cell_type is None:
            self.by_cell_type_parameters=extract_parameters(
                client,
                self.by_cell_type,
                self.by_cell_type_design,
                "cell_type")
        
    def compute_model_qc(self):
        """
        Will hang if model is not finished. 
        returns None, sets self.by_cell_qc, self.by_cre_qc to dictionaries of 
        QC information comparing nb params of each direction of the ortho 
        to the data means.

        Meant for debugging / manual inspection.
        """
        self.by_cell_qc=ortho._nb_versus_means(params=self.by_cell_type_parameters,
            design_matricies=dict_unwrap(self.by_cell_type_design),
            scMPRAdat=self.training_data)
        self.by_cre_qc=ortho._nb_versus_means(params=self.by_cre_parameters,
            design_matricies=dict_unwrap(self.by_cre_design),
            scMPRAdat=self.training_data)
    
    @staticmethod
    def _nb_versus_means(params,design_matricies,scMPRAdat):
        """
        Takes a model & design matricies corresponding to one direction of an ortho
        (A set of 'by_cre' or a 'by_cell-type') and the original training data and produces a QC dictionary
        regressing data means against nb estimates.
        Used for quality control.
        """
        
        assert sorted(params.keys)==sorted(list(design_matricies.keys())), "mismatched model"
        split=params.broken_on
        anti=anti_split(split)

        data=scMPRAdat.data

        QC={}
        for model_level in params.keys:
            #model level are the levels of the split
            #e.g. a by-cell-type model will have cell-type values for model_level
            subset=data[data[split]==model_level]
            
            data_means=subset.groupby(anti)["umis_mpra_bc"].agg("mean").sort_values()
            data_means.name="mean(umis_mpra_bc)"

            mu_estimates=params.nb[model_level].result()

            mu_summary=mu_estimates.join(data_means,how="left")
            # Fit regression
            try:
                slope, intercept, r_value, p_value, std_err = linregress(mu_summary["mean(umis_mpra_bc)"], mu_summary["mu"])

                #store regression info
                ret={'success':True,
                    'dat':mu_summary,
                    'slope':slope,
                    'intercept':intercept,
                    'r_value':r_value,
                    'p_value':p_value,
                    'std_err':std_err
                }
            except ValueError:
                print(f"regression error on {model_level}")
                ret={'success':False,
                    'dat':mu_summary,
                    'slope':None,
                    'intercept':None,
                    'r_value':None,
                    'p_value':None,
                    'std_err':None
                }
            
            QC[model_level]=ret
        return QC

class parameters:
    """
    Stores triples of parameters
    - nb (negative binomial mean), zi (zero inflation), theta (dispersion parameter)

    for 'broken_on' (by cell type, by cre, or whatever models)

    """
    def __init__(self, nb, zi, theta, broken_on):
        self.nb=nb
        self.zi=zi
        self.theta=theta

        self.broken_on = broken_on

        assert nb.keys() == zi.keys()
        assert zi.keys() == theta.keys()

        self.keys=list(nb.keys())
        

    def _unflatten_futures(self,client):
        """
        Wraps all the models in futures
        """
        for key in self.keys:
            self.nb[key]=client.submit(lambda x: x, self.nb[key])
            self.zi[key]=client.submit(lambda x: x, self.zi[key])
            self.theta[key]=client.submit(lambda x: x, self.theta[key])
        
    
    def flattened_copy(self):
        """
        Makes a copy where the members are not futures but just objects.
        """
        #make a shallow copy
        ret=parameters(
            nb={key:None for key in self.nb},
            zi={key:None for key in self.zi},
            theta={key:None for key in self.theta},
            broken_on=self.broken_on
        )
        #fill shallow copy with result data...
        for key in self.keys:
            ret.nb[key]=self.nb[key].result()
            ret.zi[key]=self.zi[key].result()
            ret.theta[key]=self.theta[key].result()
        
        return ret
    
    def save(self,path):
        """
        Saves parameters to a filepath.
        Will hang if computation is not done yet
        """
        with open(path,"wb") as f:
            pickle.dump(self.flattened_copy(),f)

    @staticmethod
    def load(client,path):
        """
        Loads parameters from a path
        & returns it. Requires a client to wrap the individual models
        in futures for use on a dask cluster. 
        """
        with open(path,"rb") as f:
            ret=pickle.load(f)
            if not ret is None:
                ret._unflatten_futures(client)
        return ret


def extract_parameters(client,model,design,split):
    
    def _extract_mu(split,model,design_matrix):
        #unpack information from the design matrix.
        X=design_matrix["nb_regressors"]
        model_type=design_matrix["model_type"]
        reference=design_matrix["reference"]

        #multiply design matrix by weights
        linear_mu= X @ model["weights"]["x_mu"]
        #undo the link function to get predictions for each cell
        mu_predictions=np.exp(linear_mu)

        #now we have one mu value for each measurement, and we want to average 
        #to get one value for each level...
        
        #Get the level names for each row in the design matrix
        #e.g. the cell type for each measurement for a by-cre model,
        #the cre name for each measurement for a by-cell-type model...
        level_type=anti_split(split)

        if model_type=="contrastable":
            #multiple levels. 
            row_labeling=undo_one_hot_encoding(X)
            row_labeling=row_labeling[f"{level_type}, contr.treatment(base='{reference}')"]
            row_labeling=row_labeling.str.removeprefix("T.")
        elif model_type=="simulation_only":
            #only one level
            row_labeling=pd.Series(np.repeat(reference,len(mu_predictions)))
        #Merge those level names onto the predictions for each row
        mu_predictions_df=pd.DataFrame({level_type:row_labeling,'mu':mu_predictions.squeeze()})

        #aggregate predictions to one per level name
        mu_summary=mu_predictions_df.groupby(level_type).agg("mean")

        return mu_summary

    def _extract_zi(model,design_matrix):
        Z=design_matrix["zi_regressors"]
        ## ZI ##
        # multiply design matrix by weights
        linear_zi=(Z.to_numpy() @ model["weights"]["x_pi"])
        zi_predictions=linear_zi=1/(1+np.exp(-linear_zi))
        
        #extract names 
        replicate_labeling=undo_one_hot_encoding(Z)["rep_id"]
        replicate_labeling=replicate_labeling.str.removeprefix("T.")

        #apply names to ZI
        zi_predictions_df=pd.DataFrame({'rep_id':replicate_labeling,'zi':zi_predictions.squeeze()})
        #aggregate
        zi_summary=zi_predictions_df.groupby("rep_id").agg("mean")

        return zi_summary
    
    def _extract_theta(model):
        theta=np.exp(model['weights']['theta'].squeeze())
        return theta

    mus={}
    zis={}
    thetas={}

    #not worth parallelizing this loop: submission time would be greater than saved time.
    for level in model.model:

        mus[level]=client.submit(_extract_mu,
                         split=split,
                         model=model.model[level],
                         design_matrix=design[level])
        
        zis[level]=client.submit(_extract_zi,
                         model=model.model[level],
                         design_matrix=design[level])
        
        thetas[level]=client.submit(_extract_theta,
                            model=model.model[level])
    
    return parameters(nb=mus,zi=zis,theta=thetas,broken_on=split)

def cast_multiindex_to_str_inplace(df):
    """Convert all levels of a MultiIndex to strings, modifying df in-place."""
    if not isinstance(df.index, pd.MultiIndex):
        raise ValueError("Index is not a MultiIndex.")

    new_index = pd.MultiIndex.from_tuples(
        [tuple(str(x) for x in tup) for tup in df.index],
        names=df.index.names
    )
    df.index = new_index  # triggers inplace update of index

def describe_parameters(parameters,dat,split):
    """
    Produces a convienient single-dataframe description of one split model. 
    'parameters' is a parameters object
    Requires that you pass original data as dat to compute cell-numbers
    Leaves non-split columns as one-hot. Returns "split" column as str categorical
    """

    #change to return a dask instead of pandas dataframe

    anti=anti_split(split)

    #count cells per group
    cell_counts=dat.groupby([split,anti,"rep_id"]).size()
    cell_counts=pd.DataFrame({"cells":cell_counts})
    
    #cast rep id to string just in case.
    cell_counts.index = cell_counts.index.set_levels(
        cell_counts.index.levels[1].astype(str), level=1
    )

    flattened_param=flatten_param_representation(parameters.flattened_copy(),split=split)
    #flattened_param["rep_id"]=flattened_param["rep_id"].astype(str)
    flattened_param=flattened_param.reset_index().set_index([split,anti,"rep_id"])

    cast_multiindex_to_str_inplace(flattened_param)
    cast_multiindex_to_str_inplace(cell_counts)
    
    lost_params=flattened_param[flattened_param.isna().any(axis=1)]
    if not lost_params.empty:
        logger.warning(f"Losing params {lost_params}")
        flattened_param=flattened_param.dropna()

    cell_counts=cell_counts.dropna()

    working=flattened_param.join(cell_counts,how="inner")
    working["r"]=working["theta"]
    working["sigmasquare"]=working["mu"]**2/working["r"]+working["mu"]
    working["p"]=working["mu"]/working["sigmasquare"]
    working=working.reset_index()


    return working

def flatten_param_representation(params, split: str):
    dfs=[]
    for key in params.keys:
        nb=params.nb[key].reset_index()
        zi=params.zi[key].reset_index()
        cartesian = nb.merge(zi, how='cross')
        cartesian[split] = np.repeat(key, len(cartesian))
        cartesian['theta'] = np.repeat(params.theta[key], len(cartesian))

        anti=anti_split(split)
        index_cols=[split,anti,"rep_id"]
        #index_cols = (
        #    params.nb[key].index.names +
        #    params.zi[key].index.names +
        #    [split]
        #)
        cartesian.set_index(index_cols, inplace=True)

        dfs.append(cartesian)
    return pd.concat(dfs)

def _flatten_param_representation(client: Client, params, split: str):
    
    #params_future = client.scatter(params, broadcast=True)
    # probably want to change this later once dat is always a future
    keys = params.keys  # assume this is list-like

    def process_key(key, params):
        nb_reset = params.nb[key].reset_index()
        zi_reset = params.zi[key].reset_index()
        cartesian = nb_reset.merge(zi_reset, how='cross')

        cartesian[split] = np.repeat(key, len(cartesian))
        cartesian['theta'] = np.repeat(params.theta[key], len(cartesian))

        index_cols = (
            params.nb[key].index.names +
            params.zi[key].index.names +
            [split]
        )
        cartesian.set_index(index_cols, inplace=True)
        return cartesian

    futures = [
        client.submit(process_key, key, params_future)
        for key in keys
    ]
    results = client.gather(futures)
    return pd.concat(results)


def anti_split(split):
    """ 
    Returns the opposite split for a given split.
    If splits are extended beyond 'cre_id' and 'cell_type',
    this function should be extended to handle those cases.
    Specifcally, all but split should be returned.
    """
    if split == "cre_id":
        return "cell_type"
    elif split == "cell_type":
        return "cre_id"
    else:
        raise ValueError(f"Unsupported split: {split}")


def get_cell_counts(client: Client, dat: pd.DataFrame, split: str):
    """
    Takes a dask client and a pandas DataFrame `dat` containing MPRA data.
    """
    # Broadcast dat to all workers once
    # probably want to change this later once dat is always a future
    dat_future = client.scatter(dat, broadcast=True)

    def process_key(key, dat):
        relevant_subset = dat[dat[split] == key]

        relevant_subset=relevant_subset.drop(columns=[split])

        anti=anti_split(split)
        formula = Formula(f"umis_mpra_bc ~ C({anti}) + C(rep_id) - 1")
        _, mat = formula.get_model_matrix(relevant_subset, output='pandas',ensure_full_rank=False)
        
        mat = mat.value_counts()

        mat = pd.DataFrame(mat)
        mat = mat.rename({'count': 'cells'}, axis=1)
        mat[split] = np.repeat(key, len(mat))

        #fix the index to include the split-on column...
        old_index=mat.index.names
        mat.reset_index(inplace=True)
        mat.set_index(old_index+[split],inplace=True)

        
        return mat

    keys = dat[split].unique()
    futures = [client.submit(process_key, key, dat_future) for key in keys]
    results = client.gather(futures)
    return pd.concat(results)

def auto_partition(pdf, target_mb_per_partition=PARTITION_SIZE_MB):
    """
    Convert a pandas DataFrame to a Dask DataFrame with automatic partition sizing.
    
    Parameters:
    - pdf: input pandas DataFrame
    - target_mb_per_partition: desired memory usage per partition (in megabytes)
    
    Returns:
    - ddf: Dask DataFrame with chosen number of partitions

    Minimum of 2!
    """
    
    pdf=pdf.reset_index()#dask doesn't like multi-indexes, so reset to single index.
    est_bytes = pdf.memory_usage(index=True, deep=True).sum()
    target_bytes = target_mb_per_partition * 1_000_000
    npartitions = max(2, int(np.ceil(est_bytes / target_bytes)))
    return dd.from_pandas(pdf, npartitions=npartitions)

def _simulate_from_description(description):
    """
    Simulate from a description dask dataframe
    """
    # Repeat rows by 'cells' count, exploding to one row per cell
    df=description
    repeated_df = df.loc[df.index.repeat(df['cells'])].reset_index(drop=True)
    #repeated_df = description.loc[description.index.repeat(description['cells'])].reset_index(drop=True)

    # Simulate NB and ZI in numpy
    r = repeated_df['r'].to_numpy()
    p = repeated_df['p'].to_numpy()
    zi = repeated_df['zi'].to_numpy()

    nb = np.random.negative_binomial(n=r, p=p)
    keep_mask = np.random.binomial(n=1, p=1 - zi)
    zinb = nb * keep_mask

    repeated_df['zinb_sample'] = zinb
    return repeated_df

def simulate_from_description(description):
    """
    Simulate from a description dask dataframe
    """

    def simulate_partition(df):
        # Repeat rows by 'cells' count, exploding to one row per cell
        repeated_df = df.loc[df.index.repeat(df['cells'])].reset_index(drop=True)

        # Simulate NB and ZI in numpy
        r = repeated_df['r'].to_numpy()
        p = repeated_df['p'].to_numpy()
        zi = repeated_df['zi'].to_numpy()

        nb = np.random.negative_binomial(n=r, p=p)
        keep_mask = np.random.binomial(n=1, p=1 - zi)
        zinb = nb * keep_mask

        repeated_df['zinb_sample'] = zinb
        return repeated_df
    # Use auto_partition to repartition the description DataFrame
    #description = auto_partition(description)
    #npartitions=
    description=description.repartition(npartitions=2)
    return description.map_partitions(simulate_partition)

class simulation_batch:
    """
    Class which takes a single <ortho> object and simulates replicates.
    Optionally, fits additional ortho objects to simulations & plots their paremeter spread
    Useful for estimating variance of an experimental setup...
    """
    
    #consolidate split and parameter validity checking

    splits=["cre_id","cell_type"]

    def __init__(self,primordial,partition_mb=50):
        #the initial 
        self.primordial=primordial
        description_primordial_by_cre=None
        description_primordial_by_cell_type=None
        
        self.simulated_from_cre=[]
        self.simulated_from_cell_type=[]

        self.ortho_simulated_cre=[]
        self.ortho_simulated_cell_type=[]

        self.partition_mb=partition_mb
    
    def describe_primordial(self):
        """Generates and saves descriptions of the primordial which are necessary for subsequent simulation"""
        self.description_primordial_by_cre=describe_parameters(parameters=self.primordial.by_cre_parameters,
                                                                   dat=self.primordial.training_data.data,
                                                                   split="cre_id")

        self.description_primordial_by_cre=auto_partition(self.description_primordial_by_cre,
                                                              self.partition_mb)

        self.description_primordial_by_cell_type=describe_parameters(parameters=self.primordial.by_cell_type_parameters,
                                                                     dat=self.primordial.training_data.data,
                                                                     split="cell_type")
        
        self.description_primordial_by_cell_type=auto_partition(self.description_primordial_by_cell_type,
                                                                    self.partition_mb)

    @unimplemented
    def clear_simulations(self):
        """Removes simulated data. Does not remove models fit to simulated data. Useful for reducing object size"""
        pass
        #del self.simulations
        #self.simulations=[]
    
    def _round_of_simu(self,client,description_primordial):
        """
        Simulates one replicate from a given description. 
        """
        #s for simulated
        s=simulate_from_description(description_primordial).compute()
        s=undo_one_hot_encoding(s)
        s=s.rename({'zinb_sample':'umis_mpra_bc'},axis=1)[['rep_id','cre_id','cell_type','umis_mpra_bc']]
        
        ret=scMPRA_data()
        ret.flag_synthetic()
        ret.data=s
        return ret
    
    def simulate_many(self,client,n):
        """
        simulates n replicates

        todo: additional parallelism
        todo: pick which set of models to create : by cre, by cell-type, or both
        currently all both
        """
        #add parallel scatter-gather here!
        for _ in range(0,n):
            self.simulated_from_cre.append(self._round_of_simu(client,self.description_primordial_by_cre))
            self.simulated_from_cell_type.append(self._round_of_simu(client,self.description_primordial_by_cell_type))

    def fit_to_simulations(self,client):
        """
        Fits ortho models to all simulated datasets.
        Depends on simulate_many having been called. 
        """
        
        #add parallel scatter-gather here!

        #from cre
        for s in self.simulated_from_cre:
            recap=ortho()
            recap.criss_cross(client=client,dat=s)
            recap.extract_params(client)
            self.ortho_simulated_cre.append(recap)

        #from cell
        for s in self.simulated_from_cell_type:
            recap=ortho()
            recap.criss_cross(client=client,dat=s)
            recap.extract_params(client)
            self.ortho_simulated_cell_type.append(recap)
        
        self._flatten_all_parameters()
    
    def _flatten_all_parameters(self):
        """Flatten all simulated ortho parameters into single dataframes for convienient access"""
        #old from_cre and from_cell_type are the orthos...

        splits=["cre_id","cell_type"]
        vars=["nb","zi","theta"]

        nbs=[]
        zi_cre=[]
        zi_cell=[]

        theta_cre=[]
        theta_cell=[]

        def flatten_results(working,split:str,var:str):
            """Internal, temporary function"""
            if split not in splits or var not in vars:
                raise ValueError("invalid chocie")
            #add check for valid working split var
            
            if var=="theta":
                if split=="cre_id":
                    working=working.by_cre_parameters.flattened_copy()
                elif split=="cell_type":
                    working=working.by_cell_type_parameters.flattened_copy()
                
                working=pd.DataFrame([working.theta]).T
                working=working.reset_index()
                working.columns=[split,"theta"]
                working['theta']=working['theta'].astype(float)
                return working

            if split=="cre_id":
                working=working.by_cre_parameters.flattened_copy()
            elif split=="cell_type":
                working=working.by_cell_type_parameters.flattened_copy()
            if var =="nb":
                working=pd.concat(working.nb)
            elif var =="zi":
                working=pd.concat(working.zi)
            working.index=working.index.set_names(split,level=0)
            working=working.reset_index()
            working=undo_one_hot_encoding(working)
            return working

        #add primordial for reference
        for split in splits:
            for var in ["nb","zi","theta"]:
                #nb
                prim_flattened=flatten_results(self.primordial,split,var)
                prim_flattened["id"]=f"primordial {split}"
                prim_flattened["rep"]="primordial"
                if var=="nb":
                    nbs.append(prim_flattened)
                elif var=="zi":
                    if split=="cre_id":
                        zi_cre.append(prim_flattened)
                    elif split=="cell_type":
                        zi_cell.append(prim_flattened)
                elif var=="theta":
                    if split=="cre_id":
                        theta_cre.append(prim_flattened)
                    elif split=="cell_type":
                        theta_cell.append(prim_flattened)
        
        #iterate over all simulation models & flatten out their parameters
        for _from in splits:
            orthos=None
            if _from=="cre_id":
                orthos=self.ortho_simulated_cre
            else:
                orthos=self.ortho_simulated_cell_type
            
            for _to in splits:
                
                id_string=f"{_from}->{_to}"
                for i, ortho in enumerate(orthos):

                    #ugly duplicated code, fix...

                    #nb
                    working=flatten_results(ortho,split=_to,var="nb")
                    working['id']=id_string
                    working['rep']=f"rep {i}"
                    nbs.append(working)

                    #zi
                    working=flatten_results(ortho,split=_to,var="zi")
                    working['id']=id_string
                    working['rep']=f"rep {i}"
                    if _to=="cre_id":
                        zi_cre.append(working)
                    else:
                        zi_cell.append(working)

                    #theta
                    working=flatten_results(ortho,split=_to,var="theta")
                    working['id']=id_string
                    working['rep']=f"rep {i}"
                    if _to=="cre_id":
                        theta_cre.append(working)
                    else:
                        theta_cell.append(working)

        # concat and add to object
        self._nbs=pd.concat(nbs)
        self._zi_cre=pd.concat(zi_cre)
        self._zi_cell=pd.concat(zi_cell)
        self._theta_cre=pd.concat(theta_cre)
        self._theta_cell=pd.concat(theta_cell)
    
    def plot_theta_spread(self, split):
        """Plots the spread of the thetas of simulated data with primordial for reference."""
        #no need to copy : no mutation
        #theta=theta.copy()
        
        theta=None
        if split=="cre_id":
            theta=self._theta_cre
        elif split=="cell_type":
            theta=self._theta_cell

        plt.figure(figsize=(12, 6))

        #basic violin plot
        sns.violinplot(
            data=theta, x=split, y='theta',
            inner=None, palette='Set1'
        )

        # Regular points
        sns.scatterplot(
            data=theta[~theta['id'].isin(['primordial cre_id', 'primordial cell_type'])],
            x=split, y='theta',
            hue=split, style='id',
            palette='Set1', edgecolor='black', s=50, alpha=0.7, legend='brief'
        )

        # Primordial points overlay (bold black Xs)
        sns.scatterplot(
            data=theta[theta['id'].isin(['primordial cre_id', 'primordial cell_type'])],
            x=split, y='theta',
            edgecolor='black', marker='X',
            s=120, zorder=10,hue='id'
        )

        plt.xlabel(f'{split}')
        plt.ylabel('theta')
        plt.title(f'theta parameters by {split}')
        plt.xticks(rotation=45, ha='right')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.show()

    def plot_zi_spread(self,split):
        """Plots the spread of the zi parameters of simulated data with primordial for reference."""
        zi=None
        if split=="cre_id":
            zi=self._zi_cre
        elif split=="cell_type":
            zi=self._zi_cell

        
        zi=zi.copy()
        zi['group']=zi[split] + " | " + zi["rep_id"]

        plt.figure(figsize=(12, 6))
        
        #basic violin plot
        sns.violinplot(
            data=zi, x='group', y='zi',
            inner=None, palette='Set1'
        )

        # Regular points
        sns.scatterplot(
            data=zi[~zi['id'].isin(['primordial cre_id', 'primordial cell_type'])],
            x='group', y='zi',
            hue=split, style='id',
            palette='Set1', edgecolor='black', s=50, alpha=0.7#, legend='brief'
        )

        # Primordial points overlay (bold black Xs)
        sns.scatterplot(
            data=zi[zi['id'].isin(['primordial cre_id', 'primordial cell_type'])],
            x='group', y='zi',
            edgecolor='black', marker='X',
            s=120, zorder=10, hue='id'
        )
        

        plt.xlabel(f'{split} | rep_id')
        plt.ylabel('zi')
        plt.title(f'zi parameters by {split} and rep_id')
        plt.xticks(rotation=45, ha='right')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.show()

    def plot_nb_spread(self, groups_per_plot=-1):
        nbs = self._nbs
        nbs = nbs.reset_index(drop=True)
        nbs['group'] = nbs['cell_type'] + " | " + nbs['cre_id']

        all_groups = nbs['group'].unique()
        num_groups = len(all_groups)

        if groups_per_plot == -1:
            groups_per_plot = num_groups

        for i in range(0, num_groups, groups_per_plot):
            group_chunk = all_groups[i:i+groups_per_plot]
            subset = nbs[nbs['group'].isin(group_chunk)]

            plt.figure(figsize=(12, 6))
            sns.violinplot(
                data=subset, x='group', y='mu',
                inner=None, palette='Set1'
            )
            sns.scatterplot(
                data=subset[~subset['id'].isin(['primordial cre_id', 'primordial cell_type'])],
                x='group', y='mu',
                hue='cell_type', style='id',
                palette='Set1', s=50, alpha=0.7
            )
            sns.scatterplot(
                data=subset[subset['id'].isin(['primordial cre_id', 'primordial cell_type'])],
                x='group', y='mu',
                edgecolor='black', marker='X',
                s=120, zorder=10, hue='id'
            )
            plt.xlabel('cell_type | cre_id')
            plt.ylabel('mu')
            plt.title(f'nb parameters (mu) by cell_type and cre_id — groups {i+1} to {min(i+groups_per_plot, num_groups)}')
            plt.xticks(rotation=45, ha='right')
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.show()

    def save(self,path,name):
        """
        Simple save of a full simulation batch.
        creates directory 'name' in 'path' (and many subdirectories).

        Caveats:
        1. Does not keep a copy of training data in simulated orthos, since 
        2. It does not save the parameters objects, since they would be redundant with the summary dataframes.
        Instead, you can use 
        - description_primordial_by_cre
        - description_primordial_by_cell_type
        to get the values of the primordial parameters, and
        - _nbs
        - _zi_cre
        - _zi_cell
        - _theta_cre
        - _theta_cell
        for the parameter values of the simulated replicates.

        (or you can have them recomputed).

        Will block & wait for results if not done computing
        """
        # Create output directory #
        root_path=Path(path)/name
        try:
            root_path.mkdir(parents=True)
        except FileExistsError:
            logger.error(f"Cowardly refusing to overwrite {root_path}.")
            raise(FileExistsError)
        # save primordial ortho if present #
        if self.primordial is not None:
            self.primordial.save(path=root_path,name="primordial_ortho")

        # save the actual simulated replicate dfs if present. #

        #loop over each set of simulated datasets
        for from_by in ["simulated_from_cre","simulated_from_cell_type"]:
        
            simulations=getattr(self,from_by)
            
            #make an output directory

            simulations_path=root_path/from_by
            simulations_path.mkdir()
            
            #if we actually have simulated datsets of this class:
            if not simulations is None:

                #dump each as a parquet
                for idx, simulated_dataset in enumerate(simulations):
                    simulated_dataset.to_parquet(simulations_path/f"{idx}.scmpra")

        # save the simulated replicate orthos if present. #

        #generally, the relationship will be : models mean there are simulated data
        #simulated data doesn't mean there are models.

        for from_by in ["ortho_simulated_cre","ortho_simulated_cell_type"]:
            ortho_simulations=getattr(self,from_by)
            
            if not ortho_simulations is None:
                orthos_path=root_path/from_by
                for idx, simulated_ortho in enumerate(ortho_simulations):
                    #We strip the data since we just saved it above. 
                    simulated_ortho.save(path=orthos_path,name=f"{idx}",strip_training_data=True)


        #save the primordial descriptions, if they exist.
        for from_by in ["description_primordial_by_cre","description_primordial_by_cell_type"]:
            if not getattr(self,from_by) is None:
                description_path=root_path/from_by
                getattr(self,from_by).to_parquet(description_path,compression="gzip")
        
        # save the parameters of the simulated reps if present
        sim_parameters_root=root_path/"simulated_rep_parameters"
        sim_parameters_root.mkdir()

        for from_by in ["_nbs","_zi_cre","_zi_cell","_theta_cre","_theta_cell"]:
            if not getattr(self,from_by) is None:
                param_path=sim_parameters_root/from_by
                pa_data_table=pa.Table.from_pandas(getattr(self,from_by),preserve_index=True)
                pq.write_table(pa_data_table,param_path,compression="gzip")

    @classmethod
    def load(cls,client,path,name):
        """
        Loads a batch saved with save_batch.
        See the caveats noted in the docstring for save_batch.
        Requires you pass a client to re-wrap futures. 
        """
        # check directory #
        root_path=Path(path)/name
        
        if not root_path.is_dir():
            logger.error(f"No directory {root_path}.")
            raise(FileNotFoundError)
        
        # create the batch object #
        ret=simulation_batch(None)
        
        # load primordial orthos if present #
        primordial_path=root_path/"primordial_ortho"
        if not root_path.is_dir():
            logger.info("No primordial found, skipping...")
        else:
            ret.primordial=ortho.load(client,path=root_path,name="primordial_ortho")

        
        # load the actual simulated replicate dfs if present. #

        for from_by in ["simulated_from_cre","simulated_from_cell_type"]:
            setattr(ret,from_by,[])
            replicates_root=root_path/from_by
            for simulated_replicate_path in replicates_root.iterdir():
                replicate=scMPRA_data.from_parquet(simulated_replicate_path)
                getattr(ret,from_by).append(replicate)

        # Load the actual simulated replicate orthos #

        for from_by in ["ortho_simulated_cre","ortho_simulated_cell_type"]:
            setattr(ret,from_by,[])
            ortho_replicates_root=root_path/from_by
            for simulated_ortho_path in ortho_replicates_root.iterdir():
                replicate_ortho=ortho.load(client,ortho_replicates_root,simulated_ortho_path)
                getattr(ret,from_by).append(replicate_ortho)

        # load the primordial descriptions, if they exist.

        for from_by in ["description_primordial_by_cre","description_primordial_by_cell_type"]:
            primordial_path=root_path/from_by
            if primordial_path.exists():
                setattr(ret,from_by,dd.read_parquet(primordial_path))
        
        # load the parameters of the simulated reps if present
        sim_parameters_root=root_path/"simulated_rep_parameters"
        if sim_parameters_root.exists():
            for from_by in ["_nbs","_zi_cre","_zi_cell","_theta_cre","_theta_cell"]:
                pa_data_table=pq.read_table(sim_parameters_root/from_by)
                pd_data_table=pa_data_table.to_pandas(types_mapper=pd.ArrowDtype)
                setattr(ret,from_by,pd_data_table)
        
        return ret


def versus_truth(ground_truth_mu:pd.DataFrame,inp_ortho:ortho):
    """
    Function takes a dataframe of ground truth values for each CRE, cell-type combination and compares to estimated parameters.

    Note that mean absolute percentage error is only reported for cases where the truth values is nonzero.

    TODO: clean up duplicate code
    """
    ret_mse=[]#mean squared error
    ret_mbe=[]#mean biased error
    ret_mape=[]#mean absolute percentage error
    
    def mse(df):
        return np.mean((df["mu"]-df["true_mu"])**2)
    def mbe(df):
        return np.mean(df["mu"]-df["true_mu"])
    def mape(df):
        df_copy=df[df["true_mu"]!=0].copy()
        return 100*np.mean(np.abs(df_copy["mu"]-df_copy["true_mu"])/df_copy["true_mu"])
    
    ## By cell_type ##
    if inp_ortho.by_cell_type_parameters is None:
        ret_mse.append(None)
        ret_mbe.append(None)
        ret_mape.append(None)
        logger.warning("Missing parameters. Maybe run ortho.extract_params(client) first.")
    else:
        by_cell_type=describe_parameters(parameters=inp_ortho.by_cell_type_parameters,
                            dat=inp_ortho.training_data.data,
                            split="cell_type")
        
        by_cell_type=by_cell_type.set_index(['cell_type','cre_id'])

        comp_by_cell_type=by_cell_type.join(ground_truth_mu,how="inner")
        comp_by_cell_type=comp_by_cell_type[["mu","true_mu"]].drop_duplicates()
        
        ret_mse.append(mse(comp_by_cell_type))
        ret_mbe.append(mbe(comp_by_cell_type))
        ret_mape.append(mape(comp_by_cell_type))

    ## by_cre ##
    if inp_ortho.by_cre_parameters is None:
        ret_mse.append(None)
        ret_mbe.append(None)
        logger.warning("Missing parameters. Maybe run ortho.extract_params(client) first.")
    else:
        by_cre=describe_parameters(parameters=inp_ortho.by_cre_parameters,
                            dat=inp_ortho.training_data.data,
                            split="cre_id")
        
        by_cre=by_cre.set_index(['cell_type','cre_id'])

        comp_by_cre=by_cre.join(ground_truth_mu,how="inner")
        comp_by_cre=comp_by_cre[["mu","true_mu"]].drop_duplicates()
        ret_mse.append(mse(comp_by_cre))
        ret_mbe.append(mbe(comp_by_cre))
        ret_mape.append(mape(comp_by_cre))
    
    return pd.DataFrame({"by":["cell_type","cre_id"],
                         "mean_squared_error":ret_mse,
                         "mean_biased_error":ret_mbe,
                         "mean_absolute_percent_error":ret_mape})

@unimplemented
def volcano(results:experiment_model):
    """
    Volcano plot of p value versus log fold change
    """
    pass
