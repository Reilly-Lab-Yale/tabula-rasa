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

import itertools

import scipy
from scipy.stats import linregress, chi2, norm, mannwhitneyu

import statsmodels.api as sm
import statsmodels.discrete.discrete_model as smd

import patsy
from tensorzinb.tensorzinb import TensorZINB
from statsmodels.stats.multitest import fdrcorrection
from formulaic import Formula

from dask.distributed import Client, Future, get_client

import dask.dataframe as dd
import dask.array as da
import dask

import os

from enum import Enum
from typing import List, Dict, Sequence, Tuple, Optional

from dataclasses import dataclass, replace
import json
import tarfile
import tempfile

from sklearn.metrics import precision_recall_curve, average_precision_score, roc_curve, roc_auc_score, confusion_matrix

import warnings

# tensorflow import for Wald test Hessian/SE computation
try:
    import tensorflow as tf
    tf.compat.v1.disable_eager_execution()  # graph mode only
except Exception as _tf_err:
    tf = None

import scipy.linalg as la
from numpy.linalg import LinAlgError
import uuid

#internal imports
from .utils import unimplemented
from .utils import bcs_to_lut
from .utils import undo_one_hot_encoding
from .utils import dict_wrap, dict_unwrap
from .utils import one_versus_all, find_treatment_column
from .utils import generate_barcodes, sample_from_library
from .utils import alpha_for_expected_groups, sample_crp_groups, _plot_test_bars
from .utils import one_library_replicate, pow_curve
logger = logging.getLogger("scMPRAforge")

def dump_df_debug(df, prefix="debug_df", outdir="."):
    uid = uuid.uuid4().hex
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    path = outdir / f"{prefix}_{uid}.tsv.gz"
    df.to_csv(path, sep="\t", index=False, compression="gzip")

    logger.info(f"[debug] dumped df to: {path}")
    return path

def dump_df_pickle_debug(df, prefix="debug_df", outdir="."):
    uid = uuid.uuid4().hex
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    path = outdir / f"{prefix}_{uid}.pkl"
    df.to_pickle(path)

    print(f"[debug] pickled df to: {path}")
    return path


def load_df_pickle_debug(path):
    path = Path(path)
    df = pd.read_pickle(path)

    print(f"[debug] loaded df from: {path}")
    return df

MIN_PTS=3
PARTITION_SIZE_MB=50
DEFAULT_SIGNIFICANCE_THRESHOLD=0.05
#If abs(a-b) < FLOATING_POINT_DIFF, we conclude a=b.
FLOATING_POINT_DIFF=1e-8

# Table schemas (centralized to avoid drift), open to packing this up into the class object if we feel that is cleaner
HYPOTHESIS_REQUIRED = {"comparison_CRE", "comparison_cell_type"}
HYPOTHESIS_OPTIONAL = {"reference_CRE", "reference_cell_type", "meta"}
HYPOTHESIS_ALL = HYPOTHESIS_REQUIRED | HYPOTHESIS_OPTIONAL

WARN_MULTI_TRANSFECTION_PERCENT=2.0

ERROR_TEST_NAN_PERCENT=5

RESULT_REQUIRED = {
    "test_type", "test_statistic", "p_value", "fold_change", "bh_p", "flattened"
}
RESULT_ALL = HYPOTHESIS_ALL | RESULT_REQUIRED

THREADS_DEFAULT = 5

#functions
@unimplemented
def always_unfinished():
    """tests unimplemented decorator."""
    pass

def clobber_mkdir(dir):
    try:
        dir.mkdir()
    except FileExistsError:
        logger.warn(f"Directory {dir} already exists, continuing.")

def helloworld():
    print("hello world!")

def _require_tensorflow():
    if tf is None:
        raise ImportError(
            "TensorFlow is required for Wald test Hessian/SE computation. "
            "Please install tensorflow>=2.x."
        )

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
    ### UPDATE: ENG refactored slightly to incorporate hypothesis table and simplify some lines
    #Performs subset checking so that extra columns are allowed.
    #Matching multiple definitions however is NOT allowed. 
    
    #If columns match no definitions the table is malformed: so this is the default. 
    cols = set(map(str, column_names))  # tolerate ArrowDtype/Index types
    ret = 'malformed'
    matches = 0

    # Read-wise MPRA (use the names used elsewhere in the codebase)
    if {'cell_bc','rep_id','cre_id','cell_type','mpra_bc','umi','reads'} <= cols:
        ret = 'mpra_readwise'; matches += 1

    # UMI-wise MPRA
    if {'rep_id','cre_id','cell_type','umis_mpra_bc'} <= cols:
        ret = 'mpra_umiwise'; matches += 1

    # Hypotheses (must have required; optional ok)
    if HYPOTHESIS_REQUIRED <= cols and cols <= (HYPOTHESIS_ALL | cols):
        # Distinguish result vs hypothesis next:
        if RESULT_REQUIRED <= cols:
            ret = 'results'
        else:
            ret = 'hypotheses'
        matches += 1

    return ret if matches == 1 else 'malformed'

import re

def _extract_square(s):
    """
    Utility function which extracts the contents of [T. ] or [ ] brackets from  string `s`
    Assumes one pair of brackets
    If no brackets are found, returns `s` untouched.
    """
    tm = re.search(r"\[T\.(.*?)\]", s)
    m = re.search(r"\[(.*?)\]", s)
    if tm:
        return tm.group(1)
    elif m:
        return m.group(1)
    else:
        return s
    

def _matricies_to_order(matricies):
    """
    Helper function. Extracts column index order from matricies for use elsewhere.
    Replaces `Intercept` with `reference`.
    """
    zi_idx=matricies["zi_regressors"].columns.to_list()
    zi_idx=[_extract_square(s) if s !="Intercept" else "reference" for s in zi_idx]

    nb_idx=matricies["nb_regressors"].columns.to_list()
    nb_idx=[_extract_square(s) if s !="Intercept" else "reference" for s in nb_idx]

    return {'zi_idx':zi_idx,'nb_idx':nb_idx}

def _mom_from_training_data(data,split,subset,indicies):
    """
    Helper function implementing warm start method of moments for parameter initalization.
    See [[Fixing zinb initialization]] and [[LFC is beta]] for math.

    `indicies` are the return of of `_matricies_to_order`
    """
    
    anti=anti_split(split)

    #clean up raw training data
    raw=data[["rep_id","cell_type","cre_id","umis_mpra_bc"]]
    raw=data[data[split]==subset]
    raw=raw.drop(columns=split)

    #collapse to summary statistics
    nb_stats = (
        raw.groupby(anti)
        .agg(
            mean_umis_mpra_bc=('umis_mpra_bc', 'mean'),
            var_umis_mpra_bc=('umis_mpra_bc', 'var')
        )
    )

    #get that set which are valid nb
    nb_stats["valid_nb"]=nb_stats["var_umis_mpra_bc"] > nb_stats["mean_umis_mpra_bc"]

    # compute rep-level counts including zeros
    zi_stats = (
        raw.groupby(['rep_id', anti])
        .agg(
            n=('umis_mpra_bc', 'count'),
            n_zero=('umis_mpra_bc', lambda x: (x == 0).sum())
        )
    )

    
    # merge global mean/var into rep-level
    zi_stats = zi_stats.reset_index().merge(nb_stats, on=anti, how='left')

    #so a row in zi_stats represents 'for replicate [rep_id] 
    #we see [cre_id] [n] times, of which [n_zero] data points are zero
    #that cre, across all reps, has a mean of [mean_umis_mpra_bc] and a 
    #variance of [var_umis_mpra_bc].' 
    
    
    
    zi_stats["gross_zero_prop"]=zi_stats["n_zero"]/zi_stats["n"]
    
    #We now want to figure out how many zeros are expected from the nb portion
    
    
    #compute nb params
    zi_stats["p"]=zi_stats["mean_umis_mpra_bc"]/zi_stats["var_umis_mpra_bc"]
    zi_stats["r"]=zi_stats["mean_umis_mpra_bc"]**2 / (zi_stats["var_umis_mpra_bc"] - zi_stats["mean_umis_mpra_bc"])
    
    zi_stats["nb_zero_prop"]=np.nan
    #fill out valid nb cases with zero proportion...
    zi_stats.loc[zi_stats["valid_nb"],"nb_zero_prop"]=zi_stats["p"]**zi_stats["r"]
    #fill out non-valid nb cases with zero portion using poisson
    zi_stats.loc[~zi_stats["valid_nb"],"nb_zero_prop"]=np.exp(-zi_stats["mean_umis_mpra_bc"])

    assert ~any(zi_stats["nb_zero_prop"].isna())


    zi_stats["zero_inflation"]=zi_stats["gross_zero_prop"]-zi_stats["nb_zero_prop"]
    zi_stats["zero_inflation"]=np.clip(zi_stats["zero_inflation"],0,1)
    
    zi_stats=zi_stats.groupby("rep_id")["zero_inflation"].mean()

    # now let's reindex
    zi_stats=zi_stats.reindex(indicies["zi_idx"])
    #logistic function
    zi_beta=1/(1 + np.exp(-zi_stats.to_numpy()))

    #now calculate nb betas
    working_nb=nb_stats.copy()
    ref=working_nb.loc["reference"]["mean_umis_mpra_bc"]
    working_nb["fc"]=working_nb["mean_umis_mpra_bc"]/ref
    working_nb["lfc"]=np.log(working_nb["fc"])
    working_nb["beta"]=working_nb["lfc"]
    working_nb["beta"].loc["reference"]=np.log(working_nb["mean_umis_mpra_bc"].loc["reference"])
    
    nb_betas=working_nb["beta"]

    
    #now calculate theta betas
    working_nb=nb_stats.copy()
    thetas=working_nb["mean_umis_mpra_bc"]**2/(working_nb["var_umis_mpra_bc"]-working_nb["mean_umis_mpra_bc"])
    thetas=thetas[working_nb["valid_nb"]]
    beta_theta=np.log(np.mean(thetas))

    #with our table constructed, we can now extract our actual parameter estimates
    #first, the betas for the nb portion

    nb_betas=nb_betas.reindex(indicies["nb_idx"])

    init={}
    init["x_mu"]=nb_betas.to_numpy().reshape(-1, 1)
    init["x_pi"]=zi_beta.reshape(-1, 1)
    init["theta"]=np.array([[beta_theta]])


    return init



@unimplemented
def skew_spread():
    """
    Creates a ground-truth dataframe of an scMPRA experiment
    that is meant to test skew
    (see readme for ground truth dataframe specification)
    """
    pass

def recombinator(primary,secondary):
    """
    All pairs of (All pairs of primary), secondary.
    two duplicate `secondary` entries in each element.
    """
    combos=itertools.combinations(primary,2)
    return [(i,j,k,k) for (i,j) in combos for k in secondary]

def activity_spread(cell_types:List[str],
    minimum:float,
    maximum:float,
    minp_value:float,
    total:int,
    frac_active:int,
    ct_specificity:float):
    """
    Creates a ground-truth dataframe of an scMPRA experiment
    with a controllable number of active CREs.
    (see readme for ground truth dataframe specification)
    Assumes experiment is interested in activity vs a known negative control,
    not skew. 
    
    This is useful to create datasets with a balance 
    of active and inactive CREs which roughly
    approx. real libraries. 

    `total` is the total number of CREs to create.
    `frac_active` is the fraction of the library that is active elements
    `ct_specificity` is how much cell-type specificity there is. Its 
    how many different values a given CRE will take across different cell-types. So 
    specificity=2 implies that, on average, there will be two different 
    means for each CRE across cell-types. 
    """
    ### first, we create the reference minP control ###
    reference=pd.DataFrame({"cell_type":cell_types})
    reference["true_mean"]=minp_value
    reference["cre_id"]="reference"

    ### second, we create a DF of all the inactive CREs. ###
    total_inactive =int(total*(1-frac_active))
    inactive_names=[f"inactive_{i}"for i in range(0,total_inactive)]
    inactive_tuples=[i for i in itertools.product(inactive_names,cell_types)]
    inactive=pd.DataFrame(inactive_tuples,columns=["cre_id","cell_type"])
    inactive["true_mean"]=minp_value
    
    ### third, we create a df of active elements. ###
    total_active=int(total*frac_active)
    active_names=[f"active_{i}"for i in range(0,total_active)]
    active=pd.DataFrame(active_names,columns=["cre_id"])
    active["cell_types"]=[cell_types for _ in range(len(active))]
    #we have a list of all cell types in each.
    #now we want to randomly decide which cell_types are the same
    #and which are different. 
    alpha = alpha_for_expected_groups(n=len(cell_types), K_target=ct_specificity)
    def _apply_crp(row):
        groups = sample_crp_groups(n=len(cell_types), alpha=alpha)
        return groups

    active["groups"]=active.apply(_apply_crp,axis=1)

    #calculate and print the percent of active CREs with any cell-type specificity
    def _different_groups(row):
        return len(set(row["groups"]))==1
    is_cell_type_specific=active.apply(_different_groups,axis=1)
    logger.info(f"{sum(is_cell_type_specific)/len(is_cell_type_specific)*100}% of active elements are not cell-type specific.")
    
    #generate the actual means for each group for each CRE. 
    def _generate_means(row):
        uniq_groups=list(set(row["groups"]))
        means=np.random.uniform(low=minimum,
            high=maximum,
            size=len(uniq_groups)
        )
        means=means.tolist()
        means_dict=dict(zip(uniq_groups,means))
        means_rep=[means_dict[i] for i in row["groups"]]
        return means_rep

    active["means"]=active.apply(_generate_means,axis=1)

    #dont need group ID anymore
    active=active.drop(columns="groups")

    #explode out means
    active=active.explode(["cell_types","means"])
    #rename to singular
    active=active.rename({"cell_types":"cell_type","means":"true_mean"},axis=1)
    
    # create final df!
    final_gt=pd.concat([reference,inactive,active])

    ### create hypothesis set ###
    #first, we want every CRE vs reference
    unique_CREs=final_gt["cre_id"].unique()
    idx=np.where(unique_CREs=="reference")
    unique_CREs=np.delete(unique_CREs,idx)
    all_cre_vs_ref=pd.DataFrame([(i,"reference",c,c) for i in unique_CREs for c in cell_types],
    columns=["comparison_CRE","reference_CRE","comparison_cell_type","reference_cell_type"])


    #now we want every cell type vs all other cell-types
    ct_tests=pd.DataFrame(recombinator(cell_types,unique_CREs),
    columns=["reference_cell_type","comparison_cell_type","comparison_CRE","reference_CRE"])

    hypothesis_set=HypothesisSet.from_dataframe(pd.concat([all_cre_vs_ref,ct_tests]))
    

    return (final_gt, hypothesis_set)

def simple_spread(cell_types:List[str],
    min:float,
    max:float,
    fineness:int=10,
    hypothesis_type:str="cartesian"):
    """
    TODO: remove hypothesis_type & make obligate cartesian.
    TODO: extract cartesian code to its own function.
    Create a ground truth dataframe tiling all cell-types.
    with synthetic CREs at a variety of strengths.

    Returns a tuple of (ground truth, hypothesis object) 

    Useful for simulation and power calculations.
    (see readme for ground truth dataframe specification)

    min, max are the min & max MPRA UMI / cell values.
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
    if hypothesis_type=="cartesian":
        unique_CREs=final_ground_truth["cre_id"].unique()
        unique_cell_types=final_ground_truth["cell_type"].unique()


        cre_comparisons=pd.DataFrame(
            recombinator(primary=unique_CREs,
                secondary=unique_cell_types),
            columns=["comparison_CRE","reference_CRE","comparison_cell_type","reference_cell_type"])

        cell_type_comparisons=pd.DataFrame(
            recombinator(primary=unique_cell_types,
                secondary=unique_CREs),
            columns=["comparison_cell_type","reference_cell_type","comparison_CRE","reference_CRE"])

        hypothesis_set=HypothesisSet.from_dataframe(
            pd.concat([cre_comparisons,cell_type_comparisons])
        )
    else:
        assert 1==2,"Unrecognized hypothesis set type!"

    
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
    
    def draw_nb(self,size, rng = np.random.default_rng()):
        """
        returns a 1d numpy vector of draws from the nb
        model of the object. 
        """
        #return rng.negative_binomial(n=self.r, p=self.p, size=size)
        return scipy.stats.nbinom.rvs(self.r,self.p,size=size,random_state=rng)
    
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

    total_uniq_mpra_bc:int=None

    by_cre_theta:float=None
    by_cell_type_theta:float=None
    theta:float=None

    zi:float=None
    by_cre_zi:float=None
    by_cell_type_zi:float=None

    reference_activity:float=None
    by_cell_type_reference_activity:float=None
    by_cre_reference_activity:float=None

    cells_per_cell_type:dict=None

    excess_tfection:float=None
    total_tfection:float=None
    total_uniq_mpra_bc:int=None

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
            type_map = {}
            # dump members as parquet files
            for name, value in self.__dict__.items():
                path = os.path.join(tmpdir, f"{name}.parquet")
                if isinstance(value, pd.DataFrame):
                    type_map[name] = "dataframe"
                    value.to_parquet(path, engine="pyarrow", index=True)
                elif isinstance(value,simple_count):
                    type_map[name] = "simplecount"
                    value.to_dataframe().to_parquet(path, engine="pyarrow", index=True)
                elif isinstance(value, pd.Series):
                    type_map[name] = "series"
                    value.to_frame(name).to_parquet(path, engine="pyarrow", index=True)
                else:
                    type_map[name] = "scalar"
                    pd.DataFrame({name: [value]}).to_parquet(path, engine="pyarrow", index=False)

            type_map_path = os.path.join(tmpdir, "types.json")
            with open(type_map_path, "w", encoding="utf-8") as f:
                json.dump(type_map, f)

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

            type_map = {}
            type_map_path = os.path.join(tmpdir, "types.json")
            if os.path.exists(type_map_path):
                with open(type_map_path, "r", encoding="utf-8") as f:
                    type_map = json.load(f)

            # load parquet members back
            for fname in os.listdir(tmpdir):
                if fname.endswith(".parquet"):
                    name = fname[:-8]
                    path = os.path.join(tmpdir, fname)
                    df = pd.read_parquet(path, engine="pyarrow")
                    val_type = type_map.get(name)
                    if val_type == "simplecount":
                        val = simple_count.from_dataframe(df)
                    elif val_type == "series":
                        val = df.iloc[:, 0]
                        val.name = df.columns[0] if len(df.columns) else name
                    elif val_type == "dataframe":
                        val = df
                    elif val_type == "scalar":
                        val = df.iloc[0, 0]
                    else:
                        if df.shape == (1, 1):
                            val = df.iloc[0, 0]
                        elif df.shape[1] == 1:
                            val = df.iloc[:, 0]
                        else:
                            val = df
                    setattr(ret, name, val)

        #re-initalize 
        if isinstance(ret.transfection_model, pd.DataFrame):
            ret.transfection_model=simple_count.from_dataframe(ret.transfection_model)
        if isinstance(ret.library_model, pd.DataFrame):
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
            
            ## reference ##
            if var=="by_cell_type_parameters":
                reference_mus=[]
                for key in getattr(inp,var).nb:
                    df=getattr(inp,var).nb[key].result()
                    assert len(df.loc["reference"])==1; "Multi reference"
                    reference_mus.append(df.loc["reference"]["mu"])
                ret.by_cell_type_reference_activity=np.mean(reference_mus)
            elif var=="by_cre_parameters":
                ret.by_cre_reference_activity=np.mean(getattr(inp,var).nb["reference"].result()["mu"].to_list())
            
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
                ret.by_cell_type_zi=zis.mean(axis=1)#.mean()
            elif var=="by_cre_parameters":
                ret.by_cre_zi=zis.mean(axis=1)#.mean()

        
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

        
        #save the preferred model direction
        ret.preferred=preferred
        
        #choose representative parameters based on preferred model direction
        if preferred=="by_cell_type":
            ret.zi=ret.by_cell_type_zi
            ret.theta=ret.by_cell_type_theta
            ret.reference_activity=ret.by_cell_type_reference_activity
        elif preferred=="by_cre":
            ret.zi=ret.by_cre_zi
            ret.theta=ret.by_cre_theta
            ret.reference_activity=by_cre_reference_activity
        else:
            assert False, "Unrecognized direction."

        #total number of MPRA barcodes
        ret.total_uniq_mpra_bc=len(inp.training_data.data["mpra_bc"].unique())
        
        #calculate post-hoc overtransfection. 
        tfection=inp.training_data.data.groupby(["rep_id","cell_bc"])["mpra_bc"].nunique().reset_index()
        tfection=tfection.rename({"mpra_bc":"unique_mpra_bc"},axis=1)
        observed=tfection["unique_mpra_bc"]
        tfection["tot_plasmid"]=np.log(1-observed/ret.total_uniq_mpra_bc)/np.log((ret.total_uniq_mpra_bc-1)/ret.total_uniq_mpra_bc)
        total_tfection=tfection['tot_plasmid'].sum()
        excess=total_tfection-tfection["unique_mpra_bc"].sum()
        logger.info(f"Computed a total of {excess} estimated collision events, out of a total of {total_tfection}, or {excess/total_tfection*100}%")
        ret.excess_tfection=excess
        ret.total_tfection=total_tfection
        
        return ret

from pathlib import Path

working_dir = Path(__file__).resolve().parent
SHENDURE_BOUNDS=Bounds.from_tgz(working_dir/"presets/shendure_bounds.tgz")
COHEN_BOUNDS=Bounds.from_tgz(working_dir/"presets/cohen_bounds.tgz")

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
    
    def flatten_overtransfection(self):
        """
        If you have simulated a dataset, it will probably have some degree of 
        overtransfection (same MPRA bc transfected into the same cell multiple times).
        This function flattens such events, as they would be observed in a real dataset.
        """
        if not self.metadata.get("synthetic"):
            raise ValueError("Can't flatten overtransfection on an emperical dataset.")
        
        if "overtransfection_flattened" in self.operations:
            logger.warning("Overtransfection flattening already performed. Skipping.")
            return
        
        groupby_columns=list(set(self.data.columns)-{"umis_mpra_bc"})
        
        self.data= self.data.groupby(groupby_columns).agg("sum").reset_index()
        self.operations.append("overtransfection_flattened")
    
    def overtransfected(self, log=True, threshold_pct=WARN_MULTI_TRANSFECTION_PERCENT):
        """
        Return True iff the overall percent of cells with >=1 multi-transfection
        (same mpra_bc observed >1 time in the same cell within a replicate)
        exceeds `threshold_pct`. Logging is optional. Also flags overtransfection 
        (or lack thereof) in metadata.

        Uses a scale-free metric: (# cells with >=1 dup) / (total # cells) * 100
        """

        if not self.metadata.get("synthetic"):
            raise ValueError("Can't directly test if an emperical dataset is overtranfected. Use a Bounds object to extract a transfection model, from which you can predict the degree of overtransfection.")
        
        if "overtransfection_flattened" in self.operations:
            raise ValueError("This dataset has already had its overtransfection flattened & so it can't be computed.")

        df = self.data
        # Count per (rep, cell, mpra_bc)
        triplet_counts = (
            df.groupby(["rep_id", "cell_bc", "mpra_bc"])
            .size()
        )

        # For each cell, did ANY barcode appear more than once?
        cell_has_dup = (
            triplet_counts.gt(1)
                        .groupby(level=["rep_id", "cell_bc"])
                        .any()
        )

        n_cells = int(cell_has_dup.size)
        n_cells_with_dup = int(cell_has_dup.sum())
        percent = (n_cells_with_dup / n_cells * 100.0) if n_cells else 0.0

        if log:
            msg = (f"{n_cells_with_dup}/{n_cells} cells "
                f"({percent:.3f}%) have ≥1 multi-transfection event.")
            if percent > threshold_pct:
                logger.warning(msg)
                logger.warning(
                    f"Multi-transfections exceed threshold of {threshold_pct:.3f}%!"
                )
            else:
                logger.info(msg)

        self.metadata["overtransfected"]=percent
        return percent > threshold_pct
    
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
        barcodes for each cre_id
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
        all_combos = self.data[['cell_type', 'cre_id']].drop_duplicates()

        # Compute dropped combos
        dropped_combos = pd.merge(all_combos, valid_combos, on=['cell_type', 'cre_id'], how='outer', indicator=True)
        dropped_combos = dropped_combos[dropped_combos['_merge'] == 'left_only'][['cell_type', 'cre_id']]

        # Warn if reference was filtered out
        ref_mask = (
            (dropped_combos["cell_type"] == "reference")
            | (dropped_combos["cre_id"] == "reference")
        )

        if ref_mask.any():
            n_ref = ref_mask.sum()
            logger.warning(
                f"ortho_filter removed {n_ref} combinations involving 'reference' "
            )

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

    X = X.astype(pd.SparseDtype("int", fill_value=0))
    
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



def _tensorzinb_fit(matricies,name,init_method="nb",init_vals=None):
    """
    Takes matricies & produces a single tensorzinb model.
    init_method takes "nb", "ones" or "pass".
    mom (method of moments) only implemented for by_cell_type models at the moment.
    """
    if init_method=="nb":
        zinbo = TensorZINB(endog=matricies["regressand"]["umis_mpra_bc"].to_numpy().squeeze(),
                        exog=matricies["nb_regressors"].to_numpy(),
                        exog_infl=matricies["zi_regressors"].to_numpy())
        result = zinbo.fit(return_history=True,init_method="nb")#reset_keras_session=True)
        del zinbo
    elif init_method=="pass":
        if not init_vals:
            raise ValueError("init_vals required for init_method=pass")

        zinbo = TensorZINB(endog=matricies["regressand"]["umis_mpra_bc"].to_numpy().squeeze(),
                    exog=matricies["nb_regressors"].to_numpy(),
                    exog_infl=matricies["zi_regressors"].to_numpy())
        
        result = zinbo.fit(return_history=True,init_weights=init_vals)

        del zinbo

    elif init_method=="ones":
        num_feat_zi = matricies["zi_regressors"].to_numpy().shape[1]
        num_feat_nb = matricies["nb_regressors"].to_numpy().shape[1]
        if matricies["regressand"]["umis_mpra_bc"].to_numpy().squeeze().ndim == 1:
            num_out = 1
        else:
            num_out = y.shape[1]
        
        ones_init = {}
        ones_init["x_mu"] = np.ones((num_feat_nb, num_out), dtype=np.float32)
        ones_init["x_pi"] = np.ones((num_feat_zi, num_out), dtype=np.float32)
        ones_init["theta"] = np.ones((1, num_out), dtype=np.float32)
        
        zinbo_ones = TensorZINB(endog=matricies["regressand"]["umis_mpra_bc"].to_numpy().squeeze(),
                    exog=matricies["nb_regressors"].to_numpy(),
                    exog_infl=matricies["zi_regressors"].to_numpy())
        
        result = zinbo_ones.fit(return_history=True,init_weights=ones_init)

        del zinbo_ones
    else:
        raise ValueError(f"Unrecognized initalization type {init_method}.")
    
    if pd.isnull(result["llf_total"]):
        logger.warning(f"Unconverged model in {name}.")
    
    return result

def standard_fit(client,data,split):
    """
    Takes an scMPRA object and produces a set of models along one axis,
    specified by split.
    """

    data=data.data
    levels=data[split].unique()

    if split=="cell_type":
        mat_resource={"CELL_DESIGN":1}
    else:
        mat_resource={"CRE_DESIGN":1}
    
    mats_futures = {
        t: client.submit(
            _smart_matrix,
            data=data[data[split]==t],
            split=split,
            resources=mat_resource
        )
        for t in levels
    }

    if split=="cell_type":
        #init_method="pass"
        #init_vals={
        #    t:client.submit(_mom_from_training_data, 
        #        data=data,
        #        split="cell_type",
        #        subset=t,
        #        indicies=client.submit(_matricies_to_order, matricies=mats_futures[t])
        #        )
        #    for t in levels
        #}
        #TEMP OVERRIDE
        init_method="nb"
        init_vals={t:None for t in levels}
    else:
        init_method="nb"
        init_vals={t:None for t in levels}

    tzinb_futures = {
        t: client.submit(
                _tensorzinb_fit,
                mats_futures[t],
                t,
                init_method=init_method,
                init_vals=init_vals[t],
                resources={'FIT': 1}
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

        self.wald_precomp = None  # WaldPrecomp or None

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

        if getattr(self, "wald_precomp", None) is None:
            with open(full_path/"wald_precomp.pkl","wb") as f:
                pickle.dump(None, f)
        else:
            self.wald_precomp.save(full_path/"wald_precomp.pkl")

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
   
        try:
            wp = WaldPrecomp.load(client, full_path/"wald_precomp.pkl")
        except Exception:
            wp = None
        setattr(ret_ortho, "wald_precomp", wp)

        return ret_ortho

    
    def _condense_dat(self,dat):
        if dat==None:#passed dat is none, look for cached training data.
            if self.training_data==None:
                raise RuntimeError("No training data supplied.")
            dat=self.training_data
        else:#dat is not none, save it
            self.training_data=dat.copy()
        return dat

    
    def fit_by_cre_models(self,client,dat=None):
        dat=self._condense_dat(dat)
        self.by_cre, self.by_cre_design=standard_fit(client,
                                                     dat,
                                                     split="cre_id")
        self.by_cre.label_regressors(client,self.by_cre_design)
        

        
    def fit_by_cell_type_models(self,client,dat=None):
        dat=self._condense_dat(dat)
        self.by_cell_type, self.by_cell_type_design=standard_fit(client,
                                                        dat,
                                                        split="cell_type")
        self.by_cell_type.label_regressors(client,self.by_cell_type_design)
        
    
    def criss_cross(self,client,dat):
        """
        Makes by_cre and by_cell_type models.
        """
        self.fit_by_cre_models(client=client,dat=dat)
        self.fit_by_cell_type_models(client=client,dat=dat)
        
        
    
    def extract_params(self,client):
        """
        Extracts parameters for all models in the object
        Silently passes either / both directions if not computed previous
        """

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

    def precompute_wald(
        self,
        client: Client,
        *,
        cov_method: str = "sandwich",
    ):
        """
        Compute and cache Wald precomputations (SEs, covariances, name maps) for
        every by-cell-type model and every by-CRE model in this ortho.
        Stores results in self.wald_precomp (as Futures); persists with save().

        Parameters
        ----------
        client : dask.distributed.Client
            Dask client for dispatching work.
        cov_method : {"sandwich", "opg"}
            - "sandwich": use H^{-1} J H^{-1} covariance (default).
            - "opg":      use OPG-only covariance J^{-1}.
        """
        if self.training_data is None:
            raise RuntimeError(
                "precompute_wald requires self.training_data to subset matrices."
            )

        if (self.by_cell_type is None) and (self.by_cre is None):
            raise RuntimeError(
                "precompute_wald requires at least one of by_cell_type, by_cre models."
            )

        # Map cov_method -> opg_only flag
        cov_method = cov_method.lower()
        if cov_method not in {"sandwich", "opg"}:
            raise ValueError("cov_method must be 'sandwich' or 'opg'")
        opg_only = (cov_method == "opg")

        by_ct = {}
        by_cr = {}
        #Very Lame retry hack due to extremely rare failures in _hessian_se
        #TODO: debug intermitant `AlreadyExistsError` properly once precompute_wald is faster.
        with dask.annotate(retries=10):
            if not self.by_cell_type is None:
                for ct in self.by_cell_type.model.keys():
                    model_f = self.by_cell_type.model[ct]
                    design_f = self.by_cell_type_design[ct]
                    df_ct = self.training_data.data[
                        self.training_data.data["cell_type"] == ct
                    ]
                    by_ct[ct] = client.submit(
                        _build_wald_precomp_for_subset,
                        model_f,
                        design_f,
                        df_ct,
                        opg_only=opg_only,   # <-- pass the flag
                    )
            else:
                by_ct=None

            if not self.by_cre is None:
                for cr in self.by_cre.model.keys():
                    model_f = self.by_cre.model[cr]
                    design_f = self.by_cre_design[cr]
                    df_cr = self.training_data.data[
                        self.training_data.data["cre_id"] == cr
                    ]
                    by_cr[cr] = client.submit(
                        _build_wald_precomp_for_subset,
                        model_f,
                        design_f,
                        df_cr,
                        opg_only=opg_only,   # <-- pass the flag
                    )
            else:
                by_cr=None

        self.wald_precomp = WaldPrecomp(by_cell_type=by_ct, by_cre=by_cr)

    def make_wald_eval_bundle(self) -> dict:
        """
        Build a small, pickle-friendly snapshot with everything the workers
        need to evaluate Wald tests. No Dask Futures inside.
        """
        if self.wald_precomp is None:
            raise RuntimeError("Call ortho.precompute_wald(...) first.")

        wp = self.wald_precomp.flattened_copy()  # resolve futures
        logger.info("flattened ortho successfully")

        by_ct = {}
        if wp.by_cell_type is None:
            by_ct=None
        else:
            for ct, entry in wp.by_cell_type.items():
                model = _to_plain(self.by_cell_type.model[ct])
                by_ct[str(ct)] = _pack_model_block(model, entry)

        by_cr = {}
        if wp.by_cre is None:
            by_cr=None
        else:
            for cr, entry in wp.by_cre.items():
                model = _to_plain(self.by_cre.model[cr])
                by_cr[str(cr)] = _pack_model_block(model, entry)

        return {"by_cell_type": by_ct, "by_cre": by_cr}

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
    dense = working.copy()
    for col in dense.columns:
        if pd.api.types.is_sparse(dense[col].dtype):
            dense[col] = dense[col].sparse.to_dense()
    working=dense
    working["r"]=working["theta"]
    working["sigmasquare"]=working["mu"]**2/working["r"]+working["mu"]
    working["p"]=working["mu"]/working["sigmasquare"]
    #handle case where mu is zero
    working.loc[working["mu"]==0.0,"p"]=0.0
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


def simulate_from_description(description):
    """
    Simulate from a description dataframe.
    Assumes input is one transfection event per row
    Removes ground-truth rows.
    """

    # Simulate NB and ZI in numpy
    r = description['r'].to_numpy(dtype=float)
    p = description['p'].to_numpy(dtype=float)
    zi = description['zi'].to_numpy(dtype=float)

    #print(description['r'].map(type).value_counts())
    #print(description['p'].map(type).value_counts())

    #print(f"type(r):{type(r)}; type(p): {type(p)}; r: {r} ; p:{p}")
    
    nb = np.random.negative_binomial(n=r, p=p)
    keep_mask = np.random.binomial(n=1, p=1 - zi)
    zinb = nb * keep_mask

    description['zinb_sample'] = zinb

    return description


class simulation_batch:
    """
    DEPRECATED
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
                raise ValueError("invalid choice")
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


class WaldPrecompEntry:
    """
    Minimal payload needed to evaluate Wald tests quickly for a single fitted model.
    All numpy; safe to pickle.
    """
    __slots__ = ("xmu_names", "se_x_mu", "cov_nb", "k_nb", "debug_msg")

    def __init__(self, xmu_names, se_x_mu, cov_nb, k_nb, debug_msg=None):
        self.xmu_names = list(map(str, xmu_names))
        self.se_x_mu   = np.asarray(se_x_mu, dtype=float)
        self.cov_nb    = np.asarray(cov_nb, dtype=float)
        self.k_nb      = int(k_nb)
        self.debug_msg = debug_msg  # string or None

    def name_to_idx(self):
        # build lazily to keep pickle small
        return {nm: j for j, nm in enumerate(self.xmu_names)}


class WaldPrecomp:
    """
    Mirrors the shape of `parameters`: dicts keyed by split level.
    Values are WaldPrecompEntry objects (or Futures thereof when live on a cluster).
    """
    def __init__(self, by_cell_type=None, by_cre=None):
        self.by_cell_type = by_cell_type or {}
        self.by_cre       = by_cre or {}

    def _unflatten_futures(self, client: Client):
        # wrap raw objects in futures, for symmetry with other classes
        for d in (self.by_cell_type, self.by_cre):
            if d is None:
                continue
            for k in list(d.keys()):
                d[k] = client.submit(lambda x: x, d[k])

    def flattened_copy(self):
        # gather futures to plain objects
        def _gather(d):
            if d is None:
                return None
            return {k: (v.result() if isinstance(v, Future) else v) for k, v in d.items()}
        return WaldPrecomp(by_cell_type=_gather(self.by_cell_type),
                           by_cre=_gather(self.by_cre))

    def save(self, path: str | Path):
        with open(path, "wb") as f:
            pickle.dump(self.flattened_copy(), f)

    @staticmethod
    def load(client: Client, path: str | Path) -> "WaldPrecomp":
        with open(path, "rb") as f:
            obj: WaldPrecomp = pickle.load(f)
        # re-wrap as futures so downstream access is uniform
        obj._unflatten_futures(client)
        return obj

def _normalize_cell_type_label(label, scmpra: scMPRA_data):
    """
    Map a user-facing cell type name to the internal label used in the data/models.
    If the dataset has a reference cell type set, that name is mapped to 'reference'.
    """
    if label is None or (isinstance(label, float) and pd.isna(label)):
        return None
    s = str(label)
    if scmpra.reference_cell_type and s == scmpra.reference_cell_type:
        return "reference"
    return s

def _normalize_cre_label(label, scmpra: scMPRA_data):
    """
    Map a user-facing CRE name to the internal label used in the data/models.
    If negative controls were flattened to 'reference', map those original names to 'reference'.
    """
    if label is None or (isinstance(label, float) and pd.isna(label)):
        return None
    s = str(label)
    if "cre_id_original" in scmpra.data.columns:
        df = scmpra.data[["cre_id", "cre_id_original"]].dropna()
        # Was this original label collapsed to 'reference'?
        if ((df["cre_id_original"] == s) & (df["cre_id"] == "reference")).any():
            return "reference"
    return s

class HypothesisSet:
    """
    A validated container for hypothesis rows.

    Columns:
      - comparison_CRE (str) [T]
      - comparison_cell_type (str) [T]
      - reference_CRE (str or NA) [F]
      - reference_cell_type (str or NA) [F]
      - meta (str or NA) [F]

    Rules:
      - If both reference_* are NA: interpret as "compare to zero" (activity vs 0).
            - maybe change this to use the references from the scMPRA data object? discuss what we want the appropriate default to be
      - If exactly one of reference_CRE / reference_cell_type is provided: malformed.
    """

    def __init__(self, df: pd.DataFrame):
        self.df = df.copy()
        self._coerce_and_validate()

    # -------- Construction / IO --------
    @classmethod
    def from_dataframe(cls, df: pd.DataFrame) -> "HypothesisSet":
        return cls(df)

    @classmethod
    def from_tsv(cls, path: str) -> "HypothesisSet":
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=True)
        return cls(df)

    @classmethod
    def from_csv(cls, path: str) -> "HypothesisSet":
        df = pd.read_csv(path, dtype=str, keep_default_na=True)
        return cls(df)

    # @classmethod
    ## Probably not needed for this class but it's here if we change our minds
    # def from_parquet(cls, path: str) -> "HypothesisSet":
    #     t = pq.read_table(path)
    #     df = t.to_pandas(types_mapper=pd.ArrowDtype)
    #     # convert ArrowDtype to pandas strings where appropriate
    #     for col in df.columns:
    #         if col in HYPOTHESIS_ALL:
    #             df[col] = df[col].astype("string")
    #     return cls(df)

    def to_dataframe(self) -> pd.DataFrame:
        return self.df.copy()

    def to_tsv(self, path: str) -> None:
        self.df.to_csv(path, sep="\t", index=False)

    def to_csv(self, path: str) -> None:
        self.df.to_csv(path, index=False)

    # def to_parquet(self, path: str) -> None:
    #     table = pa.Table.from_pandas(self.df, preserve_index=False)
    #     pq.write_table(table, path, compression="gzip")

    # -------- Validation / Utilities --------
    def _coerce_and_validate(self) -> None:
        # Ensure expected cols exist; extra cols are allowed (and kept)
        for col in HYPOTHESIS_REQUIRED | HYPOTHESIS_OPTIONAL:
            if col not in self.df.columns:
                self.df[col] = pd.Series([pd.NA] * len(self.df), dtype="string")

        # Cast to string dtype (allows NA)
        for col in HYPOTHESIS_ALL:
            self.df[col] = self.df[col].astype("string")

        # Required present
        missing_req = [c for c in HYPOTHESIS_REQUIRED if self.df[c].isna().any()]
        if missing_req:
            raise ValueError(f"Missing required entries in columns: {missing_req}")

        # Rule: either both reference_* are NA, or both are present (non-NA)
        # Debating whether this rule is actually necessary or not
        ref_cre_na = self.df["reference_CRE"].isna()
        ref_ct_na  = self.df["reference_cell_type"].isna()
        malformed_mask = (ref_cre_na ^ ref_ct_na)  # xor -> exactly one NA
        if malformed_mask.any():
            bad_rows = list(self.df.index[malformed_mask][:10])
            raise ValueError(
                "Malformed hypotheses: rows contain only one of reference_CRE / reference_cell_type. "
                f"Example indices: {bad_rows} (showing up to 10)."
            )

    def is_zero_reference(self) -> pd.Series:
        """Row-wise boolean: True if comparing against implicit zero (both reference_* are NA)."""
        return self.df["reference_CRE"].isna() & self.df["reference_cell_type"].isna()

    def add_meta(self, series_like) -> None:
        """Attach/overwrite a meta column (useful for painting plots)."""
        self.df["meta"] = pd.Series(series_like, index=self.df.index, dtype="string")

    def __len__(self):
        return len(self.df)


# ---------------------------------------------------------------------
# Convenience hypothesis builders (wrap one_versus_all)
# ---------------------------------------------------------------------


def make_by_celltype_hypotheses(
    *,
    comparison_cell_type: str,
    counts: "scMPRA_data",
    comparison_cres: "list[str] | str" = "all",
    reference_cre: str | None = "reference",
    meta: str | None = None,
) -> "HypothesisSet":
    """
    Build hypotheses that test many CREs within a single cell type
    (CRE varies; cell_type fixed). This is the natural input for the
    by-cell-type Wald test (CRE vs baseline CRE in that cell type).

    Examples:
        hs = make_by_celltype_hypotheses(
                comparison_cell_type="NeuroectodermBrain",
                counts=shendure,
                comparison_cres="all",
                reference_cre="reference",   # flattened minP/noP
                meta="emvar_screen")

    Notes:
        - We set BOTH reference columns per the table spec:
            reference_CRE = `reference_cre`
            reference_cell_type = `comparison_cell_type`
        - Passing reference_cre=None will generate the "compare‐to‐zero" flavor,
          but the current Wald code ignores zero and interprets baseline from the model.
    """
    if not isinstance(counts, scMPRA_data):
        raise TypeError("counts must be an scMPRA_data object.")

    cell_type = str(comparison_cell_type)

    # What CREs exist in this cell type?
    df = counts.data
    available = (
        df.loc[df["cell_type"] == cell_type, "cre_id"]
        .astype(str).unique().tolist()
    )

    if comparison_cres == "all":
        cand = available.copy()
        # If the negative control is labeled "reference", users usually
        # don’t want a “reference vs reference” row; drop it.
        if "reference" in cand:
            cand.remove("reference")
    else:
        cand = [str(x) for x in comparison_cres]
        missing = sorted(set(cand) - set(available))
        if missing:
            warnings.warn(f"[make_by_celltype_hypotheses] Skipping CREs not present in '{cell_type}': {missing}")
            cand = [c for c in cand if c in available]

    # Build (cre, cell_type) pairs for one_versus_all (comparison_on='cre')
    pairs = [(cre, cell_type) for cre in cand]

    hyp_df = one_versus_all(
        pairs,
        comparison_on="cre",
        reference_CRE=reference_cre,
        reference_cell_type=cell_type,
        meta=meta,
    )
    return HypothesisSet.from_dataframe(hyp_df)


def make_by_cre_hypotheses(
    *,
    comparison_cre: str,
    counts: "scMPRA_data",
    comparison_cell_types: "list[str] | str" = "all",
    reference_cell_type: str | None = None,
    meta: str | None = None,
) -> "HypothesisSet":
    """
    Build hypotheses that test many cell types for one CRE
    (cell_type varies; CRE fixed). This is the natural input for the
    by-CRE Wald test (cell_type vs baseline cell type for the same CRE).

    Examples:
        hs = make_by_cre_hypotheses(
                comparison_cre="CRE123",
                counts=shendure,
                comparison_cell_types="all",
                reference_cell_type="Pluripotent",
                meta="cell_specificity")

    Notes:
        - We set BOTH reference columns per the table spec:
            reference_CRE = `comparison_cre`
            reference_cell_type = provided (or inferred)
        - If `reference_cell_type` is not provided, we try:
            counts.reference_cell_type, then literal "reference" if present.
    """
    if not isinstance(counts, scMPRA_data):
        raise TypeError("counts must be an scMPRA_data object.")

    cre = str(comparison_cre)

    # What cell types exist for this CRE?
    df = counts.data
    available = (
        df.loc[df["cre_id"] == cre, "cell_type"]
        .astype(str).unique().tolist()
    )

    # Pick/validate reference cell type
    # Normalize the reference to internal labeling (e.g., "Pluripotent" -> "reference")
    ref_ct = reference_cell_type
    if ref_ct is None:
        if getattr(counts, "reference_cell_type", None):
            ref_ct = counts.reference_cell_type
        elif "reference" in available:
            ref_ct = "reference"
        else:
            raise ValueError(
                "reference_cell_type not provided and could not be inferred "
                "(counts.reference_cell_type not set and 'reference' not found for this CRE)."
            )
    ref_ct = _normalize_cell_type_label(ref_ct, counts)

    if comparison_cell_types == "all":
        cand = [ct for ct in available if ct != ref_ct]
    else:
        # also normalize any user-provided labels before comparing to 'available'
        cand = [_normalize_cell_type_label(x, counts) for x in comparison_cell_types]
        missing = sorted(set(cand) - set(available))
        if missing:
            warnings.warn(f"[make_by_cre_hypotheses] Skipping cell types not present for '{cre}': {missing}")
            cand = [c for c in cand if c in available]
        cand = [c for c in cand if c != ref_ct]

    # Build (cell_type, cre) pairs for one_versus_all (comparison_on='cell_type')
    pairs = [(ct, cre) for ct in cand]
    hyp_df = one_versus_all(
        pairs,
        comparison_on="cell_type",
        reference_CRE=cre,                 # same CRE = “within-CRE” contrast
        reference_cell_type=ref_ct,        # baseline cell type
        meta=meta,
    )
    return HypothesisSet.from_dataframe(hyp_df)

def make_all_by_celltype_hypotheses(
    *,
    counts: "scMPRA_data",
    reference_cre: str | None = "reference",
    meta: str | None = None,
    include_cell_types: "list[str] | None" = None,
    exclude_cell_types: "list[str] | None" = None,
) -> "HypothesisSet":
    """
    Build hypotheses for *every* cell type in the dataset, testing all CREs
    within each cell type against the provided baseline CRE (default: 'reference').

    Effectively: concat over cell types of
        make_by_celltype_hypotheses(comparison_cell_type=ct, comparison_cres="all")

    Parameters
    ----------
    counts : scMPRA_data
        Dataset.
    reference_cre : str or None
        Baseline CRE label used in the by-cell-type test (typically the flattened negative controls -> 'reference').
        If None, you get the "compare-to-zero" flavor, though the Wald code currently focuses on baseline contrasts.
    meta : str or None
        Optional label propagated to the 'meta' column.
    include_cell_types : list[str] or None
        If provided, restrict to this whitelist of cell types (user-facing labels OK; they’ll be normalized).
    exclude_cell_types : list[str] or None
        If provided, drop these cell types (applied after include).

    Returns
    -------
    HypothesisSet
        A validated hypothesis set spanning all requested cell types.
    """
    if not isinstance(counts, scMPRA_data):
        raise TypeError("counts must be an scMPRA_data object.")

    all_cts = sorted(map(str, counts.data["cell_type"].unique().tolist()))

    # Optional include/exclude
    if include_cell_types is not None:
        want = set(map(str, include_cell_types))
        all_cts = [ct for ct in all_cts if ct in want]
    if exclude_cell_types is not None:
        drop = set(map(str, exclude_cell_types))
        all_cts = [ct for ct in all_cts if ct not in drop]

    frames = []
    for ct in all_cts:
        # keep user-facing ct string; per-row normalization happens in the test functions
        hs_ct = make_by_celltype_hypotheses(
            comparison_cell_type=ct,
            counts=counts,
            comparison_cres="all",
            reference_cre=reference_cre,
            meta=meta,
        )
        frames.append(hs_ct.to_dataframe())

    if not frames:
        return HypothesisSet.from_dataframe(pd.DataFrame(columns=list(HYPOTHESIS_ALL)))

    big = pd.concat(frames, ignore_index=True)
    return HypothesisSet.from_dataframe(big)


def make_all_by_cre_hypotheses(
    *,
    counts: "scMPRA_data",
    reference_cell_type: str | None = None,
    meta: str | None = None,
    include_cres: "list[str] | None" = None,
    exclude_cres: "list[str] | None" = None,
    drop_reference_cre: bool = True,
) -> "HypothesisSet":
    """
    Build hypotheses for *every* CRE in the dataset, testing all cell types
    within each CRE against a baseline cell type (defaults to the dataset’s
    reference cell type if available, otherwise tries literal 'reference').

    Effectively: concat over CREs of
        make_by_cre_hypotheses(comparison_cre=cre, comparison_cell_types='all')

    Parameters
    ----------
    counts : scMPRA_data
        Dataset.
    reference_cell_type : str or None
        Baseline cell type label for within-CRE contrasts. If None, we try
        counts.reference_cell_type, else literal 'reference' if present. 
        (Same behavior as make_by_cre_hypotheses.)
    meta : str or None
        Optional label propagated to the 'meta' column.
    include_cres : list[str] or None
        If provided, restrict to this whitelist of CREs.
    exclude_cres : list[str] or None
        If provided, drop these CREs (applied after include).
    drop_reference_cre : bool
        If True (default), skip the collapsed negative-control CRE named 'reference'
        to avoid generating “reference vs baseline cell type” rows.

    Returns
    -------
    HypothesisSet
        A validated hypothesis set spanning all requested CREs.
    """
    if not isinstance(counts, scMPRA_data):
        raise TypeError("counts must be an scMPRA_data object.")

    all_cres = sorted(map(str, counts.data["cre_id"].unique().tolist()))
    if drop_reference_cre and "reference" in all_cres:
        all_cres.remove("reference")

    # Optional include/exclude
    if include_cres is not None:
        want = set(map(str, include_cres))
        all_cres = [c for c in all_cres if c in want]
    if exclude_cres is not None:
        drop = set(map(str, exclude_cres))
        all_cres = [c for c in all_cres if c not in drop]

    frames = []
    for cre in all_cres:
        hs_cre = make_by_cre_hypotheses(
            comparison_cre=cre,
            counts=counts,
            comparison_cell_types="all",
            reference_cell_type=reference_cell_type,
            meta=meta,
        )
        frames.append(hs_cre.to_dataframe())

    if not frames:
        return HypothesisSet.from_dataframe(pd.DataFrame(columns=list(HYPOTHESIS_ALL)))

    big = pd.concat(frames, ignore_index=True)
    return HypothesisSet.from_dataframe(big)

def make_bootstrap_activity_hypotheses(
    *,
    counts: "scMPRA_data",
    comparison_cres: "list[str] | str" = "all",
    controls: "list[str] | str | None" = None,
    meta: str | None = "bootstrap_activity",
) -> "HypothesisSet":
    """
    Build a hypothesis set for the bootstrap activity test:
      - ONE ROW PER CRE (not per cell type)
      - comparison_cell_type = "ALL" (sentinel; ignored by the test)
      - reference_cell_type = "ALL" (to satisfy the 'both present or both NA' rule)
      - reference_CRE carries the control label(s); the test will union all unique controls
        present in the HS when building its bundle.

    Params
    ------
    counts : scMPRA_data
        Your dataset.
    comparison_cres : list[str] | "all"
        Which CREs to test. "all" = every CRE in counts (we’ll drop any that are also controls).
    controls : list[str] | str | None
        Control CRE label(s). If None, we try to infer "reference" if present.
    meta : str | None
        Optional meta label.

    Returns
    -------
    HypothesisSet
    """
    if not hasattr(counts, "data"):
        raise TypeError("counts must be an scMPRA_data object with a `.data` DataFrame.")

    df = counts.data
    all_cres = sorted(map(str, df["cre_id"].unique().tolist()))

    # Controls: explicit -> as provided; else try 'reference'
    if controls is None:
        controls_list = ["reference"] if "reference" in all_cres else []
        if not controls_list:
            raise ValueError("controls not provided and 'reference' not present in counts.")
    elif isinstance(controls, str):
        controls_list = [controls]
    else:
        controls_list = list(map(str, controls))

    # Comparison CREs
    if comparison_cres == "all":
        comp = [c for c in all_cres if c not in set(controls_list)]
    else:
        comp = [str(c) for c in comparison_cres if str(c) in all_cres and str(c) not in set(controls_list)]
    if not comp:
        raise ValueError("No comparison CREs remain after excluding controls.")

    # Build rows: one row per CRE; put the controls into reference_CRE.
    # It's okay if multiple distinct control labels appear across rows; the bundle
    # will union the unique set from reference_CRE.
    rows = []
    for cre in comp:
        for ctrl in controls_list:
            rows.append(
                {
                    "comparison_CRE": cre,
                    "comparison_cell_type": "MAX",
                    "reference_CRE": ctrl,
                    "reference_cell_type": "MAX",
                    "meta": meta,
                }
            )

    return HypothesisSet.from_dataframe(pd.DataFrame(rows))


def coerce_bootstrap_activity_from_hs(hs: "HypothesisSet") -> "HypothesisSet":
    """
    If you already have a hypothesis set (e.g., from make_all_by_cre_hypotheses),
    collapse it to the bootstrap-activity shape:
      - de-duplicate to one row per comparison_CRE
      - set comparison_cell_type = reference_cell_type = "ALL"
      - keep reference_CRE as-is (we’ll union controls later)
    """
    df = hs.to_dataframe().copy()

    # prefer rows that already have a reference_CRE
    df = df.dropna(subset=["comparison_CRE"])
    # de-duplicate by CRE, keeping first
    df = df.sort_index().drop_duplicates(subset=["comparison_CRE"], keep="first")

    df["comparison_cell_type"] = "MAX"
    df["reference_cell_type"] = "MAX"

    # If reference_CRE is entirely NA here, raise (we need it to infer controls)
    if df["reference_CRE"].isna().all():
        raise ValueError("Cannot infer controls: 'reference_CRE' is NA for all rows in the provided HS.")

    return HypothesisSet.from_dataframe(df[["comparison_CRE", "comparison_cell_type", "reference_CRE", "reference_cell_type", "meta"]])
class ResultSet(HypothesisSet):
    """
    Extends HypothesisSet with result columns:
      - test_type (str)     [T]
      - test_statistic (float) [T]
      - p_value (float)     [T]
      - fold_change (float) [T]
      - bh_p (float)        [T]
      - flattened (bool)    [T]
    """

    def __init__(self, df: pd.DataFrame):
        super().__init__(df)
        self._validate_results_block()

    def _validate_results_block(self) -> None:
        missing = [c for c in RESULT_REQUIRED if c not in self.df.columns]
        if missing:
            raise ValueError(f"Results table missing required columns: {missing}")

        # Basic dtype coercions (tolerant)
        self.df["test_type"] = self.df["test_type"].astype("string")
        for c in ("test_statistic", "p_value", "fold_change", "bh_p"):
            self.df[c] = pd.to_numeric(self.df[c], errors="coerce")
        # flattened must be boolean-ish
        if self.df["flattened"].dtype != bool:
            self.df["flattened"] = self.df["flattened"].astype("bool")

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame) -> "ResultSet":
        return cls(df)

def _bh_adjust(pvals: pd.Series) -> pd.Series:
    """Benjamini–Hochberg FDR control on a 1D array-like of p-values."""
    p = pd.Series(pvals, dtype=float).fillna(1.0).clip(0, 1)
    # alpha only affects the boolean 'reject'; the adjusted p-values are
    # independent of alpha for BH, so any alpha is fine here.
    _, p_adj = fdrcorrection(p.values, alpha=0.05, is_sorted=False)
    return pd.Series(p_adj, index=p.index)
    

def canonicalize_hypotheses(hs: HypothesisSet, scmpra: scMPRA_data, inplace: bool = False) -> HypothesisSet:
    df = hs.df if inplace else hs.df.copy()

    ct_cols  = ["comparison_cell_type", "reference_cell_type"]
    cre_cols = ["comparison_CRE", "reference_CRE"]

    # 1) Cell types: reference name -> "reference" (vectorized)
    ref_ct = getattr(scmpra, "reference_cell_type", None)
    if ref_ct is not None:
        # mask is a DataFrame of booleans with same shape as df[ct_cols]
        m = df[ct_cols].eq(ref_ct)
        if m.values.any():
            df[ct_cols] = df[ct_cols].mask(m, "reference")

    # 2) CREs: any original CRE that was flattened to "reference" -> "reference" (vectorized)
    #    We only need the set of originals that ended up as 'reference'
    if "cre_id_original" in getattr(scmpra, "data", pd.DataFrame()).columns:
        collapsed = scmpra.data.loc[
            scmpra.data["cre_id"] == "reference", "cre_id_original"
        ].astype(str).unique()
        if len(collapsed) > 0:
            collapsed_set = set(collapsed)
            m = df[cre_cols].isin(collapsed_set)
            if m.values.any():
                df[cre_cols] = df[cre_cols].mask(m, "reference")

    # keep dtype consistent
    for col in ct_cols + cre_cols:
        df[col] = df[col].astype("string")

    return hs if inplace else HypothesisSet.from_dataframe(df)

def _to_plain(obj):
    # resolve dask Futures to plain objects
    from dask.distributed import Future
    return obj.result() if isinstance(obj, Future) else obj






# ---- Wald Test -------------------------------------
def _zinb_loglik_tf(params, exog, exog_infl, endog, return_per_obs=False):
    """
    TF implementation of the ZINB log-likelihood.

    params = concat([x_mu, x_pi, log_theta])

    If return_per_obs=False (default):
        returns a scalar total log-likelihood (previous existing behavior).

    If return_per_obs=True:
        returns (ll_sum, ll_vec), where:
          - ll_vec is shape (N,)   per-observation log-likelihoods (up to const.)
          - ll_sum is scalar sum(ll_vec)
    """
    N = tf.cast(tf.shape(endog)[0], tf.float64)

    num_features      = tf.shape(exog)[1]
    num_infl_features = tf.shape(exog_infl)[1]

    # Split parameter vector
    x_mu      = params[:num_features]
    x_pi      = params[num_features:num_features + num_infl_features]
    raw_log_theta = params[-1]

    # --- CLIPPED linear predictors and dispersion (same as prev. current version) ---

    # μ part: η_mu = X β_mu, clipped before exp
    eta_mu = tf.matmul(exog, tf.expand_dims(x_mu, axis=-1))
    eta_mu = tf.clip_by_value(eta_mu, -20.0, 20.0)
    mu     = tf.exp(eta_mu)

    # π logits for ZI part
    pi_logits = tf.matmul(exog_infl, tf.expand_dims(x_pi, axis=-1))
    pi_logits = tf.clip_by_value(pi_logits, -20.0, 20.0)

    # θ (dispersion) bounded away from 0/∞ in log-space
    log_theta = tf.clip_by_value(raw_log_theta, -10.0, 10.0)
    theta     = tf.exp(log_theta)

    # zero-inflation logits -> log(q0), log(q1)
    log_q0 = -tf.nn.softplus(-pi_logits)
    log_q1 = log_q0 - pi_logits

    y = tf.cast(endog, tf.float64)

    # NB log-likelihood (for y>0)
    t1 = tf.math.lgamma(y + theta)
    t2 = -tf.math.lgamma(theta)
    t3 = theta * log_theta
    t4 = y * tf.math.log(mu + 1e-8)
    ty = tf.math.log(mu + theta + 1e-8)
    t5 = -(theta + y) * ty
    nb_case = t1 + t2 + t3 + t4 + t5 + log_q1

    # Zero case
    p1        = theta * (log_theta - ty) + log_q1
    zero_case = tf.reduce_logsumexp(tf.stack([log_q0, p1], axis=0), axis=0)

    # ll per observation (up to the log-factorial constant)
    ll = tf.where(y < 1e-8, zero_case, nb_case)   # shape (N, 1)

    if return_per_obs:
        # Per-observation loglik (we don't include log(y!) since it drops out of scores)
        ll_vec = tf.reshape(ll, [-1])           # (N,)
        ll_sum = tf.reduce_sum(ll_vec)          # scalar
        return ll_sum, ll_vec

    # ---- existing “total log-likelihood” reduction for compatibility ----
    mean_neg_ll = -tf.reduce_mean(ll, axis=0)                   # (num_outputs,)
    log_fact    = tf.reduce_sum(tf.math.lgamma(endog + 1), axis=0)
    llfs        = -(mean_neg_ll * N + log_fact)                 # (num_outputs,)
    log_likelihood = tf.reduce_sum(llfs)
    return log_likelihood

def _setup_params_from_fit(zinb_model_fit):
    """
    Extract params vector and [ optional, commented out for now TF variable] (x_mu, x_pi, theta) from a single
    fitted TensorZINB result dict with labeled weights.
    """
    x_mu = zinb_model_fit['weights']['x_mu']
    x_pi = zinb_model_fit['weights']['x_pi']
    log_theta = zinb_model_fit['weights']['theta'].flatten()  # already log-space

    # N.B. x_mu and x_pi are already pandas Series with names; keep arrays here
    params = np.concatenate([np.asarray(x_mu).ravel(),
                             np.asarray(x_pi).ravel(),
                             np.asarray(log_theta).ravel()])
    # params_tensor = tf.Variable(params, dtype=tf64)
    return params.astype(np.float64, copy=False) #, params_tensor

def _pack_model_block(model_dict, entry, *, include_cov_nb=False):
    w = model_dict['weights']['x_mu']
    d = {
        "xmu_names": list(w.index),
        "xmu":      np.asarray(w).ravel(),
        "se_x_mu":  np.asarray(entry.se_x_mu),
        "k_nb":     int(entry.k_nb),
        "debug_msg": entry.debug_msg,
    }
    if include_cov_nb:
        d["cov_nb"] = np.asarray(entry.cov_nb)
    return d

def _model_matrices_for_subset(df_subset, design_dict, nb_formula, zi_formula):
    """
    Build endog/exog/exog_infl as numpy and their TF constants for a given subset.
    """
    # y, X = Formula(nb_formula).get_model_matrix(df_subset, output='pandas')
    # Z = Formula(zi_formula).get_model_matrix(df_subset, output='pandas')
    X, y, Z = design_dict['nb_regressors'], design_dict['regressand'], design_dict['zi_regressors']

    endog = y.to_numpy().reshape((-1, 1))
    exog = X.to_numpy()
    exog_infl = Z.to_numpy()

    # exog_tensor = tf.constant(exog, dtype=tf.float64)
    # endog_tensor = tf.constant(endog, dtype=tf.float64)
    # exog_infl_tensor = tf.constant(exog_infl, dtype=tf.float64)
    return (y, X, Z), (endog.astype(np.float64, copy=False), exog.astype(np.float64, copy=False), exog_infl.astype(np.float64, copy=False)) # , (endog_tensor, exog_tensor, exog_infl_tensor)

def _estimate_cov_se(
    params_np,
    exog_np,
    infl_np,
    endog_np,
    *,
    ridge_h=1e-8,
    ridge_j=1e-12,
    batch_size=4096,
    opg_only=False,
):
    """
    Compute covariance via either:
      - Sandwich: cov = H^{-1} J H^{-1}, where
          H = -∂^2 log L / ∂θ∂θᵀ  (observed information)
          J = Σ_i s_i s_iᵀ        (outer product of per-observation scores)
      - OPG-only: cov ≈ J^{-1} (Outer Product of Gradients)

    Uses TF1 graph mode; float64 throughout. Batches per-observation scores for memory safety.

    Parameters
    ----------
    params_np : np.ndarray
        1D parameter vector θ = concat([x_mu, x_pi, log_theta]).
    exog_np : np.ndarray
        NB regressors (N × K_nb).
    infl_np : np.ndarray
        ZI regressors (N × K_zi).
    endog_np : np.ndarray
        Response counts (N × 1).
    ridge_h : float
        Small ridge added to diag(H_info) for stability (sandwich mode).
    ridge_j : float
        Small ridge added to diag(J) for stability.
    batch_size : int
        Batch size for per-observation score computation.
    opg_only : bool
        If True, skip the Hessian and use OPG-only covariance (J^{-1});
        if False, use full sandwich H^{-1} J H^{-1}.
    """
    
    g = tf.Graph()
    
    with g.as_default():
        # Shapes (placeholders accept variable N along axis 0)
        P = int(params_np.size)
        params = tf.compat.v1.placeholder(tf.float64, shape=[P],                           name="params")
        exog   = tf.compat.v1.placeholder(tf.float64, shape=[None, exog_np.shape[1]],      name="exog")
        infl   = tf.compat.v1.placeholder(tf.float64, shape=[None, infl_np.shape[1]],      name="infl")
        endog  = tf.compat.v1.placeholder(tf.float64, shape=[None, 1],                     name="endog")

        # total + per-observation loglik
        ll_sum, ll_vec = _zinb_loglik_tf(params, exog, infl, endog, return_per_obs=True)

        # Hessian only if we're not in OPG-only mode
        if not opg_only:
            H = tf.hessians(ll_sum, params)[0]  # (P, P)

        # Per-observation score vectors: s_i = ∂ ll_i / ∂ params
        def _grad_one(li):
            gi = tf.gradients(li, params)[0]
            # Replace NaNs/Infs with 0 in the graph to avoid propagating nastiness
            return tf.where(tf.math.is_finite(gi), gi, tf.zeros_like(gi))

        S = tf.map_fn(_grad_one, ll_vec, dtype=tf.float64)  # shape (N, P)

    with tf.compat.v1.Session(graph=g) as sess:
        sess.run(tf.compat.v1.global_variables_initializer())

        # Hessian evaluation (sandwich mode only)
        if not opg_only:
            H_np = sess.run(
                H,
                feed_dict={params: params_np, exog: exog_np, infl: infl_np, endog: endog_np},
            )
            if not np.all(np.isfinite(H_np)):
                raise FloatingPointError("Sandwich path: non-finite Hessian")

        # Batch J accumulation
        P = params_np.size
        J = np.zeros((P, P), dtype=np.float64)
        N = endog_np.shape[0]

        for start in range(0, N, batch_size):
            stop = min(start + batch_size, N)
            S_b = sess.run(
                S,
                feed_dict={
                    params: params_np,
                    exog:   exog_np[start:stop],
                    infl:   infl_np[start:stop],
                    endog:  endog_np[start:stop],
                },
            )
            # Safety: zero out non-finite scores
            S_b[~np.isfinite(S_b)] = 0.0
            J += S_b.T @ S_b

    # HC1 correction and ridge on J
    n, p = endog_np.shape[0], P
    if n > p:
        J *= (n / (n - p))
    J[np.diag_indices_from(J)] += ridge_j

    # ---------- OPG-only path ----------
    if opg_only:
        used = "inv"
        try:
            if not np.all(np.isfinite(J)):
                raise FloatingPointError("OPG path: non-finite J")
            cov = np.linalg.inv(J)
        except (LinAlgError, FloatingPointError):
            used = "pinv"
            cov = la.pinvh(J, check_finite=False)

        cov = 0.5 * (cov + cov.T)
        se  = np.sqrt(np.clip(np.diag(cov), 0.0, np.inf))

        # Compact diagnostics for OPG
        try:
            evals_J = la.eigvalsh(J, check_finite=False)
            jmin    = float(evals_J.min())
            jmax    = float(evals_J.max())
            condJ   = (np.inf if jmin <= 0 else jmax / jmin)
            diag = (
                f"opg: Jmin={jmin:.2e},Jmax={jmax:.2e},condJ={condJ:.2e}; "
                f"sol={used}"
            )
        except Exception:
            diag = f"opg: sol={used}"

        return se, cov, diag

    # ---------- Sandwich path (uses H as well) ----------
    H_info = -H_np
    H_info = 0.5 * (H_info + H_info.T)
    H_info[np.diag_indices_from(H_info)] += ridge_h

    used = "chol"
    try:
        if not np.all(np.isfinite(H_info)):
            raise FloatingPointError("non-finite H_info")
        c, low = la.cho_factor(H_info, check_finite=False)
        X   = la.cho_solve((c, low), J, check_finite=False)
        cov = la.cho_solve((c, low), X, check_finite=False)
    except (LinAlgError, ValueError, FloatingPointError):
        used = "solve"
        try:
            X   = la.solve(H_info, J, assume_a="sym", check_finite=False)
            cov = la.solve(H_info, X, assume_a="sym", check_finite=False)
        except (LinAlgError, ValueError):
            used = "pinv"
            H_pinv = la.pinvh(H_info, check_finite=False)
            cov    = H_pinv @ J @ H_pinv

    cov = 0.5 * (cov + cov.T)
    se  = np.sqrt(np.clip(np.diag(cov), 0.0, np.inf))

    # Compact diagnostics for sandwich
    try:
        evals_H = la.eigvalsh(H_info, check_finite=False)
        evals_J = la.eigvalsh(J,      check_finite=False)
        hmin    = float(evals_H.min())
        hmax    = float(evals_H.max())
        jmin    = float(evals_J.min())
        jmax    = float(evals_J.max())
        condH   = (np.inf if hmin <= 0 else hmax / hmin)
        condJ   = (np.inf if jmin <= 0 else jmax / jmin)
        diag = (
            f"sandwich: Hmin={hmin:.2e},Hmax={hmax:.2e},condH={condH:.2e}; "
            f"Jmin={jmin:.2e},Jmax={jmax:.2e},condJ={condJ:.2e}; sol={used}"
        )
    except Exception:
        diag = f"sandwich: sol={used}"

    return se, cov, diag
def _hessian_se_graph(params_np, exog_np, infl_np, endog_np, *, ridge=1e-8):
    """
    Graph-mode only:
      - builds a self-contained graph,
      - evaluates gradients/Hessian in a local Session,
      - adds tiny ridge to stabilize inversion.
    """
    g = tf.Graph()
    with g.as_default():
        P = int(params_np.size)
        params = tf.compat.v1.placeholder(tf.float64, shape=[P], name="params")
        exog   = tf.compat.v1.placeholder(tf.float64, shape=[None, exog_np.shape[1]], name="exog")
        infl   = tf.compat.v1.placeholder(tf.float64, shape=[None, infl_np.shape[1]], name="infl")
        endog  = tf.compat.v1.placeholder(tf.float64, shape=[None, 1], name="endog")

        ll = _zinb_loglik_tf(params, exog, infl, endog)   

        grad = tf.gradients(ll, params)[0]                # (P,)
        # Build Hessian by differentiating each grad component
        hcols = [tf.gradients(grad[i], params)[0] for i in range(P)]
        hess  = tf.stack(hcols, axis=1)                   # (P, P)

        # cfg = tf.compat.v1.ConfigProto(
        #     intra_op_parallelism_threads=1,
        #     inter_op_parallelism_threads=1,
        #     allow_soft_placement=True,
        # )
        with tf.compat.v1.Session(graph=g) as sess: # , config=cfg
            sess.run(tf.compat.v1.global_variables_initializer())
            G, H = sess.run(
                [grad, hess],
                feed_dict={
                    params: params_np,
                    exog:   exog_np,
                    infl:   infl_np,
                    endog:  endog_np,
                },
            )

    # Check for NaNs in gradient or Hessian
    if np.isnan(G).any():
        uid = str(uuid.uuid4())
        for name, arr in [
            ("G", G),
            ("params", params_np),
            ("exog", exog_np),
            ("infl", infl_np),
            ("endog", endog_np),
        ]:
            np.save(f"{uid}_{name}.npy", arr)

        raise FloatingPointError(
            f"NaNs detected in gradient; dumped arrays to disk with prefix {uid}"
        )
    if np.isnan(H).any():
        raise FloatingPointError("NaNs detected in Hessian")

    # Check for non-finite values (Inf / -Inf)
    if not np.all(np.isfinite(G)):
        raise FloatingPointError("non-finite gradient")
    if not np.all(np.isfinite(H)):
        raise FloatingPointError("non-finite Hessian")

    H_info = -H
    H_info[np.diag_indices_from(H_info)] += ridge

    try:
        cov = np.linalg.inv(H_info)
    except LinAlgError:
        cov = la.pinvh(H_info)

    se = np.sqrt(np.diag(cov))
    return se, cov, H

def _slice_se(standard_errors, exog_cols, infl_cols):
    """
    Split the stacked SE vector into (se_x_mu, se_x_pi, se_theta) by actual sizes.
    """
    k_nb = len(exog_cols)
    k_zi = len(infl_cols)
    se_x_mu = standard_errors[:k_nb]
    se_x_pi = standard_errors[k_nb:k_nb+k_zi]
    se_theta = standard_errors[k_nb+k_zi:]
    return se_x_mu, se_x_pi, se_theta

def _build_wald_precomp_for_subset(
    model_dict,
    design_dict,
    df_subset,
    *,
    opg_only: bool = False,
) -> WaldPrecompEntry:
    """
    Worker-safe function: builds X/Z/y, computes covariance/SEs for the NB
    block via either sandwich or OPG, splits SE, returns WaldPrecompEntry.

    Parameters
    ----------
    opg_only : bool
        If True, use OPG-only covariance (J^{-1});
        if False, use full sandwich H^{-1} J H^{-1}.
    """
    _require_tensorflow()

    nb_formula, zi_formula = design_dict["nb_formula"], design_dict["zi_formula"]
    (_, X, Z), (endog_np, exog_np, infl_np) = _model_matrices_for_subset(
        df_subset, design_dict, nb_formula, zi_formula
    )
    params_np = _setup_params_from_fit(model_dict)

    try:
        # se_all, cov = _hessian_se_graph(params_np, exog_np, infl_np, endog_np, ridge=1e-8)
        se_all, cov, diag = _estimate_cov_se(
            params_np,
            exog_np,
            infl_np,
            endog_np,
            ridge_h=1e-8,
            ridge_j=1e-12,
            opg_only=opg_only,
        )
        se_x_mu, _, _ = _slice_se(se_all, X.columns, Z.columns)
        k_nb = len(X.columns)
        cov_nb = cov[:k_nb, :k_nb]

        # Diagnostics
        zero_var_cols = [
            c
            for c in X.columns
            if np.allclose(exog_np[:, X.columns.get_loc(c)], 0)
        ]

        cov_has_nan = np.isnan(cov_nb).any()
        cov_has_inf = np.isinf(cov_nb).any()
        se_has_nan = np.isnan(se_all).any()
        se_has_inf = np.isinf(se_all).any()

        if not (cov_has_nan or cov_has_inf):
            cond_nb = np.linalg.cond(cov_nb)
        else:
            cond_nb = np.inf

        issues = []
        if se_has_nan:
            issues.append("NaNs in SEs")
        if se_has_inf:
            issues.append("Infs in SEs")
        if cov_has_nan:
            issues.append("NaNs in cov_nb")
        if cov_has_inf:
            issues.append("Infs in cov_nb")
        if len(zero_var_cols) > 0:
            issues.append(
                f"zero-var in X: {zero_var_cols[:3]}"
                f"{'...' if len(zero_var_cols) > 3 else ''}"
            )
        if cond_nb > 1e12 and not (cov_has_nan or cov_has_inf):
            issues.append(f"ill-conditioned cov_nb (cond={cond_nb:.2e})")

        # `diag` already encodes "sandwich: ..." or "opg: ..."
        if issues:
            debug_msg = diag + "; " + "; ".join(issues)
        else:
            debug_msg = diag
    except Exception as e:
        # Any failure: fill NaNs and report which path failed
        se_x_mu = np.full(len(X.columns), np.nan)
        cov_nb = np.full((len(X.columns), len(X.columns)), np.nan)
        mode = "opg" if opg_only else "sandwich"
        debug_msg = f"{mode} failure: {type(e).__name__}: {e}"

    # Collinearity check
    try:
        rank_X = np.linalg.matrix_rank(X)
        rank_Z = np.linalg.matrix_rank(Z)
        if rank_X < X.shape[1]:
            debug_msg += f"; colinearity in X (rank {rank_X}/{X.shape[1]})"
        if rank_Z < Z.shape[1]:
            debug_msg += f"; colinearity in Z (rank {rank_Z}/{Z.shape[1]})"
    except Exception as e:
        debug_msg += f"; colinearity check failed: {e.__class__.__name__}"

    xmu_names = list(model_dict["weights"]["x_mu"].index)
    return WaldPrecompEntry(
        xmu_names=xmu_names,
        se_x_mu=se_x_mu,
        cov_nb=cov_nb,
        k_nb=len(X.columns),
        debug_msg=debug_msg,
    )

def _wald_by_celltype_row(row: dict, bundle: dict):
    ct  = row["comparison_cell_type"]
    cre = row["comparison_CRE"]
    blk = bundle["by_cell_type"].get(ct)

    if blk is None:
        return {
            "test_statistic": np.nan,
            "p_value": np.nan,
            "fold_change": np.nan,
            "flattened": False,
            "wald_debug": "no precomp block for cell_type",
        }

    debug_base = blk.get("debug_msg", "ok")

    if cre == "reference":
        return {
            "test_statistic": 0.0,
            "p_value": 1.0,
            "fold_change": 1.0,
            "flattened": False,
            "wald_debug": debug_base,
        }

    col = find_treatment_column(blk["xmu_names"], "cre_id", cre)
    if col is None:
        return {
            "test_statistic": np.nan,
            "p_value": np.nan,
            "fold_change": np.nan,
            "flattened": False,
            "wald_debug": f"no contrast term for {cre}",
        }

    j    = blk["xmu_names"].index(col)
    beta = float(blk["xmu"][j])
    se   = float(blk["se_x_mu"][j])

    # sanity / failure modes
    if not np.isfinite(se):
        return {
            "test_statistic": np.nan,
            "p_value": np.nan,
            "fold_change": np.nan,
            "flattened": False,
            "wald_debug": f"{debug_base}; non-finite se for {cre}",
        }

    if se == 0.0:
        return {
            "test_statistic": np.nan,
            "p_value": np.nan,
            "fold_change": 1.0,
            "flattened": False,
            "wald_debug": f"{debug_base}; se==0 for {cre}",
        }


    z = beta / se
    eps = np.finfo(float).tiny  # ~1e-308
    p = max(chi2.sf(z*z, 1) , eps) # survival function instead of cdf should be more stable with regard to tiny values

    return {
        "test_statistic": z,
        "p_value": p,
        "fold_change": float(np.exp(beta)),
        "flattened": False,
        "wald_debug": debug_base,
    }

def _wald_by_cre_row(row: dict, bundle: dict):
    cre = row["comparison_CRE"]
    ct  = row["comparison_cell_type"]
    blk = bundle["by_cre"].get(cre)

    if blk is None:
        return {
            "test_statistic": np.nan,
            "p_value": np.nan,
            "fold_change": np.nan,
            "flattened": False,
            "wald_debug": "no precomp block for cre",
        }

    debug_base = blk.get("debug_msg", "ok")

    if ct == "reference":
        return {
            "test_statistic": 0.0,
            "p_value": 1.0,
            "fold_change": 1.0,
            "flattened": False,
            "wald_debug": debug_base,
        }

    col = find_treatment_column(blk["xmu_names"], "cell_type", ct)
    if col is None:
        return {
            "test_statistic": np.nan,
            "p_value": np.nan,
            "fold_change": np.nan,
            "flattened": False,
            "wald_debug": f"no contrast term for {ct}",
        }

    j    = blk["xmu_names"].index(col)
    beta = float(blk["xmu"][j])
    se   = float(blk["se_x_mu"][j])

    if not np.isfinite(se):
        return {
            "test_statistic": np.nan,
            "p_value": np.nan,
            "fold_change": np.nan,
            "flattened": False,
            "wald_debug": f"{debug_base}; non-finite se for {ct}",
        }

    if se == 0.0:
        return {
            "test_statistic": np.nan,
            "p_value": np.nan,
            "fold_change": np.nan,
            "flattened": False,
            "wald_debug": f"{debug_base}; se==0 for {ct}",
        }

    z = beta / se
    eps = np.finfo(float).tiny  # ~1e-308
    p = max(chi2.sf(z*z, 1) , eps) # survival function instead of cdf should be more stable with regard to tiny values

    return {
        "test_statistic": z,
        "p_value": p,
        "fold_change": float(np.exp(beta)),
        "flattened": False,
        "wald_debug": debug_base,
    }

def _wald_row_fn(row: dict, bundle: dict):
    ref_ct  = row.get("reference_cell_type")
    ref_cre = row.get("reference_CRE")
    comp_ct  = row["comparison_cell_type"]
    comp_cre = row["comparison_CRE"]

    # by-cell-type if reference cell type missing or equals the comparison
    if pd.isna(ref_ct) or (ref_ct == comp_ct):
        return _wald_by_celltype_row(row, bundle)

    # by-CRE if reference CRE missing or equals the comparison
    if pd.isna(ref_cre) or (ref_cre == comp_cre):
        return _wald_by_cre_row(row, bundle)

    # crossed case not supported
    return {
        "test_statistic": np.nan,
        "p_value": np.nan,
        "fold_change": np.nan,
        "flattened": False,
        "wald_debug": "crossed comparison not supported",
    }

def _wald_make_bundle(hypotheses, models_or_counts, **kw):
    # expects an ortho object with make_wald_eval_bundle()
    return models_or_counts.make_wald_eval_bundle()
# ---- Mann–Whitney U / Wilcoxon rank-sum -------------------------------------

def _mwu_make_bundle(hypotheses, models_or_counts, **kw):
    """
    Build a lookup table of UMI counts keyed by (cell_type, cre_id).

    This lets the row_fn do O(1) dict lookups instead of re-filtering the
    full counts dataframe for every hypothesis row, which drastically
    reduces overhead when running under Dask.
    """
    if hasattr(models_or_counts, "data"):
        df = models_or_counts.data[["cell_type", "cre_id", "umis_mpra_bc"]].copy()
    elif hasattr(models_or_counts, "training_data"):
        df = models_or_counts.training_data.data[["cell_type", "cre_id", "umis_mpra_bc"]].copy()
    else:
        raise TypeError("MWU requires a scMPRA_data object (UMI-wise) or training data attribute to be in ortho object.")

    df["cell_type"] = df["cell_type"].astype(str)
    df["cre_id"] = df["cre_id"].astype(str)

    # Group once: (cell_type, cre_id) → np.array of counts
    grouped = (
        df.groupby(["cell_type", "cre_id"])["umis_mpra_bc"]
          .apply(lambda s: s.to_numpy())
    )
    counts_dict = {key: arr for key, arr in grouped.items()}

    return {"counts": counts_dict}

def _mwu_row_fn(
    row,
    bundle,
    *,
    method="auto",
    alternative="two-sided",
    pseudocount=0.01,
):
    """
    Compute Mann–Whitney U p-value for the appropriate within-cell / within-CRE
    comparison, and a descriptive fold change based on a log1p-mean summary.

    Assumes `bundle["counts"]` is a dict:
        (cell_type, cre_id) -> np.ndarray of umis_mpra_bc
    built once in _mwu_make_bundle.
    """
    counts = bundle["counts"]

    comp_ct  = str(row["comparison_cell_type"])
    comp_cre = str(row["comparison_CRE"])
    ref_ct   = row.get("reference_cell_type")
    ref_cre  = row.get("reference_CRE")

    # Normalize possible NaNs to None for easier checks
    if pd.isna(ref_ct):
        ref_ct = None
    if pd.isna(ref_cre):
        ref_cre = None

    # helper to grab a group safely
    empty = np.empty(0, dtype=float)

    def get_group(ct, cre):
        return counts.get((str(ct), str(cre)), empty)

    # --- Decide comparison axis: mirror the Wald dispatcher’s rules ---

    # by-cell-type => compare CREs within the same cell type
    if (ref_ct is None) or (str(ref_ct) == comp_ct):
        base_cre = "reference" if ref_cre is None else str(ref_cre)
        g1 = get_group(comp_ct, comp_cre)
        g0 = get_group(comp_ct, base_cre)

    # by-CRE => compare cell types within the same CRE
    elif (ref_cre is None) or (str(ref_cre) == comp_cre):
        base_ct = "reference" if ref_ct is None else str(ref_ct)
        g1 = get_group(comp_ct, comp_cre)
        g0 = get_group(base_ct, comp_cre)

    else:
        # crossed case not supported in this simple MWU
        return {
            "test_statistic": np.nan,
            "p_value": np.nan,
            "fold_change": np.nan,
            "flattened": False,
            "ref_mean": np.nan,
            "comp_mean": np.nan,
        }

    # If either side has no observations, MWU is undefined; return NA-ish.
    if (g1.size == 0) or (g0.size == 0):
        return {
            "test_statistic": np.nan,
            "p_value": np.nan,
            "fold_change": np.nan,
            "flattened": False,
            "ref_mean": np.nan,
            "comp_mean": np.nan,
        }

    # Mann–Whitney U (SciPy >=1.7 supports method="auto")
    try:
        stat, p = mannwhitneyu(g1, g0, alternative=alternative, method="auto")
    except TypeError:
        # for older SciPy, fall back without 'method'
        stat, p = mannwhitneyu(g1, g0, alternative=alternative)

    # Descriptive summary: log1p mean (≈ geometric-ish on counts) with pseudocount FC
    s1 = float(np.exp(np.mean(np.log1p(g1))) - 1.0)
    s0 = float(np.exp(np.mean(np.log1p(g0))) - 1.0)

    fc = (s1 + pseudocount) / (s0 + pseudocount)

    return {
        "test_statistic": float(stat),
        "p_value": float(p),
        "fold_change": float(fc),
        "flattened": False,
        "ref_mean": s0,
        "comp_mean": s1,
    }
# ---- Bootstrap activity measurement ------------------------------------------
# ---- BOOTSTRAP ACTIVITY (empirical p vs controls) ----------------------------
# ==============================
# Bootstrap flavor (activity)
# ==============================
# --- bootstrap support structs ---

@dataclass
class _BootRepGroupCT:
    # one biological replicate worth of data
    cell_type:  np.ndarray              # per-row cell_type (string)
    norm_umis:  np.ndarray              # per-row normalized_umis_mpra_bc (float)
    idx_by_cre_ct: dict[tuple[str,str], np.ndarray]   # (cre, ct) -> row indices
    idx_ctrl_by_ct: dict[str, np.ndarray]             # ct -> union of control rows
    n_int_by_cre_ct: dict[tuple[str,str], int]        # observed #integrations
    median_int_nonctrl: int | None

@dataclass
class _BootBundleCT:
    by_rep: dict[str, _BootRepGroupCT]
    control_cres: tuple[str, ...]
    n_int_strategy: str                 # "as_observed" | "median_non_reference"
# ---- worker-side helpers ----
def _bootstrap_build_bundle(
    hypotheses: "HypothesisSet",
    models_or_counts,
    client=None,
    *,
    n_bootstraps: int = 10_000,
    n_int_strategy: str = "match_cre",           # {"match_cre", "median_controls"}
    pseudocount: float = 1e-8,
    rng_seed: int | None = None,
    rep_to_biol: "dict[str,str] | None" = None,  # optional mapping if you want to elevate rep_id -> biol_rep
    **kw,
):
    """
    Build a compact bundle with per-integration values we can ship once to workers.

    Safeguards:
      - Infers control CREs from hypotheses.reference_CRE (non-null uniques).
      - Falls back to 'rep_id' when 'biol_rep' is absent (or allows a mapping).
      - Chooses metric column automatically.
      - Drops any controls absent from the counts (warns).
      - Filters to only CREs referenced by the hypothesis set (comparison + control).
    """
    # ---- Validate counts object ----
    counts = models_or_counts
    if not hasattr(counts, "data"):
        raise TypeError("bootstrap_activity expects a scMPRA_data-like object with a `.data` DataFrame.")

    df = counts.data
    needed = {"cell_type", "cre_id", "rep_id", "cell_bc", "transfection_bc"}
    missing = sorted(needed - set(df.columns))
    if missing:
        raise ValueError(f"Counts table is missing required columns: {missing}")

    # ---- Choose metric column ----
    if "normalized_umis_mpra_bc" in df.columns:
        metric_col = "normalized_umis_mpra_bc"
    elif "umis_mpra_bc" in df.columns:
        metric_col = "umis_mpra_bc"
        warnings.warn(
            "[bootstrap_activity] Using raw 'umis_mpra_bc' because "
            "'normalized_umis_mpra_bc' was not found."
        )
    else:
        raise ValueError("Neither 'normalized_umis_mpra_bc' nor 'umis_mpra_bc' present in counts table.")

    # ---- Derive/ensure biol_rep ----
    if "biol_rep" in df.columns:
        biol = df["biol_rep"].astype(str)
    elif rep_to_biol is not None:
        # Mapping provided
        biol = df["rep_id"].astype(str).map(rep_to_biol).fillna(df["rep_id"].astype(str))
    else:
        # Fall back: treat rep_id as biological replicate
        biol = df["rep_id"].astype(str)
    df = df.assign(biol_rep=biol)

    # ---- Controls from hypotheses ----
    hdf = hypotheses.to_dataframe()
    controls_from_hs = sorted(set(hdf["reference_CRE"].dropna().astype(str).unique()))
    if not controls_from_hs:
        raise ValueError("No control CREs found in hypothesis set (column 'reference_CRE').")

    # ---- Comparison CREs present in HS ----
    compare_cres = sorted(set(hdf["comparison_CRE"].dropna().astype(str).unique()))
    all_wanted_cres = set(controls_from_hs).union(compare_cres)

    # ---- Filter counts to only needed CREs (comp + controls) ----
    df = df[df["cre_id"].astype(str).isin(all_wanted_cres)].copy()
    if df.empty:
        raise ValueError("After filtering to hypothesis CREs + controls, no rows remain in counts.")

    # ---- Validate controls presence ----
    present_controls = sorted(set(df.loc[df["cre_id"].astype(str).isin(controls_from_hs), "cre_id"].astype(str).unique()))
    missing_controls = sorted(set(controls_from_hs) - set(present_controls))
    if missing_controls:
        warnings.warn(
            f"[bootstrap_activity] The following control CREs are not present in counts and will be ignored: {missing_controls}"
        )
    if not present_controls:
        raise ValueError("[bootstrap_activity] No control CREs from the hypothesis set are present in counts.")

    # ---- Build compact per-integration table ----
    # Integration unit = unique (cell_bc, transfection_bc) within (biol_rep, cell_type, cre_id)
    # Value for each integration = mean(metric_col) across those rows (normally 1:1, but defensively aggregate)
    keys = ["biol_rep", "cell_type", "cre_id", "cell_bc", "transfection_bc"]
    g = (
        df.assign(_val=df[metric_col].astype(float))
          .groupby(keys, as_index=False)["_val"].mean()
          .rename(columns={"_val": "value"})
    )
    # We only keep the minimal table needed for worker-side bootstraps
    integrations = g  # columns: biol_rep, cell_type, cre_id, cell_bc, transfection_bc, value

    # ---- Bundle ----
    bundle = {
        "integrations": integrations,
        "metric_col": metric_col,
        "controls": present_controls,         # control labels (strings)
        "n_bootstraps": int(n_bootstraps),
        "n_int_strategy": str(n_int_strategy),
        "pseudocount": float(pseudocount),
        "rng_seed": None if rng_seed is None else int(rng_seed),
    }
    return bundle

def _bootstrap_row_fn(row: dict, bundle: dict, **kw):
    """
    Bootstrap activity test.
    If max_over_celltypes=True (default), reproduces the old pipeline behavior:
      - for each biological replicate, resample integrations for the CRE-of-interest
        and for controls across ALL cell types,
      - within each bootstrap draw, compute per-CT means and take the MAX across CTs,
      - compare the CRE bootstrap distribution of MAX(CT mean) to the controls'.
    Returns:
      test_statistic = median of CRE bootstrap MAX(CT mean)
      p_value        = empirical two-sided p (|null| >= |obs|)
      fold_change    = median(CRE bootstrap MAX) / median(Control bootstrap MAX) (+ε)
    """


    # options / defaults
    B        = int(bundle["n_bootstraps"])
    strategy = bundle["n_int_strategy"]          # "match_cre" or "median_controls" (kept for compat)
    pc       = float(bundle["pseudocount"])
    seed     = bundle.get("rng_seed", None)
    max_over_ct = kw.get("max_over_celltypes", True)

    comp_cre = str(row["comparison_CRE"])

    # Early neutral return for control-vs-control
    if comp_cre.lower() == "reference":
        return {"test_statistic": 0.0, "p_value": 1.0, "fold_change": 1.0, "flattened": False}

    integ    = bundle["integrations"]             # columns: biol_rep, cell_type, cre_id, cell_bc, transfection_bc, value
    controls = set(bundle["controls"])

    # We ignore comparison_cell_type when max_over_celltypes=True (match old behavior)
    # If user explicitly disabled max_over_celltypes, we fall back to  previous per-CT logic.
    if not max_over_ct:
        # fall back to  prior per-CT implementation 
        return _bootstrap_row_fn_per_ct(row, bundle, **kw)  # can keep previous function under this name

    # Prepare RNG stable per CRE
    rng = np.random.default_rng(None if seed is None else (hash(("BOOTMAX", seed, comp_cre)) % (2**32 - 1)))

    # Biological replicates present for this CRE or controls
    biols = sorted(integ["biol_rep"].astype(str).unique().tolist())
    if not biols:
        warnings.warn("[bootstrap_activity] No biological replicates present.")
        return {"test_statistic": np.nan, "p_value": np.nan, "fold_change": np.nan, "flattened": False}

    # Collect bootstrap MAX(CT mean) across reps (we’ll average across reps per draw)
    A_boot_max_per_rep = []   # list of arrays, each shape (B,)
    C_boot_max_per_rep = []   # list of arrays, each shape (B,)
    NULL_boot_diffs_per_rep = []

    have_any = False

    for bi in biols:
        sub = integ.loc[integ["biol_rep"] == bi]
        A_df = sub.loc[sub["cre_id"] == comp_cre, ["cell_type", "value"]]
        C_df = sub.loc[sub["cre_id"].isin(controls), ["cell_type", "value"]]

        if A_df.empty or C_df.empty:
            continue

        have_any = True

        # encode CTs as codes for fast grouping
        ct_levels = pd.Categorical(pd.concat([A_df["cell_type"], C_df["cell_type"]]).astype(str))
        # We must re-map separately to keep consistent codes for both A and C
        # Create a code map on the union levels
        ct_union = ct_levels.categories
        ct_map = {ct: i for i, ct in enumerate(ct_union)}

        A_ct = A_df["cell_type"].astype(str).map(ct_map).to_numpy(dtype=np.int64)
        A_v  = A_df["value"].to_numpy(dtype=float)
        C_ct = C_df["cell_type"].astype(str).map(ct_map).to_numpy(dtype=np.int64)
        C_v  = C_df["value"].to_numpy(dtype=float)

        nA = A_v.size
        nC_all = C_v.size
        if nA == 0 or nC_all == 0:
            continue

        # choose nC per strategy
        if strategy in ("match_cre", "observed"):
            nC = nA
        elif strategy in ("median_controls", "median_non_reference"):
            nC = int(max(1, np.median([nC_all])))
        else:
            warnings.warn(f"[bootstrap_activity] Unknown n_int_strategy '{strategy}', using 'match_cre'.")
            nC = nA

        # helper to compute MAX over CT means for a bootstrap sample
        def max_ct_mean(values, ct_codes, idx):
            # values[idx] are the resampled per-integration values
            vv  = values[idx]
            cc  = ct_codes[idx]
            # accumulate sum and count per CT
            k   = len(ct_union)
            sums   = np.bincount(cc, weights=vv, minlength=k)
            counts = np.bincount(cc, minlength=k)
            # avoid div0: set means only where count>0
            means = np.zeros(k, dtype=float)
            nz = counts > 0
            means[nz] = sums[nz] / counts[nz]
            return means.max() if nz.any() else 0.0

        # Pre-sample indices: (B, nA) and (B, nC)
        A_idx = rng.integers(0, nA, size=(B, nA), endpoint=False)
        C_idx = rng.integers(0, nC_all, size=(B, nC), endpoint=False)

        # Compute bootstrap MAX(CT mean) for CRE and Controls
        A_max = np.empty(B, dtype=float)
        C_max = np.empty(B, dtype=float)
        for b in range(B):
            A_max[b] = max_ct_mean(A_v, A_ct, A_idx[b])
            C_max[b] = max_ct_mean(C_v, C_ct, C_idx[b])

        A_boot_max_per_rep.append(A_max)
        C_boot_max_per_rep.append(C_max)

                # NEW ⬇︎ Null bootstrap: pool A and C, resample two groups of sizes (nA, nC)
        P_v  = np.concatenate([A_v, C_v])
        P_ct = np.concatenate([A_ct, C_ct])
        P_n  = P_v.size
        P_idx_A = rng.integers(0, P_n, size=(B, nA), endpoint=False)
        P_idx_C = rng.integers(0, P_n, size=(B, nC), endpoint=False)

        null_diffs = np.empty(B, dtype=float)
        for b in range(B):
            maxA = max_ct_mean(P_v, P_ct, P_idx_A[b])
            maxC = max_ct_mean(P_v, P_ct, P_idx_C[b])
            null_diffs[b] = maxA - maxC
        NULL_boot_diffs_per_rep.append(null_diffs)

    if not have_any:
        warnings.warn(f"[bootstrap_activity] No usable integrations for CRE='{comp_cre}' across biological replicates.")
        return {"test_statistic": np.nan, "p_value": np.nan, "fold_change": np.nan, "flattened": False}

    # Combine across reps: average the MAX(CT mean) across reps per bootstrap draw
    A_mat = np.vstack(A_boot_max_per_rep)      # (R, B)
    C_mat = np.vstack(C_boot_max_per_rep)      # (R, B)
    A_bar = A_mat.mean(axis=0)                 # (B,)
    C_bar = C_mat.mean(axis=0)                 # (B,)

    # Observed difference = medians of CRE and Control boot max (matches summary centering)
    obs_diff = float(np.median(A_bar) - np.median(C_bar))
    
    # NEW ⬇︎ Build combined NULL distribution (average per-rep null diffs across reps)
    NULL_mat = np.vstack(NULL_boot_diffs_per_rep)  # (R, B)
    NULL_bar = NULL_mat.mean(axis=0)               # (B,)

    # NEW ⬇︎ Two-sided empirical p-value against null
    p = float((np.sum(np.abs(NULL_bar) >= np.abs(obs_diff)) + 1) / (NULL_bar.size + 1))

    # Reportables (matching old summaries)
    test_statistic = float(np.median(A_bar))                       # q50 of CRE bootstrap max across CTs
    fold_change    = float((np.median(A_bar) + pc) / (np.median(C_bar) + pc))

    return {
        "test_statistic": test_statistic,
        "p_value": p,
        "fold_change": fold_change,
        "flattened": False,
    }

def _bootstrap_row_fn_per_ct(row: dict, bundle: dict, **kw):
    """
    Evaluate one hypothesis row via bootstrap resampling.

    Interpretation:
      - We assume the hypothesis fixes a (comparison_cell_type, comparison_CRE),
        and its `reference_CRE` in the hypotheses are the controls (already folded
        into the bundle).
      - For each biological replicate present, we:
          * determine the # of integrations for the CRE-of-interest (N_A)
          * determine the control sample size (N_C), either N_A ("match_cre") or
            median(#integrations) among control integrations for that replicate
          * draw with replacement N_A from CRE-of-interest integrations and N_C
            from pooled controls (both within the same cell_type & biol_rep)
          * compute mean difference (CRE - Control)
        We average the mean differences across replicates to obtain a combined
        statistic. Repeating this B times yields a bootstrap distribution; the
        empirical two-sided p-value is based on |diff| >= |obs_diff|.
      - Fold change is computed on the replicate-averaged means
        with a small pseudocount.

    Safeguards return NaNs with clear conditions:
      - missing cell_type or CRE rows
      - CRE is literally the control label (returns neutral stats)
      - no control integrations present for that cell_type
      - zero integration counts in either group
    """
    try:
        comp_ct  = str(row["comparison_cell_type"])
        comp_cre = str(row["comparison_CRE"])

        # Early neutral return if user asked for the control vs control
        if comp_cre.lower() == "reference":
            return {"test_statistic": 0.0, "p_value": 1.0, "fold_change": 1.0, "flattened": False}

        integ = bundle["integrations"]
        controls = set(bundle["controls"])
        B = int(bundle["n_bootstraps"])
        strategy = bundle["n_int_strategy"]
        pc = float(bundle["pseudocount"])
        seed = bundle.get("rng_seed", None)

        # Subset to this cell type once
        sub_ct = integ.loc[integ["cell_type"] == comp_ct]
        if sub_ct.empty:
            warnings.warn(f"[bootstrap_activity] No integrations for cell_type='{comp_ct}'.")
            return {"test_statistic": np.nan, "p_value": np.nan, "fold_change": np.nan, "flattened": False}

        # Identify biological replicates present for this CT and CRE/control
        biols = sorted(sub_ct["biol_rep"].astype(str).unique().tolist())
        if not biols:
            warnings.warn(f"[bootstrap_activity] No biological replicates found after subsetting to cell_type='{comp_ct}'.")
            return {"test_statistic": np.nan, "p_value": np.nan, "fold_change": np.nan, "flattened": False}

        # Make RNG (per-row) stable but different if seed provided
        rng = np.random.default_rng(None if seed is None else (hash((seed, comp_ct, comp_cre)) % (2**32 - 1)))

        # For each biol_rep, prepare arrays for bootstrapping
        per_rep_obs_diff = []
        per_rep_boot_diffs = []
        per_rep_A_means    = []   # NEW: store bootstrap CRE means per rep (shape: (B,))
        per_rep_C_means    = []   # NEW: store bootstrap Control means per rep (shape: (B,))

        have_any = False

        for bi in biols:
            sub = sub_ct.loc[sub_ct["biol_rep"] == bi]
            A = sub.loc[sub["cre_id"] == comp_cre, "value"].to_numpy(dtype=float)
            C = sub.loc[sub["cre_id"].isin(controls), "value"].to_numpy(dtype=float)
           
            if A.size == 0 or C.size == 0:
                            continue
            have_any = True

            # choose N for CRE and Controls
            nA = A.size
            if strategy in ("match_cre", "observed"):
                nC = nA
            elif strategy in ("median_controls", "median_non_reference"):
                nC = int(max(1, np.median([C.size])))
            else:
                warnings.warn(f"[bootstrap_activity] Unknown n_int_strategy '{strategy}', using 'match_cre'.")
                nC = nA
            
            # Observed per-rep mean difference
            obs_diff_rep = (A.mean() - C.mean())
            per_rep_obs_diff.append(obs_diff_rep)

            # Bootstrap: sample means for each group and difference
            # Pre-sample indices for speed
            # A_boot: shape (B, nA), C_boot: shape (B, nC)
            A_idx = rng.integers(0, nA, size=(B, nA), endpoint=False)
            C_idx = rng.integers(0, C.size, size=(B, nC), endpoint=False)
            A_means = A[A_idx].mean(axis=1)
            C_means = C[C_idx].mean(axis=1)
            diffs = A_means - C_means   # shape (B,)
            per_rep_boot_diffs.append(diffs)
            per_rep_A_means.append(A_means)  # NEW
            per_rep_C_means.append(C_means) 
        if not have_any:
            # nothing to compare for this row
            warnings.warn(
                f"[bootstrap_activity] No usable integrations for CRE='{comp_cre}' at cell_type='{comp_ct}' "
                "(missing CRE or control integrations across biological replicates)."
            )
            return {"test_statistic": np.nan, "p_value": np.nan, "fold_change": np.nan, "flattened": False}

        # Combine across biological reps
        # Observed = mean of per-rep observed differences
        obs_diff = float(np.mean(per_rep_obs_diff))

        # Combined bootstrap distribution = average per bootstrap across reps
        # (stack to array (R, B) then take mean over axis=0 -> (B,))
        boot_mat = np.vstack(per_rep_boot_diffs)   # (R, B)
        boot_combined = boot_mat.mean(axis=0)      # (B,)

        # Empirical two-sided p-value
        p = float((np.sum(np.abs(boot_combined) >= np.abs(obs_diff)) + 1) / (boot_combined.size + 1))

        # NEW: combine bootstrap CRE/Control means across reps (average per draw),
        # then take medians for the summary stat and FC.
        A_boot_mat = np.vstack(per_rep_A_means)         # (R, B)
        C_boot_mat = np.vstack(per_rep_C_means)         # (R, B)
        A_bar_boot = A_boot_mat.mean(axis=0)            # (B,)
        C_bar_boot = C_boot_mat.mean(axis=0)            # (B,)

        test_statistic = float(np.median(A_bar_boot)) if A_bar_boot.size else np.nan
        if A_bar_boot.size and C_bar_boot.size:
            fold_change = float((np.median(A_bar_boot) + pc) / (np.median(C_bar_boot) + pc))
        else:
            fold_change = np.nan

        return {
            "test_statistic": test_statistic,
            "p_value": p,
            "fold_change": fold_change,
            "flattened": False,
        }

    except Exception as e:
        warnings.warn(f"[bootstrap_activity] row failed with error: {e}")
        return {"test_statistic": np.nan, "p_value": np.nan, "fold_change": np.nan, "flattened": False}


# ---- The tiny switchboard ----------------------------------------------------
#switchboard to hole types of hypothesis tests that have been implemented so far
TESTS = {
    "wald": {"make_bundle": _wald_make_bundle, "row_fn": _wald_row_fn, "defaults": {}},
    "mwu":  {"make_bundle": _mwu_make_bundle,  "row_fn": _mwu_row_fn,  "defaults": {"method": "auto", "alternative": "two-sided"}},
    "bootstrap_activity": {"make_bundle": _bootstrap_build_bundle, "row_fn": _bootstrap_row_fn, "defaults": {
            "n_int_strategy": "median_non_reference",   # or "observed"
            "n_bootstraps": 10_000,
            "rng_seed": 42,
            "max_over_celltypes": True,
        }
        },
}

class HypothesisTester:
    """
    Orchestrates running a test function on each hypothesis row.
    You supply `test_fn` that implements a single-row comparison and returns:
        dict(test_type, test_statistic, p_value, fold_change, flattened)
    The runner adds BH (`bh_p`) and merges back with the hypothesis columns to return a ResultSet.
    """

    def __init__(self, test_type: str, **overrides):
        if test_type not in TESTS:
            raise ValueError(f"Unknown test_type '{test_type}'. Choose one of {sorted(TESTS)}.")
        self.test_type = test_type
        self._make_bundle = TESTS[test_type]["make_bundle"]
        self._row_fn      = TESTS[test_type]["row_fn"]
        self.kw = dict(TESTS[test_type].get("defaults", {}))
        self.kw.update(overrides or {})

    def run(self, hypotheses: HypothesisSet, models_or_counts: "scMPRA_data | ortho", client=None) -> ResultSet:
        # canonicalize if we can see the counts
        df_h = hypotheses.to_dataframe().copy()
        scmpra = getattr(models_or_counts, "training_data", None) or getattr(models_or_counts, "data", None)
        if scmpra is not None and hasattr(scmpra, "data"):
            df_h = canonicalize_hypotheses(HypothesisSet.from_dataframe(df_h),
                                           scmpra if hasattr(scmpra, "data") else models_or_counts).to_dataframe()

        # build the bundle once
        bundle = self._make_bundle(hypotheses=HypothesisSet.from_dataframe(df_h),
                                   models_or_counts=models_or_counts, **self.kw)

        recs = df_h.to_dict(orient="records")
        if client is not None:
            bundle_f = client.scatter(bundle, broadcast=True)
            # futures = client.map(self._row_fn, recs, repeat(bundle_f), pure=True)
            futures = [client.submit(self._row_fn, r, bundle_f, **self.kw) for r in recs]
            results = client.gather(futures)
        else:
            results = [self._row_fn(r, bundle, **self.kw) for r in recs]

        out = pd.concat([df_h.reset_index(drop=True), pd.DataFrame(results)], axis=1)
        out["test_type"] = self.test_type
        out["bh_p"] = _bh_adjust(out["p_value"])
        if "flattened" in out:
            out["flattened"] = out["flattened"].astype(bool)
        return ResultSet.from_dataframe(out)

def load_hypothesis_set(filepath):
    """
    Loads a hypothesis (or result) table from disk and returns a HypothesisSet or ResultSet.
    Supports .tsv/.csv/.parquet by extension.
    """
    path = Path(filepath)
    ext = path.suffix.lower()
    if ext == ".tsv":
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=True)
    elif ext == ".csv":
        df = pd.read_csv(path, dtype=str, keep_default_na=True)
    # elif ext in {".parquet", ".pq"}:
    #     df = pq.read_table(path).to_pandas(types_mapper=pd.ArrowDtype)
    else:
        raise ValueError(f"Unsupported file extension: {ext}")

    ttype = table_type(df.columns)
    if ttype == "hypotheses":
        return HypothesisSet.from_dataframe(df)
    elif ttype == "results":
        return ResultSet.from_dataframe(df)
    else:
        raise ValueError(f"Loaded table does not look like a hypothesis/results table (got type '{ttype}').")

def hypothesis_tester(scmpra_models_or_data, hypotheses: HypothesisSet, flavor="wald", test_fn=None):
    """
    Backward-compatible facade.

    Provide either:
      - test_fn: a per-row callable used by HypothesisTester (preferred while migrating)
      - or later, we can route based on `flavor` to built-in tests.

    Returns a ResultSet.
    """
    if test_fn is None:
        raise NotImplementedError("Please provide `test_fn` for now. Built-in tests will be added to route by `flavor`.")
    runner = HypothesisTester(test_fn=test_fn, test_type_name=flavor)
    return runner.run(hypotheses)


class de_novo_simulation:
    """
    Class for simulating new datasets.

    A single instance of this object should be used to 
    represent n simulated replicates of one experimental setup.
    For different parameters, initalize multiple objects.

    Takes MPRA libraries in standard form. To generate a synthetic libraries, create
    a ground_truth and use simulate_library.

    Initalized over a directory, where it reads and writes.
    Directory is specififed by 'location', 'name' pair (pointing to directory
    location/name). 'location', 'name' and derived paths are stored 
    in self object but never saved to disc, so you can rename or move
    the simulation as you please.

    The general workflow is:
    - __init__; sets the initial simulation parameters
        OR loads a simulation previously initalized.
    - gamut; automatically runs:
        - _simulate_transfection
        - _simulate_transcription
        - _realize_simulations
    Then optionally
    - fit_orthos
    - precompute_wald (OPG or sandwich)
    Then, for each hypothesis test / method
    - test
    Summarization can be performed & saved with
    - PRC_aggregated
        - Also saves AUPRC and other metrics to state
    - volcano
    - p_value_calibration
    - gt_vs_p
    Once you are done you can 
    - save
    NOTE: You really want to call `save` when you are done!

    Note that saving is blocking, but generally the object is written
    to dump data to disc as it goes, avoiding keeping too much data in memory. So 
    if you are running many simulations, you can queue them all up, then
    just call save on all of them in a loop or similar.
    
    Plots are lightweight to generate and so are not saved.

    Efficiency could be improved by avoiding saving information multiple times,
    but we sacrifice a little disc space to enhance reproducability & debugging
    (e.g. saving twice allows us to keep both the filtered and unfiltered 
    versions of the datasets, e.g. saving the ground-truth locally instead of a reference
    to another location keeps the object self-contained...) I try to keep things 
    gzip compressed to offset this.
    """
    
    def __init__(self,location,name,client,
        libraries:list[pd.DataFrame]=None,
        library_mapping:str | list[int]=None,
        n_sims:int=None,
        experiment_bounds:Bounds=None,
        ground_truth:pd.DataFrame=None,
        flatten_overtransfection=None,
        negative_controls=["reference"],
        reference_cell_type="reference"
        ):
        """
        'location', 'name' and 'client' are always mandatory. 
        (client *can* be set to None, but only for summary functions: graphing PRC, computing confusion matricies, etc...)
        The rest of the parameters are optional iff loading a previously initalized simulation batch.
        Otherwise they are also mandatory.

        'libraries' takes a list of 1+ MPRA libraries in standard format.

        'experiment_bounds' takes a 'Bounds' object describing the experiment. You may wish to use one of the
        presets extracted from emperical data...
        
        'library_mapping' takes one of
        - the string 'one_library': all replicates are simulated from one library
            - e.g. simulating from a real library
        - the string 'corresponding': all replicates are simulated from a corresponding library
            - e.g. prospectively evaluating theoretical performace under some experimental setup
        - A list of ints: this is just passed through. The list should be of length of the
          number of simulations, and each element is an int pointing to a library.

        'n_sims'
        - number of sims.
        """
        self.location = Path(location)
        self.name = Path(name)
        self.client = client
        self.testqueue = []
        
        fullp=(self.location/name)
        fullp.mkdir(parents=True, exist_ok=True)
        statep=fullp/"state.parquet"
        
        #save in object so other functions can access
        self.fullp=fullp
        self.statep=statep
        self.descripd=fullp/"descriptions"
        libp=fullp/"libraries"
        self.libp=libp
        self.scmpradatp=fullp/"simulated_scmpra"
        self.orthod=fullp/"orthos"
        self.sandwichd=fullp/"sandwich"
        self.opgd=fullp/"opg"
        self.testd=fullp/"tests"

        if statep.exists():
            logger.info(f"'state.parquet' found for '{name}', loading.")

            #load state
            self.state=pd.read_parquet(statep)

            #load ground_truth
            self.ground_truth= pd.read_csv(fullp/"ground_truth.tsv.gz",sep="\t",compression="gzip",index_col=0)

            #load futures tracker
            rawfut=pd.read_csv(fullp/"futures.tsv.gz",sep="\t",compression="gzip", index_col=0)
            flat=rawfut.to_numpy().ravel()

            futures = [client.submit(lambda x: x, v, pure=False) for v in flat]
            arr = pd.array(futures, dtype="object").to_numpy().reshape(rawfut.shape)
            
            self.futures=pd.DataFrame(arr,
                                      index=rawfut.index,
                                      columns=rawfut.columns)
        else:
            logger.info(f"No 'state.parquet' found for '{name}'. Initalizing new object.")
            
            self.state = pd.DataFrame()
            
            self.futures = pd.DataFrame(index=range(n_sims))
            
            required = ("libraries", 
                        "library_mapping",
                        "n_sims",
                        "experiment_bounds",
                        "ground_truth",
                        "flatten_overtransfection")
                        
            
            #lightweight, just save.
            self.set_state_field("n_sims",n_sims)
            self.set_state_field("reference_cell_type",reference_cell_type)
            self.set_state_field("negative_controls",negative_controls)
            self.set_state_field("alpha",DEFAULT_SIGNIFICANCE_THRESHOLD)
            self.set_state_field("flatten_overtransfection",flatten_overtransfection)
            
            if any(x is None for x in [libraries,library_mapping,n_sims,experiment_bounds,ground_truth]):
                raise ValueError(f"When initalizing a new object, required params include all of: {required}.")

            #if `libraries` is a single dataframe, put it in a one-element list. 
            if isinstance(libraries, pd.DataFrame):
                libraries=[libraries]

            #sanity checks on, and setting of library mapping 

            if library_mapping=="one_library":
                if not len(libraries)==1: raise ValueError(f"Specified one_library, but passed {len(libraries)} libraries.")
                self.set_state_field("library_mapping",[0 for i in range(0,n_sims)])
            
            elif library_mapping == "corresponding":
                if not len(libraries)==n_sims: raise ValueError(f"Corresponding library mapping requires n_sims to be the same as the number of libraries.")
                self.set_state_field("library_mapping",[i for i in range(0,n_sims)])

            elif isinstance(library_mapping, list):
                if all(isinstance(x, int) for x in library_mapping):
                    #check if length is correct
                    if len(library_mapping)!=n_sims: raise ValueError(f"n_sims={n_sims}, =/= len(library_mapping)={len(library_mapping)}")
                    #check if any ints fall outside the appropriate range 
                    for i in library_mapping:
                        if i<0 or i>len(libraries)+1: raise ValueError("At least one of the ints in your library mapping refers to a library you do not have")
                    #all the checks passed!
                    self.set_state_field("library_mapping",library_mapping)
                else:
                    raise ValueError(f"Library mapping {library_mapping}, does not contain ints.")
            else:
                raise ValueError(f"Unrecognized library mapping {library_mapping}.")
            
            #save the passed libraries
            libp.mkdir(parents=True, exist_ok=True)

            for idx,lib in enumerate(libraries):
                lib.to_csv(libp/f"{idx}.tsv.gz",sep="\t",compression='gzip')
            
            #saved the passed bounds object
            experiment_bounds.to_tgz(fullp/"experiment_bounds.tgz")

            #save the ground truth
            self.ground_truth=cast_string_keys(ground_truth, ["cell_type", "cre_id"])
            ground_truth.to_csv(fullp/"ground_truth.tsv.gz",sep="\t",compression="gzip")

            #save the state
            self.save()
    
    def _simulate_transfection(self):
        """
        Simulates transfection. Most of the logic
        offloaded to non-method function _simulate_transfection.
        """
        n_sims=self.get_state_field("n_sims")

        #ground truth
        ground_truth=self.ground_truth

        #load the experiment bounds
        experiment_bounds=Bounds.from_tgz(self.fullp/"experiment_bounds.tgz")


        #create an output directory for the descriptions
        self.descripd.mkdir(parents=True, exist_ok=True)
        
        def _simulate_transfection_helper(experiment_bounds,
                    ground_truth,
                    library,
                    pth):
            """
            Little wrappper for _simulate_transfection
            which handles output to disc and returns a `True` boolean when writing is done
            """
            transfected=_simulate_transfection(
                    experiment_bounds=experiment_bounds,
                    ground_truth=ground_truth,
                    library=lib)
            
            transfected.to_csv(pth,
                    sep="\t",
                    compression="gzip")
            
            return True

        transfection_tracker=[]
        
        for idx in range(0,n_sims):
            #load the corresponding library
            lib=pd.read_csv(self.libp/f"{idx}.tsv.gz",sep="\t",compression='gzip',index_col=0)

            #submit a job to simulate transfection using the helper function
            ret=self.client.submit(_simulate_transfection_helper,
                experiment_bounds=experiment_bounds,
                ground_truth=ground_truth,
                library=lib,
                pth=self.descripd/f"{idx}.tsv.gz")
            
            transfection_tracker.append(ret)
        
        self.futures["transfection"]=transfection_tracker

    def set_state_field(self,field,value):
        self.state[field]=[value]
        self.save()
    
    def get_state_field(self,field):
        return self.state[field][0]
    
    def gamut(self):
        self._simulate_transfection()
        self._simulate_transcription()
    
    def _merge_in_ground_truth(self,hypothesis_set_name,test_type,index):
        """
        Returns a dataframe which a merge of the data of a results object
        (of a test and particular index)
        and the ground-truth. 
        Used as part of test evaluations.

        Note that function is a collector and WILL hang if ANY tests are not completed.
        """
        self._block_until_all_tests_are_done()

        #init & load relevant hypothesis set...
        hypod=self.testd/hypothesis_set_name

        testd=hypod/test_type

        if not testd.is_dir():
            raise FileNotFoundError(f"Could not find test {test_type} for hypothesis set \'{name}\'.")

        index=int(index)
        if index<0 or index>self.get_state_field("n_sims")-1:
            raise RuntimeError(f"Invalid index {index}, out of range [{0},{self.get_state_field('n_sims')-1}]")

        results=ResultSet.from_tsv(testd/f"{index}_results.tsv").df

        #merge in ground truth for `comparison` 
        merged=results.merge(self.ground_truth,
            left_on=["comparison_CRE","comparison_cell_type"],
            right_on=["cre_id","cell_type"]
        )
        
        merged=merged.drop(columns=["cre_id","cell_type"])
        #rename 2x given that  
        merged=merged.rename({"true_mean":"comparison_truth",
                            "mu":"comparison_truth"},axis=1)

        #merge in ground truth for `reference` 
        merged=merged.merge(self.ground_truth,
            left_on=["reference_CRE","reference_cell_type"],
            right_on=["cre_id","cell_type"]
        )
        merged=merged.drop(columns=["cre_id","cell_type"])
        merged=merged.rename({"true_mean":"reference_truth",
                            "mu":"reference_truth"},axis=1)

        #ground truth effect size
        merged["gt_effect_size"]=merged["comparison_truth"]/merged["reference_truth"]
        #ground truth null hypothesis that the CREs are the same : true or false?
        merged["gt_null"]=abs(merged["gt_effect_size"]-1)<FLOATING_POINT_DIFF

        merged["reject_null"]=merged["bh_p"]<self.get_state_field("alpha")

        return merged

    def _classifier_summary(self,hypothesis_set_name,test_type):
        """
        Treating the given test as a classifier, function computes AUROC, AUPRC for all replicates. 
        """
        self._block_until_all_tests_are_done()

        ret={'replicate':[],'auroc':[],'auprc':[]}

        for idx in range(0,self.get_state_field("n_sims")):
            df=self._merge_in_ground_truth(hypothesis_set_name,
                    test_type,
                    index=idx)

            df["meta"] = df["meta"].fillna(0)
            
            #df drop nans
            nona=df.copy()
            nona=nona.dropna()

            if len(nona)!=len(df):
                percent=(len(df)-len(nona))/len(df)*100
                logger.info(f"Dropped {len(df)-len(nona)} or {percent:.1f}% of tests with NA values for test:{test_type}, rep:{idx}.")
                if percent>ERROR_TEST_NAN_PERCENT:
                    raise RuntimeError("Percent of dropped tests greater than threshold, aborting.")

            
            df=nona
            
            y_true = (~df["gt_null"]).astype(int)
            y_score = 1.0 - df["p_value"]

            auroc = roc_auc_score(y_true, y_score)
            auprc = average_precision_score(y_true, y_score)

            ret['replicate'].append(idx)
            ret['auroc'].append(auroc)
            ret['auprc'].append(auprc)
        
        return pd.DataFrame(ret)
           
    def _all_classifier_summary(self,hypothesis_set_name):
        """
        Function computes AUROC, AUPRC, for all tests.
        """
        tests=self.list_tests()
        if not hypothesis_set_name in tests.keys():
            raise RuntimeError(f"Unable to find hypothesis set '{hypothesis_set_name}.'")
        
        tests=tests[hypothesis_set_name]

        ret=[]
        
        for test in tests:
            df=self._classifier_summary(hypothesis_set_name,test_type=test)
            df["test"]=test
            ret.append(df)

        return pd.concat(ret)
    
    def median_performance_curve(self,hypothesis_set_name,performance_type,test_types=None,include_alpha=True):
        """
        Plots the ROC curve for the replicate with the median auROC or the PRC curve for the replicate with the median auPRC.
        - hypothesis_set_name: string of hypothesis set name.
        - test_types: list of strings of tests you want to plot
        - performance_type: ROC or PRC
        
        Note that function is a collector and WILL hang if ANY tests are not completed.
        """
        if performance_type not in ["ROC","PRC"]:
            raise ValueError("Invalid choice for performance_type. Valid choices are 'ROC' and 'PRC'")

        if performance_type=="ROC":
            relevant_au="auroc"
        else:
            relevant_au="auprc"

        #For each test, find the median replicate
        def replicate_with_median(group):
            med = group[relevant_au].median()
            idx = (group[relevant_au] - med).abs().idxmin()
            return group.loc[idx, "replicate"]

        all_summary=self._all_classifier_summary(hypothesis_set_name=hypothesis_set_name)
        if not test_types is None:
            all_summary = all_summary[all_summary["test"].isin(test_types)]

        median_reps = (
            all_summary.groupby("test", group_keys=False)
            .apply(replicate_with_median)
            .to_dict()
        )
        
        #makes a dict of test:median_rep_index

        def plot_curves(curves, kind="roc", ax=None, title=None, p_value=None, point_kwargs=None):
            """
            Plot multiple ROC or PR curves on the same axes, optionally marking a point
            corresponding to a fixed p-value cutoff.

            Assumes y_score = 1.0 - p_value.

            Parameters
            ----------
            curves : list of dicts
                Each dict must have:
                    {
                    "name": "Test A",
                    "y_true": array-like (0/1),
                    "y_score": array-like (1 - p_value)
                    }
            kind : {"roc", "prc"}
            ax : matplotlib Axes, optional
            title : str, optional
            p_value : float or None
                If provided, plot a point corresponding to p_value <= this cutoff.
                If None, no point is plotted.
            point_kwargs : dict, optional
                Passed to ax.scatter for the threshold point.
            """
            kind = kind.lower()
            if kind not in {"roc", "prc"}:
                raise ValueError("kind must be 'roc' or 'prc'")

            if ax is None:
                fig, ax = plt.subplots()

            if point_kwargs is None:
                point_kwargs = dict(s=60, zorder=5)

            score_threshold = None
            if p_value is not None:
                score_threshold = 1.0 - p_value

            for c in curves:
                y_true = np.asarray(c["y_true"])
                y_score = np.asarray(c["y_score"])
                name = c.get("name", "model")

                if kind == "roc":
                    fpr, tpr, _ = roc_curve(y_true, y_score)
                    auc = roc_auc_score(y_true, y_score)
                    ax.plot(fpr, tpr, label=f"{name} (AUC={auc:.3f})")

                else:  # PRC
                    precision, recall, _ = precision_recall_curve(y_true, y_score)
                    ap = average_precision_score(y_true, y_score)
                    ax.plot(recall, precision, label=f"{name} (AP={ap:.3f})")

                # ----- threshold point -----
                if score_threshold is not None:
                    y_pred = (y_score >= score_threshold).astype(int)

                    n_pos = np.sum(y_pred)
                    n_total = len(y_pred)

                    print(f"Called positive: {n_pos} / {n_total} ({n_pos/n_total:.4%})")

                    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()

                    if kind == "roc":
                        fpr_pt = fp / (fp + tn) if (fp + tn) > 0 else 0.0
                        tpr_pt = tp / (tp + fn) if (tp + fn) > 0 else 0.0
                        ax.scatter(fpr_pt, tpr_pt, **point_kwargs)

                    else:  # PRC
                        recall_pt = tp / (tp + fn) if (tp + fn) > 0 else 0.0
                        precision_pt = tp / (tp + fp) if (tp + fp) > 0 else 1.0
                        ax.scatter(recall_pt, precision_pt, **point_kwargs)

            # ----- baselines + labels -----
            if kind == "roc":
                ax.plot([0, 1], [0, 1], linestyle="--", linewidth=1, label="chance")
                ax.set_xlabel("False Positive Rate")
                ax.set_ylabel("True Positive Rate")
                ax.set_title(title or "ROC Curve")
            else:
                base = np.mean(np.asarray(curves[0]["y_true"]))
                ax.hlines(base, 0, 1, linestyles="--", linewidth=1, label=f"baseline={base:.3f}")
                ax.set_xlabel("Recall")
                ax.set_ylabel("Precision")
                ax.set_title(title or "Precision–Recall Curve")

            ax.legend()
            ax.grid(True, alpha=0.3)
            return ax
        

        curves=[]
        #Now, for each test/median replicate...
        for test in median_reps.keys():
            eval_df=self._merge_in_ground_truth(hypothesis_set_name=hypothesis_set_name,
                                    test_type=test,
                                    index=median_reps[test])
            eval_df["meta"] = eval_df["meta"].fillna(0)
            eval_df=eval_df.dropna()

            y_true = np.asarray((~eval_df["gt_null"]).astype(int))
            y_score = np.asarray(1.0 - eval_df["p_value"])

            curves.append({"name":test,"y_true":y_true,"y_score":y_score})
        
        p_value=None
        if include_alpha:
            p_value=self.get_state_field("alpha")

        plot_curves(curves,kind=performance_type,p_value=p_value)

    def performance_barchart(self,hypothesis_set_name,metric):
        """
        metric is one of auroc, auprc
        """
        self._block_until_all_tests_are_done()
        df=self._all_classifier_summary(hypothesis_set_name)
        _plot_test_bars(df,metric)

    def _switch_cov_method(self,cov_method):
        """
        Helper function which validates passed cov method.
        """
        if cov_method =="opg":
            precompd=self.opgd
        elif cov_method =="sandwich":
            precompd=self.sandwichd
        else:
            raise ValueError(f"Unrecognized cov_method \'{cov_method}\'. Valid options are \'sandwich\' and \'opg\'.")
        return precompd
    
    def precompute_wald(self,cov_method="sandwich"):
        """
        cov_method : {"sandwich", "opg"}
        - See ortho.precompute_wald for more information.
        """
        
        if "ortho" not in self.futures.columns:
            raise RuntimeError("Cannot precompute wald if orthos have not been fit. You want to call `fit_orthos`.")

        #pull out the futures tracking 
        ortho_futures=self.futures["ortho"].to_list()

        #pick output dir destination based on covariance matrix estimation procedure
        precompd=self._switch_cov_method(cov_method)
        
        #make output directory
        precompd.mkdir(exist_ok=True)

        #function to submit
        def _precomp_wald_helper(ortho_future, path_input, path_output, name):
            client=get_client()
            #load the ortho
            orth=ortho.load(client,
                path=path_input,
                name=name
            )
            #perform precompute
            orth.precompute_wald(client,cov_method=cov_method)
            #save to disc
            orth.save(path=path_output,
                    name=name)
            
            return True

        precomp_tracker=[]

        for idx in range(0,self.get_state_field("n_sims")):
            
            r=self.client.submit(_precomp_wald_helper,
                                ortho_future=ortho_futures[idx],
                                path_input=self.orthod,
                                path_output=precompd,
                                name=str(idx))
            precomp_tracker.append(r)
        
        self.futures[cov_method]=precomp_tracker

    def fit_orthos(self, direction="both", serial_orthos: bool = False):
        """
        Applies ortho filtering fits & saves orthos for all simulated replicates.
        Note that this can spawn some very heavy functions!

        direction can be 'both', 'by_cre' or 'by_cell_type'

        Orchestration is intentionally performed on the client (rather than by
        submitting a long-running helper task per replicate). This avoids keeping
        idle helper tasks resident on workers while they wait on upstream futures.

        If serial_orthos is True, each replicate will be fit strictly in index
        order. Otherwise the method still orchestrates from the client, but does
        not enforce a dependency chain between replicate fits.
        """
        valid_directions = ["both", "by_cre", "by_cell_type"]
        if direction not in valid_directions:
            raise ValueError(f"Invalid direction {direction}, valid directions are {valid_directions}")

        # Check to make sure previous step has at least been queued.
        if "transcription" not in self.futures.columns:
            raise RuntimeError(
                "Tried to fit orthos, but transcription has not yet been simulated! "
                "You probably want to run 'gamut' first."
            )

        # Pull out the futures tracking
        tscription_futures = self.futures["transcription"].to_list()

        # Make output directory
        self.orthod.mkdir(exist_ok=True)

        def _fit_one_ortho(idx: int):
            data = scMPRA_data.from_parquet(self.scmpradatp / f"{idx}.scmpra")
            data.ortho_filter()

            primordial = ortho()
            if direction == "both":
                primordial.criss_cross(client=self.client, dat=data)
            elif direction == "by_cre":
                primordial.fit_by_cre_models(client=self.client, dat=data)
            elif direction == "by_cell_type":
                primordial.fit_by_cell_type_models(client=self.client, dat=data)

            primordial.extract_params(self.client)
            primordial.save(path=self.orthod, name=str(idx))

            # Track completion with an already-resolved future so downstream
            # steps that expect futures continue to work unchanged.
            return self.client.scatter(True, hash=False)

        n_sims = self.get_state_field("n_sims")
        ortho_tracker = [None] * n_sims

        if serial_orthos:
            for idx in range(0, n_sims):
                tscription_futures[idx].result()
                ortho_tracker[idx] = _fit_one_ortho(idx)
        else:
            future_to_idx = {fut.key: idx for idx, fut in enumerate(tscription_futures)}
            for fut in dask.distributed.as_completed(tscription_futures):
                idx = future_to_idx[fut.key]
                ortho_tracker[idx] = _fit_one_ortho(idx)

        self.futures["ortho"] = ortho_tracker
    
    def _simulate_transcription(self):
        """
        Simulates transcription.
        Most of the logic is offloaded to the non-method simulate_from_description.
        """
        #check to make sure previous step has at least been queued...
        if "transfection" not in self.futures.columns:
                raise RuntimeError("Tried to simulate transcription, but transfection has not yet been simulated!")
        
        tfection_futures=self.futures["transfection"].to_list()

        #make output directory
        self.scmpradatp.mkdir(exist_ok=True)
        
        n_sims=self.get_state_field("n_sims")

        #helper function to be submitted to the cluster...
        def _simulate_transcription_helper(tfection_fut,description_path,path,negative_controls,reference_cell_type,flatten):
            
            #load a description
            description=pd.read_csv(description_path,
                                    sep="\t",
                                    compression="gzip")

            working=simulate_from_description(description)
            working=working.rename(columns={'zinb_sample':'umis_mpra_bc'})
            scd=scMPRA_data()
    
            scd.data=working
            
            scd.set_negative_controls(negative_controls)
            scd.set_reference_cell(reference_cell_type)
            scd.flag_synthetic()
            scd.overtransfected()
            if flatten:
                scd.flatten_overtransfection()

            scd.to_parquet(path)

            return True
            
        transcription_tracker=[]
        #main loop
        for idx in range(0,n_sims):
            #compute the output path
            path=self.scmpradatp/f"{idx}.scmpra"
            #submit the job
            r=self.client.submit(_simulate_transcription_helper,
                               tfection_fut=tfection_futures[idx],
                               description_path=self.descripd/f"{idx}.tsv.gz",
                               path=path,
                               negative_controls=self.get_state_field("negative_controls"),
                               reference_cell_type=self.get_state_field("reference_cell_type"),
                               flatten=self.get_state_field("flatten_overtransfection"))
            #append the 
            transcription_tracker.append(r)
        
        self.futures["transcription"]=transcription_tracker

    def _block_until_all_tests_are_done(self):
        dummy = [f.result() for f in self.testqueue]
    
    def save(self):
        state_path=self.location/self.name/Path("state.parquet")
        self.state.to_parquet(state_path)

        self._block_until_all_tests_are_done()
        
        
        #block until all queued steps are done
        fut=self.futures.to_numpy().ravel()
        dask.distributed.wait(fut)
        
        results = [f.result() for f in fut]
        
        arr = pd.array(results, dtype="bool").to_numpy().reshape(self.futures.shape)
        resolved = pd.DataFrame(arr, index=self.futures.index, columns=self.futures.columns)

        resolved.to_csv(self.fullp/"futures.tsv.gz",
                sep="\t",
                compression="gzip")

    def mwu(self,name):
        """
        Takes a name of a hypotheses set added previously with 
        `add_hypothesis_set` and runs mann whitney u tests.
        NOTE: Tests are performed directly on simulated data.
        This means that ortho filtering is not applied!
        """
        #init & load relevant hypothesis set...
        hypod=self.testd/name
        hypof=hypod/"hypotheses.tsv"
        if not hypof.is_file():
            raise FileNotFoundError(f"Could not find hypothesis set \'{name}\'.")
        
        hypotheses=HypothesisSet.from_tsv(hypof)

        testd=hypod/"mwu"

        testd.mkdir()
        
        #check to make sure previous step has at least been queued.
        if "transcription" not in self.futures.columns:
            raise RuntimeError("Tried to fit perfrom mwu testing, but transcription has not yet been simulated! You probably want to run 'gamut' first.")
        
        #pull out the futures tracking 
        tscription_futures=self.futures["transcription"].to_list()
        
        #little helper function to shuttle to workers
        def _mwu_helper(tscription_future,
                        path_scmpradat,
                        path_output,
                        hypothesis_set):
            #load the data
            dat = scMPRA_data.from_parquet(path_scmpradat)
            #test
            tester = HypothesisTester("mwu")
            results = tester.run(hypothesis_set, dat)
            #save
            results.to_tsv(path_output)
            
            return True
        
        #submit jobs
        for idx in range(0,self.get_state_field("n_sims")):
            r=self.client.submit(_mwu_helper,
                    tscription_future=tscription_futures[idx],
                    path_scmpradat=self.scmpradatp/f"{idx}.scmpra",
                    path_output=testd/f"{idx}_results.tsv",
                    hypothesis_set=hypotheses)
            
            self.testqueue.append(r)
            
    def wald(self, name, cov_method="sandwich", serial_orthos: bool = True):
        """
        Takes a name of a hypotheses set added previously with
        `add_hypothesis_set` and runs wald tests.
        Requires that `precompute_wald` have been computed previously.

        cov_method : {"sandwich", "opg"}
        - See ortho.precompute_wald for more information.

        If serial_orthos is True, only one _wald_helper will run at a time by
        chaining each submission to depend on the previous wald future.
        """
        # pick input dir destination based on covariance matrix estimation procedure
        precompd = self._switch_cov_method(cov_method)

        # check to make sure previous step has at least been queued.
        if cov_method not in self.futures.columns:
            raise RuntimeError(
                f"wald precompute with coariance estimation method '{cov_method}' not yet run."
            )

        # extract precomp futures
        precomp_futures = self.futures[cov_method]

        # init & load relevant hypothesis set...
        hypod = self.testd / name
        hypof = hypod / "hypotheses.tsv"
        if not hypof.is_file():
            raise FileNotFoundError(f"Could not find hypothesis set '{name}'.")

        hypotheses = HypothesisSet.from_tsv(hypof)

        testd = hypod / f"wald_{cov_method}"
        testd.mkdir()

        # define helper function to submit
        def _wald_helper(
            precomp_future,
            input_dir,
            name,
            path_output,
            hypothesis_set,
            _prev_wald=None,
        ):
            # NOTE: precomp_future and _prev_wald are passed to enforce dependencies.
            # We intentionally do NOT call .result() on them here to avoid blocking a worker.

            client = get_client()

            # load the relevant ortho
            test_ortho = ortho.load(client=client, path=input_dir, name=name)

            # run the tests
            tester = HypothesisTester("wald")
            results = tester.run(hypothesis_set, test_ortho, client=None)

            # save
            results.to_tsv(path_output)

            return True

        # submit jobs
        prev_wald_future = None
        n_sims = self.get_state_field("n_sims")

        for idx in range(0, n_sims):
            kwargs = dict(
                precomp_future=precomp_futures[idx],
                input_dir=precompd,
                name=str(idx),
                path_output=testd / f"{idx}_results.tsv",
                hypothesis_set=hypotheses,
            )

            if serial_orthos:
                kwargs["_prev_wald"] = prev_wald_future

            r = self.client.submit(_wald_helper, **kwargs)

            self.testqueue.append(r)
            prev_wald_future = r
        
    def add_hypothesis_set(self,name:str,hypotheses:HypothesisSet):
        """
        Function adds a hypothesis set to the object for testing.
        'name' will be the name of a directory, so usual caveats apply based on your filesystem
        (avoid special characters, spaces...)
        To run the actual tests, call one of the associated functions (test_wald, test_mwu).
        """
        hypod=self.testd/name
        
        try:
            hypod.mkdir(parents=True)
        except FileExistsError as e:
            raise FileExistsError(
                    f"Hypothesis test with name {name} already exists"
                ) from e
        
        hypotheses.to_tsv(hypod/"hypotheses.tsv")
        
    def list_tests(self):
        """
        Lists all hypothesis sets and the tests which have been run on them.
        """
        root = self.testd
        tree = {}

        for level1 in sorted(p for p in root.iterdir() if p.is_dir()):
            tree[level1.name] = sorted(
                p.name for p in level1.iterdir() if p.is_dir()
            )

        return tree


def cast_string_keys(df, keys):
    df = df.copy()
    for k in keys:
        df[k] = (
            df[k]
            .astype("string")   # pandas StringDtype, not object
            .str.strip()        # kill whitespace
        )
    return df

def _simulate_transfection(experiment_bounds:Bounds,
                        ground_truth:pd.DataFrame,
                        library:pd.DataFrame):
    """
    Simulates transfection, producing a description dataframe
    from which transcription can be simulated...

    See README spec for details on ground truth dataframe.
    You can easially create one with the helper function `simple_spread`. 
    """

    ground_truth=cast_string_keys(ground_truth,["cell_type", "cre_id"])

    #known before you start : "to be optimized":
    #  cells per cell-type is a fixed parameter
    #  barcodes per CRE is a fixed parameter (library)

    def _simulate_single_replicate_transfection():
        #get cells_per_cell_type
        cells_per_cell_type=experiment_bounds.cells_per_cell_type
        cells_per_cell_type.name= 'cells_per_cell_type'
        cells_per_cell_type.index.name= 'cell_type'    
        #copy cells df out
        cells_df=pd.DataFrame(cells_per_cell_type).reset_index()
        #repeat to create multiple cells for each cell type
        cells_df=cells_df.loc[cells_df.index.repeat(cells_df["cells_per_cell_type"])]
        #drop number of cells
        cells_df=cells_df.drop(columns=["cells_per_cell_type"]).reset_index(drop=True)
        cells_df["cell_bc"]=generate_barcodes(length=20,count=len(cells_df))
        #now get "how many MPRA constructs transfected into each cell"
        #since this can draw a zero, some cells may effectively drop out at this
        cells_df["num_transfected"]=experiment_bounds.transfection_model.draw_nb(len(cells_df))
        #duplicate so we have one row for each transfection event
        cells_df=cells_df.loc[cells_df.index.repeat(cells_df["num_transfected"])].reset_index(drop=True)
        #drop num transfected, since it is no longer required
        cells_df=cells_df.drop(columns=["num_transfected"])

        #for each transfection event, sample an MPRA barcode
        drawn_library=sample_from_library(library=library,
                                size=len(cells_df))
        
        #logger.info(f"A: drawn_library cols: {drawn_library.columns}, types: {drawn_library.dtypes}")
        #logger.info(f"A: cells_df cols: {cells_df.columns}, types: {cells_df.dtypes}")
        
        #merge dataframes
        
        cells_df=cells_df.merge(drawn_library,
                                left_index=True,
                                right_index=True,
                                validate="one_to_one")
        
        #logger.info(f"1 {cells_df.isna().sum().sum()}")
        #logger.info(f"B: cells_df cols: {cells_df.columns}, types: {cells_df.dtypes}")
        
        keys = ["cell_type", "cre_id"]
        
        # drop library abundance, we don't care anymore.
        cells_df=cells_df.drop(columns=["abundance"])
        cells_df = cast_string_keys(cells_df, ["cell_type", "cre_id"])
        
        #logger.info(f"C: cells_df cols: {cells_df.columns}, types: {cells_df.dtypes}")
        #logger.info(f"C.5: cells_df cols: {ground_truth.columns}, types: {ground_truth.dtypes}")

        # merge in ground truth

        #dump_df_pickle_debug(cells_df,prefix="permerge_cells_df")
        #dump_df_pickle_debug(ground_truth,prefix="premerge_gt")
        
        # left: maybe by chance an MPRA bc was never transfected.
        cells_df=cells_df.merge(ground_truth,
                                on=["cell_type","cre_id"],
                                validate="many_to_one",
                                how="left")
        
        #dump_df_pickle_debug(cells_df,prefix="postmerge")

        #logger.info(f"D: cells_df cols: {cells_df.columns}, types: {cells_df.dtypes}")
        
        #logger.info(f"D {cells_df.isna().sum().sum()}")

        #check to make sure there were no NAs introduced
        #assert not cells_df.isna().any().any(), "DataFrame contains NA values! Check to make sure cell type & CRE names in all parameters match."

        #dump_df_debug(cells_df)
        
        #TEMP DEBUG OVERRIDE
        #return cells_df
        
        bad = cells_df[cells_df.isna().any(axis=1)]
        assert bad.empty, (
            "DataFrame contains NA values in these rows:\n"
            f"{bad}\n\n"
            "Check that cell type & CRE names match across parameters."
        )
        
        return cells_df

    #loop & generate one round of transfection per replicate
    #I think _simulate_single_replicate_transfection is too small to be profitably parallelized
    #that is, the overhead will probably be more expensive than gains from parallelization
    #so we won't bother.
    all_rep_cells_df=[]
    for rep_id in experiment_bounds.zi.index:
        working=_simulate_single_replicate_transfection()
        working["rep_id"]=rep_id
        working["zi"]=experiment_bounds.zi.at[rep_id]
        all_rep_cells_df.append(working)
    
    all_rep_cells_df=pd.concat(all_rep_cells_df)
    all_rep_cells_df["theta"]=experiment_bounds.theta
    all_rep_cells_df=all_rep_cells_df.rename({"true_mean":"mu"},axis=1)
    

    ret=all_rep_cells_df
    
    #compute alternate parametrization
    #redundant code with `def describe_parameters()`
    ret["r"]=ret["theta"]
    ret["sigmasquare"]=ret["mu"]**2/ret["r"]+ret["mu"]
    ret["p"]=ret["mu"]/ret["sigmasquare"]
    #handle case where mu is zero
    ret.loc[ret["mu"]==0.0,"p"]=0.0

    return ret


def volcano(results: "ResultSet", title = None, bh_thresh=0.05, fc_thresh=1.0):
    """
    Volcano plot using BH-corrected p-values (bh_p) versus log2 fold change.

    Parameters
    ----------
    results : ResultSet
        Must contain columns: 'fold_change', 'bh_p'.
    bh_thresh : float
        FDR threshold for significance (default 0.05).
    fc_thresh : float
        Absolute log2 fold change threshold for vertical lines.
    """
    df = results.to_dataframe().copy()

    # Calculate log2 fold change and -log10(BH p-value)
    df["log2FC"] = np.log2(df["fold_change"].replace(0, np.nan))
    df["neg_log10_bh"] = -np.log10(df["bh_p"])

    # Define significance mask
    sig_mask = (df["bh_p"] < bh_thresh) & (df["log2FC"].abs() > fc_thresh)

    # Map significance to descriptive labels
    sig_labels = {
        True: f"Significant (BH p<{bh_thresh}, |log2FC|>{fc_thresh})",
        False: "Not significant"
    }
    df["sig_label"] = sig_mask.map(sig_labels)

    # Plot
    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        data=df,
        x="log2FC",
        y="neg_log10_bh",
        hue="sig_label",
        palette={
            f"Significant (BH p<{bh_thresh}, |log2FC|>{fc_thresh})": "red",
            "Not significant": "grey"
        },
        alpha=0.7,
        edgecolor=None
    )

    # Threshold lines
    plt.axhline(-np.log10(bh_thresh), color="black", linestyle="--", lw=1)
    plt.axvline(fc_thresh, color="black", linestyle="--", lw=1)
    plt.axvline(-fc_thresh, color="black", linestyle="--", lw=1)

    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10(BH-adjusted p-value)")
    plt.title(title)
    plt.tight_layout()
    plt.show()
