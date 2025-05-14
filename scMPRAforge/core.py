#        1         2         3         4         5         6         7         8
#2345678901234567890123456789012345678901234567890123456789012345678901234567890

#all the main functions.

#external imports
import seaborn as sns
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np

import logging

import statsmodels.discrete.count_model as smdc
import patsy
from tensorzinb.tensorzinb import TensorZINB
from formulaic import Formula

from enum import Enum

#internal imports
from .utils import unimplemented
from .utils import bcs_to_lut
logger = logging.getLogger("scMPRAforge")

#functions
@unimplemented
def always_unfinished():
    """tests unimplemented decorator."""
    pass


def helloworld():
    print("hello world!")
    pass


#actually probably easier to do an 'inflate' function to avoid too much mem usage...?
@unimplemented
def mpra_unstack():
    pass

@unimplemented
def mpra_stack():
    pass

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
    
    (We could extend to type-checking as well, but that seems a tad draconian.)
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




def load_scMPRA_data(filepath):
    """
    Arguments
        filepath <str>
    Returns
        <pd.DataFrame>

    Loads tsv scMPRA data from `filepath`.
    
    """
    tab=pd.read_csv(filepath,sep="\t")
    tabtype=table_type(tab.columns)
    assert tabtype=="mpra_readwise" or tabtype=="mpra_umiwise", "Malformed table."
    return tab



def graph_chimeric(scmpra_data, *args, **kwargs):
    """
    Arguments
        scmpra_data <pd.DataFrame>
        *args
        **kwargs
    Returns
        <matplotlib.axes._axes.Axes>

    Takes `scmpra_data`, a pandas dataframe of read-wise MPRA data (see docs) 
    and plots a histogram of frequency of reads per UMI using seaborn.histplot. 

    All other arguments are passed to the histplot call to allow graph 
    customization. Particular useful are `bins`, `binrange`, and `log_scale`
    """
    assert table_type(scmpra_data.columns) == "mpra_readwise"
    
    sns.histplot(scmpra_data['reads'], *args, **kwargs)

    plt.xlabel('Reads')
    plt.ylabel('Frequency')
    plt.title('Histogram of Reads')
    plt.show()


def cut_chimeric_reads(scmpra_data, threshold):
    """
    Arguments
        scmpra_data : <pandas.DataFrame> of read-wise scMPRA data 
        threshold : <int>
    
    Returns
        <pandas.DataFrame> of read-wise MPRA data
    
    subsets to those UMIs which lie ABOVE the number-of-reads threshold, 
    removing chimeric reads. 
    """
    assert table_type(scmpra_data.columns) == "mpra_readwise"
    assert threshold >=0, "threshold must be greater than zero."
    
    #Trim
    ret=scmpra_data[scmpra_data["reads"]>threshold]

    original_umi_count=len(scmpra_data["umi"].unique())
    cut_umi_count=len(ret["umi"].unique())

    logger.info(f"Original={original_umi_count} UMIs, Cut={cut_umi_count} UMIs, Lost={original_umi_count-cut_umi_count} UMIs.")

    return ret

#        1         2         3         4         5         6         7         8
#2345678901234567890123456789012345678901234567890123456789012345678901234567890

def read_wise_to_umi_wise(scmpra_data,keep_reads=False,bypass_consistency_check=False):
    """
    Arguments
        scmpra_data : <pandas.DataFrame> of read-wise scMPRA data
        keep_reads : <bool>
    Returns
        <pandas.DataFrame> of umi-wise scMPRA data

    Converts read-wise to UMI-wise table (see readme for spec).
    """
    if not bypass_consistency_check:
        assert table_type(scmpra_data.columns)=="mpra_readwise","Malformed table."

    grouping_columns = [col for col in scmpra_data.columns if col not in ['umi', 'reads']]


    aggregations = {
        'umis': ('umi', 'nunique')  # Count unique UMIs
    }

    # Conditionally include 'reads' sum
    if keep_reads:
        aggregations['reads'] = ('reads', 'sum')

    return scmpra_data.groupby(grouping_columns).agg(**aggregations).reset_index()
    

def flatten_barcode_errors(df, barcode_column,*args,**kwargs):
    """
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

def flatten_mpra_barcodes(scmpra_data):
    """
    Arguments
        scmpra_data : <pd.DataFrame> of umi-wise, full MPRA
    returns
        <pd.DataFrame> of umi-wise, flattened MPRA
    """
    
    assert table_type(scmpra_data.columns)=="mpra_umiwise","Malformed table."
    
    # sum umi counts across MPRA barcodes
    # - if reads present, sum that too

    grouping_columns = [col for col in scmpra_data.columns if col not in ['reads', 'umis','mpra_bc']]


    aggregations={
        'umis':('umis','sum'),
        'mpra_bcs':('mpra_bc','nunique')
    }

    if "reads" in scmpra_data.columns:
        aggregations['reads']=('reads','sum')
    
    return scmpra_data.groupby(grouping_columns).agg(**aggregations).reset_index()


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

"""

#suggested formulas. 
#SUGGESTED_NB=['reads_mpra_bc ~ C(cell_type)*C(cre_id)',
#    'reads_mpra_bc ~ C(cell_type)',
#    'reads_mpra_bc ~ umis_transfection_bc:C(cell_type) + umis_transfection_bc:C(cre_id) + umis_transfection_bc:C(cell_type):C(cre_id) -1']
#SUGGESTED_ZI=['C(replicate)']
#SUGGESTED_BREAKBY=['']



def fit(client,
    nb_formula:str,
    zi_formula:str,
    broken_on:str,
    round_down_threshold:int=4,
    dry:bool=False):
    """
    dry = dry run: don't actually fit anything, just return an experiment_model 
    object initalized
    """
    ret=experiment_model(nb_formula=nb_formula,zi_formula=zi_formula,broken_on=broken_on,round_down_threshold=round_down_threshold)
    if dry:
        return ret
    
    


#        1         2         3         4         5         6         7         8
#2345678901234567890123456789012345678901234567890123456789012345678901234567890
class experiment_model:
    """
    """

    def __init__(self,nb_formula:str,
            zi_formula:str,
            broken_on:str,
            round_down_threshold:int=4):

        self.round_down_threshold=round_down_threshold
        self.model=None
        #(creating from scratch so no model until after fit is called)
        

@unimplemented
def volcano(results:experiment_model):
    """
    Volcano plot of p value versus log fold change
    """
    pass

@unimplemented
def extract_parameters():
    pass