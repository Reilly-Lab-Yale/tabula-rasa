#        1         2         3         4         5         6         7         8
#2345678901234567890123456789012345678901234567890123456789012345678901234567890

#all the main functions.

#external imports
import seaborn as sns
import pandas as pd
import numpy as np

#internal imports
from .utils import unimplemented

@unimplemented
def always_unfinished():
    """tests unimplemented decorator."""
    pass


def helloworld():
    print("hello world!")
    pass


@unimplemented
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
    pass

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
    #assert table_type(table.columns)==
    pass



@unimplemented
def load_scMPRA_data(filepath):
    """
    Arguments
        filepath <str>
    Returns
        <pd.DataFrame>

    Loads scMPRA data from `filepath`. Determines whether it is umi or read-wise.
    
    """
    pass

#        1         2         3         4         5         6         7         8
#2345678901234567890123456789012345678901234567890123456789012345678901234567890

@unimplemented
def graph_qc_metrics(scmpra_data, *args, **kwargs):
    """
    Arguments
        scmpra_data <pd.DataFrame>
        *args
        **kwargs
    Returns
        <matplotlib.axes._axes.Axes>

    Takes `scmpra_data`, a pandas dataframe of read-wise MPRA data (see docs) 
    and returns a histogram of frequency of reads per UMI using seaborn.histplot. 
    All other arguments are passed to the histplot call to allow graph 
    customization.
    """
    assert table_type(scmpra_data.columns) == "mpra_readwise"
    pass

@unimplemented
def cut_chimeric_reads(scmpra_data, threshold):
    """
    Arguments
        `scmpra_data` <pandas.DataFrame> of read-wise scMPRA data 
        `threshold` <int>.
    
    Returns
        a pandas dataframe of umi-wise MPRA data
    
    subsets to those UMIs which lie above the number-of-reads threshold, 
    removing chimeric reads. 
    """
    assert table_type(scmpra_data.columns) == "mpra_readwise"
    pass

@unimplemented
def hypothesis_tester(scmpra_data, hypotheses, flavor="wald"):
    """
    Arguments
        scmpra_data: 
        hypotheses :
    `flavor`, a string with values 

    a pandas dataframe of umi-wise scMPRA data
    """
    assert table_type(scmpra_data.columns) == "mpra_umiwise"
    assert table_type(hypotheses.columns) == "hypotheses"
    pass

@unimplemented
def fit(scmpra_data,round_down_threshold=4):
    """
        requires collapsed-on-barcodes, umi-wise
        impl. round-down for very low CREs...
    """
    pass

@unimplemented
def extract_parameters():
    pass