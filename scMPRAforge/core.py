#        1         2         3         4         5         6         7         8
#2345678901234567890123456789012345678901234567890123456789012345678901234567890

#all the main functions.

#external imports
import seaborn as sns
import pandas as pd
import numpy as np
import warnings

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
    #assert table_type(table.columns)=="hypothesis" or table_type(table.columns)=="results"
    #return table
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



@unimplemented
def graph_qc_metrics(scmpra_data, actually_plot=True *args, **kwargs):
    """
    Arguments
        scmpra_data <pd.DataFrame>
        actually_plot <bool>
        *args
        **kwargs
    Returns
        <matplotlib.axes._axes.Axes>

    Takes `scmpra_data`, a pandas dataframe of read-wise MPRA data (see docs) 
    and returns a histogram of frequency of reads per UMI using seaborn.histplot. 
    if actually_plot, will do the plotting. Otherwise will simply return the plot

    All other arguments are passed to the histplot call to allow graph 
    customization.
    """
    assert table_type(scmpra_data.columns) == "mpra_readwise"
    #if actually_plot
    #   plot in-place to avoid making the user import seaborn if they don't want to
    pass

@unimplemented
def cut_chimeric_reads(scmpra_data, threshold):
    """
    Arguments
        scmpra_data <pandas.DataFrame> of read-wise scMPRA data 
        threshold <int>.
    
    Returns
        a pandas dataframe of umi-wise MPRA data
    
    subsets to those UMIs which lie above the number-of-reads threshold, 
    removing chimeric reads. 
    """
    assert table_type(scmpra_data.columns) == "mpra_readwise"
    assert threshold >=0
    pass

#        1         2         3         4         5         6         7         8
#2345678901234567890123456789012345678901234567890123456789012345678901234567890

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

@unimplemented
def fit(scmpra_data,round_down_threshold=4):
    """
    Arguments
        scmpra_data
    Returns
        {
            'model':<statsmodels.discrete.count_model.ZeroInflatedNegativeBinomialResultsWrapper>,
            'flattened': <pd.DataFrame>
        }

    """
    # for all (cre, cell-type) combos with less than 4 umis:
    #   add index to list `too_low`

    #select out the cre-celltype combos we have to ditch because we have no UMIs making fitting impossible. 
    # flattened=scmpra_data.groupby()[too_low,["cre_id","cell_type"]]

    #drop those flattened barcodes from the original

    #scmpra_data=scmpra_data.drop(too_low).copy()

    #if nrow(scmpra_data) ==0:
    #   raise error("")
    
    #!!do the fitting...!!
    #model.fit will produce an object of type <statsmodels.discrete.count_model.ZeroInflatedNegativeBinomialResultsWrapper>
    
    #if !converged:
    #   print("error, model failed to converge.")

    #if r^2 <= 0.6
    #warnings.warn("be careful, model fit is pretty bad : pseudo-r^2 is only f{r^2}")

    pass

@unimplemented
def volcano(results):
    """
    Volcano plot of p value versus log fold change
    """
    pass

@unimplemented
def extract_parameters():
    pass