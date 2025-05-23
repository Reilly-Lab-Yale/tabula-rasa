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

MIN_PTS=3

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

extract_triples()

triples

"""

#suggested formulas. 
#SUGGESTED_NB=['reads_mpra_bc ~ C(cell_type)*C(cre_id)',
#    'reads_mpra_bc ~ C(cell_type)',
#    'reads_mpra_bc ~ umis_transfection_bc:C(cell_type) + umis_transfection_bc:C(cre_id) + umis_transfection_bc:C(cell_type):C(cre_id) -1']
#SUGGESTED_ZI=['C(replicate)']
#SUGGESTED_BREAKBY=['']



def create_matricies(nb_formula,zi_formula,data):
    y, X=Formula(nb_formula).get_model_matrix(data,output='pandas')
    Z=Formula(zi_formula).get_model_matrix(data,output='pandas')
    return(X, y, Z)

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


def toosmall(y):
    if sum(y["umis_mpra_bc"].to_numpy()>0)<MIN_PTS:
        return True
    else:
        return False

def tensorzinb_fit(p):
    X,y,Z=p
    if toosmall(y):
        return "too_small"

    from tensorzinb.tensorzinb import TensorZINB

    zinbo=TensorZINB(y["umis_mpra_bc"].to_numpy().reshape((-1,1)),X,exog_infl=Z.to_numpy())#,same_dispersion=True
    zinb_result=None
    try:
        zinb_result=zinbo.fit(init_method="nb")
    except InvalidIndexError:
        return "index_error"


    return zinb_result


def fit(client,
    data,
    nb_formula:str,
    zi_formula:str,
    broken_on:str,
    round_down_threshold:int=4,
    dry:bool=False):
    """
    dry : dry run: don't actually fit anything, just return an experiment_model 
    broken_on : "unified" for one model, or put the name of a column 
    """
    ret=experiment_model(model={},uniq_predictor={},nb_formula=nb_formula,zi_formula=zi_formula,broken_on=broken_on,round_down_threshold=round_down_threshold)
    if dry:
        return ret
    
    ###create design matricies###
    types=None
    mats_futures=None
    if broken_on == "unified":
        #unified model
        types=["unified"]
        mats_futures = {'unified':client.submit(create_matricies,
            data=data,
            zi_formula=zi_formula,
            nb_formula=nb_formula
        )}
    else:
        #not a unified model : let's get all the proper types
        #that is : the name of each value for the column we split on
        types_future = client.submit(lambda df: df[broken_on].unique(), data)
        #wait & grab result (small enough for master node)
        types = types_future.result()
        abort_on_failure(types_future,client)

        #now that we have the types, submit jobs to create design matrices from each data slice
        #we do the slicing in a little sumbitted lambda
        mats_futures = {
            t: client.submit(
                create_matricies,
                data=client.submit(lambda df, t=t: df[df[broken_on] == t], data, t),
                zi_formula=zi_formula,
                nb_formula=nb_formula
            )
            for t in types
        }


    ### create uniq_predictor ### 
    uniq_predictor_futures = {
        t: (
            client.submit(lambda tup: tup[0].drop_duplicates(), mats_futures[t]),
            client.submit(lambda tup: tup[2].drop_duplicates(), mats_futures[t])
        )
        for t in types
    }
    #recall the order is X,y,Z. We only want predictors, which are X,Z, hence 0,2

    ###create statsmodel futures###
    tzinb_futures = {
        t: client.submit(
                tensorzinb_fit,
                mats_futures[t]
            )
        for t in types
    }

    ###put data in proper class and return###
    model_future=client.submit(
        experiment_model,
        model=tzinb_futures,
        uniq_predictor=uniq_predictor_futures,
        nb_formula=nb_formula,
        zi_formula=zi_formula,
        broken_on=broken_on,
        round_down_threshold=round_down_threshold
    )

    return model_future





#        1         2         3         4         5         6         7         8
#2345678901234567890123456789012345678901234567890123456789012345678901234567890
class experiment_model:
    """
    uniq_predictor stores...
    """

    def __init__(self,
            model,
            uniq_predictor,
            nb_formula:str,
            zi_formula:str,
            broken_on:str,
            round_down_threshold:int=4):

        self.nb_formula=nb_formula
        self.zi_formula=zi_formula
        self.broken_on=broken_on
        self.round_down_threshold=round_down_threshold
        self.model=model
        self.uniq_predictor=uniq_predictor
        
        

@unimplemented
def volcano(results:experiment_model):
    """
    Volcano plot of p value versus log fold change
    """
    pass

@unimplemented
def extract_parameters():
    pass