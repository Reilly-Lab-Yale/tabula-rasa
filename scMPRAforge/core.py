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

import statsmodels.discrete.count_model as smdc
import patsy
from tensorzinb.tensorzinb import TensorZINB
from formulaic import Formula

from dask.distributed import Client
import dask.dataframe as dd
import dask.array as da

from enum import Enum

#internal imports
from .utils import unimplemented
from .utils import bcs_to_lut
from .utils import undo_one_hot_encoding
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


def filter_low_umi_count(mpra_data: pd.DataFrame, 
                         group_by: list[str] = ['cre_id', 'cell_type']):
    """
    Takes a umiwise MPRA dataframe and drops rows which would correspond to models
    with not enough cells (observations) to model.

    Returns:
        filtered: the filtered DataFrame
        dropped_groups: a dict mapping each column in group_by to a list of dropped group values
    """
    tabtype = table_type(mpra_data.columns)
    assert tabtype == "mpra_umiwise", "Malformed table."

    filtered = mpra_data.copy()
    dropped_groups = {}

    for col in group_by:
        # Find groups with too few nonzero UMIs
        low_count_mask = (
            filtered.groupby(col)["umis_mpra_bc"]
            .apply(lambda col: (col.to_numpy() > 0).sum() < MIN_PTS)
        )
        dropped = low_count_mask[low_count_mask].index.tolist()
        dropped_groups[col] = dropped

        # Keep only valid rows
        filtered = filtered[~filtered[col].isin(dropped)]

    return filtered, dropped_groups

def toosmall(y):
    if sum(y["umis_mpra_bc"].to_numpy()>0)<MIN_PTS:
        return True
    else:
        return False

def tensorzinb_fit(p):
    X, y, Z = p
    if toosmall(y):
        raise ValueError("One of the models doesn't have enough data for a good fit. Did you run `filter_low_umi_count`?")
        # return "too_small"

    from tensorzinb.tensorzinb import TensorZINB

    zinbo = TensorZINB(y["umis_mpra_bc"].to_numpy().reshape((-1, 1)), X, exog_infl=Z.to_numpy())  # ,same_dispersion=True
    zinb_result = None
    try:
        zinb_result = zinbo.fit(init_method="nb")
    except InvalidIndexError:
        return "index_error"

    return zinb_result

# Timed tensorzinb_fit wrapper

def timed_tensorzinb_fit(p, label=None):
    start = time.time()
    try:
        result = tensorzinb_fit(p)
    except Exception as e:
        result = e
    end = time.time()
    duration = end - start
    logger.warning(f"[timed_fit] {label or 'unknown'} took {duration:.2f} seconds")
    return result


def fit(client,
    data,
    nb_formula:str,
    zi_formula:str,
    broken_on:str,
    round_down_threshold:int=4,
    dry:bool=False,
    return_design_matricies=False):
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
                timed_tensorzinb_fit,
                mats_futures[t],
                label=t
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
    if return_design_matricies:
        return(model_future,mats_futures)
    else:
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



class ortho:
    """
    Stores multiple models of the same data.
    Not to be used with multiple datasets. 
    stores hypothesis sets & coresp. models
    """
    def __init__(self):
        self.by_cre=None
        self.by_cre_parameters=None
        
        self.by_cell_type=None
        self.by_cell_type_parameters=None

    
    def criss_cross(self,client,dat,retain_design_matricies=False):
        """
        Note: a little computationally intensive...
        """
        if retain_design_matricies:
            #mostly for debugging
            self.by_cre, self.by_cre_design=fit(client,dat,nb_formula="umis_mpra_bc ~ C(cell_type)-1",zi_formula="C(rep_id)-1",broken_on="cre_id",return_design_matricies=True)
            self.by_cell_type, self.by_cell_type_design=fit(client,dat,nb_formula="umis_mpra_bc ~ C(cre_id)-1",zi_formula="C(rep_id)-1",broken_on="cell_type",return_design_matricies=True)
        else:
            self.by_cre=fit(client,dat,nb_formula="umis_mpra_bc ~ C(cell_type)-1",zi_formula="C(rep_id)-1",broken_on="cre_id")
            self.by_cell_type=fit(client,dat,nb_formula="umis_mpra_bc ~ C(cre_id)-1",zi_formula="C(rep_id)-1",broken_on="cell_type")

    def extract_params(self,client):
        """Extracts parameters for all models in the object"""
        if not self.by_cre is None:
            self.by_cre_parameters=client.submit(model_to_parameters,self.by_cre)
        
        if not self.by_cell_type is None:
            self.by_cell_type_parameters=client.submit(model_to_parameters,self.by_cell_type)

class parameters:
    def __init__(self, nb, zi, theta):
        self.nb=nb
        self.zi=zi
        self.theta=theta

        assert nb.keys() == zi.keys()
        assert zi.keys() == theta.keys()

        self.keys=nb.keys()


def model_to_parameters(model):
    """
    A function to extract model parameter triples from simple model.
    Currently pretty bespoke: be careful with more complicated models... 
    """
    #should rework to use multi-indexes instead of sorting & assuming they are the same
    
    #currently assuming a single theta per model (equal variances within a model).
    #also does not intelligently sum for interaction effects : so just use for non-interaction for now...
    #Need also to test with incomplete matricies (non-cartesian product) : applying labels will probably break

    #recall the order within each tuple: (X,y,Z), but here we skip y
    #so it is just X, Z
    
    zi={}
    nb={}
    theta={}

    for key in model.model:
        ## theta ##
        theta[key]=model.model[key]['weights']['theta'].squeeze()
        
        ## ZI ##
        #extract min design matrix & sort
        #We need to sort to match the order in the tensorzinb model...
        zi_df=model.uniq_predictor[key][1]
        zi_df=zi_df.sort_values(zi_df.columns.tolist(),ascending=False)

        #multiply out & undo link
        result=(zi_df.values @ model.model[key]["weights"]["x_pi"]).squeeze()
        result=1/(1+np.exp(-result))
        #save multi-hot column names
        col=zi_df.columns.to_list()
        #add nb parameter reconstructions to multi-hot
        zi_df["zi"]=result
        #set all columns other than zi to be a multi-index indexing the zi estimates
        zi_df.set_index(col,inplace=True)
        #add to dictionary...
        zi[key]=zi_df

        ## NB ##
        #extract min design matrix & sort
        nb_df=model.uniq_predictor[key][0]
        nb_df=nb_df.sort_values(nb_df.columns.tolist(),ascending=False)

        #multiply out & undo link
        result=(nb_df.values @ model.model[key]['weights']['x_mu']).squeeze()
        result=np.exp(result)
        #save multi-hot column names
        col=nb_df.columns.to_list()
        #add zi parameter reconstruction to multi-hot
        nb_df["nb"]=result
        #set all columns other than nb to be a multi-index indexing the nb estimates
        nb_df.set_index(col,inplace=True)
        #add to dictionary
        nb[key]=nb_df


    
    return parameters(nb=nb, zi=zi, theta=theta)


def describe_parameters(client,parameters,dat,split):
    """
    Produces a convienient single-dataframe description of one split model. 
    'parameters' is a parameters object
    Requires that you pass original data as dat to compute cell-numbers
    Leaves non-split columns as one-hot. Returns "split" column as str categorical
    """

    #change to return a dask instead of pandas dataframe

    cell_counts=get_cell_counts(client,dat,split=split)
    flattened_param=flatten_param_representation(client,parameters,split=split)
    working=cell_counts.join(flattened_param)
    working["r"]=np.exp(working["theta"])
    working["sigmasquare"]=working["nb"]**2/working["r"]+working["nb"]
    working["p"]=working["nb"]/working["sigmasquare"]


    return working

def flatten_param_representation(client: Client, params, split: str):
    
    params_future = client.scatter(params, broadcast=True)
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
    Simulate from a description dask dataframe using map_partitions instead of da.random for speed.
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

    return description.map_partitions(simulate_partition)

class simulation_batch:
    """
    Class which takes a single <ortho> object and simulates replicates.
    Optionally, fits additional ortho objects to simulations & plots their paremeter spread
    Useful for estimating variance of an experimental setup...
    """
    #several functions modify state in a way necessary for subsq functions: add checks to make sure prev. has been called.
    
    #consolidate split and parameter validity checking

    splits=["cre_id","cell_type"]

    def __init__(self,primordial,partition_mb=50):
        #the initial 
        self.primordial=primordial
        
        self.simulated_from_cre=[]
        self.simulated_from_cell_type=[]

        self.ortho_simulated_cre=[]
        self.ortho_simulated_cell_type=[]

        self.partition_mb=partition_mb
    
    def describe_primordial(self,client,dat):
        """Generates and saves descriptions of the primordial which are necessary for subsequent simulation"""
        self.description_primordial_by_cre=describe_parameters(client,
                                                                   parameters=self.primordial.by_cre_parameters.result(),
                                                                   dat=dat,
                                                                   split="cre_id")

        self.description_primordial_by_cre=auto_partition(self.description_primordial_by_cre,
                                                              self.partition_mb)

        self.description_primordial_by_cell_type=describe_parameters(client,
                                                                         parameters=self.primordial.by_cell_type_parameters.result(),
                                                                         dat=dat,
                                                                         split="cell_type")
        
        self.description_primordial_by_cell_type=auto_partition(self.description_primordial_by_cell_type,
                                                                    self.partition_mb)

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
        return s
    
    def simulate_many(self,client,n):
        """
        simulates n replicates
        Depends on 
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
    
    def _flatten_all_parameters(self):
        """Flatten all parameters in preparation for plotting"""
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
                    working=working.by_cre_parameters.result()
                elif split=="cell_type":
                    working=working.by_cell_type_parameters.result()
                
                working=pd.DataFrame([working.theta]).T
                working=working.reset_index()
                working.columns=[split,"theta"]
                working['theta']=working['theta'].astype(float)
                return working

            if split=="cre_id":
                working=working.by_cre_parameters.result()
            elif split=="cell_type":
                working=working.by_cell_type_parameters.result()
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
            color='black', marker='X',
            s=120, label='primordial', zorder=10
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
            palette='Set1', edgecolor='black', s=50, alpha=0.7, legend='brief'
        )

        # Primordial points overlay (bold black Xs)
        sns.scatterplot(
            data=zi[zi['id'].isin(['primordial cre_id', 'primordial cell_type'])],
            x='group', y='zi',
            color='black', marker='X',
            s=120, label='primordial', zorder=10
        )
        

        plt.xlabel(f'{split} | rep_id')
        plt.ylabel('zi')
        plt.title(f'zi parameters by {split} and rep_id')
        plt.xticks(rotation=45, ha='right')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.show()

    def plot_nb_spread(self):
        nbs=self._nbs
        nbs['group'] = nbs['cell_type'] + " | " + nbs['cre_id']

        plt.figure(figsize=(12, 6))

        # Violin plot
        sns.violinplot(
            data=nbs, x='group', y='nb',
            inner=None, palette='Set1'
        )

        # Regular points
        sns.scatterplot(
            data=nbs[~nbs['id'].isin(['primordial cre_id', 'primordial cell_type'])],
            x='group', y='nb',
            hue='cell_type', style='id',
            palette='Set1', edgecolor='black', s=50, alpha=0.7, legend='brief'
        )

        # Primordial points overlay (bold black Xs)
        sns.scatterplot(
            data=nbs[nbs['id'].isin(['primordial cre_id', 'primordial cell_type'])],
            x='group', y='nb',
            color='black', marker='X',
            s=120, label='primordial', zorder=10
        )

        plt.xlabel('cell_type | cre_id')
        plt.ylabel('nb')
        plt.title('nb parameters by cell_type and cre_id')
        plt.xticks(rotation=45, ha='right')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.show()

@unimplemented
def volcano(results:experiment_model):
    """
    Volcano plot of p value versus log fold change
    """
    pass
