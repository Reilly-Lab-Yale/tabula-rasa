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

from scipy.stats import linregress
#import statsmodels.discrete.count_model as smdc
import patsy
from tensorzinb.tensorzinb import TensorZINB
from formulaic import Formula

from dask.distributed import Client
from dask.distributed import Future

import dask.dataframe as dd
import dask.array as da

from enum import Enum

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

    TODO: move to the inside of the scMPRA object
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
    
    @unimplemented
    @classmethod
    def from_json(cls,filepath):
        """
        Returns a <scMPRA_data> object with data loaded from `filepath`.
        """
        pass

    @unimplemented
    def to_json(self,filepath):
        """
        Dump object to filepath
        """
        pass
    
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

    def save(self,path,name):
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
        simple_write(self.training_data,"training_data.pkl")

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
    def __init__(self, nb, zi, theta, broken_on):
        self.nb=nb
        self.zi=zi
        self.theta=theta

        self.broken_on = broken_on

        assert nb.keys() == zi.keys()
        assert zi.keys() == theta.keys()

        self.keys=list(nb.keys())

    #def _flatten_out_futures(self):
    #    """
    #    For internal use only. 
    #    Will bork downstream operations on the object.
    #    """
    #    for key in self.nb:
    #        self.nb[key]=self.nb[key].result()
    #        self.zi[key]=self.zi[key].result()
    #        self.theta[key]=self.theta[key].result()
        

    def _unflatten_futures(self,client):
        """
        For internal use only
        wraps all the models in futures
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


def describe_parameters(client,parameters,dat,split):
    """
    Produces a convienient single-dataframe description of one split model. 
    'parameters' is a parameters object
    Requires that you pass original data as dat to compute cell-numbers
    Leaves non-split columns as one-hot. Returns "split" column as str categorical
    """

    #change to return a dask instead of pandas dataframe


    #count cells per group
    cell_counts=dat.groupby([split,"rep_id"])["umis_mpra_bc"].agg("sum")
    cell_counts=pd.DataFrame({"cells":cell_counts})
    
    #cast rep id to string just in case.
    cell_counts.index = cell_counts.index.set_levels(
        cell_counts.index.levels[1].astype(str), level=1
    )

    flattened_param=flatten_param_representation(client,parameters,split=split)
    flattened_param=flattened_param.reset_index().set_index([split,"rep_id"])
    
    working=cell_counts.join(flattened_param,how="left")
    working["r"]=working["theta"]
    working["sigmasquare"]=working["mu"]**2/working["r"]+working["mu"]
    working["p"]=working["mu"]/working["sigmasquare"]
    working=working.reset_index()


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
    description=description.repartition(partition_size="100KB")#5KB
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
    
    def describe_primordial(self,client):
        """Generates and saves descriptions of the primordial which are necessary for subsequent simulation"""
        self.description_primordial_by_cre=describe_parameters(client,
                                                                   parameters=self.primordial.by_cre_parameters.result(),
                                                                   dat=self.primordial.training_data.data,
                                                                   split="cre_id")

        self.description_primordial_by_cre=auto_partition(self.description_primordial_by_cre,
                                                              self.partition_mb)

        self.description_primordial_by_cell_type=describe_parameters(client,
                                                                         parameters=self.primordial.by_cell_type_parameters.result(),
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
