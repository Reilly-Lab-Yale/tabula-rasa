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



#        1         2         3         4         5         6         7         8
#2345678901234567890123456789012345678901234567890123456789012345678901234567890
class experiment_model:
    """
        type 
    """
    
    #low-tech, maybe replace w/ enum
    valid_formula_types=['auto','tfection_umi','custom_formula']

    def __init__(self,scmpra_data,formula_type='auto',round_down_threshold=4,formula=None,broken_on=None):
        #formula and broken_on allow creation of custom formulas.
        #broken_on allows breaking into submodels on the basis of some column values

        #saving relevant info        
        self.scmpra_data=scmpra_data
        self.round_down_threshold=round_down_threshold
        self.model=None
        #(creating from scratch so no model until after fit is called)
        

        if formula_type not in self.valid_formula_types:
            assert False; 'Invalid formula type.'
            
        if formula_type=='auto':
            #will try to guess formula type based on 
            #table_type(scmpra_data)
            assert False; 'unimplemented'

            #formula_type='whatever we just detected'

        if formula_type=='dna':
            #using plasmid DNA count as DNA count
            assert False; 'unimplemented'
        elif formula_type=='tfection_umi':
            #using umis_transfection_bc as our proxy for DNA count...
            self.formula="reads_mpra_bc ~ umis_transfection_bc:C(cell_type) + \
                umis_transfection_bc:C(cre_id) + umis_transfection_bc:C(cell_type):C(cre_id) -1"
        elif formula_type=='num_mpra_bc':
            #using the number of unique MPRA barcodes as our proxy for DNA count
            assert False; 'unimplemented'

        self.formula_type=formula_type

    @unimplemented
    @classmethod
    def from_pickle(filepath,get_data):
        #create instance by loading previously created model from disc
        pass

    @unimplemented
    @classmethod
    def to_pickle(self):
        pass

    def fit(self,library="tensor"):
        # for all (cre, cell-type) combos with less than 4 umis:
        #   add index to list `too_low`

        #select out the cre-celltype combos we have to ditch because we have no UMIs making fitting impossible. 
        # flattened=scmpra_data.groupby()[too_low,["cre_id","cell_type"]]

        #drop those flattened barcodes from the original

        #scmpra_data=scmpra_data.drop(too_low).copy()

        #if nrow(scmpra_data) ==0:
        #   raise error("")
        
        #!!do the fitting...!!

        print("creating matrices...")
        #old patsy code
        #y, X = patsy.dmatrices(self.formula,
        #                        self.scmpra_data, return_type='dataframe')
        #Z = patsy.dmatrix("C(rep_id)", self.scmpra_data, return_type='dataframe')

        y, X=Formula(self.formula).get_model_matrix(self.scmpra_data,output='pandas')
        Z=Formula('C(rep_id)').get_model_matrix(self.scmpra_data,output='pandas')

        if library =="statsmodels":
            print("fitting with statsmodels")

            zinb_model = smdc.ZeroInflatedNegativeBinomialP(y, X, exog_infl=Z)

            n_count_params = zinb_model.exog.shape[1]      # Count model parameters
            n_infl_params = zinb_model.exog_infl.shape[1]    # Inflation model parameters
            n_total = n_count_params + n_infl_params + 1 # adding 1 for alpha
            start_params = np.full(n_total, 0.1)

            self.model = zinb_model.fit(start_params=start_params,maxiter=1000)
        elif library=="tensor":
            print("fitting with tensorzinb")
            zinbo=TensorZINB(y['reads_mpra_bc'].to_numpy().reshape((-1,1)),X.to_numpy(),exog_infl=Z.to_numpy())#,same_dispersion=True
            self.model=zinbo.fit(init_method="nb")
        else:
            assert False; "Invalid library value."


        
        
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