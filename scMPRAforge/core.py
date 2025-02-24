#        1         2         3         4         5         6         7         8
#2345678901234567890123456789012345678901234567890123456789012345678901234567890

#all the main functions.

#external imports
import seaborn as sns

#internal imports
from .utils import unimplemented

@unimplemented
def always_unfinished():
    """tests unimplemented decorator"""
    pass


def helloworld():
    print("hello world!")
    pass

@unimplemented
def load_hypothesis_set():
    pass

@unimplemented
def load_scMPRA_data(filepath,wise,collapsed):
    """
        (for development, assume read-wise, uncollapsed)
        Takes 
        
    """
    pass

@unimplemented
def graph_qc_metrics(scmpra_data, *args, **kwargs):
    """
        Takes `scmpra_data`, a pandas dataframe of read-wise MPRA data 
        (see docs) and plots a histogram of frequency of reads per UMI
        using seaborn.histplot. All other arguments are passed to
        the histplot call to allow graph customization.
    """
    pass

@unimplemented
def cut_chimeric_reads(scmpra_data, threshold):
    """
        Takes `scmpra_data`, a pandas dataframe of read-wise scMPRA data 
        (see docs) and `threshold` <int>.
        Returns a pandas dataframe of umi-wise MPRA data (see docs)
        subset to those UMIs which lie above the number-of-reads
        threshold.
    """
    pass

@unimplemented
def hypothesis_tester(scmpra_data, hypotheses, flavor="wald"):
    """
        Takes scmpra_data, a pandas dataframe of umi-wise scMPRA data
    """
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