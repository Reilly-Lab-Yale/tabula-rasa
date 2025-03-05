#utility functions

#external imports
import functools
import logging
import umi_tools

logger = logging.getLogger("scMPRAforge")

unimplemented_functions = []  # List to track unimplemented functions

def unimplemented(func):
    """
    Decorator to mark functions as unimplemented.
    Adds them to a global tracking list and breaks when called. 
    """
    global unimplemented_functions
    unimplemented_functions.append(func.__name__)

    note = "\n\n**Note:** This function is not yet implemented."
    if func.__doc__:
        func.__doc__ += note
    else:
        func.__doc__ = note.strip()

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        raise NotImplementedError
        return None  # Placeholder return value

    return wrapper

def list_unimplemented():
    """Returns the list of unimplemented functions."""
    return unimplemented_functions


#        1         2         3         4         5         6         7         8
#2345678901234567890123456789012345678901234567890123456789012345678901234567890



def bcs_to_lut(bc,threshold=1,encoding="utf-8",*args,**kwargs):
    """
    Arguments
        bc <dict> of string keys and integer occuracnce count values
        threshold <int> edit distance
        encoding <str> string encoding

    Returns
        <dict>

    A simple wrapper for umi_tools.UMIClusterer(). `threshold` is the edit 
    distance passed to . args and kwargs are passed to umi_tools.UMIClusterer
    constructor. 
    Produces a lookup table (lut) where the keys are erronious & correct barcodes
    and the values are all the corrected barcodes.
    """

    #convert the strings in bc to bytes
    
    bc={key.encode(encoding):value for key, value in bc.items()}
    
    #Why is this written as an object..? Why not just a function..? such a strange decision.
    clusterer=umi_tools.UMIClusterer(*args,**kwargs)
    fixed=clusterer(bc,threshold=1)
    #the first barcode in each sub-list is the one supported by the most counts, 
    #so we will consider those as the 'correct' values. 
    
    #could maybe vectorize to speed up, but only part of pre-processing so not priority. 
    ret={}
    for cluster in fixed:
        correct_bc=cluster[0]
        for bc in cluster:
            ret[bc]=correct_bc

    # convert the bytes in the dictionary back to strings
    ret={key.decode(encoding):value.decode(encoding) for key, value in ret.items()}

    return ret

## tools for easy generation of hypotheses
@unimplemented
def one_versus_all():
    """
    Provide one negative control, and a list of others to compare against, and 
    this function will generate a hypothesis list comparing all vs it...
    Useful for a quick "which elements are expressed". 
    """
    pass