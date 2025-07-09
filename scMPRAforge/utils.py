#utility functions

#external imports
import functools
import logging
import umi_tools
import re
from dask.distributed import Future

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

import time


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
    
    start_time = time.time()
    fixed=clusterer(bc,threshold=1)
    end_time = time.time()
    print(f"Clustering execution time: {end_time - start_time:.6f} seconds")

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

def undo_one_hot_encoding(df, reference_label="reference"):
    """Should work for Dask and pandas DataFrames. Converts one-hot back to categorical.
    Rows with all zeros (intercept-only rows) are labeled with `reference_label`.
    """
    df = df.copy()
    pattern = re.compile(r'^C\((.+?)\)\[(.+?)\]$')
    one_hot_groups = {}

    # Step 1: Group one-hot columns by prefix
    for col in df.columns:
        match = pattern.match(col)
        if match:
            var, level = match.groups()
            one_hot_groups.setdefault(var, []).append((col, level))

    # Step 2: Undo one-hot encoding per group
    for var, cols_levels in one_hot_groups.items():
        cols, levels = zip(*cols_levels)
        subdf = df[list(cols)]

        row_sums = subdf.sum(axis=1)

        # Check for multi-hot rows (invalid)
        if (row_sums > 1).any():
            raise ValueError(
                f"Multi-hot encoding detected in variable '{var}'. Each row must have only one active level or be intercept-only."
            )

        # Determine the level for each row
        inferred = subdf.idxmax(axis=1)
        level_map = {col: level for col, level in cols_levels}
        result = inferred.map(level_map)

        # Assign reference label to rows where all one-hot cols are zero
        result[row_sums == 0] = reference_label

        df[var] = result
        df = df.drop(columns=list(cols))

    return df



def find_dask_future_paths(obj, seen=None, path=""):
    """
    Recursively find paths to all Dask Future objects in a structure.
    Doesn't explore futures themselves, so check those yourself.

    Returns a list of string paths like ['foo.bar', "baz[0]['x']"]
    """
    if seen is None:
        seen = set()

    results = []

    obj_id = id(obj)
    if obj_id in seen:
        return results
    seen.add(obj_id)

    if isinstance(obj, Future):
        results.append(path or "<root>")
        return results

    if isinstance(obj, dict):
        for k, v in obj.items():
            key_str = repr(k)
            subpath = f"{path}[{key_str}]" if path else f"[{key_str}]"
            results.extend(find_dask_future_paths(v, seen, subpath))
        return results

    if isinstance(obj, (list, tuple, set)):
        for i, v in enumerate(obj):
            subpath = f"{path}[{i}]" if path else f"[{i}]"
            results.extend(find_dask_future_paths(v, seen, subpath))
        return results

    if hasattr(obj, '__dict__'):
        for attr, val in vars(obj).items():
            subpath = f"{path}.{attr}" if path else attr
            results.extend(find_dask_future_paths(val, seen, subpath))
        return results

    return results

## tools for easy generation of hypotheses
@unimplemented
def one_versus_all():
    """
    Provide one negative control, and a list of others to compare against, and 
    this function will generate a hypothesis list comparing all vs it...
    Useful for a quick "which elements are expressed". 
    """
    pass
