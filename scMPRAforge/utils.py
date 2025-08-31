#utility functions

#external imports
import functools
import logging
import umi_tools
import re
import numpy as np
import pandas as pd
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



def dict_wrap(client,dic):
    """
    Takes a dictionary and wraps all of its values in dask futures using the provided client.
    """
    ret={}
    for key in dic:
        ret[key]=client.submit(lambda x: x,dic[key])
    return ret

def dict_unwrap(dic):
    """
    Takes a dict of dask futures & gets their results. 
    Should replace with client.gather!
    Note: Will hang if computations not done
    Note: Not recursive! 
    Note: Since it pulls all the data to the control process, can take a lot of memory!
    """
    ret={}
    for key in dic:
        ret[key]=dic[key].result()
    return ret

def find_treatment_column(xmu_names: list[str], factor: str, level: str) -> str | None:
    """
    Find the design column name for a treatment-coded factor/level, being tolerant to
    the presence/absence of 'contr.treatment(...)' in the name.

    Returns the matching column name or None if not found.
    """
    # simplest form
    simple = f"C({factor})[T.{level}]"
    if simple in xmu_names:
        return simple

    # tolerant scan, e.g. "C(cell_type, contr.treatment(base='reference'))[T.Level]"
    suffix = f"))[T.{level}]"
    prefix = f"C({factor},"
    for nm in xmu_names:
        if nm.startswith(prefix) and nm.endswith(suffix):
            return nm
    return None

def generate_barcodes(length, count):
    if 4**length < count:
        raise ValueError("Not enough unique barcodes of given length.")

    digit_to_base = ['A', 'C', 'G', 'T']
    barcodes = []

    for i in range(count):
        barcode = []
        n = i
        for _ in range(length):
            barcode.append(digit_to_base[n % 4])
            n //= 4
        # If length not fully filled, pad with 'A's
        while len(barcode) < length:
            barcode.append('A')
        barcodes.append(''.join(reversed(barcode)))

    return barcodes

def simulate_library(CREs,library_model):
    """
    In the event that you do not already have an MPRA library cloned,
    This function takes a np string array of CRE names and a library model from a bounds object
    and produces a table mapping each CRE to a set of random 20-mer MPRA barcodes. 
    """
    CREs=CREs.unique()
    mpra_barcodes_per_CRE=library_model.draw_nb(len(CREs))

    
    ret=pd.DataFrame({'cre_id':CREs,'n_barcodes':mpra_barcodes_per_CRE})

    #repeat each row the number of times equal to the number of barcodes
    ret=ret.loc[ret.index.repeat(ret["n_barcodes"])].reset_index(drop=True)
    #no longer need barcode per cre count, drop it
    ret=ret.drop(columns=["n_barcodes"])
    #add barcodes
    ret["mpra_bc"]=generate_barcodes(length=20,count=len(ret))
    ret["abundance"] = np.abs(np.random.randn(len(ret)))
    ret["abundance"]=ret["abundance"]/sum(ret["abundance"])

    return ret

def sample_from_library(library,size):
    """
    Takes an MPRA library in standard form and samples "size" rows. 
    Uses abundance to sample using inverse transform sampling. 
    TODO: add table type check. 
    """
    assert sum(library["abundance"])-1.0<0.0001; "Abundance must sum to 1."
    
    library=library.reset_index(drop=True)
    library["cum_abundance"]=library["abundance"].cumsum()

    random_vector=np.random.rand(size)
    indices = np.searchsorted(library["cum_abundance"].values, random_vector)
    sample = library.iloc[indices].reset_index(drop=True)
    sample=sample.drop(columns=["cum_abundance"])

    return sample

## tools for easy generation of hypotheses


def one_versus_all(
    comparisons,
    *,
    comparison_on: str,
    reference_CRE: str | None = None,
    reference_cell_type: str | None = None,
    meta: str | None = None
) -> pd.DataFrame:
    """
    Provide one negative control, and a list of others to compare against, and 
    this function will generate a hypothesis list comparing all vs it...
    Useful for a quick "which elements are expressed". 

    Generate a hypothesis table comparing each member of `comparisons` against a shared reference.

    Parameters
    ----------
    comparisons : iterable of str
        The values to place in the *comparison* column specified by `comparison_on`.
        If `comparison_on == "cre"`, you must pass corresponding `comparison_cell_type` via tuples.
        If `comparison_on == "cell_type"`, you must pass corresponding `comparison_CRE` via tuples.
        To keep API simple, we accept:
            - comparison_on == "cre":  comparisons = [(cre_id, cell_type), ...]
            - comparison_on == "cell_type": comparisons = [(cell_type, cre_id), ...]

    reference_CRE : str or None
    reference_cell_type : str or None
        If both None -> implicit zero comparison (activity vs 0).
        If exactly one provided -> ValueError (malformed by spec).

    meta : str or None
        Optional label for 'meta' column.

    Returns
    -------
    pd.DataFrame with columns:
        comparison_CRE, comparison_cell_type, reference_CRE, reference_cell_type, meta
    """
    rows = []
    if (reference_CRE is None) ^ (reference_cell_type is None):
        raise ValueError("Provide both reference_CRE and reference_cell_type, or neither (zero reference).")

    if comparison_on not in {"cre", "cell_type"}:
        raise ValueError("comparison_on must be 'cre' or 'cell_type'.")

    for pair in comparisons:
        if comparison_on == "cre":
            try:
                cre, ct = pair
            except Exception:
                raise ValueError("For comparison_on='cre', pass comparisons as iterable of (cre_id, cell_type).")
            rows.append({
                "comparison_CRE": cre,
                "comparison_cell_type": ct,
                "reference_CRE": reference_CRE,
                "reference_cell_type": reference_cell_type,
                "meta": meta
            })
        else:
            try:
                ct, cre = pair
            except Exception:
                raise ValueError("For comparison_on='cell_type', pass comparisons as iterable of (cell_type, cre_id).")
            rows.append({
                "comparison_CRE": cre,
                "comparison_cell_type": ct,
                "reference_CRE": reference_CRE,
                "reference_cell_type": reference_cell_type,
                "meta": meta
            })

    return pd.DataFrame.from_records(rows, columns=[
        "comparison_CRE","comparison_cell_type","reference_CRE","reference_cell_type","meta"
    ])

