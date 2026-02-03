#utility functions

#external imports
import functools
import logging
import umi_tools
import re
import numpy as np
import pandas as pd
from dask.distributed import Future
import itertools
import time
import seaborn as sns
import matplotlib.pyplot as plt

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

from pathlib import Path
def chkdir(path):
    p = Path(path)
    if not p.is_dir():
        raise FileNotFoundError(f"Directory not found: {p}")

def chkfile(path):
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"File not found: {p}")

def zero_pad_ground_truth(gt_df):
    """
    Takes a ground truth dataframe (see specification in readme) 
    and fills in missing `cell_type`, `cre_id` combinations with 
    `true_mean`=zero.
    """
    #make a df of all zeroes
    #all combos of cell_type and cre_id
    zeroes = itertools.product(gt_df["cell_type"].unique(),gt_df["cre_id"].unique())
    #cast to df
    zeroes = pd.DataFrame([i for i in zeroes],columns=["cell_type","cre_id"])
    #get those combos not present in the gt df
    missing_combos=gt_df.merge(zeroes,how="outer",on=["cell_type","cre_id"],indicator=True)
    missing_combos=missing_combos[missing_combos["_merge"]=="right_only"]
    missing_combos=missing_combos[["cell_type","cre_id"]]
    #init zero true mean
    missing_combos["true_mean"]=0.0
    #stack
    return pd.concat([gt_df,missing_combos])

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
    """
    Should work for Dask and pandas DataFrames, sparse & dense. Converts one-hot back to categorical.
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

        # Sparse-friendly row sum
        row_sums = subdf.sum(axis=1)

        # Check for multi-hot rows - convert to dense only for validation
        row_sums_check = row_sums.values if hasattr(row_sums, 'sparse') else row_sums
        if (row_sums_check > 1).any():
            raise ValueError(
                f"Multi-hot encoding detected in variable '{var}'. Each row must have only one active level or be intercept-only."
            )

        # Initialize result with reference label
        result = pd.Series([reference_label] * len(subdf), index=subdf.index)
        
        # Iterate through columns and assign levels where active
        # This avoids idxmax which can be expensive for sparse data
        for col, level in cols_levels:
            active_mask = subdf[col] > 0  # Works efficiently with sparse data
            result[active_mask] = level

        df[var] = pd.Categorical(result)
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
    If object is none, simply returns none
    """
    if dic is None:
        return None
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

import math
def alpha_for_expected_groups(n, K_target):
    """
    Helper for sample_crp_groups
    Choose alpha so E[K_n] ~= K_target using H_n ≈ log n + gamma.
    Technically only valid for large n, but good enough for our purposes
    """
    gamma = 0.5772156649015329
    return max(1e-12, K_target / (math.log(n) + gamma))

def sample_crp_groups(n, alpha, rng=None):
    """
    Chinese Restaurant Process partition of n items.
    Returns: np.array of length n with group ids in 0..K-1.
    """
    if rng is None:
        rng = np.random.default_rng()
    groups = np.full(n, -1, dtype=int)
    # Track current group sizes
    sizes = []  # list of counts per existing group
    for i in range(n):
        # Prob of joining existing group k is sizes[k] / (alpha + i)
        # Prob of creating new group is alpha / (alpha + i)
        total = alpha + i
        if len(sizes) == 0 or rng.random() < alpha / total:
            # new group
            sizes.append(1)
            groups[i] = len(sizes) - 1
        else:
            # join existing: pick proportional to sizes
            k = rng.choice(len(sizes), p=np.array(sizes) / (total - alpha))
            sizes[k] += 1
            groups[i] = k
    return groups


def _plot_test_bars(
    df,
    metric="auroc",
    *,
    test_col="test",
    replicate_col="replicate",
    figsize=(8, 4.5),
    estimator="mean",     # "mean" or "median"
    err="sd",             # "sd" or "se"
    show_points=True,
    point_jitter=0.18,
    point_size=4,
    order="median_desc",  # "median_desc", "mean_desc", or a list of test names
    ylims="auto",         # "auto" or (ymin, ymax)
    ax=None,
    title=None,
):
    """
    Helper function for de_novo_sim.

    Bar chart per test with replicate points overlaid and error bars.

    - Bars show mean/median across replicates.
    - Error bars show SD or SE across replicates.
    - Y-axis is fixed to a reasonable range per metric (AUROC/AUPRC) by default.
    """
    required = {test_col, replicate_col, metric}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    plot_df = df[[test_col, replicate_col, metric]].copy()
    plot_df[metric] = plot_df[metric].astype(float)

    sns.set_theme(style="whitegrid", context="talk")

    # Choose aggregation
    if estimator not in {"mean", "median"}:
        raise ValueError("estimator must be 'mean' or 'median'")
    if err not in {"sd", "se"}:
        raise ValueError("err must be 'sd' or 'se'")

    # Decide order
    tests = plot_df[test_col].unique().tolist()
    if isinstance(order, list):
        order_list = order
    else:
        if order == "median_desc":
            order_list = (
                plot_df.groupby(test_col)[metric].median().sort_values(ascending=False).index.tolist()
            )
        elif order == "mean_desc":
            order_list = (
                plot_df.groupby(test_col)[metric].mean().sort_values(ascending=False).index.tolist()
            )
        else:
            order_list = sorted(tests)

    # Fixed y-lims defaults that are usually sensible for these metrics
    if ylims == "auto":
        m = metric.lower()
        if m == "auroc":
            y_min, y_max = 0.5, 1.0
        elif m == "auprc":
            y_min, y_max = 0.0, 1.0
        else:
            # generic: pad around observed
            lo, hi = plot_df[metric].min(), plot_df[metric].max()
            pad = 0.05 * (hi - lo if hi > lo else 1.0)
            y_min, y_max = lo - pad, hi + pad
    else:
        y_min, y_max = ylims

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # Build barplot with desired estimator + error
    agg = "mean" if estimator == "mean" else "median"
    # seaborn expects callable for estimator
    est_fn = {"mean": lambda x: x.mean(), "median": lambda x: x.median()}[agg]

    # errorbars:
    # seaborn 0.12+: errorbar can be "sd" or "se"
    sns.barplot(
        data=plot_df,
        x=test_col,
        y=metric,
        order=order_list,
        estimator=est_fn,
        errorbar=err,     # "sd" or "se"
        capsize=0.15,
        ax=ax,
    )

    if show_points:
        sns.stripplot(
            data=plot_df,
            x=test_col,
            y=metric,
            order=order_list,
            jitter=point_jitter,
            size=point_size,
            alpha=0.8,
            ax=ax,
        )

    ax.set_xlabel("")
    ax.set_ylabel(metric.upper())
    ax.set_ylim(y_min, y_max)
    ax.set_title(title or f"{metric.upper()} by test ({estimator} ± {err})")
    ax.tick_params(axis="x", rotation=25)

    fig.tight_layout()
    return fig, ax

def one_library_replicate(root,min,max,reps,client,flatten_overtransfection,bound):
    """
    Notebook helper function.
    Creates a de_novo_simulation in root, with a random name.
    Assumes corresponding libraries. 
    """
    #create ground truth dataframe
    rng = np.random.default_rng()
    cre_gt=rng.uniform(min,max,size=n_cres-1)
    cre_gt=np.append(cre_gt,minP)
    names=[f"synthcre_{i}" for i in range(0,n_cres-1)]+["reference"]

    gt_df=pd.DataFrame({"cre_id":names,"mu":cre_gt})
    gt_df["cell_type"]="reference"

    # simulate libraries
    libraries=[scm.simulate_library(CREs=gt_df["cre_id"],
                 library_model=scm.SHENDURE_BOUNDS.library_model)
                 for i in range(n_sims)]
    
    #initalize the simulated replicate
    name=uuid.uuid4().hex[:8]
    
    sim=scm.de_novo_simulation(location=root,
                            name=f"sim_{name}",
                            client=client,
                            libraries=libraries,
                            library_mapping="corresponding",
                            flatten_overtransfection=False,
                            n_sims=n_sims,
                            experiment_bounds=bound,
                            ground_truth=gt_df)
    sim.gamut()
    
    return name, sim

def pow_curve(df,n_bins=100):
    """
    Takes a df of bool reject_null and fc (fold change)
    such as is produced by `sum_pow` and makes a power plot.
    You can choose the number of bins with `n_bins`
    """
    # bin comp_mean
    df = df.copy()
    df["fc"] = pd.cut(df["fc"], bins=n_bins)

    # aggregate: fraction of True in each bin
    binned = (
        df.groupby("fc", observed=True)["reject_null"]
        .mean()
        .reset_index(name="reject_frac")
    )

    # bin centers for plotting
    binned["bin_center"] = binned["fc"].apply(lambda x: x.mid)

    # plot
    sns.lineplot(
        data=binned,
        x="bin_center",
        y="reject_frac",
        marker="o"
    )

    plt.xlabel("fc (binned)")
    plt.ylabel("Power")
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.show()