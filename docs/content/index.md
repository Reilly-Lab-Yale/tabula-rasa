<a id="scMPRAforge.utils"></a>

# scMPRAforge.utils

<a id="scMPRAforge.utils.unimplemented_functions"></a>

#### unimplemented\_functions

List to track unimplemented functions

<a id="scMPRAforge.utils.unimplemented"></a>

#### unimplemented

```python
def unimplemented(func)
```

Decorator to mark functions as unimplemented.
Adds them to a global tracking list and breaks when called.

<a id="scMPRAforge.utils.list_unimplemented"></a>

#### list\_unimplemented

```python
def list_unimplemented()
```

Returns the list of unimplemented functions.

<a id="scMPRAforge.utils.bcs_to_lut"></a>

#### bcs\_to\_lut

```python
def bcs_to_lut(bc, threshold=1, encoding="utf-8", *args, **kwargs)
```

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

<a id="scMPRAforge.utils.undo_one_hot_encoding"></a>

#### undo\_one\_hot\_encoding

```python
def undo_one_hot_encoding(df, reference_label="reference")
```

Should work for Dask and pandas DataFrames, sparse & dense. Converts one-hot back to categorical.
Rows with all zeros (intercept-only rows) are labeled with `reference_label`.

<a id="scMPRAforge.utils.find_dask_future_paths"></a>

#### find\_dask\_future\_paths

```python
def find_dask_future_paths(obj, seen=None, path="")
```

Recursively find paths to all Dask Future objects in a structure.
Doesn't explore futures themselves, so check those yourself.

Returns a list of string paths like ['foo.bar', "baz[0]['x']"]

<a id="scMPRAforge.utils.dict_wrap"></a>

#### dict\_wrap

```python
def dict_wrap(client, dic)
```

Takes a dictionary and wraps all of its values in dask futures using the provided client.

<a id="scMPRAforge.utils.dict_unwrap"></a>

#### dict\_unwrap

```python
def dict_unwrap(dic)
```

Takes a dict of dask futures & gets their results.
Should replace with client.gather!
Note: Will hang if computations not done
Note: Not recursive!
Note: Since it pulls all the data to the control process, can take a lot of memory!

<a id="scMPRAforge.utils.find_treatment_column"></a>

#### find\_treatment\_column

```python
def find_treatment_column(xmu_names: list[str], factor: str,
                          level: str) -> str | None
```

Find the design column name for a treatment-coded factor/level, being tolerant to
the presence/absence of 'contr.treatment(...)' in the name.

Returns the matching column name or None if not found.

<a id="scMPRAforge.utils.simulate_library"></a>

#### simulate\_library

```python
def simulate_library(CREs, library_model)
```

In the event that you do not already have an MPRA library cloned,
This function takes a np string array of CRE names and a library model from a bounds object
and produces a table mapping each CRE to a set of random 20-mer MPRA barcodes.

<a id="scMPRAforge.utils.sample_from_library"></a>

#### sample\_from\_library

```python
def sample_from_library(library, size)
```

Takes an MPRA library in standard form and samples "size" rows. 
Uses abundance to sample using inverse transform sampling. 
TODO: add table type check.

<a id="scMPRAforge.utils.one_versus_all"></a>

#### one\_versus\_all

```python
def one_versus_all(comparisons,
                   *,
                   comparison_on: str,
                   reference_CRE: str | None = None,
                   reference_cell_type: str | None = None,
                   meta: str | None = None) -> pd.DataFrame
```

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

<a id="scMPRAforge.utils.alpha_for_expected_groups"></a>

#### alpha\_for\_expected\_groups

```python
def alpha_for_expected_groups(n, K_target)
```

Helper for sample_crp_groups
Choose alpha so E[K_n] ~= K_target using H_n ≈ log n + gamma.
Technically only valid for large n, but good enough for our purposes

<a id="scMPRAforge.utils.sample_crp_groups"></a>

#### sample\_crp\_groups

```python
def sample_crp_groups(n, alpha, rng=None)
```

Chinese Restaurant Process partition of n items.
Returns: np.array of length n with group ids in 0..K-1.

<a id="scMPRAforge.presets.presets"></a>

# scMPRAforge.presets.presets

<a id="scMPRAforge.presets"></a>

# scMPRAforge.presets

<a id="scMPRAforge.core"></a>

# scMPRAforge.core

<a id="scMPRAforge.core.always_unfinished"></a>

#### always\_unfinished

```python
@unimplemented
def always_unfinished()
```

tests unimplemented decorator.

<a id="scMPRAforge.core.table_type"></a>

#### table\_type

```python
def table_type(column_names)
```

**Arguments**:

  column_names <pd.Index>

**Returns**:

  <str>
  
  Returns the putative table type, one of:
  - mpra_umiwise
  - mpra_readwise
  - hypotheses
  - results
  - malformed
  
  Is kind with respect to extra columns & optional columns.
  
  (We could extend to type-checking as well, but that seems a tad draconian / unpythonic.)

<a id="scMPRAforge.core.skew_spread"></a>

#### skew\_spread

```python
@unimplemented
def skew_spread()
```

Creates a ground-truth dataframe of an scMPRA experiment
that is meant to test skew
(see readme for ground truth dataframe specification)

<a id="scMPRAforge.core.recombinator"></a>

#### recombinator

```python
def recombinator(primary, secondary)
```

All pairs of (All pairs of primary), secondary.
two duplicate `secondary` entries in each element.

<a id="scMPRAforge.core.activity_spread"></a>

#### activity\_spread

```python
def activity_spread(cell_types: List[str], minimum: float, maximum: float,
                    minp_value: float, total: int, frac_active: int,
                    ct_specificity: float)
```

Creates a ground-truth dataframe of an scMPRA experiment
with a controllable number of active CREs.
(see readme for ground truth dataframe specification)
Assumes experiment is interested in activity vs a known negative control,
not skew. 

This is useful to create datasets with a balance 
of active and inactive CREs which roughly
approx. real libraries. 

`total` is the total number of CREs to create.
`frac_active` is the fraction of the library that is active elements
`ct_specificity` is how much cell-type specificity there is. Its 
how many different values a given CRE will take across different cell-types. So 
specificity=2 implies that, on average, there will be two different 
means for each CRE across cell-types.

<a id="scMPRAforge.core.simple_spread"></a>

#### simple\_spread

```python
def simple_spread(cell_types: List[str],
                  min: float,
                  max: float,
                  fineness: int = 10,
                  hypothesis_type: str = "cartesian")
```

TODO: remove hypothesis_type & make obligate cartesian.
TODO: extract cartesian code to its own function.
Create a ground truth dataframe tiling all cell-types.
with synthetic CREs at a variety of strengths.

Returns a tuple of (ground truth, hypothesis object) 

Useful for simulation and power calculations.
(see readme for ground truth dataframe specification)

min, max are the min & max MPRA UMI / cell values.

<a id="scMPRAforge.core.load_hypothesis_set"></a>

#### load\_hypothesis\_set

```python
@unimplemented
def load_hypothesis_set(filepath)
```

Arguments
    filepath <str>
Returns
    <pd.DataFrame>

Loads a hypothesis or hypothesis+results set from disc.

<a id="scMPRAforge.core.simple_count"></a>

## simple\_count Objects

```python
class simple_count()
```

This class stores information pertaining to a simple negative binomial model.
It is low performance and NOT used for primary modeling of RNA-sequencing data.
Instead, it us used for small, discrete modeling tasks whch need flexibility but not 
performance.

<a id="scMPRAforge.core.simple_count.from_data"></a>

#### from\_data

```python
def from_data(data)
```

Initializes the object from count (not frequency) data.
Computes poisson and nb using statsmodels.

<a id="scMPRAforge.core.simple_count.draw_nb"></a>

#### draw\_nb

```python
def draw_nb(size)
```

returns a 1d numpy vector of draws from the nb
model of the object.

<a id="scMPRAforge.core.simple_count.plot"></a>

#### plot

```python
def plot(max_bins: int = 25, binwidth: int | None = None)
```

Plot a histogram of the data with *integer-aligned coarse bins* and overlay
NB and Poisson fits aggregated to the same bins.

Parameters
----------
max_bins : int
    Target maximum number of bars shown (used only if `binwidth` is None).
    The method will choose an integer bin width w >= 1 so that the span of
    the data is displayed in ~max_bins bars.
binwidth : int | None
    If provided, forces that integer bin width (w). If None, a width is
    computed from `max_bins`.

<a id="scMPRAforge.core.simple_count.to_dataframe"></a>

#### to\_dataframe

```python
def to_dataframe() -> pd.DataFrame
```

Represent the object as a single-row DataFrame.
All non-hidden attributes are included.

<a id="scMPRAforge.core.simple_count.from_dataframe"></a>

#### from\_dataframe

```python
@classmethod
def from_dataframe(cls, df: pd.DataFrame) -> "simple_count"
```

Recreate an object from a single-row DataFrame.

<a id="scMPRAforge.core.Bounds"></a>

## Bounds Objects

```python
@dataclass()
class Bounds()
```

Describes the boundaries of an experiment.
Specifically the UMI portion.

Bounds objects just average zero inflation parameters, as there
is no reason to complicate simulations with multiple values.
(If multiple ZI are expected for some reason, just adjust geometrically
or run multiple simulation batches.)

by_cell_type_theta and by_cell_type_zi are probably 
more accurate for downstream simulation than their 
by_cre counterparts, assuming that there are more 
cres than cell types, which should be the case for most
datasets.

<a id="scMPRAforge.core.Bounds.from_ortho"></a>

#### from\_ortho

```python
@classmethod
def from_ortho(cls, inp, preferred="by_cell_type")
```

Takes an ortho object and abstracts out its bounds.

This is an aggregation function and will hang if ortho is not done fitting yet. 

This function requires that the ortho still have its training data 
so we can extract things like "number of cells per cell-type" and 
"number of MPRA barcodes per cell".

Note that this function is totally replicate-agnostic. It averages estimated zero 
inflation across replicates.

<a id="scMPRAforge.core.scMPRA_data"></a>

## scMPRA\_data Objects

```python
class scMPRA_data()
```

Wrapper around a pandas dataframe of MPRA data. 
The primary purpose of the object is to record what operations have been performed on the data
(Pandas does not support metadata)

Could possibly replace with an anndata object.
Alternatively. also allow pass-through of pandas operations & record them... 
Alternatively, just implement a couple common operations (subsetting & friends) manually

<a id="scMPRAforge.core.scMPRA_data.flatten_overtransfection"></a>

#### flatten\_overtransfection

```python
def flatten_overtransfection()
```

If you have simulated a dataset, it will probably have some degree of 
overtransfection (same MPRA bc transfected into the same cell multiple times).
This function flattens such events, as they would be observed in a real dataset.

<a id="scMPRAforge.core.scMPRA_data.overtransfected"></a>

#### overtransfected

```python
def overtransfected(log=True, threshold_pct=WARN_MULTI_TRANSFECTION_PERCENT)
```

Return True iff the overall percent of cells with >=1 multi-transfection
(same mpra_bc observed >1 time in the same cell within a replicate)
exceeds `threshold_pct`. Logging is optional. Also flags overtransfection 
(or lack thereof) in metadata.

Uses a scale-free metric: (# cells with >=1 dup) / (total # cells) * 100

<a id="scMPRAforge.core.scMPRA_data.describe_transfection"></a>

#### describe\_transfection

```python
def describe_transfection()
```

Returns a simple_count object describing the number of transfections per cell, 
as proxied by number of unique MPRA barcodes per cell.

Returns a simple_count object
Drawing from one of the distributions (or a similar distribution shifted to a different MOI) 
can be used to simulate transfection.

<a id="scMPRAforge.core.scMPRA_data.describe_library"></a>

#### describe\_library

```python
def describe_library()
```

Returns a simple_count object describing the number of unique MPRA
barcodes for each

<a id="scMPRAforge.core.scMPRA_data.set_negative_controls"></a>

#### set\_negative\_controls

```python
def set_negative_controls(negative_controls: list[str])
```

Takes a list of CRE names that we consider to be negative controls and give them all the name "negative_control", lumping all their data together.

<a id="scMPRAforge.core.scMPRA_data.copy"></a>

#### copy

```python
def copy(exclude=())
```

Return a deepcopy of the object, optionally excluding fields.

<a id="scMPRAforge.core.scMPRA_data.from_tsv"></a>

#### from\_tsv

```python
@classmethod
def from_tsv(cls, filepath)
```

Returns a <scMPRA_data> object with data loaded from `filepath`.

<a id="scMPRAforge.core.scMPRA_data.from_parquet"></a>

#### from\_parquet

```python
@classmethod
def from_parquet(cls, path)
```

Returns a <scMPRA_data> object with data loaded from `path`.
Takes full path, /path/to/data.scmpra.

<a id="scMPRAforge.core.scMPRA_data.to_parquet"></a>

#### to\_parquet

```python
def to_parquet(path: str)
```

Saves to a parquet file using gzip compression.
Takes full path, /path/to/data.scmpra
WILL clobber existing files with the same path.

<a id="scMPRAforge.core.scMPRA_data.graph_chimeric"></a>

#### graph\_chimeric

```python
def graph_chimeric(*args, **kwargs)
```

TODO: test again now that its moved to scMPRA data obj

Arguments
    self
    *args
    **kwargs

Takes `scmpra_data`, a pandas dataframe of read-wise MPRA data (see docs) 
and plots a histogram of frequency of reads per UMI using seaborn.histplot. 

All other arguments are passed to the histplot call to allow graph 
customization. Particular useful are `bins`, `binrange`, and `log_scale`

<a id="scMPRAforge.core.scMPRA_data.read_wise_to_umi_wise"></a>

#### read\_wise\_to\_umi\_wise

```python
def read_wise_to_umi_wise(keep_reads=False)
```

Converts read-wise to UMI-wise (see readme for spec).

TODO: test again now that its moved to scMPRA data obj

<a id="scMPRAforge.core.scMPRA_data.cut_chimeric_reads"></a>

#### cut\_chimeric\_reads

```python
def cut_chimeric_reads(threshold)
```

Arguments
    self
    threshold : <int>

subsets to those UMIs which lie ABOVE the number-of-reads threshold, 
removing chimeric reads.

<a id="scMPRAforge.core.scMPRA_data.ortho_filter"></a>

#### ortho\_filter

```python
def ortho_filter()
```

Removes combinations of cre_id, cell_type which have less than MIN_PTS non-zero observations. 
This is much stricter than filter_low_umi_count

<a id="scMPRAforge.core.scMPRA_data.round_down_zeroes"></a>

#### round\_down\_zeroes

```python
@unimplemented
def round_down_zeroes()
```

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

<a id="scMPRAforge.core.flatten_barcode_errors"></a>

#### flatten\_barcode\_errors

```python
@unimplemented
def flatten_barcode_errors(df, barcode_column, *args, **kwargs)
```

Need to re-work to work with scMPRA data object
Arguments
    df <pandas.DataFrame>
    barcode_column <str>
Returns
    <pandas.DataFrame>

Uses umitools to flatten different barcodes which are likely only different
due to sequencing errors. Passes *args,**kwargs upstream to bcs_to_lut.

<a id="scMPRAforge.core.apply_deseq"></a>

#### apply\_deseq

```python
@unimplemented
def apply_deseq()
```

R quarantine zone.

<a id="scMPRAforge.core.abort_on_failure"></a>

#### abort\_on\_failure

```python
def abort_on_failure(future, client)
```

Call this function with a completed future that is strictly necessary for the task at hand.
If it didn't work, we'll crash.

<a id="scMPRAforge.core.standard_fit"></a>

#### standard\_fit

```python
def standard_fit(client, data, split)
```

Takes an scMPRA object and produces a set of models along one axis,
specified by split.

<a id="scMPRAforge.core.experiment_model"></a>

## experiment\_model Objects

```python
class experiment_model()
```

Contrains one full model of a dataset.
.model is a dict
- keys are whatever the levels of "split" are.
- values are futures of tensorzinb return dictionaries.

<a id="scMPRAforge.core.experiment_model.label_regressors"></a>

#### label\_regressors

```python
def label_regressors(client, design_matricies)
```

Takes design matricies used to generate the model &
modifies self in-place to have

<a id="scMPRAforge.core.experiment_model.flattened_copy"></a>

#### flattened\_copy

```python
def flattened_copy()
```

Makes a copy where the members are not futures but just objects

<a id="scMPRAforge.core.experiment_model.save"></a>

#### save

```python
def save(path)
```

Saves experimentmodel to a filepath.
Will hang if computation is not done yet

<a id="scMPRAforge.core.experiment_model.load"></a>

#### load

```python
@staticmethod
def load(client, path)
```

Loads an experimentmodel from a path
& returns it. Requires a client to wrap the individual models
in futures for use on a dask cluster.

<a id="scMPRAforge.core.ortho"></a>

## ortho Objects

```python
class ortho()
```

Stores multiple models of the same data
one set of by_cre models, and one set of by cell type models
Not to be used with multiple datasets.

<a id="scMPRAforge.core.ortho.save"></a>

#### save

```python
def save(path, name, strip_training_data=False)
```

Simple pickle save.

Will block & wait for results if not done computing

creates directory 'name' in 'path'

<a id="scMPRAforge.core.ortho.load"></a>

#### load

```python
@classmethod
def load(cls, client, path, name)
```

loads from a filepath, wrapping in futures on the provided cluster where appropriate

<a id="scMPRAforge.core.ortho.clean"></a>

#### clean

```python
@unimplemented
def clean(kill_list="auto")
```

Deletes intermediate values to save space. 

`kill_list` is any or all of "training_data", "design_matricies", "models", "parameters"

alternatively, "auto" is equivalent to ["training_data", "design_matricies"]

<a id="scMPRAforge.core.ortho.criss_cross"></a>

#### criss\_cross

```python
def criss_cross(client, dat)
```

Makes by_cre and by_cell_type models.

Note: a little computationally intensive...
retain_metadata will keep some information 'dat' in self.training_data
The actual MPRA data will be stripped to save space, but metadata will be retained

<a id="scMPRAforge.core.ortho.annotate_models"></a>

#### annotate\_models

```python
def annotate_models(client)
```

Adds regressor names to each model

<a id="scMPRAforge.core.ortho.extract_params"></a>

#### extract\_params

```python
def extract_params(client)
```

Extracts parameters for all models in the object

<a id="scMPRAforge.core.ortho.compute_model_qc"></a>

#### compute\_model\_qc

```python
def compute_model_qc()
```

Will hang if model is not finished. 
returns None, sets self.by_cell_qc, self.by_cre_qc to dictionaries of 
QC information comparing nb params of each direction of the ortho 
to the data means.

Meant for debugging / manual inspection.

<a id="scMPRAforge.core.ortho.precompute_wald"></a>

#### precompute\_wald

```python
def precompute_wald(client: Client)
```

Compute and cache Wald precomputations (SEs, covariances, name maps) for
every by-cell-type model and every by-CRE model in this ortho.
Stores results in self.wald_precomp (as Futures); persists with save().

<a id="scMPRAforge.core.ortho.make_wald_eval_bundle"></a>

#### make\_wald\_eval\_bundle

```python
def make_wald_eval_bundle() -> dict
```

Build a small, pickle-friendly snapshot with everything the workers
need to evaluate Wald tests. No Dask Futures inside.

<a id="scMPRAforge.core.parameters"></a>

## parameters Objects

```python
class parameters()
```

Stores triples of parameters
- nb (negative binomial mean), zi (zero inflation), theta (dispersion parameter)

for 'broken_on' (by cell type, by cre, or whatever models)

<a id="scMPRAforge.core.parameters.flattened_copy"></a>

#### flattened\_copy

```python
def flattened_copy()
```

Makes a copy where the members are not futures but just objects.

<a id="scMPRAforge.core.parameters.save"></a>

#### save

```python
def save(path)
```

Saves parameters to a filepath.
Will hang if computation is not done yet

<a id="scMPRAforge.core.parameters.load"></a>

#### load

```python
@staticmethod
def load(client, path)
```

Loads parameters from a path
& returns it. Requires a client to wrap the individual models
in futures for use on a dask cluster.

<a id="scMPRAforge.core.cast_multiindex_to_str_inplace"></a>

#### cast\_multiindex\_to\_str\_inplace

```python
def cast_multiindex_to_str_inplace(df)
```

Convert all levels of a MultiIndex to strings, modifying df in-place.

<a id="scMPRAforge.core.describe_parameters"></a>

#### describe\_parameters

```python
def describe_parameters(parameters, dat, split)
```

Produces a convienient single-dataframe description of one split model. 
'parameters' is a parameters object
Requires that you pass original data as dat to compute cell-numbers
Leaves non-split columns as one-hot. Returns "split" column as str categorical

<a id="scMPRAforge.core.anti_split"></a>

#### anti\_split

```python
def anti_split(split)
```

Returns the opposite split for a given split.
If splits are extended beyond 'cre_id' and 'cell_type',
this function should be extended to handle those cases.
Specifcally, all but split should be returned.

<a id="scMPRAforge.core.get_cell_counts"></a>

#### get\_cell\_counts

```python
def get_cell_counts(client: Client, dat: pd.DataFrame, split: str)
```

Takes a dask client and a pandas DataFrame `dat` containing MPRA data.

<a id="scMPRAforge.core.auto_partition"></a>

#### auto\_partition

```python
def auto_partition(pdf, target_mb_per_partition=PARTITION_SIZE_MB)
```

Convert a pandas DataFrame to a Dask DataFrame with automatic partition sizing.

**Arguments**:

  - pdf: input pandas DataFrame
  - target_mb_per_partition: desired memory usage per partition (in megabytes)
  

**Returns**:

  - ddf: Dask DataFrame with chosen number of partitions
  
  Minimum of 2!

<a id="scMPRAforge.core.simulate_from_description"></a>

#### simulate\_from\_description

```python
def simulate_from_description(description)
```

Simulate from a description dataframe.
Assumes input is one transfection event per row
Removes ground-truth rows.

<a id="scMPRAforge.core.simulation_batch"></a>

## simulation\_batch Objects

```python
class simulation_batch()
```

DEPRECATED
Class which takes a single <ortho> object and simulates replicates.
Optionally, fits additional ortho objects to simulations & plots their paremeter spread
Useful for estimating variance of an experimental setup...

<a id="scMPRAforge.core.simulation_batch.describe_primordial"></a>

#### describe\_primordial

```python
def describe_primordial()
```

Generates and saves descriptions of the primordial which are necessary for subsequent simulation

<a id="scMPRAforge.core.simulation_batch.clear_simulations"></a>

#### clear\_simulations

```python
@unimplemented
def clear_simulations()
```

Removes simulated data. Does not remove models fit to simulated data. Useful for reducing object size

<a id="scMPRAforge.core.simulation_batch.simulate_many"></a>

#### simulate\_many

```python
def simulate_many(client, n)
```

simulates n replicates

todo: additional parallelism
todo: pick which set of models to create : by cre, by cell-type, or both
currently all both

<a id="scMPRAforge.core.simulation_batch.fit_to_simulations"></a>

#### fit\_to\_simulations

```python
def fit_to_simulations(client)
```

Fits ortho models to all simulated datasets.
Depends on simulate_many having been called.

<a id="scMPRAforge.core.simulation_batch.plot_theta_spread"></a>

#### plot\_theta\_spread

```python
def plot_theta_spread(split)
```

Plots the spread of the thetas of simulated data with primordial for reference.

<a id="scMPRAforge.core.simulation_batch.plot_zi_spread"></a>

#### plot\_zi\_spread

```python
def plot_zi_spread(split)
```

Plots the spread of the zi parameters of simulated data with primordial for reference.

<a id="scMPRAforge.core.simulation_batch.save"></a>

#### save

```python
def save(path, name)
```

Simple save of a full simulation batch.
creates directory 'name' in 'path' (and many subdirectories).

Caveats:
1. Does not keep a copy of training data in simulated orthos, since 
2. It does not save the parameters objects, since they would be redundant with the summary dataframes.
Instead, you can use 
- description_primordial_by_cre
- description_primordial_by_cell_type
to get the values of the primordial parameters, and
- _nbs
- _zi_cre
- _zi_cell
- _theta_cre
- _theta_cell
for the parameter values of the simulated replicates.

(or you can have them recomputed).

Will block & wait for results if not done computing

<a id="scMPRAforge.core.simulation_batch.load"></a>

#### load

```python
@classmethod
def load(cls, client, path, name)
```

Loads a batch saved with save_batch.
See the caveats noted in the docstring for save_batch.
Requires you pass a client to re-wrap futures.

<a id="scMPRAforge.core.versus_truth"></a>

#### versus\_truth

```python
def versus_truth(ground_truth_mu: pd.DataFrame, inp_ortho: ortho)
```

Function takes a dataframe of ground truth values for each CRE, cell-type combination and compares to estimated parameters.

Note that mean absolute percentage error is only reported for cases where the truth values is nonzero.

TODO: clean up duplicate code

<a id="scMPRAforge.core.WaldPrecompEntry"></a>

## WaldPrecompEntry Objects

```python
class WaldPrecompEntry()
```

Minimal payload needed to evaluate Wald tests quickly for a single fitted model.
All numpy; safe to pickle.

<a id="scMPRAforge.core.WaldPrecomp"></a>

## WaldPrecomp Objects

```python
class WaldPrecomp()
```

Mirrors the shape of `parameters`: dicts keyed by split level.
Values are WaldPrecompEntry objects (or Futures thereof when live on a cluster).

<a id="scMPRAforge.core.HypothesisSet"></a>

## HypothesisSet Objects

```python
class HypothesisSet()
```

A validated container for hypothesis rows.

Columns:
  - comparison_CRE (str) [T]
  - comparison_cell_type (str) [T]
  - reference_CRE (str or NA) [F]
  - reference_cell_type (str or NA) [F]
  - meta (str or NA) [F]

Rules:
  - If both reference_* are NA: interpret as "compare to zero" (activity vs 0).
        - maybe change this to use the references from the scMPRA data object? discuss what we want the appropriate default to be
  - If exactly one of reference_CRE / reference_cell_type is provided: malformed.

<a id="scMPRAforge.core.HypothesisSet.is_zero_reference"></a>

#### is\_zero\_reference

```python
def is_zero_reference() -> pd.Series
```

Row-wise boolean: True if comparing against implicit zero (both reference_* are NA).

<a id="scMPRAforge.core.HypothesisSet.add_meta"></a>

#### add\_meta

```python
def add_meta(series_like) -> None
```

Attach/overwrite a meta column (useful for painting plots).

<a id="scMPRAforge.core.make_by_celltype_hypotheses"></a>

#### make\_by\_celltype\_hypotheses

```python
def make_by_celltype_hypotheses(*,
                                comparison_cell_type: str,
                                counts: "scMPRA_data",
                                comparison_cres: "list[str] | str" = "all",
                                reference_cre: str | None = "reference",
                                meta: str | None = None) -> "HypothesisSet"
```

Build hypotheses that test many CREs within a single cell type
(CRE varies; cell_type fixed). This is the natural input for the
by-cell-type Wald test (CRE vs baseline CRE in that cell type).

**Examples**:

  hs = make_by_celltype_hypotheses(
  comparison_cell_type="NeuroectodermBrain",
  counts=shendure,
  comparison_cres="all",
  reference_cre="reference",   # your flattened minP/noP
  meta="emvar_screen")
  

**Notes**:

  - We set BOTH reference columns per the table spec:
  reference_CRE = `reference_cre`
  reference_cell_type = `comparison_cell_type`
  - Passing reference_cre=None will generate the "compare‐to‐zero" flavor,
  but the current Wald code ignores zero and interprets baseline from the model.

<a id="scMPRAforge.core.make_by_cre_hypotheses"></a>

#### make\_by\_cre\_hypotheses

```python
def make_by_cre_hypotheses(*,
                           comparison_cre: str,
                           counts: "scMPRA_data",
                           comparison_cell_types: "list[str] | str" = "all",
                           reference_cell_type: str | None = None,
                           meta: str | None = None) -> "HypothesisSet"
```

Build hypotheses that test many cell types for one CRE
(cell_type varies; CRE fixed). This is the natural input for the
by-CRE Wald test (cell_type vs baseline cell type for the same CRE).

**Examples**:

  hs = make_by_cre_hypotheses(
  comparison_cre="CRE123",
  counts=shendure,
  comparison_cell_types="all",
  reference_cell_type="Pluripotent",
  meta="cell_specificity")
  

**Notes**:

  - We set BOTH reference columns per the table spec:
  reference_CRE = `comparison_cre`
  reference_cell_type = provided (or inferred)
  - If `reference_cell_type` is not provided, we try:
  counts.reference_cell_type, then literal "reference" if present.

<a id="scMPRAforge.core.make_all_by_celltype_hypotheses"></a>

#### make\_all\_by\_celltype\_hypotheses

```python
def make_all_by_celltype_hypotheses(
        *,
        counts: "scMPRA_data",
        reference_cre: str | None = "reference",
        meta: str | None = None,
        include_cell_types: "list[str] | None" = None,
        exclude_cell_types: "list[str] | None" = None) -> "HypothesisSet"
```

Build hypotheses for *every* cell type in the dataset, testing all CREs
within each cell type against the provided baseline CRE (default: 'reference').

Effectively: concat over cell types of
    make_by_celltype_hypotheses(comparison_cell_type=ct, comparison_cres="all")

Parameters
----------
counts : scMPRA_data
    Dataset.
reference_cre : str or None
    Baseline CRE label used in the by-cell-type test (typically the flattened negative controls -> 'reference').
    If None, you get the "compare-to-zero" flavor, though the Wald code currently focuses on baseline contrasts.
meta : str or None
    Optional label propagated to the 'meta' column.
include_cell_types : list[str] or None
    If provided, restrict to this whitelist of cell types (user-facing labels OK; they’ll be normalized).
exclude_cell_types : list[str] or None
    If provided, drop these cell types (applied after include).

Returns
-------
HypothesisSet
    A validated hypothesis set spanning all requested cell types.

<a id="scMPRAforge.core.make_all_by_cre_hypotheses"></a>

#### make\_all\_by\_cre\_hypotheses

```python
def make_all_by_cre_hypotheses(
        *,
        counts: "scMPRA_data",
        reference_cell_type: str | None = None,
        meta: str | None = None,
        include_cres: "list[str] | None" = None,
        exclude_cres: "list[str] | None" = None,
        drop_reference_cre: bool = True) -> "HypothesisSet"
```

Build hypotheses for *every* CRE in the dataset, testing all cell types
within each CRE against a baseline cell type (defaults to the dataset’s
reference cell type if available, otherwise tries literal 'reference').

Effectively: concat over CREs of
    make_by_cre_hypotheses(comparison_cre=cre, comparison_cell_types='all')

Parameters
----------
counts : scMPRA_data
    Dataset.
reference_cell_type : str or None
    Baseline cell type label for within-CRE contrasts. If None, we try
    counts.reference_cell_type, else literal 'reference' if present. 
    (Same behavior as make_by_cre_hypotheses.)
meta : str or None
    Optional label propagated to the 'meta' column.
include_cres : list[str] or None
    If provided, restrict to this whitelist of CREs.
exclude_cres : list[str] or None
    If provided, drop these CREs (applied after include).
drop_reference_cre : bool
    If True (default), skip the collapsed negative-control CRE named 'reference'
    to avoid generating “reference vs baseline cell type” rows.

Returns
-------
HypothesisSet
    A validated hypothesis set spanning all requested CREs.

<a id="scMPRAforge.core.make_bootstrap_activity_hypotheses"></a>

#### make\_bootstrap\_activity\_hypotheses

```python
def make_bootstrap_activity_hypotheses(
        *,
        counts: "scMPRA_data",
        comparison_cres: "list[str] | str" = "all",
        controls: "list[str] | str | None" = None,
        meta: str | None = "bootstrap_activity") -> "HypothesisSet"
```

Build a hypothesis set for the bootstrap activity test:
  - ONE ROW PER CRE (not per cell type)
  - comparison_cell_type = "ALL" (sentinel; ignored by the test)
  - reference_cell_type = "ALL" (to satisfy the 'both present or both NA' rule)
  - reference_CRE carries the control label(s); the test will union all unique controls
    present in the HS when building its bundle.

Params
------
counts : scMPRA_data
    Your dataset.
comparison_cres : list[str] | "all"
    Which CREs to test. "all" = every CRE in counts (we’ll drop any that are also controls).
controls : list[str] | str | None
    Control CRE label(s). If None, we try to infer "reference" if present.
meta : str | None
    Optional meta label.

Returns
-------
HypothesisSet

<a id="scMPRAforge.core.coerce_bootstrap_activity_from_hs"></a>

#### coerce\_bootstrap\_activity\_from\_hs

```python
def coerce_bootstrap_activity_from_hs(hs: "HypothesisSet") -> "HypothesisSet"
```

If you already have a hypothesis set (e.g., from make_all_by_cre_hypotheses),
collapse it to the bootstrap-activity shape:
  - de-duplicate to one row per comparison_CRE
  - set comparison_cell_type = reference_cell_type = "ALL"
  - keep reference_CRE as-is (we’ll union controls later)

<a id="scMPRAforge.core.ResultSet"></a>

## ResultSet Objects

```python
class ResultSet(HypothesisSet)
```

Extends HypothesisSet with result columns:
  - test_type (str)     [T]
  - test_statistic (float) [T]
  - p_value (float)     [T]
  - fold_change (float) [T]
  - bh_p (float)        [T]
  - flattened (bool)    [T]

<a id="scMPRAforge.core._BootRepGroupCT"></a>

## \_BootRepGroupCT Objects

```python
@dataclass
class _BootRepGroupCT()
```

<a id="scMPRAforge.core._BootRepGroupCT.cell_type"></a>

#### cell\_type

per-row cell_type (string)

<a id="scMPRAforge.core._BootRepGroupCT.norm_umis"></a>

#### norm\_umis

per-row normalized_umis_mpra_bc (float)

<a id="scMPRAforge.core._BootRepGroupCT.idx_by_cre_ct"></a>

#### idx\_by\_cre\_ct

(cre, ct) -> row indices

<a id="scMPRAforge.core._BootRepGroupCT.idx_ctrl_by_ct"></a>

#### idx\_ctrl\_by\_ct

ct -> union of control rows

<a id="scMPRAforge.core._BootRepGroupCT.n_int_by_cre_ct"></a>

#### n\_int\_by\_cre\_ct

observed `integrations`

<a id="scMPRAforge.core._BootBundleCT"></a>

## \_BootBundleCT Objects

```python
@dataclass
class _BootBundleCT()
```

<a id="scMPRAforge.core._BootBundleCT.n_int_strategy"></a>

#### n\_int\_strategy

"as_observed" | "median_non_reference"

<a id="scMPRAforge.core.HypothesisTester"></a>

## HypothesisTester Objects

```python
class HypothesisTester()
```

Orchestrates running a test function on each hypothesis row.
You supply `test_fn` that implements a single-row comparison and returns:
    dict(test_type, test_statistic, p_value, fold_change, flattened)
The runner adds BH (`bh_p`) and merges back with the hypothesis columns to return a ResultSet.

<a id="scMPRAforge.core.load_hypothesis_set"></a>

#### load\_hypothesis\_set

```python
def load_hypothesis_set(filepath)
```

Loads a hypothesis (or result) table from disk and returns a HypothesisSet or ResultSet.
Supports .tsv/.csv/.parquet by extension.

<a id="scMPRAforge.core.hypothesis_tester"></a>

#### hypothesis\_tester

```python
def hypothesis_tester(scmpra_models_or_data,
                      hypotheses: HypothesisSet,
                      flavor="wald",
                      test_fn=None)
```

Backward-compatible facade.

Provide either:
  - test_fn: a per-row callable used by HypothesisTester (preferred while migrating)
  - or later, we can route based on `flavor` to built-in tests.

Returns a ResultSet.

<a id="scMPRAforge.core.de_novo_simulation"></a>

## de\_novo\_simulation Objects

```python
class de_novo_simulation()
```

Class for simulating datasets anew.

<a id="scMPRAforge.core.de_novo_simulation.__init__"></a>

#### \_\_init\_\_

```python
def __init__(simulation_replicates: int,
             experiment_bounds: Bounds,
             ground_truth: pd.DataFrame,
             library: pd.DataFrame,
             negative_controls: list[str] = ["reference"],
             reference_cell_type: str = "reference")
```

See readme for formatting of ground_truth & library dataframes

<a id="scMPRAforge.core.de_novo_simulation.gamut"></a>

#### gamut

```python
def gamut(client)
```

Run the full simulation: transfection->transcription->realization

Probably the only method you will ever need to call on de_novo_simulation

<a id="scMPRAforge.core.de_novo_simulation.save"></a>

#### save

```python
def save(path, name)
```

Note that this function saves self.simulated_scMPRA
but not self.simulated, since the latter is intermediate.

<a id="scMPRAforge.core.de_novo_simulation.crosstab"></a>

#### crosstab

```python
def crosstab(test, index)
```

COLLECTOR FUNCTION
Evaluates the performance of a test as a classifier, using 
BH corrected p-value cutoff.

<a id="scMPRAforge.core.de_novo_simulation.prc"></a>

#### prc

```python
def prc(test, index)
```

COLLECTOR

<a id="scMPRAforge.core.de_novo_simulation.roc"></a>

#### roc

```python
def roc(test, index)
```

COLLECTOR

<a id="scMPRAforge.core.volcano"></a>

#### volcano

```python
def volcano(results: "ResultSet", title=None, bh_thresh=0.05, fc_thresh=1.0)
```

Volcano plot using BH-corrected p-values (bh_p) versus log2 fold change.

Parameters
----------
results : ResultSet
    Must contain columns: 'fold_change', 'bh_p'.
bh_thresh : float
    FDR threshold for significance (default 0.05).
fc_thresh : float
    Absolute log2 fold change threshold for vertical lines.

