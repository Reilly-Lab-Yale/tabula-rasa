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
  
- `TODO` - move to the inside of the scMPRA object

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

<a id="scMPRAforge.core.scMPRA_data.from_json"></a>

#### from\_json

```python
@unimplemented
@classmethod
def from_json(cls, filepath)
```

Returns a <scMPRA_data> object with data loaded from `filepath`.

<a id="scMPRAforge.core.scMPRA_data.to_json"></a>

#### to\_json

```python
@unimplemented
def to_json(filepath)
```

Dump object to filepath

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

<a id="scMPRAforge.core.hypothesis_tester"></a>

#### hypothesis\_tester

```python
@unimplemented
def hypothesis_tester(scmpra_models, hypotheses, flavor="wald")
```

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

<a id="scMPRAforge.core.abort_on_failure"></a>

#### abort\_on\_failure

```python
def abort_on_failure(future, client)
```

Call this function with a completed future that is strictly necessary for the task at hand.
If it didn't work, we'll crash.

<a id="scMPRAforge.core.nb_versus_means"></a>

#### nb\_versus\_means

```python
def nb_versus_means(params, design_matricies, scMPRAdat)
```

Takes a model & design matricies corresponding to one direction of an ortho
(A set of 'by_cre' or a 'by_cell-type') and the original training data and produces a QC dictionary
regressing data means against nb estimates.
Used for quality control.

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
def save(path, name)
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

<a id="scMPRAforge.core.describe_parameters"></a>

#### describe\_parameters

```python
def describe_parameters(client, parameters, dat, split)
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

Simulate from a description dask dataframe

<a id="scMPRAforge.core.simulation_batch"></a>

## simulation\_batch Objects

```python
class simulation_batch()
```

Class which takes a single <ortho> object and simulates replicates.
Optionally, fits additional ortho objects to simulations & plots their paremeter spread
Useful for estimating variance of an experimental setup...

<a id="scMPRAforge.core.simulation_batch.describe_primordial"></a>

#### describe\_primordial

```python
def describe_primordial(client, dat)
```

Generates and saves descriptions of the primordial which are necessary for subsequent simulation

<a id="scMPRAforge.core.simulation_batch.clear_simulations"></a>

#### clear\_simulations

```python
def clear_simulations()
```

Removes simulated data. Does not remove models fit to simulated data. Useful for reducing object size

<a id="scMPRAforge.core.simulation_batch.simulate_many"></a>

#### simulate\_many

```python
def simulate_many(client, n)
```

simulates n replicates

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

<a id="scMPRAforge.core.volcano"></a>

#### volcano

```python
@unimplemented
def volcano(results: experiment_model)
```

Volcano plot of p value versus log fold change

<a id="scMPRAforge.utils"></a>

# scMPRAforge.utils

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

Should work for Dask and pandas DataFrames. Converts one-hot back to categorical.
Rows with all zeros (intercept-only rows) are labeled with `reference_label`.

<a id="scMPRAforge.utils.find_dask_future_paths"></a>

#### find\_dask\_future\_paths

```python
def find_dask_future_paths(obj, seen=None, path="")
```

Recursively find paths to all Dask Future objects in a structure.
Doesn't explore futures themselves, so check those yourself.

Returns a list of string paths like ['foo.bar', "baz[0]['x']"]

<a id="scMPRAforge.utils.make_present_dict"></a>

#### make\_present\_dict

```python
def make_present_dict(futures)
```

Takes a dict of dask futures & gets their results.

Note: Will hang if computations not done
Note: Not recursive!
Note: Since it pulls all the data to the control process, can take a lot of memory!

<a id="scMPRAforge.utils.one_versus_all"></a>

#### one\_versus\_all

```python
@unimplemented
def one_versus_all()
```

Provide one negative control, and a list of others to compare against, and 
this function will generate a hypothesis list comparing all vs it...
Useful for a quick "which elements are expressed".

