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
  
  (We could extend to type-checking as well, but that seems a tad draconian.)

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

<a id="scMPRAforge.core.load_scMPRA_data"></a>

#### load\_scMPRA\_data

```python
def load_scMPRA_data(filepath)
```

Arguments
    filepath <str>
Returns
    <pd.DataFrame>

Loads tsv scMPRA data from `filepath`.

<a id="scMPRAforge.core.graph_chimeric"></a>

#### graph\_chimeric

```python
def graph_chimeric(scmpra_data, *args, **kwargs)
```

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

<a id="scMPRAforge.core.cut_chimeric_reads"></a>

#### cut\_chimeric\_reads

```python
def cut_chimeric_reads(scmpra_data, threshold)
```

Arguments
    scmpra_data : <pandas.DataFrame> of read-wise scMPRA data 
    threshold : <int>

Returns
    <pandas.DataFrame> of read-wise MPRA data

subsets to those UMIs which lie ABOVE the number-of-reads threshold, 
removing chimeric reads.

<a id="scMPRAforge.core.read_wise_to_umi_wise"></a>

#### read\_wise\_to\_umi\_wise

```python
def read_wise_to_umi_wise(scmpra_data, keep_reads=False, bypass_consistency_check=False)
```

Arguments
    scmpra_data : <pandas.DataFrame> of read-wise scMPRA data
    keep_reads : <bool>
Returns
    <pandas.DataFrame> of umi-wise scMPRA data

Converts read-wise to UMI-wise table (see readme for spec).

<a id="scMPRAforge.core.flatten_barcode_errors"></a>

#### flatten\_barcode\_errors

```python
def flatten_barcode_errors(df, barcode_column, *args, **kwargs)
```

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

<a id="scMPRAforge.core.fit"></a>

#### fit

```python
def fit(client, data, nb_formula: str, zi_formula: str, broken_on: str, round_down_threshold: int = 4, dry: bool = False, return_design_matricies=False)
```

dry : dry run: don't actually fit anything, just return an experiment_model 
broken_on : "unified" for one model, or put the name of a column

<a id="scMPRAforge.core.experiment_model"></a>

## experiment\_model Objects

```python
class experiment_model()
```

uniq_predictor stores...

<a id="scMPRAforge.core.ortho"></a>

## ortho Objects

```python
class ortho()
```

Stores multiple models of the same data.
Not to be used with multiple datasets. 
stores hypothesis sets & coresp. models

<a id="scMPRAforge.core.ortho.criss_cross"></a>

#### criss\_cross

```python
def criss_cross(client, dat, retain_design_matricies=False)
```

Note: a little computationally intensive...

<a id="scMPRAforge.core.ortho.extract_params"></a>

#### extract\_params

```python
def extract_params(client)
```

Extracts parameters for all models in the object

<a id="scMPRAforge.core.model_to_parameters"></a>

#### model\_to\_parameters

```python
def model_to_parameters(model)
```

A function to extract model parameter triples from simple model.
Currently pretty bespoke: be careful with more complicated models...

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

Simulate from a description dask dataframe.

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
Depends on

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
def undo_one_hot_encoding(df)
```

should work for dask and pandas dataframes.

<a id="scMPRAforge.utils.one_versus_all"></a>

#### one\_versus\_all

```python
@unimplemented
def one_versus_all()
```

Provide one negative control, and a list of others to compare against, and 
this function will generate a hypothesis list comparing all vs it...
Useful for a quick "which elements are expressed".

