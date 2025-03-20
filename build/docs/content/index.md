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

<a id="scMPRAforge.utils.one_versus_all"></a>

#### one\_versus\_all

```python
@unimplemented
def one_versus_all()
```

Provide one negative control, and a list of others to compare against, and 
this function will generate a hypothesis list comparing all vs it...
Useful for a quick "which elements are expressed".

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
def read_wise_to_umi_wise(scmpra_data, keep_reads=False)
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

<a id="scMPRAforge.core.flatten_mpra_barcodes"></a>

#### flatten\_mpra\_barcodes

```python
def flatten_mpra_barcodes(scmpra_data)
```

Arguments
    scmpra_data : <pd.DataFrame> of umi-wise, full MPRA
returns
    <pd.DataFrame> of umi-wise, flattened MPRA

<a id="scMPRAforge.core.fit"></a>

#### fit

```python
@unimplemented
def fit(scmpra_data, round_down_threshold=4)
```

Arguments
    scmpra_data
Returns
    {
        model : <statsmodels.discrete.count_model.ZeroInflatedNegativeBinomialResultsWrapper>,
        flattened : <pd.DataFrame>
    }

<a id="scMPRAforge.core.volcano"></a>

#### volcano

```python
@unimplemented
def volcano(results)
```

Volcano plot of p value versus log fold change

