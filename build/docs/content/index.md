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

<a id="scMPRAforge.core"></a>

# scMPRAforge.core

<a id="scMPRAforge.core.always_unfinished"></a>

#### always\_unfinished

```python
@unimplemented
def always_unfinished()
```

tests unimplemented decorator

<a id="scMPRAforge.core.load_scMPRA_data"></a>

#### load\_scMPRA\_data

```python
@unimplemented
def load_scMPRA_data(filepath, wise, collapsed)
```

(for development, assume read-wise, uncollapsed)
Takes

<a id="scMPRAforge.core.graph_qc_metrics"></a>

#### graph\_qc\_metrics

```python
@unimplemented
def graph_qc_metrics(scmpra_data, *args, **kwargs)
```

Takes `scmpra_data`, a pandas dataframe of read-wise MPRA data 
(see docs) and plots a histogram of frequency of reads per UMI
using seaborn.histplot. All other arguments are passed to
the histplot call to allow graph customization.

<a id="scMPRAforge.core.cut_chimeric_reads"></a>

#### cut\_chimeric\_reads

```python
@unimplemented
def cut_chimeric_reads(scmpra_data, threshold)
```

Takes `scmpra_data`, a pandas dataframe of read-wise scMPRA data 
(see docs) and `threshold` <int>.
Returns a pandas dataframe of umi-wise MPRA data (see docs)
subset to those UMIs which lie above the number-of-reads
threshold.

<a id="scMPRAforge.core.hypothesis_tester"></a>

#### hypothesis\_tester

```python
@unimplemented
def hypothesis_tester(scmpra_data, hypotheses, flavor="wald")
```

Takes scmpra_data, a pandas dataframe of umi-wise scMPRA data

<a id="scMPRAforge.core.fit"></a>

#### fit

```python
@unimplemented
def fit(scmpra_data, round_down_threshold=4)
```

requires collapsed-on-barcodes, umi-wise
impl. round-down for very low CREs...

