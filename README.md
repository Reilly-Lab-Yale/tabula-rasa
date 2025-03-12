(Low priority) reminder : change asserts to thrown errors. 

Assumptions:
1. Barcode flattening is mandatory. Though ideally we would model each barcode as per MPRAmodel for greater statistical power, the sparse nature of scMPRA data means that we won't be able to do this. 
2. Removal of "false zeroes" through clonotype analysis and transfection reporters is part of pre-processing. 
3. We are interested in changes in CRE activity within and between cell-types
4. I've run ahead with the "one big model" approach (instead of the "one model per CRE" approach) but we can change this later based on some exploratory analysis. 
5. We could potentially add MPRA barcode variability as another mixture-model component : however I expect this to be well-captured in a single set of per-CRE-cell-type ZINB paramater sets (should be tested). (Ought to be tested separately...)

(We may eventually wish to turn the ad-hoc pre-processing steps into common tools.)

Run `pydoc-markdown` in the repo root to update the docs. (Could automate this with a github action). 

# MPRA data formatting

All on-disc data should be tsv or (some high-performace format to be chosen later).

Note that most of our code does not *require* that any of the barcode sequences should be actual nucleotide sequences. So you can replace them with, for example, numerical combinatorial barcoding sub-barcode ID strings with no consequence. The exception is barcode-deduplication. 

However, given python's weak typing, I'd recommend not putting something that could be obviously misinterpreted as another type (e.g. replicate a single int as this could easially be misinterpreted to suggest an ordinal position).

Data can be umi-wise or read-wise.

(Collapsing cells, replicates, and MPRA barcodes is not meaningful for us, since this would collapse what we believe to be true biological samples. However we may frequently *summarize* (e.g. compute mean UMIs per cell, or model a distribution where all UMI count originating from one CRE in one cell-type (regardless of MPRA barcode) are considered to have come from one triplicate of zinb parameters).

In memory, dataframes use named columns with dummy row indicies. 

**Read-wise**
| Column name | Type             | Description                 | Mandatory? |
| ----------- | ---------------- | --------------------------- | ---------- |
| cell_bc     | str (nucleotide) | cell barcode                | T          |
| rep_id      | str              | replicate id                | T          |
| cre_id      | str              | CRE id or name              | T          |
| cell_type   | str              | cell-type                   | T          |
| mpra_bc     | str (nucleotide) | MPRA reporter barcode       | T          |
| umi         | str (nucleotide) | unique molecular identifier | T          |
| reads       | int              | number of reads             | T          |


**Umi-wise, full MPRA**
| Column name | Type             | Description                             | Mandatory? |
| ----------- | ---------------- | --------------------------------------- | ---------- |
| cell_bc     | str (nucleotide) | cell barcode                            | T          |
| rep_id      | str              | replicate id                            | T          |
| cre_id      | str              | CRE id or name                          | T          |
| cell_type   | str              | cell-type                               | T          |
| mpra_bc     | str (nucleotide) | MPRA reporter barcode                   | T          |
| umis        | int              | number of unique molecular identifiers  | T          |
| reads       | int              | number of reads, summed across all UMIs | F          |

**Umi-wise, flattened MPRA**
| Column name | Type             | Description                             | Mandatory? |
| ----------- | ---------------- | --------------------------------------- | ---------- |
| cell_bc     | str (nucleotide) | cell barcode                            | T          |
| rep_id      | str              | replicate id                            | T          |
| cre_id      | str              | CRE id or name                          | T          |
| cell_type   | str              | cell-type                               | T          |
| umis        | int              | number of unique molecular identifiers  | T          |
| mpra_bcs     | int | Numer of MPRA barcodes                   | F          |
| reads       | int              | number of reads, summed across all UMIs | F          |


Also note that all these strings are really factors / categorical data, and will be treated as such. 

 # Hypothesis formatting
 A hypothesis set is a table with the following format:

| Column name          | Type | Description    | Mandatory? |
| -------------------- | ---- | -------------- | ---------- |
| comparison_CRE       | str  | CRE id or name | T          |
| comparison_cell_type | str  | cell-type      | T          |
| reference_CRE        | str  | CRE id or name | F          |
| reference_cell_type  | str  | cell-type      | F          |
| meta                 | str  | metadata       | F          |

- If no reference CRE is provided, the package will assume that we are comparing to zero (looking for any activity at all). 
- A row that contains only one of reference_CRE, reference_cell_type is considered malformed. 
- For assessing variant effects, we recommend the convention that reference columns refer to the reference allele, and comparison columns refer to the alternate allele.
- The metadata column can also be used to paint plots, so a categorical like "negative_control", "positive_control", "emvar", "putative-brain-specific" or similar would work well. No strict requirements: modify names as suits your experimental design.

A result table is the same as a hypothesis table with the following additional columns:

| Column name    | Type  | Description                                            | Mandatory? |
| -------------- | ----- | ------------------------------------------------------ | ---------- |
| test_type      | str   | which test was performed                               | T          |
| test_statistic | float | the test-statistic for that test                       | T          |
| p_value        | float | type 1 error probability                               | T          |
| fold_change    | float | between ref and comparison                             | T          |
| bh_p           | float | benjamini hochberg corrected p-value                   | T          |
| flattened      | bool  | whether the CRE was flattened due to insufficient UMIs | T          |

# Model formatting

(If we change our mind and instead desire a per-CRE model approach, we could replace this with a pandas dataframe). 
 
Model can be broken into a cohort parameters + a parameter table (but the reverse is not trivial).

# Cohort parameters formatting

A dictionary with keys
- alpha: alpha value for the dataset (constant term in mean-variance relationship). 

# Parameter table formatting

Note that this object is closely associated with a gross

# Preprocessing tools

scMPRA is a new technique and there are as of yet no standard flows for the initial analysis. These preprocessing tools are only for common-to but also-specific-to scMPRA tasks and are not a substituite for standard scRNA-seq analysis. 

 1. cell-filterer (kill from a list of low-quality cells)
 2. MPRABC-CRE association (assign CREs to UMIs on the basis of some lookup table generated from DNA-sequencing). 