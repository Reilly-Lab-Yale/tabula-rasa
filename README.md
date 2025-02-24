Assumptions:
1. Barcode flattening is mandatory. Though ideally we would model each barcode as per MPRAmodel for greater statistical power, the sparse nature of scMPRA data means that we won't be able to do this. 
2. Removal of "false zeroes" through clonotype analysis and transfection reporters is part of pre-processing. 
3. We are interested in changes in CRE activity within and between cell-types
4. I've run ahead with the "one big model" approach (instead of the "one model per CRE" approach) but we can change this later based on some exploratory analysis. 

(We may eventually wish to turn the ad-hoc pre-processing steps into common tools.)


# MPRA data formatting

All on-disc data should be tsv or (some high-performace format to be chosen later).

Note that nothing in our code should *require* that any of the barcode sequences should be actual nucleotide sequences. So you can replace them with, for example, numerical combinatorial barcoding sub-barcode ID strings with no consequence.

However, given python's weak typing, I'd recommend not putting something that could be obviously misinterpreted as another type (e.g. replicate a single int as this could easially be misinterpreted to suggest an ordinal position).

Data can be umi-wise or read-wise, and with or without collapsed barcodes, for a total of four possible formats. Explicit column names and types are given below for each combination.

(Collapsing cells or replicates is not meaningful for us, since this would collapse what we believe to be true biological replicates. We may *summarize* (e.g. mean UMIs per cell).

**Read-wise UNcollapsed-on-barcodes**
| Column name | Type             | Description                 | Mandatory? |
| ----------- | ---------------- | --------------------------- | ---------- |
| cell_bc     | str (nucleotide) | cell barcode                | T          |
| rep_id      | str              | replicate id                | T          |
| cre_id      | str              | CRE id or name              | T          |
| cell_type   | str              | cell-type                   | T          |
| mpra_bc     | str (nucleotide) | MPRA reporter barcode       | T          |
| umi         | str (nucleotide) | unique molecular identifier | T          |
| reads       | int              | number of reads             | T          |

**Umi-wise UNcollapsed-on-barcodes**
| Column name | Type             | Description                             | Mandatory? |
| ----------- | ---------------- | --------------------------------------- | ---------- |
| cell_bc     | str (nucleotide) | cell barcode                            | T          |
| rep_id      | str              | replicate id                            | T          |
| cre_id      | str              | CRE id or name                          | T          |
| cell_type   | str              | cell-type                               | T          |
| mpra_bc     | str (nucleotide) | MPRA reporter barcode                   | T          |
| umis        | int              | number of unique molecular identifiers  | T          |
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
- For assessing variant effects, we recommend the convention that reference columns refer to the reference allele, and comparison columns refer to the alternate allele. You can explicitly specify this in the metadata column. 

 # Preprocessing tools

scMPRA is a new technique and there are as of yet no standard flows for the initial analysis. These preprocessing tools are only for common-to but also-specific-to scMPRA tasks and are not a substituite for standard scRNA-seq analysis. 

 1. cell-filterer (kill from a list of low-quality cells)
 2. MPRABC-CRE association (assign CREs to UMIs on the basis of some lookup table generated from DNA-sequencing). 