# introdution

(take from abstract)

- Find the preprint HERE
- Find the data zenodo HERE (dump of /gpfs/gibbs/pi/reilly/tabula_data)
- Find the github archival zenodo HERE

Skip down to "install" to get started!

This repository is organzied as follows:
```

├── docs					# documentation
├── notebooks				# jupyter notebooks
│   ├── object_creation		# creates objects from scMPRA data
│   ├── results				# figures and numerical results
│   └── preprocessing		# preprocesses emperical scMPRA data
└── scMPRAforge				# the package proper
    ├── vingettes			# notebooks demonstrating the use of the package
    ├── presets				# presets to be distributed with the package
    └── tests				# tests
```

"scMPRAforge" contains the actual pacakge

"Preprocessing" processes emperical data for "object_creation".

"object_creation" creates and saves several python objects (`ortho`, `scMPRA_data`, `Bounds` and similar) useful in multiple analyses run in "results". Intermediate values are not distributed. Simulated datasets & processed emperical datasets are distributed separately, as `.scmpra` files on a data zenodo **LINK**. Bounds objects are stored in presets directory in main package.

"results" contains notebooks generating results in the paper. This includes figure panels and numbers. 

# key assumptions
(move to paper)
1. Though ideally we would model each barcode (like the approach in MPRAmodel) for greater statistical power, the sparse nature of scMPRA data means that we won't be able to do this. Not that the barcodes are retained, just not modeled separaately. 
2. Removal of "false zeroes" through clonotype analysis and transfection reporters is part of pre-processing. 
3. We are interested in changes in CRE activity within and between cell-types


# install
(temporary : replace with pip install scMPRAforge once we have actually )
```
conda create -n env_tensorzinb python=3.10
conda activate env_tensorzinb
echo "numpy <2" >> $CONDA_PREFIX/conda-meta/pinned
#navigate to tensor zinb
#install apple silicon deps as below if necessary
python setup.py install
```

For apple silicon, may need
```
conda install -c apple tensorflow-deps
python -m pip install tensorflow-macos==2.9.2
python -m pip install tensorflow-metal==0.5.1
```

```
conda install statsmodels 
```

Additional required packages
```
conda install matplotlib seaborn bioconda::umi_tools formulaic dask dask-jobqueue
```

may need to install `tensorflow keras` if you didn't above

At this point, it may work fine. However some environments (such as mccleary gpu nodes) may require some more fnagling to get GPU accel working. Perform the following:
```
conda install cudatoolkit=11.2 cudnn=8.1.0
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_PREFIX/lib/' > $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
```

latter not necessary on bouchet

then deactivate, activate, and test with
```
python3 -c "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"
```

If you see a device, you're good to go.

(adapted from [mccleary docs](https://docs.ycrc.yale.edu/clusters-at-yale/guides/tensorflow/))

# supplementary conda environments

The above is sufficient for you to use the package, but recapitulating some of the preprocessing steps & results from the paper require other conda environments. 
Fill in 
- [ ] biopython (used for some data carpentry and single cell processing steps)
- [ ] parse_scanpy (used for scanpy : particular scanpy version particularly important to get same leiden clustering)
- [ ] speedracer (used for some rapid text-processing steps)

# Data formats

This section describes the various standard tabular formats used internally.

## MPRA data formatting

Note that most of our code does not *require* that any of the barcode sequences should be actual nucleotide sequences. So you can replace them with, for example, numerical combinatorial barcoding sub-barcode ID strings with no consequence. The exception is barcode-deduplication. 



(Collapsing cells, replicates, and MPRA barcodes is not meaningful for us, since this would collapse what we believe to be true biological samples. However we may frequently *summarize* (e.g. compute mean UMIs per cell, or model a distribution where all UMI count originating from one CRE in one cell-type (regardless of MPRA barcode) are considered to have come from one triplicate of zinb parameters).

In memory, dataframes use named columns with dummy row-indicies. 

scMPRA data can be umi-wise or read-wise.

UMI = unique molecular identifier

**Read-wise**
| Column name           | Type             | Description                          | Mandatory? |
| --------------------- | ---------------- | ------------------------------------ | ---------- |
| cell_bc               | str (nucleotide) | cell barcode                         | T          |
| rep_id                | str              | replicate id                         | T          |
| cre_id                | str              | CRE id or name                       | T          |
| cell_type             | str              | cell-type                            | T          |
| mpra_bc               | str (nucleotide) | MPRA reporter barcode                | T          |
| mpra_umi              | str (nucleotide) | MPRA umi                             | T          |
| reads_mpra            | int              | number mpra reads                    | T          |
| transfection_bc       | str (nucleotide) | transfection barcode                 | F          |
| transfection_umi      | str (nucleotide) | transfection reporter umi            | F          |
| reads_transfection_bc | int              | number transfection reporter reads   | F          |
| reads_DNA             | int              | DNA library reads (for MPRA barcode) | F          |

**Umi-wise**
| Column name           | Type             | Description                                                      | Mandatory? |
| --------------------- | ---------------- | ---------------------------------------------------------------- | ---------- |
| cell_bc               | str (nucleotide) | cell barcode                                                     | ~T          |
| rep_id                | str              | replicate id                                                     | T          |
| cre_id                | str              | CRE id or name                                                   | T          |
| cell_type             | str              | cell-type                                                        | T          |
| mpra_bc               | str (nucleotide) | MPRA reporter barcode                                            | T          |
| umis_mpra_bc          | int              | number of MPRA UMIs                                              | T          |
| reads_mpra_bc         | int              | number of MPRA reads, summed across all UMIs                     | F          |
| transfection_bc       | str (nucleotide) | transfection barcode                                             | F          |
| umis_transfection_bc  | int              | number transfection reporter UMIs                                | F          |
| reads_transfection_bc | int              | number transfection reporter reads, summed across all UMIs       | F          |
| reads_DNA             | int              | DNA library reads (for MPRA barcode)                             | F          |

(Though MPRA barcodes cannot be modeled individually, its undersierable to sum across to remove them : since this would inflate strength estimates of CREs with more MPRA barcodes)

Simulated UMI-wise data doesn't currently require cell barcode. 

Note that DNA is treated as constant for a given CRE : it's never reduced on account of certain MPRA barcodes being present in certain cells.
This keeps it as a totally exogenous source of information, far from the vicissitudes of single-cell sequencing...

Also note that all these strings are really factors / categorical data, and will be treated as such. 

 ## Hypothesis formatting
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

- NOTE: for tiny pvalues where floating point sizes require rounding p, the exact value of this epsilon will depend on the computer hardware/architecture that you are operating on. 
  ```
   eps = np.finfo(float).tiny  # ~1e-308
   p = max(1.0 - chi2.cdf(z*z, 1), eps)
  ```

## Ground-truth formatting
These tables store made-up ground truth, for the purposes of simulation. 
| Column name | Type  | Description                           |
|-------------|-------|---------------------------------------|
| cell_type   | str   | cell-type                             |
| cre_id      | str   | CRE id or name                        |
| true_mean   | float | Ground-truth # MPRA barcode UMIs/cell |

## MPRA library table

| Column name | Type                | Description                       | Mandatory? |
|-------------|---------------------|-----------------------------------|------------|
| cre_id      | string              | CRE id or name                    | T          |
| mpra_bc     | string (nucleotide) | MPRA reporter barcode             | T          |
| abundance   | float               | Relative abundance in DNA library | T          |

Abundance column must sum to 1. It can be defined in different ways, but in most cases will be (reads MPRA barcode)/(total MPRA barcode reads)




