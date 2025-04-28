#!/bin/bash
set -euo pipefail

LOOKUPS_FILE=$1
PERMUTATION_SCRIPT=$2
N_PERMUTATION=$3
DATE_STR=$4
COUNTS_FILE=$5
OUTPUT_DIR=$6
BIO_REP_MAP=$7

module load miniconda
conda activate scmpra

LOOKUP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LOOKUPS_FILE")

python "$PERMUTATION_SCRIPT" --lookup_str "$LOOKUP" --n_permutations "$N_PERMUTATION" --date_str "$DATE_STR" --counts_file "$COUNTS_FILE" --output_dir "$OUTPUT_DIR" --biol_rep_map "$BIO_REP_MAP"