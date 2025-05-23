#!/bin/bash
#SBATCH --job-name=permutation_array
#SBATCH --output=/gpfs/gibbs/project/reilly/eng26/scmpra/bin/tabula-rasa/stdout/permutation_array_%A_%a.out
#SBATCH --error=/gpfs/gibbs/project/reilly/eng26/scmpra/bin/tabula-rasa/stdout/permutation_array_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --partition=ycga

# set -euo pipefail

LOOKUP_FILE=$1
PERMUTATION_SCRIPT=$2
N_PERMUTATION=$3
DATE_STR=$4
COUNTS_FILE=$5
OUTPUT_DIR=$6
BIO_REP_MAP=$7

LOOKUP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LOOKUP_FILE")

if [[ -z "$LOOKUP" ]]; then
  echo "Error: LOOKUP is empty! Check your LOOKUP_FILE and SLURM_ARRAY_TASK_ID."
  exit 1
fi

echo "Processing CRE: $LOOKUP"

python "$PERMUTATION_SCRIPT" --lookup_str "$LOOKUP" --n_permutations "$N_PERMUTATION" --date_str "$DATE_STR" --counts_file "$COUNTS_FILE" --output_dir "$OUTPUT_DIR" --biol_rep_map "$BIO_REP_MAP"