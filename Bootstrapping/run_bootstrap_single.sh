#!/bin/bash
#SBATCH --job-name=bootstrap_array
#SBATCH --output=/gpfs/gibbs/project/reilly/eng26/scmpra/bin/tabula-rasa/stdout/bootstrap_array_%A_%a.out
#SBATCH --error=/gpfs/gibbs/project/reilly/eng26/scmpra/bin/tabula-rasa/stdout/bootstrap_array_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=ycga

set -euo pipefail

module load miniconda
conda activate scmpra

LOOKUP_FILE="$1"
BOOTSTRAP_SCRIPT="$2"
N_BOOTSTRAP="$3"
DATE_STR="$4"
COUNTS_FILE="$5"
OUTPUT_DIR="$6"
BIO_REP_MAP="$7"

LOOKUP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LOOKUP_FILE")

if [[ -z "$LOOKUP" ]]; then
  echo "Error: LOOKUP is empty! Check your LOOKUP_FILE and SLURM_ARRAY_TASK_ID."
  exit 1
fi

echo "Processing CRE: $LOOKUP"

python "$BOOTSTRAP_SCRIPT" \
    --lookup_str "$LOOKUP" \
    --n_bootstrap "$N_BOOTSTRAP" \
    --date_str "$DATE_STR" \
    --counts_file "$COUNTS_FILE" \
    --output_dir "$OUTPUT_DIR" \
    --biol_rep_map "$BIO_REP_MAP"