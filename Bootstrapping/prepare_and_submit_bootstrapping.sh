#!/bin/bash
set -euo pipefail

# === Global Paths ===
DATA_PATH="/gpfs/gibbs/pi/reilly/tabula_data/shendure/bootstrapping"         # <-- all your data here
BIN_PATH="/gpfs/gibbs/project/reilly/eng26/scmpra/bin/tabula-rasa/Bootstrapping"       # <-- all your code here

# === Specific files/folders ===
COUNTS_FILE="$DATA_PATH/shendure_counts_grouped.txt"
LOOKUPS_FILE="$DATA_PATH/lookups.txt"
BOOTSTRAP_SCRIPT="$BIN_PATH/bootstrap_subsampling.py"
LOOKUPS_SCRIPT="$BIN_PATH/generate_lookups.py"
BIO_REP_MAP="$DATA_PATH/biol_rep_mapping.txt"
OUTPUT_DIR="$DATA_PATH/bootstrap_outputs"
N_BOOTSTRAP=100
DATE_STR=$(date +%Y%m%d)

# === Step 1: Generate lookups.txt ===
echo "Generating lookups file..."
python "$LOOKUPS_SCRIPT" "$COUNTS_FILE" "$LOOKUPS_FILE"

# === Step 2: Count number of CREs ===
NUM_LOOKUPS=$(wc -l < "$LOOKUPS_FILE")
echo "Number of lookups to process: $NUM_LOOKUPS"

# === Step 3: Submit array job ===
echo "Submitting array job..."

sbatch --array=1-"$NUM_LOOKUPS" <<EOT
#!/bin/bash
#SBATCH --job-name=bootstrap_array
#SBATCH --output=$BIN_PATH/../stdout/bootstrap_array_%A_%a.out
#SBATCH --error=$BIN_PATH/../stdout/bootstrap_array_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=ycga

set -euo pipefail

# Get the lookup corresponding to this array index
LOOKUP=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "$LOOKUPS_FILE")

echo "Processing CRE: \$LOOKUP"

# Run the bootstrap script
python "$BOOTSTRAP_SCRIPT" \
    "\$LOOKUP" \
    "$N_BOOTSTRAP" \
    "$DATE_STR" \
    "$DATA_PATH" \
    "$OUTPUT_DIR" \
    --biol_rep_map "$BIO_REP_MAP"
EOT