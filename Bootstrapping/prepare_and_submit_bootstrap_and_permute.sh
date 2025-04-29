#!/bin/bash
set -euo pipefail

module load miniconda
conda activate scmpra

# === Global Paths ===
DATA_PATH="/gpfs/gibbs/pi/reilly/tabula_data/shendure"
BIN_PATH="/gpfs/gibbs/project/reilly/eng26/scmpra/bin/tabula-rasa/Bootstrapping"

# === Specific files/folders ===
COUNTS_FILE="$DATA_PATH/shendure_counts_grouped.txt"
LOOKUPS_FILE="$DATA_PATH/bootstrapping/lookups.txt"
BOOTSTRAP_SCRIPT="$BIN_PATH/bootstrap_subsampling.py"
PERMUTATION_SCRIPT="$BIN_PATH/permutation_subsampling.py"   # <-- NEW
LOOKUPS_SCRIPT="$BIN_PATH/generate_lookups.py"
RUN_SINGLE_BOOTSTRAP_SCRIPT="$BIN_PATH/run_bootstrap_single.sh"
RUN_SINGLE_PERMUTATION_SCRIPT="$BIN_PATH/run_permutation_single.sh"  # <-- NEW
BIO_REP_MAP="$DATA_PATH/biol_rep_mapping.txt"
OUTPUT_DIR_BOOTSTRAP="$DATA_PATH/bootstrapping"
OUTPUT_DIR_PERMUTATION="$DATA_PATH/permutation"   # <-- NEW
N_BOOTSTRAP=$((10**4)) ## update to 10**4 after troubleshooting $((10**4))
N_PERMUTATION=$((10**4))  # you can change separately if you want $((10**4))
DATE_STR=$(date +%Y%m%d)

# === Step 1: Generate lookups.txt ===
echo "Generating lookups file..."
python "$LOOKUPS_SCRIPT" --counts_file "$COUNTS_FILE" --output_file "$LOOKUPS_FILE"

# === Step 2: Count number of CREs ===
NUM_LOOKUPS=$(wc -l < "$LOOKUPS_FILE")
echo "Number of lookups to process: $NUM_LOOKUPS"

# # === Step 3: Submit array job for bootstrapping ===
# echo "Submitting array job for bootstrapping..."
# #"$NUM_LOOKUPS"
# # 1-"$NUM_LOOKUPS"
# array_jobid_bootstrap=$(sbatch --parsable --array=211-211 "$RUN_SINGLE_BOOTSTRAP_SCRIPT" \
#     "$LOOKUPS_FILE" \
#     "$BOOTSTRAP_SCRIPT" \
#     "$N_BOOTSTRAP" \
#     "$DATE_STR" \
#     "$COUNTS_FILE" \
#     "$OUTPUT_DIR_BOOTSTRAP" \
#     "$BIO_REP_MAP"
# )
# echo "Bootstrap array job submitted with Job ID: $array_jobid_bootstrap"

echo "Skipping bootstrap step, submitting dummy job for dependency handling..."
array_jobid_bootstrap=$(sbatch --parsable --job-name=dummy_bootstrap_done --wrap="sleep 1")

# # === Step 4: Submit array job for permutation ===
# echo "Submitting array job for permutation..."
# # 
# array_jobid_permutation=$(sbatch --parsable --array=1-"$NUM_LOOKUPS" "$RUN_SINGLE_PERMUTATION_SCRIPT" \
#     "$LOOKUPS_FILE" \
#     "$PERMUTATION_SCRIPT" \
#     "$N_PERMUTATION" \
#     "$DATE_STR" \
#     "$COUNTS_FILE" \
#     "$OUTPUT_DIR_PERMUTATION" \
#     "$BIO_REP_MAP"
# )
# echo "Permutation array job submitted with Job ID: $array_jobid_permutation"
echo "Skipping permutation step, submitting dummy job for dependency handling..."
array_jobid_permutation=$(sbatch --parsable --job-name=dummy_permutation_done --wrap="sleep 1")

# === Step 5: Submit combination jobs (optional) ===
echo "Submitting combine jobs after bootstrap and permutation finish..."

COMBINE_BOOTSTRAP_SCRIPT="$BIN_PATH/combine_bootstraps.py"
COMBINE_PERMUTATION_SCRIPT="$BIN_PATH/combine_permutations.py"  

combine_jobid_bootstrap=$(sbatch --parsable --partition=ycga --time=08:00:00 --dependency=afterok:$array_jobid_bootstrap --job-name=combine_bootstrap_outputs --wrap="\
module load miniconda; \
conda activate scmpra; \
python $COMBINE_BOOTSTRAP_SCRIPT --input_dir $OUTPUT_DIR_BOOTSTRAP --output_file $OUTPUT_DIR_BOOTSTRAP/combined_bootstrap_allCREs_${DATE_STR}.tsv")

combine_jobid_permutation=$(sbatch --parsable --partition=ycga --time=08:00:00 --dependency=afterok:$array_jobid_permutation --job-name=combine_permutation_outputs --wrap="\
module load miniconda; \
conda activate scmpra; \
python $COMBINE_PERMUTATION_SCRIPT --input_dir $OUTPUT_DIR_PERMUTATION --output_file $OUTPUT_DIR_PERMUTATION/combined_permutation_allCREs_${DATE_STR}.tsv")

echo "Combine bootstrap job submitted with Job ID: $combine_jobid_bootstrap"
echo "Combine permutation job submitted with Job ID: $combine_jobid_permutation"

# === Step 6: Submit empirical p-value calculation job ===
echo "Submitting empirical p-value job..."

EMP_PVAL_SCRIPT="$BIN_PATH/compute_empirical_pvalues.py"
EMP_PVAL_OUTPUT_activity="$DATA_PATH/empirical_pvals_activity_${DATE_STR}.tsv"
EMP_PVAL_OUTPUT_specificity="$DATA_PATH/empirical_pvals_activity_${DATE_STR}.tsv"

mkdir -p "$DATA_PATH/empirical_pvalues"  # make output dir if needed

empirical_pval_jobid=$(sbatch --parsable --time=08:00:00 --partition=ycga --dependency=afterok:$combine_jobid_bootstrap,$combine_jobid_permutation --job-name=empirical_pvals --wrap="\
module load miniconda; \
conda activate scmpra; \
python $EMP_PVAL_SCRIPT --bootstrap_file $OUTPUT_DIR_BOOTSTRAP/combined_bootstrap_allCREs_${DATE_STR}.tsv \
--permutation_file $OUTPUT_DIR_PERMUTATION/combined_permutation_allCREs_${DATE_STR}.tsv \
--output_activity_file $EMP_PVAL_OUTPUT_activity \
--output_specificity_file $EMP_PVAL_OUTPUT_specificity") \


echo "Empirical p-value job submitted with Job ID: $empirical_pval_jobid"