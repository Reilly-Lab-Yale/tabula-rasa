#! /bin/bash -l

#SBATCH -J 12_11_2024
#SBATCH --mem=32G
#SBATCH --time=20:00:00
#SBATCH -o /home/eng26/project/scmpra/bin/tabula-rasa/stdout/testing_round5_%A_%a.out
#SBATCH -e /home/eng26/project/scmpra/bin/tabula-rasa/stdout/testing_round5_%A_%a.err
#SBATCH --array 1-8


module load miniconda
conda activate scmpra_r




dataset=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $1}' r_params.txt)
model_choice=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $2}' r_params.txt)
out_file=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $3}' r_params.txt)

Rscript basic_test.R --dataset  $dataset --model_choice $model_choice --out_file $out_file
