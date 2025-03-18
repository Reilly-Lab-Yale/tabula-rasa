#! /bin/bash -l

#SBATCH -J combo_model_fitting
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH -o /home/eng26/project/scmpra/bin/tabula-rasa/stdout/testing_round4_%A_%a.out
#SBATCH -e /home/eng26/project/scmpra/bin/tabula-rasa/stdout/testing_round4_%A_%a.err
#SBATCH --array 1


module load miniconda
conda activate scmpra



datadir=/gpfs/gibbs/pi/reilly/tabula_data
params_file=shendure_combo_params.txt
split_level=COMBO

id=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $1}' params_file)
model=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $2}' params_file)
formula=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $3}' params_file)
zi_param=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $4}' params_file)
maxiter=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $5}' params_file)
split=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $6}' params_file)
counts_name=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $7}' params_file)
counts=${datadir}/${counts_name}

python basic_test.py --scmpra_counts_file $counts --model_choice $model --formula $formula --maxiter $maxiter --zi_param $zi_param --temp_dir $tempdir --out_file ${id}_${split_level}_${split}_${model}

