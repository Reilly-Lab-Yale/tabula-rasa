#! /bin/bash -l

#SBATCH -J 10_03_2024
#SBATCH --mem=64G
#SBATCH --time=20:00:00
#SBATCH -N 4
#SBATCH -o /home/eng26/project/scmpra/bin/tabula-rasa/stdout/testing_round1_%A_%a.out
#SBATCH -e /home/eng26/project/scmpra/bin/tabula-rasa/stdout/testing_round1_%A_%a.err
#SBATCH --array 1-72


module load miniconda
conda activate scmpra



tempdir=/home/eng26/palmer_scratch/scmpra_temp
counts=/home/eng26/project/scmpra/data/shendure_mpra_counts_GSE217686.parq
id=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $1}' test_params.txt)
formula=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $4}' test_params.txt)
maxiter=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $5}' test_params.txt)
regfit=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $2}' test_params.txt)
model=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $3}' test_params.txt)

python basic_test.py --scmpra_counts_file $counts --model_choice $model --formula $formula --maxiter $maxiter --regularized_fit $regfit --temp_dir $tempdir --out_file ($id_$model$maxiter_$regfit_$formula)
