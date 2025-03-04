#! /bin/bash -l

#SBATCH -J 10_03_2024
#SBATCH --mem=64G
#SBATCH --time=20:00:00
#SBATCH -o /home/eng26/project/scmpra/bin/tabula-rasa/stdout/testing_round4_%A_%a.out
#SBATCH -e /home/eng26/project/scmpra/bin/tabula-rasa/stdout/testing_round4_%A_%a.err
#SBATCH --array 1-212


module load miniconda
conda activate scmpra



tempdir=/home/eng26/palmer_scratch/scmpra_temp
#counts=/home/eng26/project/scmpra/data/shendure_mpra_counts_grouped_GSE217686.parq
#counts=/home/eng26/project/scmpra/data/shendure_counts_grouped_cardio.txt
id=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $1}' CRE_indi_test_params.txt)
model=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $2}' CRE_indi_test_params.txt)
formula=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $3}' CRE_indi_test_params.txt)
zi_param=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $4}' CRE_indi_test_params.txt)
maxiter=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $5}' CRE_indi_test_params.txt)
cre=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $6}' CRE_indi_test_params.txt)
counts=/home/eng26/project/scmpra/data/shendure_mpra_counts_${cre}.parq

python basic_test.py --scmpra_counts_file $counts --model_choice $model --formula $formula --maxiter $maxiter --zi_param $zi_param --temp_dir $tempdir --out_file ${id}_${cre}_GROUPED_${model}${maxiter}_${formula}_${zi_param}

