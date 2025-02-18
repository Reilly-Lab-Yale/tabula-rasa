#! /bin/bash -l

#SBATCH -J 10_03_2024
#SBATCH --mem=32G
#SBATCH --time=20:00:00
#SBATCH -o /home/eng26/project/scmpra/bin/tabula-rasa/stdout/testing_round4_%A_%a.out
#SBATCH -e /home/eng26/project/scmpra/bin/tabula-rasa/stdout/testing_round4_%A_%a.err
#SBATCH --array 1-12


module load miniconda
conda activate scmpra



tempdir=/home/eng26/palmer_scratch/scmpra_temp
counts=/home/eng26/project/scmpra/data/shendure_mpra_counts_GSE217686.parq
#counts=/home/eng26/project/scmpra/data/shendure_counts_grouped_cardio.txt
id=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $1}' test_params.txt)
formula=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $4}' test_params.txt)
maxiter=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $5}' test_params.txt)
#regfit=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $2}' test_params.txt)
model=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $3}' test_params.txt)

#echo 'regfit: ' $regfit

python basic_test.py --scmpra_counts_file $counts --model_choice $model --formula $formula --maxiter $maxiter --regularized_fit $regfit --temp_dir $tempdir --out_file ${id}GROUPED_${model}${maxiter}_${regfit}_${formula}

#python basic_test_bambi.py --scmpra_counts_file $counts --model_choice $model --formula $formula --maxiter $maxiter --temp_dir $tempdir --out_file ${id}_cardio_${model}${maxiter}_${formula}
