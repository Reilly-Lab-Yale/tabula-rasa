#! /bin/bash -l

#SBATCH -J model_fitting
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH -o /home/eng26/project/scmpra/bin/tabula-rasa/stdout/%x_%A_%a.out
#SBATCH -e /home/eng26/project/scmpra/bin/tabula-rasa/stdout/%x_%A_%a.err
#SBATCH --array 1

###### usage ########
# sbatch --partition ycga --array 1-[num_lines_in_params] -J [name_job] --mem=128G --cpus-per-task 8 launch_model_fit.sh [params_file]
# sbatch --partition ycga --array 4 -J seelig_combo_479G --mem=479G --cpus-per-task 4 launch_model_fit.sh seelig_combo_model_params.txt 
module load miniconda
conda activate scmpra



datadir=/gpfs/gibbs/pi/reilly/tabula_data
tempdir=/home/eng26/palmer_scratch/scmpra_temp

# params_file=shendure_combo_params.txt
params_file=$1
echo $params_file

id=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $1}' ${params_file})
model=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $2}' ${params_file})
formula=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $3}' ${params_file})
zi_param=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $4}' ${params_file})
maxiter=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $5}' ${params_file})
split_level=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $6}' ${params_file})
split=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $7}' ${params_file})
counts_name=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR == row {print $8}' ${params_file})
counts=${datadir}/${counts_name}

echo 'launch python'
python basic_test.py --scmpra_counts_file $counts --model_choice $model --formula $formula --maxiter $maxiter --zi_param $zi_param --temp_dir $tempdir --out_file ${id}_${split_level}_${split}_${model}

