#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=0-7
#SBATCH -t 4:00:00

module load miniconda
conda activate env_preseq

root=/home/mcn26/palmer_scratch/raw_recap/cohen_retina/uniq_counts
files=($(ls $root))
current_file=${files[$SLURM_ARRAY_TASK_ID]}
echo "---FILE ${current_file} ---"
echo "---LC_EXTRAP---"
preseq lc_extrap -H ${root}/${current_file} -s 10000000 -n 1000
echo "---BOUND_POP---"
preseq bound_pop -H ${root}/${current_file} -s 10000000 -n 1000