#!/bin/bash
#SBATCH -c 1
#SBATCH --array=0-7
#SBATCH --mem-per-cpu=1G
#SBATCH -t 6:00:00
#SBATCH -p ycga

module load miniconda
conda activate speedracer

inroot="/home/mcn26/palmer_scratch/raw_recap/cohen_retina/filtered_fastq"
files=($(ls $inroot))
file="${files[$SLURM_ARRAY_TASK_ID]}"
outroot="/home/mcn26/palmer_scratch/raw_recap/cohen_retina/length_split_fastq"

cat $inroot/$file | pypy3 split_by_size.py $file $outroot