#!/bin/bash
#SBATCH -c 6
#SBATCH --array=0-5
#SBATCH --mem-per-cpu=1G
#SBATCH -t 6:00:00
#SBATCH -p ycga

module load miniconda
conda activate biopython

inroot="/home/mcn26/palmer_scratch/raw_recap/cohen_retina/raw_fastq/"
files=($(ls $inroot))
file="${files[$SLURM_ARRAY_TASK_ID]}"
echo $file
zcat $inroot/$file | python get_stats.py > ${file}.stats

