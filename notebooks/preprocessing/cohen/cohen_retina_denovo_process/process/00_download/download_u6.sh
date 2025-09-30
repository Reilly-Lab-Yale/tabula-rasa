#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -p ycga
#SBATCH -c 3
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-2

module load SRA-Toolkit
SRA_LIST="sra_runs.txt"
SRA_RUN=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SRA_LIST)
echo "Downloading SRA run: $SRA_RUN"
prefetch $SRA_RUN
echo "Extracting FASTQ files for SRA run: $SRA_RUN"
fastq-dump --split-files --gzip $SRA_RUN
echo "Completed processing SRA run: $SRA_RUN"
