#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=0-7
#SBATCH -t 4:00:00

root=/home/mcn26/palmer_scratch/raw_recap/cohen_retina
inp=${root}/uniq_counts
oup=${root}/histograms
files=($(ls $inp))
current_file=${files[$SLURM_ARRAY_TASK_ID]}
name=$(basename "$current_file" ".txt")

echo "Processing ${current_file}"
cat ${inp}/${current_file} | cut -f1 | sort | uniq -c | awk '{print $2 "\t" $1}' > ${oup}/${name}.hist
# we reverse the order because preseq expects 