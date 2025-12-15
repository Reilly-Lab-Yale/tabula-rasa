#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=0-7
#SBATCH -t 4:00:00

#previous step counts occurances of reads, producing
#5 AATGATGA
#3 TAGATAGATATGAG
#etc...
#we call those ocurances read counts
#this script makes a "histogram": it counts 
#how many times we see a read count of n

root=/home/mcn26/palmer_scratch/raw_recap/cohen_retina
inp=${root}/uniq_counts
oup=${root}/histograms
files=($(ls $inp))
current_file=${files[$SLURM_ARRAY_TASK_ID]}
name=$(basename "$current_file" ".txt")

echo "Processing ${current_file}"
cat ${inp}/${current_file} | cut -f1 | sort | uniq -c | awk '{print $2 "\t" $1}' | sort -g > ${oup}/${name}.hist
# we reverse the order (print 2 then 1) because preseq expects "read count" then "frequency"
# we employ one final sort because preseq is a diva and wants the data pre-sorted by read count