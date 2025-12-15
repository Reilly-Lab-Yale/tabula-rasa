#!/bin/bash
#SBATCH --job-name=uniq_counts
#SBATCH --array=1-8
#SBATCH --output=uniq_%A_%a.out
#SBATCH --error=uniq_%A_%a.err
#SBATCH --time=05:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=2

FILE_LIST="fastq_list.txt"   # list of .fastq files, one per line
OUTDIR="/home/mcn26/palmer_scratch/raw_recap/cohen_retina/uniq_counts"
INDIR="/home/mcn26/palmer_scratch/raw_recap/cohen_retina/raw_fastq"
mkdir -p "$OUTDIR"

fqname=$(sed -n "$((SLURM_ARRAY_TASK_ID))p" "$FILE_LIST")
FASTQ="${INDIR}/${fqname}"
BASENAME=$(basename "$FASTQ" .fastq)

echo "Processing $FASTQ ..."

sed -n '2~4p' "$FASTQ" | sort | uniq -c | awk '{print $1 "\t" $2}'  > "$OUTDIR/${BASENAME}_uniq.txt"

echo "Done: $OUTDIR/${BASENAME}_uniq.txt"
