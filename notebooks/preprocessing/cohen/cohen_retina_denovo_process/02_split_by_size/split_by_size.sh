#!/bin/bash
#SBATCH -J lenbucket
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 7:00:00
#SBATCH -o slurm-%A_%a.out
#SBATCH --array=1-7

module load pigz/2.7-GCCcore-12.2.0

shopt -s nullglob

# Set these at submit time or let them default to CWD.
indir="/home/mcn26/palmer_scratch/raw_recap/cohen_retina/raw_fastq"
outdir="/home/mcn26/palmer_scratch/raw_recap/cohen_retina/length_split_fastq"
mkdir -p "$outdir"

# Build file list and pick the one for this array index
mapfile -t files < <(printf '%s\n' "$indir"/*.fastq)

idx="${SLURM_ARRAY_TASK_ID:-0}"
file="${files[$idx]}"

bn="${file##*/}"
bn="${bn%.fastq.gz}"

# Stream, select every 4th line starting at line 2, bucket by length
cat -- "$file" | awk -v base="$bn" -v outdir="$outdir" '
  NR % 4 == 2 {
    len = length($0);
    out = outdir "/" base "_len" len ".txt";
    print $0 >> out;
    close(out);  # avoid too many open files
  }
'

echo "Done: $file -> $outdir/${bn}_len*.txt"