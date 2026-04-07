#!/bin/bash
#SBATCH -J q20_filter
#SBATCH -p ycga
#SBATCH -c 3
#SBATCH --mem=8G
#SBATCH -t 04:00:00
#SBATCH --array=0-3

module reset 
module load cutadapt/5.0-GCCcore-12.2.0

shopt -s nullglob

INDIR="/home/mcn26/palmer_scratch/raw_recap/cohen_retina/raw_fastq"
OUTDIR="/home/mcn26/palmer_scratch/raw_recap/cohen_retina/filtered_fastq"

# Collect R1s; pairing assumes _1.fastq / _2.fastq
mapfile -t R1S < <(printf '%s\n' "$INDIR"/*_1.fastq | sort)
idx="${SLURM_ARRAY_TASK_ID:-0}"
R1="${R1S[$idx]}"
[[ -n "${R1:-}" ]] || { echo "No R1 file for index $idx"; exit 1; }

R2="${R1/_1.fastq/_2.fastq}"
[[ -f "$R2" ]] || { echo "Missing mate for $R1"; exit 2; }

base="$(basename "$R1" _1.fastq)"
O1="$OUTDIR/${base}_1.filtered.fastq"
O2="$OUTDIR/${base}_2.filtered.fastq"

# Q20 ~= average error rate ≤ 0.01
# Cutadapt discards pairs together by default; -j uses threads.
cutadapt \
  --max-aer 0.01 \
  -j "${SLURM_CPUS_PER_TASK:-4}" \
  -o "$O1" -p "$O2" \
  "$R1" "$R2" \
  --report=minimal
