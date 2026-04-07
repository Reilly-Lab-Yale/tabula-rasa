#!/bin/bash
#SBATCH -p ycga
#SBATCH -t 12:00:00
#SBATCH -c 1
#SBATCH --mem-per-cpu=12G
#SBATCH --array=0-1
#SBATCH -J bc_extract
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

module reset
module load miniconda
conda activate biopython

root_inp="/home/mcn26/palmer_scratch/raw_recap/cohen_retina/filtered_fastq"
root_out="/home/mcn26/palmer_scratch/raw_recap/cohen_retina/main_bc_extract"

# Each array index handles one pair
SAMPLES=(SRR21787462 SRR32774353)
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

R1="$root_inp/${SAMPLE}_1.filtered.fastq"
R2="$root_inp/${SAMPLE}_2.filtered.fastq"

outdir="$root_out/$SAMPLE"
mkdir -p "$outdir"

echo "Beginning barcode extraction for $SAMPLE"

umi_tools extract \
  --extract-method=regex \
  --bc-pattern='^(?P<cell_1>.{16})(?P<umi_1>.{12})' \
  --bc-pattern2='^(?P<discard_1>GAGCTGTACAAGTAACTGCGAT){s<=3}(?P<umi_2>.{8})(?P<discard_2>AAGGAACCCGCGCTATACCGGTATCGC){s<5}(?P<umi_3>.{24})' \
  --stdin="$R1" \
  --read2-in="$R2" \
  --stdout="$outdir/${SAMPLE}.processed_R1.fastq" \
  --log="$outdir/${SAMPLE}.extract.log" \
  --read2-out="$outdir/${SAMPLE}.processed_R2.fastq"

echo "Done: $SAMPLE"
