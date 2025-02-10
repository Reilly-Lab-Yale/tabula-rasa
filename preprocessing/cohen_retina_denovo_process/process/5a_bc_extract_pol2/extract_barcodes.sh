#!/bin/bash
#SBATCH -p ycga
#SBATCH -t 12:00:00
#SBATCH -c 1
#SBATCH --mem-per-cpu=12G

module unload Python
module load miniconda
conda activate biopython

root_inp="/home/mcn26/palmer_scratch/tabula_data/raw_recap/cohen_retina/filtered_fastq"
root_out="/home/mcn26/palmer_scratch/tabula_data/raw_recap/cohen_retina/main_bc_extract"

#echo "unzipping R1..."
#gunzip $root_inp/SRR21787462_1.fastq.gz_R1_paired.fastq.gz
#echo "unzipping R2..."
#gunzip $root_inp/SRR21787462_1.fastq.gz_R2_paired.fastq.gz

echo "beginning barcode extraction..."

umi_tools extract \
	--extract-method=regex \
	--bc-pattern='^(?P<cell_1>.{16})(?P<umi_1>.{12})' \
	--bc-pattern2='^(?P<discard_1>GAGCTGTACAAGTAACTGCGAT){s<=3}(?P<umi_2>.{8})(?P<discard_2>AAGGAACCCGCGCTATACCGGTATCGC){s<5}(?P<umi_3>.{24})' \
	--stdin=$root_inp/SRR21787462_1.fastq.gz_R1_paired.fastq \
	--read2-in=$root_inp/SRR21787462_1.fastq.gz_R2_paired.fastq \
	--stdout=$root_out/processed.fastq \
	--log=$root_out/extract.log \
	--read2-out=processed_R2.fastq

#Manually renamed
# mv processed.fastq SRR21787462.fast
