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

rep_run_name="SRR21787461"

#echo "unzipping R1..."
#gunzip $root_inp/${rep_run_name}_1.fastq.gz_R1_paired.fastq.gz
#echo "unzipping R2..."
#gunzip $root_inp/${rep_run_name}_1.fastq.gz_R2_paired.fastq.gz

echo "beginning barcode extraction..."


umi_tools extract \
	--extract-method=regex \
	--bc-pattern='^(?P<cell_1>.{16})(?P<umi_1>.{12})' \
	--bc-pattern2='(?P<discard_1>SCGGTAAGCTCCCGGGAGCTTGT){s<=3}(?P<umi_2>.{8})' \
	--stdin=$root_inp/${rep_run_name}_1.fastq.gz_R1_paired.fastq \
	--stdout=${root_out}/${rep_run_name}_processed_R1.fastq \
	--read2-in=$root_inp/${rep_run_name}_1.fastq.gz_R2_paired.fastq \
	--read2-out=${root_out}/${rep_run_name}_processed_R2.fastq
	--log=$root_out/${rep_run_name}_extract.log \
