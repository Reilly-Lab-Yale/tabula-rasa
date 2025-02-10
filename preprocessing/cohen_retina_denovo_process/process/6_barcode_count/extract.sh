#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --array=0-2
#SBATCH -c 10
#SBATCH -p ycga
#SBATCH --mem-per-cpu=1G

#After the last step, read-names now have extracted barcodes
#here I extract & count unique barcode combos

filenames=("SRR21787460_processed_R1.fastq" "SRR21787461_processed_R1.fastq" "SRR21787462.fastq")

inp_root="/home/mcn26/palmer_scratch/tabula_data/raw_recap/cohen_retina/main_bc_extract"
out_root="/home/mcn26/palmer_scratch/tabula_data/raw_recap/cohen_retina/counts"

idx=0
filename=${filenames[${SLURM_ARRAY_TASK_ID}]}

echo "Processing ${filename}"

cat ${inp_root}/${filename} \
	| grep -E "^@" \
	| cut -d'_' -f2-3 \
	| cut -d' '  -f1 \
	| sort \
	| uniq -c \
	> ${out_root}/${filename}.count

