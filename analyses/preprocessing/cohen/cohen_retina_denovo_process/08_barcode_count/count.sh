#!/bin/bash
#SBATCH -t 4:00:00
#SBATCH --array=0-3
#SBATCH -c 1
#SBATCH -p ycga
#SBATCH -J counts
#SBATCH --mem-per-cpu=5G

#After the last step, read-names now have extracted barcodes
#here I extract & count unique barcode combos. 
#Note that umitools puts the bacrodes in the field-1 (lines that start with @) in both R1 & 2 output
#so we only need to process one read of the pair for each. I have arbitrarially chosen R1

paths=("main_bc_extract/SRR21787462/SRR21787462.processed_R1.fastq" "main_bc_extract/SRR32774353/SRR32774353.processed_R1.fastq" "u6_bc_extract/SRR21787461/SRR21787461.processed_R1.fastq" "u6_bc_extract/SRR21787460/SRR21787460.processed_R1.fastq" )

inp_root="/home/mcn26/palmer_scratch/raw_recap/cohen_retina/"
out_root="/home/mcn26/palmer_scratch/raw_recap/cohen_retina/counts"


filename=${paths[${SLURM_ARRAY_TASK_ID}]}

echo "Processing ${filename}"

#at this step we process both u6 and non-u6 since umitools smushes all umis together, making
#some processing identical for both

#examples:
#normal mpra looks like
#@SRR21787462.71_TTGCATTAGTGTACAA_AGCCGTGATATAATTGGTAGCGAAAATTGAAACTAGACCGGATT NB501065:669:HJL2TBGXL:1:11101:11858:1071 length=28
#then splitting `cut -d'_' -f2-3`
#TTGCATTAGTGTACAA, AGCCGTGATATAATTGGTAGCGAAAATTGAAACTAGACCGGATT NB501065:669:HJL2TBGXL:1:11101:11858:1071 length=28
#then splitting `| cut -d' '  -f1 `
#TTGCATTAGTGTACAA, AGCCGTGATATAATTGGTAGCGAAAATTGAAACTAGACCGGATT
#where TTGCATTAGTGTACAA is 16bp <cell_1>, cbc, AGCCGTGATATA is 12bp, <umi_1> umi, ATTGGTAG is 8bp <umi_2> pBC, CGAAAATTGAAACTAGACCGGATT is 24 bp, <umi_3> rBC

#U6 looks like
#@SRR21787461.241_ATCCGTCGTGGTGAGT_TTGTCAAAGACTACGTATGG NB501801:559:H3VVTBGXM:1:11101:10557:1087 length=28
#processing similarly, we get
#ATCCGTCGTGGTGAGT, TTGTCAAAGACTACGTATGG
#where ATCCGTCGTGGTGAGT is 16bp <cell_1>, TTGTCAAAGACT is 12bp actual umi, ACGTATGG is 8bp pbc

base=$(basename ${filename})

cat ${inp_root}/${filename} \
	| grep -E "^@" \
	| cut -d'_' -f2-3 \
	| cut -d' '  -f1 \
	| sort \
	| uniq -c \
	> ${out_root}/${base}.count

