#!/bin/bash
#SBATCH -p ycga
#SBATCH -t 24:00:00
#SBATCH -J srDump
module load SRA-Toolkit
fastq-dump --split-files --gzip --outdir . SRR21787462
