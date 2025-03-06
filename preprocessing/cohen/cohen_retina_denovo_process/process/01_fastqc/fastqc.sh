#!/bin/bash
#SBATCH -c 10
#SBATCH --mem-per-cpu=1G
#SBATCH -t 12:00:00
#SBATCH -p ycga

module load FastQC/0.12.1-Java-11
fastqc -v

# Directory containing the .fastq.gz files
input_dir="/home/mcn26/palmer_scratch/raw_recap/cohen_retina/raw_fastq"
output_root="/home/mcn26/palmer_scratch/raw_recap/cohen_retina/fastqc"

# Loop through all .fastq.gz files in the directory
for file in "$input_dir"/*.fastq.gz; do
    # Get the base name of the file (without the directory and extension)
    base_name=$(basename "$file" .fastq.gz)
    
    # Create a directory for the output named after the file
    output_dir="$output_root/$base_name"
    mkdir -p "$output_dir"
    
    # Run FastQC on the file and place the output in the created directory
    fastqc --noextract -t 10 -o "$output_dir" "$file"
done
