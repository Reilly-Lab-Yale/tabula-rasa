#!/bin/bash
#SBATCH -p ycga
#SBATCH -t 6:00:00
#SBATCH -c 9

module load Trimmomatic/0.39-Java-11

# Input and output directories
INPUT_DIR="/home/mcn26/palmer_scratch/tabula_data/raw_recap/cohen_retina/raw_fastq"
OUTPUT_DIR="/home/mcn26/palmer_scratch/tabula_data/raw_recap/cohen_retina/filtered_fastq"

# Path to Trimmomatic JAR (update if necessary)
TRIMMOJAR="$EBROOTTRIMMOMATIC/trimmomatic-0.39.jar"

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop through all *_R1.fastq.gz files and find their matching *_R2.fastq.gz files
for R1 in "$INPUT_DIR"/*_1.fastq.gz; do
    # Derive the corresponding R2 filename
    R2="${R1/_1.fastq.gz/_2.fastq.gz}"
    
    # Ensure the corresponding R2 file exists before proceeding
    if [[ ! -f "$R2" ]]; then
        echo "Warning: Matching R2 file for $R1 not found. Skipping..."
        continue
    fi

    # Extract base name for output files
    SAMPLE_NAME=$(basename "$R1" _R1.fastq.gz)
    
    # Define output file paths
    P1="$OUTPUT_DIR/${SAMPLE_NAME}_R1_paired.fastq.gz"
    U1="$OUTPUT_DIR/${SAMPLE_NAME}_R1_unpaired.fastq.gz"
    P2="$OUTPUT_DIR/${SAMPLE_NAME}_R2_paired.fastq.gz"
    U2="$OUTPUT_DIR/${SAMPLE_NAME}_R2_unpaired.fastq.gz"

    # Run Trimmomatic for the sample
    echo "Processing: $SAMPLE_NAME"
    java -jar "$TRIMMOJAR" PE -threads 8 \
        "$R1" "$R2" \
        "$P1" "$U1" "$P2" "$U2" \
        LEADING:30 TRAILING:30 \
        SLIDINGWINDOW:4:30

    echo "Finished: $SAMPLE_NAME"
done

echo "All files processed. Trimmed reads are in $OUTPUT_DIR."
