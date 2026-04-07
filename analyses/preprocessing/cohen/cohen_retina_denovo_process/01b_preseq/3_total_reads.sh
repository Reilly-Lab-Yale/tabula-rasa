#!/bin/bash
# computes total reads

OUTFILE="total_reads.tsv"
echo -e "filename\tsum" > "$OUTFILE"

for f in /home/mcn26/palmer_scratch/raw_recap/cohen_retina/uniq_counts/*.txt; do
    if [[ -f "$f" ]]; then
        fname=$(basename "$f")
        sum=$(awk -F'\t' '{s+=$1} END {print s}' "$f")
        echo -e "${fname}\t${sum:-0}" >> "$OUTFILE"
    fi
done