import sys
from Bio import SeqIO
import numpy as np

def calculate_nucleotide_percentage(fastq_stream):
    # Initialize a list to hold the count of nucleotides at each position
    count = []
    total_sequences = 0

    # Process each record from the stream
    for record in SeqIO.parse(fastq_stream, "fastq"):
        total_sequences += 1
        seq = str(record.seq)

        # Extend count list size if the current sequence is longer than count list
        if len(count) < len(seq):
            count.extend([{'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0} for _ in range(len(seq) - len(count))])

        # Count each nucleotide at each position
        for i, nucleotide in enumerate(seq):
            if nucleotide in count[i]:
                count[i][nucleotide] += 1

    # Calculate percentages
    percentages = []
    for pos_counts in count:
        pos_percentages = {nuc: (pos_counts[nuc] / total_sequences * 100) for nuc in pos_counts}
        percentages.append(pos_percentages)

    return percentages

# Usage: Pipe in a FASTQ file via standard input
# Example: cat example.fastq | python this_script.py
percentages = calculate_nucleotide_percentage(sys.stdin)

# Print or process the percentages as needed
for i, p in enumerate(percentages):
    print(f"Position {i+1}: {p}")
