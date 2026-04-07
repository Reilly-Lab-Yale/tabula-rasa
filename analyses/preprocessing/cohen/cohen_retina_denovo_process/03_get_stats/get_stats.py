# get_stats.py
import sys
from collections import Counter
import fileinput

def main():
    total_lines=0
    counts = None
    for lineno, line in enumerate(fileinput.input(files=sys.argv[1:] or ('-',)), start=1):
        line = line.strip()
        #skip empty lines
        if not line:
            continue

        total_lines +=1

        #init code
        if counts is None:
            counts = [Counter() for _ in range(len(line))]
        elif len(line) != len(counts):
            raise ValueError(f"different read length at line {lineno}: got {len(line)}, expected {len(counts)}")

        #increment counters
        for i, nt in enumerate(line):
            counts[i][nt] += 1

    if not counts:
        return
    
    for c in counts:
        #cast to dict from counter, change to percent from abolute count. 
        c=dict(c)
        for nt in c:
            c[nt]=c[nt]/float(total_lines)
        print(c)

if __name__ == "__main__":
    main()