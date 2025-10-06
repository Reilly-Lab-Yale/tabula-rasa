#run with pypy3 and pipe file input

import sys

def main():
    line_number=-1
    counts=[]
    init=True
    for line in sys.stdin.buffer:
        #fastq format is 
        #0: header
        #1: nucleotides
        #2: header
        #3: quality
        #repeat
        #we only want nucleotides

        line_number=line_number+1
        
        if line_number % 4 !=1:
            continue
        
        line=line.decode()
        #print(f"Line: {line}")
        
        #initalize length. Assuming all reads are the same length.
        if init:
            counts=[{}]*len(line)
            init=False
        elif len(line)!=len(counts):
            raise ValueError(f"different read length at line {line_number}")
        
        for idx, nt in enumerate(line):
            counts[idx][nt] = counts[idx].get(nt,0)+1
        
        
    
    for i in counts:
        print(i)
		
		
if __name__=="__main__":
	main()