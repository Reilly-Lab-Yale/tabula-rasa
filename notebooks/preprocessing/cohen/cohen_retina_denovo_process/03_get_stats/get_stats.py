#run with pypy3 and pipe file input

import sys

def main():
    counts=[]
    init=True
    for line in sys.stdin.buffer:
        
        line=line.decode()
        #print(f"Line: {line}")
        
        #initalize length. Assuming all reads are the same length.
        if init:
            counts=[{}]*len(line)
            init=False
        elif len(line)!=len(counts):
            raise ValueError(f"different read length")
        
        for idx, nt in enumerate(line):
            counts[idx][nt] = counts[idx].get(nt,0)+1
    
    for i in counts:
        print(i)
		
		
if __name__=="__main__":
	main()