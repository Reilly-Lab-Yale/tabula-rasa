#template for python scripts which process files line-by-line
#run with pypy3 and pipe file input


import sys

def main():
	filename = sys.argv[1]
	output_root = sys.argv[2]
	print(filename)
	print(output_root)
	line_number=-1
	output_files={}
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

		line=line.decode().strip()

		if len(line) not in output_files.keys():
			file=open(f"{output_root}/{filename}_{len(line)}nt.txt","w")
			output_files[len(line)]=file
		
		output_files[len(line)].write(line+"\n")
	
	for key in output_files.keys():
		output_files[key].close()

		
if __name__=="__main__":
	main()