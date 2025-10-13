Let's examine the consensus sequences we just generated from the cohen retina scMPRA sequencing data and see if we can understand them and so extract the relevant barcodes...

# reference materials

Recall that 'pBC' and 'cBC' are interchangable 

![fig4ab](../images/fig4ab.png)
Chromium diagram:
![chromium_diagram](../images/chromium.png)

TSOs from [10x docs](https://kb.10xgenomics.com/hc/en-us/articles/360001493051-What-is-a-template-switch-oligo-TSO). 
```
CCCATGTACTCTGCGTTGATACCACTGCTT
```

```
TTTCTTATATGGG
```



https://github.com/szhao045/scMPRA_parsing/blob/main/readprocessor.go
```go
	cellBC := read1[0:16]
	// get the UMI
	umi := read1[16:28]
	// Check the Q30 score of the read
	output.CellBC = cellBC
	output.Umi = umi
	// Fuzzy match for thr tripBC
	beforePBC := "AGTAACTGCGAT"
	afterPBC := "AAGGAACCCG"
	beforeRBC := "TACCGGTATCGC"
	afterRBC := "GGCCGCTAAG"
		coordinates, err := FuzzyMatch(beforePBC, afterPBC, beforeRBC, afterRBC, read2)
	if err != nil {
		output := Trios{}
		return output, err
	}
	// et the TBC
	tbc := read2[coordinates[0]:coordinates[1]]
	rbc := read2[coordinates[2]:coordinates[3]]
	output.RBC = rbc
	output.TBC = tbc
	// create the output struct
	return output, nil
```

# Normal MPRA

Examining the crispr screening regular dual index library (top diagram) & fig 4b left, I expect the read structure to be

```
5' [10x cbc] [UMI] [poly dT VN] [rBC] [cBC] 3'
```

**All forward reads:** 
- SRR21787462_1.filtered.fastq_28nt.txt.stats
- SRR32774353_1.filtered.fastq_28nt.txt.stats
```
NNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNN
```

So that's `[10x cbc] [UMI]`, and it's exactly 28 bp (16 bp 10x Barcode + 10 bp UMI, as per docs), as expected!

**All reverse reads:**
- SRR32774353_2.filtered.fastq_100nt.txt.stats
- SRR21787462_2.filtered.fastq_100nt.txt.stats

```
GAGCTGTACAAGTAACTGCGATMNNNNNNNAAGGAACCCGCGCTATACCGGTATCGCNNNNNNNNNNNNNNNNNNNNNNNNGGCCGCTAAGATACATTGA
GAGCTGTACAAGTAACTGCGATMNNNNNNNAAGGAACCCGCGCTATACCGGTATCGCNNNNNNNNNNNNNNNNNNNNNNNNGGCCGCTAAGATACATTGA
```

Let's break it down, annotating with the sequences from the GO snippet above. 
```
GAGCTGTACAAGTAACTGCGATMNNNNNNNAAGGAACCCGCGCTATACCGGTATCGCNNNNNNNNNNNNNNNNNNNNNNNNGGCCGCTAAGATACATTGA
     ?    |before pBC|pBC    |after pBC|  ?  |before rbc| rBC                    |afterrbc|    ?
```

8bp cBC, 24bp rBC. Makes sense. 

# U6

Examining the crispr screening dual index library (bottom diagram) & fig 4b right, I expect the read structure to be

```
5' [10x cbc] [UMI] [Capture sequence] [cBC] [feature barcode]  [tso] 3'
```

**All forward reads:**
- SRR21787460_1.filtered.fastq_28nt.txt.stats
- SRR21787461_1.filtered.fastq_28nt.txt.stats

```
NNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNN
```

So that's `[10x cbc] [UMI]`, and it's exactly 28 bp (16 bp 10x Barcode + 10 bp UMI, as per docs), as expected!

**All reverse reads:**
- SRR21787461_2.filtered.fastq_44nt.txt.stats
- SRR21787461_2.filtered.fastq_100nt.txt.stats
- SRR21787460_2.filtered.fastq_100nt.txt.stats
- SRR21787460_2.filtered.fastq_44nt.txt.stats
- SRR21787461_2.filtered.fastq_150nt.txt.stats
- SRR21787460_2.filtered.fastq_150nt.txt.stats

```
CCGGTAAGCTCCCGGGAGCTTGTMNNNNNNNGCTTTAAGGCCGG
CCGGTAAGCTCCCGGGAGCTTGTMNNNNNNNGCTTTAAGGCCGGTCCTAGCAANNNNNNNNNNNNNNNNNNNNNNNNNNNNCTGTCTCTTATACACATCT
CCGGTAAGCTCCCGGGAGCTTGTMNNNNNNNGCTTTAAGGCCGGTCCTAGCAANNNNNNNNNNNNNNNNNNNNNNNNNNNNCTGTCTCTTATACACATCT
CCGGTAAGCTCCCGGGAGCTTGTMNNNNNNNGCTTTAAGGCCGG
CCGGTAAGCTCCCGGGAGCTTGTMNNNNNNNGCTTTAAGGCCGGTCCTAGCAANNNNNNNNNNNNNNNNNNNNNNNNNNNNCTGTCTCTTATACACATCTGACGCTGCCGACGATATAATCCGAGTGTAGATCTCGGTGGTCGCCGTATC
CCGGTAAGCTCCCGGGAGCTTGTMNNNNNNNGCTTTAAGGCCGGTCCTAGCAANNNNNNNNNNNNNNNNNNNNNNNNNNNNCTGTCTCTTATACACATCTGACGCTGCCGACGATTGGCACTCGGTGTAGATCTCGGTGGTCGCCGTATC
```

Interesting. If we line them all up they match. Why exactly were there multiple read lengths ...?

Let's break down the longest one:
```
CCGGTAAGCTCCCGGGAGCTTGTMNNNNNNNGCTTTAAGGCCGGTCCTAGCAANNNNNNNNNNNNNNNNNNNNNNNNNNNNCTGTCTCTTATACACATCTGACGCTGCCGACGATTGGCACTCGGTGTAGATCTCGGTGGTCGCCGTATC
         ?             |8mer  |         C1           |    suspicious 28mer      |     <-----nextera R1          |           | P5
```

This is really really interesting. Since we are running off the ends of the reads, let's add the illumina portion to our predicted transcript

```
5' [p5] [i5] [nextera R1] [10x cbc] [UMI] [Capture sequence] [cBC] [feature barcode] [tso] [truseq r2] [i7] [p7] 3'
```


Ok, so the suspicious 28mer is probably the `[10x cbc]+[UMI]`. The 8mer is probably the cBC. Is the sequence in between the capture sequence?
Yes. After consulting [Guide RNA Specifications Compatible with Feature Barcoding technology for CRISPR Screening : CG000197 RevA](https://cdn.10xgenomics.com/image/upload/v1660261286/support-documents/CG000197_GuideRNA_SpecificationsCompatible_withFeatureBarcodingtechnology_forCRISPRScreening_Rev-A.pdf) it's perfect revcomp of C1. 

The ? sequence is some junk between the u6 promoter and the 8mer payload. 

We COULD probably flash some proportion of the reads together, but we can't for all. 