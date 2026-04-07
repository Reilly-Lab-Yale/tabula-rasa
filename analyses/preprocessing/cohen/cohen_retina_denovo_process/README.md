This directory contains code to download and preprocess cohen's ([A single-cell massively parallel reporter assay detects cell-type-specific gene regulation](https://www.nature.com/articles/s41588-022-01278-7)) ex-vivo retina expreiment into a properly formatted scMPRA data table. 
We do this because the provided preprocessed tables miss information we need. We cannot use the provided code on github as it is very incomplete. It is missing essential steps like U6 data processing and it has no documentation.

Files we process:

#GSM6614201 	retina barcodes rep1 [scMPPRA]
# SRR21787462
#GSM6614202 	retina barcodes rep2 [scMPPRA]
# SRR32774353

#GSM6614203 	retina u6 barcodes rep1 [scMPPRA]
# SRR21787461
#GSM6614204 	retina u6 barcodes rep2 [scMPPRA]
# SRR21787460

We will not process 
#GSM6614199 	retina single-cell transcriptome rep1 [scRNA-seq]
#GSM6614200 	retina single-cell transcriptome rep2 [scRNA-seq]
Since the preprocessed data are sufficient for our purposes here. 

~~"Reverse-engineer cohen scmpra" in obsidian~~ outdated