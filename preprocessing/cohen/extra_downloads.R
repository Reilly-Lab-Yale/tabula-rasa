
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188639

[GSE188639_retina_pBC_exp_rep1.csv.gz](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE188639&format=file&file=GSE188639%5Fretina%5FpBC%5Fexp%5Frep1%2Ecsv%2Egz)
[GSE188639_retina_pBC_exp_rep2.csv.gz](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE188639&format=file&file=GSE188639%5Fretina%5FpBC%5Fexp%5Frep2%2Ecsv%2Egz)

GSM6614199 	retina single-cell transcriptome rep1 [scRNA-seq]
GSM6614200 	retina single-cell transcriptome rep2 [scRNA-seq]
GSM6614201 	retina barcodes rep1 [scMPPRA]
GSM6614202 	retina barcodes rep2 [scMPPRA]
GSM6614203 	retina u6 barcodes rep1 [scMPPRA]
GSM6614204 	retina u6 barcodes rep2 [scMPPRA]



```{r}
rep1_path=paste(data_root,"/GSE188639_retina_pBC_exp_rep1.csv.gz",sep="")
rep2_path=paste(data_root,"/GSE188639_retina_pBC_exp_rep2.csv.gz",sep="")

if(download_data){
  download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE188639&format=file&file=GSE188639%5Fretina%5FpBC%5Fexp%5Frep1%2Ecsv%2Egz",rep1_path)
  download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE188639&format=file&file=GSE188639%5Fretina%5FpBC%5Fexp%5Frep2%2Ecsv%2Egz",rep2_path)
}
```

```{r}
rep1=read_csv(rep1_path)
rep2=read_csv(rep2_path)
```
```{r}
rep2
```

"pBC" I think stands for "plasmid barcode". IIRC, this is a synomym of cBC. 
We can check if this is the case by comparing the number of unique pBc to
promoter names. 

```{r}
distinct_names<-rep1 %>% distinct(name)
distinct_pBC<-rep1 %>% distinct(pBC)

stopifnot(nrow(distinct_pBC)==nrow(distinct_names))
```

Ok, it would seem so...

What is that's funny number in the first column..? Doesn't look like an index?
  
  I think these data are plasmid barcode counts. I say this because there's no 
MPRABC info, and because:
GSM6614203 	retina u6 barcodes rep1 [scMPPRA]
and
GSM6614204 	retina u6 barcodes rep2 [scMPPRA]
say
"Processed data are available on Series record".


Let's download
GSM6614199 	retina single-cell transcriptome rep1 [scRNA-seq]

[GSM6614199_retina_rep1_barcodes.tsv.gz](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6614199&format=file&file=GSM6614199%5Fretina%5Frep1%5Fbarcodes%2Etsv%2Egz)
[GSM6614199_retina_rep1_features.tsv.gz](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6614199&format=file&file=GSM6614199%5Fretina%5Frep1%5Ffeatures%2Etsv%2Egz)
[GSM6614199_retina_rep1_matrix.mtx.gz](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6614199&format=file&file=GSM6614199%5Fretina%5Frep1%5Fmatrix%2Emtx%2Egz)


GSM6614200 	retina single-cell transcriptome rep2 [scRNA-seq]
...

```{r}
retina_rep1_barcodes=paste(data_root,"/GSM6614199_retina_rep1_barcodes.tsv.gz",sep="")
retina_rep1_features=paste(data_root,"/GSM6614199_retina_rep1_features.tsv.gz",sep="")
retina_rep1_matrix=paste(data_root,"/GSM6614199_retina_rep1_matrix.mtx.gz",sep="")

if(download_data){
  download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6614199&format=file&file=GSM6614199%5Fretina%5Frep1%5Fbarcodes%2Etsv%2Egz",retina_rep1_barcodes)
  download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6614199&format=file&file=GSM6614199%5Fretina%5Frep1%5Ffeatures%2Etsv%2Egz",retina_rep1_features)
  download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6614199&format=file&file=GSM6614199%5Fretina%5Frep1%5Fmatrix%2Emtx%2Egz",retina_rep1_matrix)
}
```

```{r}
retina_rep1 <- ReadMtx(
  mtx=retina_rep1_matrix,
  features=retina_rep1_features,
  cells=retina_rep1_barcodes
)
```
Let's check out some of the features:
```{r}
head(dimnames(retina_rep1)[[1]], 10)
```

Clearly endogenous transcriptome. Does it have barcode information?

```{r}
system(paste("zcat",retina_rep1_features,"|","cut -f3 | uniq"))
```
Looks like all are labeled "gene expression". 

```{r}
system(paste("zcat",retina_rep1_features,"| cut -f1 | grep -v \"ENSMU\" "))
```
Other than EGFP and DSred all are mouse genes. So this is totally just gene expression. 

Unfortunately, 
GSM6614201 	retina barcodes rep1 [scMPPRA]
and
GSM6614202 	retina barcodes rep2 [scMPPRA]
both have no suppl files. 
They do *say* "Processed data provided as supplementary file". So let's look...



---
  
  Data could be re-constructed from the zenodo, but it looks like a better (more 
                                                                            processed) file is avilable from here:
  
  [GSE188639_mixed_cell_bulk_exp.csv.gz](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE188639&format=file&file=GSE188639%5Fmixed%5Fcell%5Fbulk%5Fexp%2Ecsv%2Egz)

```{r}
mixed_cell_path=paste(data_root,"GSE188639_mixed_cell_bulk_exp.csv.gz",sep="/")
if(download_data){
  download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE188639&format=file&file=GSE188639%5Fmixed%5Fcell%5Fbulk%5Fexp%2Ecsv%2Egz",mixed_cell_path)
}
```

```{r}
mixed_cell_data=read_csv(mixed_cell_path)
```



```{r}
mixed_cell_data
```

ok actually that's bulked and so useless'