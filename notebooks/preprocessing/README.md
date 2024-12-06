These scripts preprocess scMPRA from a couple sources into convienient eriniform tables. 

All the actual data are in `/home/mcn26/palmer_scratch/tabula_data/formatted`
The notebooks here are written to download and process all the files automatically (just set download_data to TRUE). 

Tables have columns:

cellBC <str> cell barcode
rep_id <str> replicate id
CRE_id <str> name of the CRE
cell_type_annotation <str> cell-type

MPRA_BC <string> nucleotide barcode of the MPRA reporter barcode (called mBC in the shendure paper and rBC in the cohen paper)
reads_MPRA_BC <int>
UMIs_MPRA_BC <int>

Optional:

transfection_BC <str> nucleotide barcode of the transfection reporter (called the oBC in the shendure paper, and alternately the cBC or pBC in the cohen paper)
reads_transfection_BC <int>
UMIs_transfection_BC <int>


# Table-specific notes:

in COHEN_RETINA.tsv
- MPRA_BC="0" I think indicates untransfected but pBC (transfection reporter) detected
- Consequently, I set MPRA_Bc and UMIs_MPRA_BC for these lines to NA. 
- To zero deflate untransfected, simply filter out these lines. 