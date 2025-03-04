(add dna)

These scripts preprocess scMPRA from a couple sources into convienient eriniform tables. 

Uses conda envieonment "rzone" for both r and python scripts. 

All the actual data are in `/home/mcn26/palmer_scratch/tabula_data/formatted`
The notebooks here are written to download and process all the files automatically (just set download_data to TRUE and change the input/output root path variables to your own palmer_scratch). 

Tables have columns:
- cellBC <str> cell barcode
- rep_id <str> replicate id
- CRE_id <str> name of the CRE
- cell_type_annotation <str> cell-type

- MPRA_BC <string> nucleotide barcode of the MPRA reporter barcode (called mBC in the shendure paper and rBC in the cohen paper)
- reads_MPRA_BC <int>
- UMIs_MPRA_BC <int>

Optional:

- transfection_BC <str> nucleotide barcode of the transfection reporter (called the oBC in the shendure paper, and alternately the cBC or pBC in the cohen paper)
- reads_transfection_BC <int>


# Table-specific notes:

COHEN_RETINA.tsv
- I've merged in transfected unexpressed.
- reads_transfection_BC has duplicated entries : multiple MPRA bc coresp. to same transfection BC in same cell.
- So don't sum it

COHEN_MIXED_CELL.tsv
- Data were already summed across MPRA barcodes
- (I can re-run their analysis and extract this information, but would be nontrivial amount of work).