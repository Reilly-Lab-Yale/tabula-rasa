# Coarse transfection reporter table spec (PROPOSED)

For datasets with a CRE-level transfection reporter (e.g. Cohen U6): the
reporter confirms CRE presence in a cell, so unobserved barcodes of that CRE
are true zeros. These zeros are stored in a separate reporter table and
phantom-compressed at fit time, keeping the main observation table nonzero-only.

## Table format

One row per positive (cell, CRE) detection:

| Column | Type | Description |
|--------|------|-------------|
| rep_id | string | Replicate identifier |
| cell_bc | string | Cell barcode |
| cre_id | string | CRE detected by reporter |

Optional: `umis_transfection_bc` (reporter UMI count, for QC).

Constraints: positive detections only, (rep_id, cell_bc, cre_id) unique,
identifiers match main observation table.

## Zero-expansion modes

```
observation table (nonzero only)
    +-- no reporter, no CM    --> fit on nonzero data only
    +-- coarse reporter       --> obs: phantom zeros for unobserved
    |                              barcodes of reporter-detected CREs
    +-- consider_missing      --> CM: phantom zeros for full Cartesian
                                   (reporter zeros are a strict subset)
```
