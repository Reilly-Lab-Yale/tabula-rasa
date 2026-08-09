**MPRAbase v4.9.3** -- SQLite, 400 MB, CC-BY-4.0.
Zenodo record 10920747, DOI 10.5281/zenodo.10920747.

    curl -L -o analysis/data/mprabase_v4_9.3.db \
      https://zenodo.org/api/records/10920747/files/mprabase_v4_9.3.db/content

Tables: `sample` (130 experiments), `datasets` (51 studies), `designed_library`,
`library_sequence`, `element_score`, `element_rep_score`.

### Joining samples to studies

The documented path is `sample.library_id` -> `designed_library.datasets_id` ->
`datasets`, but `designed_library` is missing two libraries (`DS0030-LID02`,
`DS0135-LID01`), so that join silently drops two experiments. Both studies exist
in `datasets`.

`sample_id` encodes the study as its first six characters (`DS0030-SID02` ->
`DS0030`), which recovers all 130. The two routes agree on all 128 rows where
both resolve; `cell_types_per_study.py` asserts this before using the prefix.