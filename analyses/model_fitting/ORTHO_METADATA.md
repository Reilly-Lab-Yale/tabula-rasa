# Ortho metadata

Every fitted ortho carries an `ortho_meta.json` sidecar describing what was
modelled and how it was computed.

Fits write their own record: `fit_by_cre_models` and `fit_by_cell_type_models`
classify their settings, and `ortho.save()` writes the result. `ortho.load()`
reads it back onto `ortho.meta`. The schema, its classification rule and its
validation live in `scMPRAforge/ortho_meta.py`, so there is one definition.

`backfill_ortho_meta.py` reconstructs the record for the 19 orthos fit before
the sidecar existed, and is the only thing that sets `backfilled: true`. It
also fills the two fields no fit can know about itself --- `canonical` and
`source_has_zero_rows` --- on records fits wrote, verifying rather than
overwriting everything else those records state. Re-running it is safe: it
never downgrades a first-hand record to a reconstruction, and a disagreement
between a fit's account of itself and its artifacts is an assertion failure
rather than a silent overwrite.

## Why a sidecar, and why these fields

`fit_mode` conflates two different things: which unobserved combinations became
zeros (a modelling choice) and whether those zeros were rows in the input or
design-row weights built at fit time (an implementation detail). A fit whose
zeros arrived in the input reads as `standard`, so its modelling choice
survives only in the directory name. The schema separates them.

A sidecar rather than an addition to the design dicts: back-fillable without
touching fitted artifacts, readable without unpickling, and unaffected by the
two design-dict schemas having diverged.

## Schema

```json
{
  "zero_expansion": "preexisting|per_delivery|per_barcode|all_combinations|all_combinations_moi",
  "expansion_stage": "data_prep|fit_time",
  "model_family": "nb|zinb",
  "reporter_resolution": "barcode|element|none",
  "reporter_source": "in_table|separate_table|absent",
  "zero_storage": "materialized|phantom_compressed",
  "dataset": "shendure|cohen|seelig|takeshi",
  "source_table": "<count table the fit read>",
  "source_has_zero_rows": true,
  "canonical": true,
  "counterfactual": true,
  "fit_script": "analyses/model_fitting/fits/<name>/fit.py",
  "fitted_at": "2026-04-03",
  "fitted_at_basis": "run_stats filename",
  "code_version": null,
  "backfilled": true,
  "backfill_basis": "..."
}
```

`zero_expansion` values: `preexisting` = zeros already rows in the input, the
fit adds none; `per_delivery` = one zero per (cell, element) detection;
`per_barcode` = one zero per barcode of each detected element;
`all_combinations` = every observable combination; `all_combinations_moi` =
same, each weighted by P(transfected).

`counterfactual` appears only when set. `canonical` and `source_has_zero_rows`
are null on a record a fit has just written and no metadata pass has seen yet:
null reads as "not yet determined", not as "no". `code_version` names the code that ran the
fit, in a self-describing form: `git:<sha>` or `git:<sha>-dirty` from a
checkout, `version:<v>` from an installed package. It is null on every
back-filled record --- none of those fits recorded the repo state they ran
against, and a sha inferred from the date would be a guess dressed as
provenance.
`backfill_basis` appears only on reconstructed records.

## Derivation

Every field comes from an artifact, never typed in:

- the fit script, for the call that produced the ortho (`set_coarse_reporter`,
  `set_consider_missing`, `phantom_compress`, `reporter_expansion`, `nb_only`,
  `moi_correct_cm`, and the source table `from_tsv`/`from_parquet` read);
- `core.py`'s `<DS>_BOUNDS = <ALIAS>_BOUNDS` aliases, for which fit each
  dataset treats as canonical;
- the design dicts and the source count table, as independent cross-checks;
- the `run_stats_*.txt` filename for the date, falling back to the ortho
  directory mtime with `fitted_at_basis` recording which was used.

Cross-checks that must agree: the directory name's `_nb_`/`_zinb_` against
`nb_only` in the fit script; the design dict's `fit_mode` and
`reporter_expansion` against the derived expansion; and, for `preexisting`,
that the source count table actually contains zero rows.

The last is one-directional. A table carrying zeros can still be refit under
consider-missing, which is how the reporter-free counterfactuals for Lalanne
et al. are constructed.

## Validation rules

Asserted on write and on read:

- `reporter_resolution == "none"` implies `zero_expansion` in
  {`all_combinations`, `all_combinations_moi`} --- consider-missing
  definitionally disregards reporter information, so it never carries a
  reporter expansion.
- `zero_expansion == "preexisting"` implies `expansion_stage == "data_prep"`
  and `zero_storage == "materialized"`.
- `zero_storage == "phantom_compressed"` implies `expansion_stage == "fit_time"`.
- at most one `canonical: true` per dataset; a dataset with none is reported.
- `zero_expansion == "per_barcode"` with `reporter_resolution == "element"` is
  legal but marked `counterfactual` --- this is the Zhao et al. case Results
  2.1 argues against, and it should be visible rather than silent.

## The 19 orthos

Written 2026-08-14. `*` canonical, `!` counterfactual.

| ortho | zero_expansion | stage | storage | reporter |
|---|---|---|---|---|
| `shendure_obs_nb` * | preexisting | data_prep | materialized | barcode / in_table |
| `shendure_obs_zinb` | preexisting | data_prep | materialized | barcode / in_table |
| `shendure_cm_{nb,zinb}_phantom` | all_combinations | fit_time | phantom | none / absent |
| `shendure_cm_moib_nb_phantom` | all_combinations_moi | fit_time | phantom | none / absent |
| `cohen_obsingle_nb_phantom` * | per_delivery | fit_time | phantom | element / separate_table |
| `cohen_obsingle_zinb_phantom` | per_delivery | fit_time | phantom | element / separate_table |
| `cohen_obs_{nb,zinb}_phantom_20260401` ! | per_barcode | fit_time | phantom | element / separate_table |
| `cohen_cm_{nb,zinb}_phantom_20260401` | all_combinations | fit_time | phantom | none / absent |
| `seelig_cm_moib_nb_phantom` * | all_combinations_moi | fit_time | phantom | none / absent |
| `seelig_cm_moib_zinb_phantom` | all_combinations_moi | fit_time | phantom | none / absent |
| `seelig_cm_{nb,zinb}_phantom` | all_combinations | fit_time | phantom | none / absent |
| `takeshi_obs_{nb,zinb}_phantom` | preexisting | data_prep | materialized | barcode / in_table |
| `takeshi_cm_{nb,zinb}_phantom` | all_combinations | fit_time | phantom | none / absent |

Takeshi has no canonical fit: it is held for v2 and `TAKESHI_BOUNDS` is
commented out in `core.py`.

Bouchet holds all 19. The local `/nfs` tree is a partial copy: the four
`cohen_*_20260401` orthos and all of takeshi are dangling symlinks there, and
the script reports them as skipped rather than passing over them silently.

## Facts worth not re-deriving

- **Lalanne (shendure) needs no reporter table.** Its transfection reporter is
  barcode-level and its zeros are already rows in `shendure_processed.tsv` ---
  748,777 of 877,064 (85.4%). The canonical fit calls `fit_by_cre_models` /
  `fit_by_cell_type_models` with no `phantom_compress`, no
  `set_coarse_reporter`, no consider-missing. "obs, fine" for shendure means
  *do nothing*, and that is correct. Takeshi's oBC design is the same shape.
- **Do not re-fit shendure through the phantom path.** `set_coarse_reporter`
  takes only `(rep_id, cell_bc, cre_id)` and dedupes to one row per
  (rep, cell, CRE); it has no barcode column, so it cannot represent what the
  oBC recorded. Neither expansion reproduces the real zero set: `single` gives
  778,248 rows against the data's 877,064 (13% fewer, because barcodes per
  detection is median 1 / mean 1.13 / max 11), and `per_barcode` would invent
  ~106M rows (~121x). The materialised fit uses better information than either.
- **Phantom compression is an optimisation, not a formula change.** It must not
  move fit values; it only stores zeros as design-row weights.
- Enabling consider-missing takes the phantom-compressed CM path in
  `standard_fit` regardless of `phantom_compress`, so every CM ortho is
  `fit_time` / `phantom_compressed`.
- `reporter_expansion="coarse"` produces the per-barcode rule; `"single"` one
  zero per delivery. Added 2026-04-06 in `248f5df` ("+single zero for coarse"),
  with `"coarse"` made the default to preserve the prior behaviour. `fit_mode`
  arrived 2026-04-05 (`1c3fb7b`), `consider_missing` 2026-03-10 (`e1d0ed1`).
  `shendure_obs_nb` was fit 2026-04-03, before any of the tags --- which
  is why its design dicts carry `fit_mode=None`.
- `fit_mode` cannot distinguish `obs` from `obsingle` (both `obs_phantom`);
  only `reporter_expansion` does, and it is absent from the CM design dicts.
- `core.py:3580-3587` infers `standard` / `obs_phantom` / `cm_phantom` for
  untagged orthos from training-data properties. The inference is correct; it
  just describes implementation, not the modelling choice.

## Naming

An ortho's name carries its dataset, its expansion and its count family.
`_phantom` records that the zeros were compressed into design-row weights, so
it belongs only on orthos whose `zero_storage` is `phantom_compressed`. The
four fits whose zeros arrived in the input --- `shendure_obs_{nb,zinb}` and
`takeshi_obs_{nb,zinb}` --- carry no such suffix.
