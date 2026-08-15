# Ortho metadata schema — implementation spec

Written 2026-08-14, before implementation. Approved in principle; build it as
described unless the user redirects.

## Why

Everything confusing about the current metadata comes from one conflation:
**what was modelled** versus **how it was computed**. `fit_mode` mixes them,
so `standard` reads as "no zero expansion" when it actually means "the
expansion happened upstream, in data prep". Nothing currently records *where*
the zeros came from, so for a materialised fit the modelling choice survives
only in the directory name.

## Facts established (do not re-derive)

- **Lalanne (shendure) needs no reporter table.** Its transfection reporter is
  barcode-level and its zeros are already rows in
  `shared/tabula_data_new/shendure/shendure_processed.tsv` — 748,777 zero rows
  of 877,064 (85.4%). The canonical fit calls `fit_by_cre_models` /
  `fit_by_cell_type_models` with **no** `phantom_compress`, **no**
  `set_coarse_reporter`, **no** consider-missing. "obs / fine" for shendure
  means *do nothing*. This is correct and intended.
- **Do NOT re-fit shendure through the phantom path.** `set_coarse_reporter`
  takes only `(rep_id, cell_bc, cre_id)` and dedupes to one row per
  (rep, cell, CRE) — it has no barcode column, so it cannot represent what the
  oBC actually recorded. Neither expansion reproduces the real zero set:
  `single` gives 778,248 rows against the data's 877,064 (13% fewer, because
  barcodes per detection is median 1 / mean 1.13 / max 11), and `per_barcode`
  would invent ~106M rows (~121x). The materialised fit uses better information
  than either.
- **Zhao (cohen)** attaches `unjoined/u6.tsv` via `set_coarse_reporter` and
  fits with `phantom_compress=True, reporter_expansion="single"`.
- **Phantom compression is an optimisation, not a formula change.** It must not
  move fit values; it only stores zeros as design-row weights instead of
  materialised rows.
- `reporter_expansion="coarse"` produces the per-barcode rule; `"single"`
  produces one zero per delivery. Added 2026-04-06 in `248f5df`
  ("+single zero for coarse"); before that, per-barcode was the only behaviour
  and `"coarse"` was made the default to preserve it. `fit_mode` arrived
  2026-04-05 (`1c3fb7b`); `consider_missing` 2026-03-10 (`e1d0ed1`).
  `shendure_obs_nb_phantom` was fit **2026-04-03**, before any of the tags.
- `fit_mode` cannot distinguish `obs` from `obsingle` (both `obs_phantom`);
  only `reporter_expansion` does, and it is absent from the CM design dicts.
- `core.py:3580-3587` already infers `standard` / `obs_phantom` / `cm_phantom`
  for untagged orthos from training-data properties. That inference is correct;
  it just describes implementation, not the modelling choice.

## Schema: `ortho_meta.json`, one per ortho directory

Sidecar, not an addition to the design dicts: back-fillable without touching
fitted artifacts, readable without unpickling, and unaffected by the two
design-dict schemas having diverged.

```json
{
  "zero_expansion": "preexisting|per_delivery|per_barcode|all_combinations|all_combinations_moi",
  "expansion_stage": "data_prep|fit_time",
  "model_family": "nb|zinb",
  "reporter_resolution": "barcode|element|none",
  "reporter_source": "in_table|separate_table|absent",
  "zero_storage": "materialized|phantom_compressed",
  "dataset": "shendure|cohen|seelig|takeshi",
  "canonical": true,
  "fit_script": "analyses/model_fitting/fits/<name>/fit.py",
  "fitted_at": "2026-04-03",
  "code_version": "<git sha at fit time, or null>",
  "backfilled": true,
  "backfill_basis": "derived from fit.py + design dicts + run_stats"
}
```

Meanings: `preexisting` = zeros already rows in the input, fit adds none;
`per_delivery` = one zero per (cell, element) detection; `per_barcode` = one
zero per barcode of each detected element; `all_combinations` = every
observable combination; `all_combinations_moi` = same, each weighted by
P(transfected).

## Validation rules (assert on write and on read)

- `reporter_resolution == "none"` implies `zero_expansion` in
  {`all_combinations`, `all_combinations_moi`} — consider-missing definitionally
  disregards reporter information, so it never carries a reporter expansion.
- `zero_expansion == "preexisting"` implies `expansion_stage == "data_prep"`
  and `zero_storage == "materialized"`.
- `zero_storage == "phantom_compressed"` implies `expansion_stage == "fit_time"`.
- `zero_expansion == "per_barcode"` with `reporter_resolution == "element"` is
  legal but flagged **counterfactual** — this is the Zhao case Results 2.1
  argues against, and it should be visible rather than silent.
- exactly one `canonical: true` per dataset.

## Back-fill

19 orthos under `/nfs/roberts/project/pi_skr2/shared/tabula_data_new/<ds>/`,
each with a `by_cell_type_design/` directory. Fit scripts are at
`/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/analyses/model_fitting/fits/<name>/fit.py`.

**Derive, never type.** Read `fit.py` for `set_coarse_reporter` /
`phantom_compress` / `reporter_expansion`, the design dicts for `fit_mode` /
`reporter_expansion`, and `run_stats_*.txt` for the date. Assert the derived
record satisfies the rules; refuse to write anything underdetermined and print
what was missing.

| ortho | zero_expansion | stage | storage | reporter |
|---|---|---|---|---|
| shendure_obs_nb_phantom **(canonical)** | preexisting | data_prep | materialized | barcode / in_table |
| shendure_obs_zinb_phantom | preexisting | data_prep | materialized | barcode / in_table |
| shendure_cm_{nb,zinb}_phantom | all_combinations | fit_time | phantom | n/a |
| shendure_cm_moib_nb_phantom | all_combinations_moi | fit_time | phantom | n/a |
| cohen_obsingle_nb_phantom **(canonical)** | per_delivery | fit_time | phantom | element / separate_table |
| cohen_obsingle_zinb_phantom | per_delivery | fit_time | phantom | element / separate_table |
| cohen_obs_{nb,zinb}_phantom_20260401 | per_barcode *(counterfactual)* | fit_time | phantom | element / separate_table |
| cohen_cm_{nb,zinb}_phantom_20260401 | all_combinations | fit_time | phantom | n/a |
| seelig_cm_moib_nb_phantom **(canonical)** | all_combinations_moi | fit_time | phantom | none / absent |
| seelig_cm_moib_zinb_phantom | all_combinations_moi | fit_time | phantom | none / absent |
| seelig_cm_{nb,zinb}_phantom | all_combinations | fit_time | phantom | none / absent |
| takeshi_obs_{nb,zinb}_phantom | **verify** | fit_time? | ? | **verify** |
| takeshi_cm_{nb,zinb}_phantom | all_combinations | fit_time | phantom | n/a |

The four takeshi entries are unverified — their fit scripts have not been read.
Takeshi is held for v2 and `TAKESHI_BOUNDS` is commented out in `core.py`.

## Consequences to handle after back-fill

- `shendure_obs_nb_phantom` becomes provably `zero_storage: materialized`,
  making `_phantom` in its name false. Renaming is a small migration, not a
  `mv`: the name appears in the manuscript Methods as a canonical preset name,
  in `scMPRAforge/presets/*.tgz` filenames, and across analysis scripts. Cheap
  before 1.0, expensive after. User has not decided.
- Consider whether `zero_expansion` should also be written into new fits going
  forward (a `standard_fit` parameter), so `backfilled` is only ever true for
  these 19.

## Also outstanding (unrelated to this spec, do not lose)

1. **The manuscript build is broken**: `bibtex needed too many passes`, stale
   aux state from interrupted runs. Fix with
   `rm -f manuscript.aux manuscript.bbl manuscript.blg manuscript.fdb_latexmk manuscript.fls manuscript.out manuscript.log`
   then `make pdf`. This is the second time; the cause is concurrent builds.
2. **The fine/coarse terminology in the draft is currently WRONG** and must be
   corrected. The user's ruling: *both* Lalanne and Zhao canonical fits are
   **fine** treatments. Lalanne is naturally fine (fine reporter, zeros already
   in the data). Zhao has a coarse reporter, but full per-barcode expansion is
   the problematic case, so `obsingle` (one zero per delivery) treats it as
   fine. The **coarse expansion** is the per-barcode rule applied to a
   coarse reporter — the Zhao counterfactual.
   - `manuscript.tex` currently says "Under fine expansion, a detection enters
     one zero for every barcode..." — backwards, fix it.
   - `plot_nb_vs_zinb_bars.py` LAYOUT currently labels
     `"Lalanne: observed, fine"` / `"Zhao: observed, coarse"` / `"Zhao:
     observed, fine"` — also backwards; Zhao's canonical (`obsingle`) is the
     fine one. Final wording not settled with the user.
   - The Fig 1 caption gloss in `sections/figures_parts.tex` was rewritten to
     match the wrong version and needs the same correction.
3. Figure re-plot rules are in memory (`figure-replot-rules.md`); only re-plot
   when asked. Fig 1D is done. A rescale-factor worklist for every figure was
   computed — Fig 4 and S4 panels are worst (0.23-0.33).
