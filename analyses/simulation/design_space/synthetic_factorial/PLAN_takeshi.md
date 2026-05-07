# Plan: include takeshi as a third empirical anchor (toggleable)

## Context

Takeshi is the 4th dataset in the analysis pipeline (the others being shendure,
cohen, seelig). Its bounds preset is canonical: `TAKESHI_BOUNDS =
takeshi_obs_nb_phantom`. Per `paper_plan.md`, the model-selection table lists
takeshi alongside the others; the dataset is unpublished as of 2026-05-06.

Goal: extend the synthetic factorial figures (marginals + attribution) to
*optionally* include takeshi as a third empirical anchor (alongside cohen-Rod
and shendure-Pluripotent), with two parallel sets of outputs:

- **default / no-takeshi**: cohen + shendure overlays only. Use these in
  the methods paper.
- **with-takeshi**: cohen + shendure + takeshi overlays. Use in internal
  reports, supplemental figures, and any analysis that benefits from a
  3-anchor comparison.

Code, data, and configs related to takeshi can sit in the repo without
restriction; the only thing that matters is which **figure file** ends up
in the paper.

## Takeshi empirical axis values (extracted 2026-05-06)

Canonical CT: **`reference` (= HepG2)**, the same convention shendure/cohen
use for their default anchor.

| Axis | takeshi-HepG2 | shendure-Pluripotent | cohen-Rod | LHS+topup combined coverage | takeshi inside? |
|---|---:|---:|---:|---|---|
| n_cells          | 1690    | 8,201    | 18,633   | [500, 50,000] | yes |
| n_cres           | 149     | 207      | 116      | [50, 2,000]   | yes |
| bcs_per_cre      | 394     | 136      | 17,244   | [3, 50,000]   | yes |
| moi (effective)  | **266** | 18       | 149      | [0.5, 200]    | **NO (above)** |
| lib_alpha_nb     | 1.72    | 0.22     | 1.39     | [0.02, 2.0]   | yes (close to upper) |
| minP             | **0.0116** | 0.041 | 0.94     | [0.02, 2.0]   | **NO (below)** |
| activity_max_mult (p95(mu)/minP) | **37** | 99 | 1.11 | [2, 8]    | **NO (above)** |

Three axes are outside the LHS+topup combined coverage for takeshi:
`moi` (266 vs max 200), `minP` (0.0116 vs min 0.02), and `activity_max_mult`
(37 vs max 8). Same clamp-and-flag handling as the existing cohen/shendure
out-of-range cases applies. Takeshi K562 and SKNSH have somewhat narrower
activity ranges (p95/minP ~9) but similar issues elsewhere.

The `reference` HepG2 anchor is the cleanest pick for parity with
shendure-Pluripotent and cohen-Rod.

## Design

### Approach: parameterize anchor set, generate both figure variants

1. **EMPIRICAL stays a 3-anchor dict.** Add takeshi-HepG2 to
   `synthetic_factorial.py::EMPIRICAL` with the values above. Sample-id
   prefix conventions, simulation logic etc are unaffected.

2. **Refactor plotting/attribution helpers to accept an optional
   `anchors` list.** Default is `["shendure-Pluripotent", "cohen-Rod"]`
   for methods-paper compatibility. Pass
   `["shendure-Pluripotent", "cohen-Rod", "takeshi-HepG2"]` for the
   3-anchor variant.

   Concretely:
   - `_plot_marginals_for_metric(df, metric, ..., anchors=DEFAULT_ANCHORS)`
   - `attribution.py` becomes parameterized: top-level `main(anchors=...)`
     produces filenames suffixed with `_no_takeshi` (default) or
     `_with_takeshi`.

3. **Output naming convention**: every figure that previously was called
   `marginals_*_combined.svg` becomes `marginals_*_combined_no_takeshi.svg`
   (the methods-paper version) and a sibling
   `marginals_*_combined_with_takeshi.svg` is also produced. Same for
   `attribution_bar.svg` -> `attribution_bar_no_takeshi.svg` /
   `attribution_bar_with_takeshi.svg` and the marginals-with-arrows variant.

   Keep the legacy `marginals_combined.svg` (no-takeshi) symlinked or
   regenerated for backwards compatibility with existing references.

4. **Attribution mechanics with 3 anchors.** With three anchors, "delta"
   isn't well-defined as a single number. Two options:

   - **(a) Pairwise** (recommended): for each pair (cohen-shendure,
     takeshi-shendure, takeshi-cohen), produce one bar chart of per-axis
     deltas. Three small bar charts per row, plus per-axis annotated
     marginals showing all three anchors with arrows between.
   - **(b) Per-anchor LOESS prediction**: just plot per-axis predicted
     power for each anchor, no deltas. Simpler but less directly answers
     the original "which axis explains the gap" question.

   Default to (a). The existing 2-anchor attribution stays as a special
   case of (a) when only cohen/shendure are passed.

5. **Color convention.** Extend `empirical_colors`:
   - shendure -> darkorange (existing)
   - cohen -> purple (existing)
   - takeshi -> teal (new)

## Files

- **Edit**:
  - `synthetic_factorial.py`:
    - Add `takeshi-HepG2` entry to `EMPIRICAL`.
    - Add `DEFAULT_ANCHORS = ["shendure-Pluripotent", "cohen-Rod"]` constant.
    - Add `empirical_colors` dict at module scope (currently inline in
      `_plot_marginals_for_metric`); add takeshi -> teal.
    - Refactor `_plot_marginals_for_metric` to accept `anchors` kwarg.
    - Refactor `phase_plot` to call `_plot_marginals_for_metric` twice per
      metric: once with `DEFAULT_ANCHORS` (filename suffix
      `_no_takeshi`) and once with all three (filename suffix
      `_with_takeshi`).
  - `attribution.py`:
    - Top-level `main(anchors_label, anchors_list)`. Default invocation
      runs both `("no_takeshi", DEFAULT_ANCHORS)` and
      `("with_takeshi", DEFAULT_ANCHORS + ["takeshi-HepG2"])`.
    - For the 3-anchor case, generate three pairwise bar charts in a
      single figure (subplots) and a single annotated-marginals figure
      with three points per panel + three arrows (or pick one
      "primary pair" -- TBD when we look at the result).
- **Edit** `PLOTTING_IDEAS.md`: append a "anchor versioning" section
  documenting the no_takeshi / with_takeshi convention.
- **Read / reuse**: `samples_power_combined.parquet` (no resim required;
  takeshi is a new overlay on the existing surface).

## Verification

1. `python attribution.py` produces 4 SVGs (2 versions x 2 plot types):
   - `attribution_bar_no_takeshi.svg`
   - `attribution_bar_with_takeshi.svg`
   - `attribution_marginals_no_takeshi.svg`
   - `attribution_marginals_with_takeshi.svg`
2. `python synthetic_factorial.py plot full` produces 10 SVGs (5 metrics x
   2 versions); the no-takeshi versions match what's currently in `output/`
   modulo the filename suffix.
3. Spot-check the with_takeshi marginals: for the 3 axes where takeshi is
   out of LHS coverage (moi, minP, activity_max_mult), the takeshi marker
   should appear at the LHS edge with a `*` flag.
4. Pairwise attribution: cohen-vs-takeshi delta on `bcs_per_cre` should be
   strongly positive (cohen 17244 vs takeshi 394, ~44x). Cohen-vs-takeshi
   on `moi` clamps both to LHS edge (cohen and takeshi both at 200) so
   delta is near zero -- expected and informative caveat.

## Phase 2: extension sweep (`topup3`) to cover takeshi + close
remaining shendure/cohen gaps on activity_max_mult

User decision (2026-05-06): extend sims to cover takeshi's three out-of-
range axes. Same opportunity to close the activity_max_mult gap on
shendure (99) and cohen (1.11) which has been outstanding since
2026-05-05.

### TOPUP3_AXIS_BOUNDS

Extend three axes; leave the other four at their original LHS bracket.
Combined coverage (full + topup + topup3) will then bracket all three
empirical anchors on every axis.

| Axis              | original LHS | topup    | topup3 (proposed)    | covers |
|-------------------|--------------|----------|----------------------|--------|
| n_cells           | [500, 5e4]   | [500, 5e4] | [500, 5e4]         | -- |
| n_cres            | [50, 2000]   | [50, 2000] | [50, 2000]         | -- |
| bcs_per_cre       | [3, 500]     | [500, 5e4] | [3, 5e4]           | -- |
| **moi**           | [0.5, 30]    | [30, 200]  | **[200, 350]**     | takeshi=266 |
| lib_alpha_nb      | [0.02, 2.0]  | [0.02, 2.0] | [0.02, 2.0]       | -- |
| **minP**          | [0.02, 2.0]  | [0.02, 2.0] | **[0.003, 0.02]** | takeshi=0.0116 |
| **activity_max_mult** | [2, 8]   | [2, 8]   | **[1, 120]**         | takeshi=37, shendure=99, cohen=1.11 |

Note `activity_max_mult` is widened both directions (1 covers cohen-low,
120 covers shendure-high); LHS will sample throughout [1, 120] log-scale.
The non-extended axes use the *full* combined range (e.g. bcs_per_cre
[3, 50000]) so each topup3 sample spans the full design space on those
axes. Other axes that already cover all anchors stay at original LHS.

### topup3 mode

Add to MODES:
```python
"topup3": dict(n_lhs=120, n_library_reps=5, n_sims=5, n_workers=50,
               worker_mem="48G", n_slices=10),
```

120 LHS samples = a touch more than topup (100) since we're extending three
axes simultaneously and want decent density across the new corners. Same
cost shape: ~10-15 hr wall on 10 parallel slices given that takeshi-empirical
work proxy (~2.6e10) is between cohen-empirical (~4.8e10) and the heavy
topup samples (~1.4e13). Driver mem stays at 64 GB (the topup OOM lesson).

### Sample-id prefix

`u` prefix for topup3 samples (`s` = full, `t` = topup, `u` = topup3 -- next
letter). Update `get_samples()` to dispatch on `CURRENT_MODE`.

### Combined parquet

`samples_power_combined.parquet` becomes the union of
`samples_power.parquet` + `samples_power_topup.parquet` +
`samples_power_topup3.parquet`. Update `phase_plot` mode=="topup3" branch
to load all three and concat. Total samples: 1000 + 99 + 120 = 1219.

### Worker / disk budget

Estimated disk: topup at ~9 GB/sample, takeshi-corner samples likely
similar (high moi but smaller ncells than topup's heavy tail). 120 * 9
= ~1.1 TB. Combined with existing 1.4 TB on scratch = ~2.5 TB total,
within group quota (4.3 TB free).

### Why now and not after Phase 1 figures

Two reasons:
1. Sim wall time (~10-15 hr) overlaps cleanly with the rest of the
   work; we lose nothing by launching it in parallel with the Phase 1
   refactor.
2. Phase 1 figures will look more authoritative if all three anchors
   sit comfortably inside the LHS for at least one figure variant.
   The methods-paper no_takeshi variant won't change but the
   with_takeshi version benefits.

## Rough effort

- Phase 1 (figure versioning + takeshi entry): ~1 hr code + 10 min replot
- Phase 2 (topup3 sim): ~30 min code + 10-15 hr cluster wall (no
  human attention required after launch)
- Final combined-plot + attribution rerun: ~30 min

Phase 1 and Phase 2 can run in parallel.
