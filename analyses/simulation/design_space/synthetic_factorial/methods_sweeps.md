# Synthetic-factorial sweep history

The synthetic-factorial design space spans seven axes that jointly determine
the statistical setup of an scMPRA experiment: `n_cells`, `n_cres`,
`bcs_per_cre`, `moi`, `lib_alpha_nb` (NB overdispersion), `minP` (per-CRE
floor activity), and `activity_max_mult` (dynamic range, p95 over minP).

Detection power is characterized from a single 5000-point Latin Hypercube
sample (the **union** sweep) over the full design box. The three earlier
sweeps (full, topup, topup3) were exploratory: they calibrated where the box
boundaries needed to be so the empirical datasets (cohen, shendure, takeshi)
land inside it. Their simulations have since been removed; the union sweep
covers a superset of their ranges on every axis, so nothing is lost in
coverage. This file records that provenance, since the leftover `s`/`t`/`u`
sample-id prefixes and the box-calibration history may otherwise be opaque.

## Sweeps

| Sweep | n_pts | LHS seed | sid prefix | Box | Role |
|---|---|---|---|---|---|
| full   | 1000 | 20260505 | `s` | `AXIS_BOUNDS`        | Exploratory: initial box, calibrated to expected empirical bounds |
| topup  | 100  | 20260506 | `t` | `TOPUP_AXIS_BOUNDS`  | Exploratory: cohen-corner extension (high `bcs_per_cre`, high `moi`) |
| topup3 | 120  | 20260507 | `u` | `TOPUP3_AXIS_BOUNDS` | Exploratory: wider `moi`, `minP`, `activity_max_mult` to reach takeshi |
| **union** | **5000** | 20260511 | `v` | `UNION_AXIS_BOUNDS` | **Canonical: independent LHS over the full union box** |

The union box is the superset of all three exploratory boxes on every axis
(`bcs_per_cre` 3-50k, `moi` 0.5-350, `minP` 0.003-2.0, `activity_max_mult`
1-120; the other three axes were unchanged throughout). So union coverage
strictly contains the earlier sweeps.

## Why a single independent union sweep (not a pooled set)

Each individual sweep is a proper Latin Hypercube (`scipy.stats.qmc`), but
concatenating sweeps is **not** itself an LHS. The exploratory sweeps each
occupy a different sub-box, so a pooled dataset acquires structural
cross-axis correlation - e.g. a high-`moi` sample could only have come from
topup/topup3, which also biased its `bcs_per_cre` and `minP`. That confounds
per-axis attribution (a LOESS marginal over a correlated design implicitly
absorbs the co-moving variables).

| Dataset | n | mean &#124;r&#124; | max &#124;r&#124; |
|---|---|---|---|
| full only            | 1000 | 0.033 | 0.104 |
| combined patchwork (full+topup+topup3) | 1220 | 0.106 | 0.444 |
| **union (canonical)** | **5000** | **0.009** | **0.022** |
| pooled (all four)    | 6220 | 0.027 | 0.094 |

Worst-pair entanglement in the patchwork: `bcs_per_cre x moi` r = +0.44,
`moi x minP` r = -0.41. The union sweep, drawn from scratch over the full
box, has cross-axis &#124;r&#124; at the level expected from sampling noise
alone (max 0.022). A pooled set (union + the three exploratory sweeps) would
still be acceptable (max 0.094, since union is 80% of it), but it trades the
cleanest design for ~24% more points whose only added density is in the
corners - density that a `frac=0.4` LOESS smoother barely uses. We therefore
analyze the union sweep alone.

## Figures

All figures derive from the union sweep's aggregated power summary,
`output/samples_power_union.parquet` (5000 rows, MWU power):

- `marginals*_union*.svg` - per-axis marginals, 5 metrics x {no/with takeshi}
  x {log, linear x}.
- `pairwise_heatmaps_union.svg` - top-3-axis pairwise power heatmaps.
- `attribution_{bar,marginals}_<pair>.svg` - per-axis attribution for the
  three anchor pairs (`cohen_vs_shendure`, `takeshi_vs_shendure`,
  `takeshi_vs_cohen`).

The raw per-rep results behind the summary (5000 samples x 5 reps x 5 sims,
64.8M rows) are consolidated into a single parquet under the shared data area
for any future re-aggregation with different metrics:
`/nfs/roberts/project/pi_skr2/shared/tabula_data/simulated/synthetic_factorial_union_2026-05-16/union_cached_results.parquet`
(338M; see the README there for schema + provenance).

## Reproducing

```bash
# Marginals + pairwise heatmap from the committed power summary (no re-sim):
sbatch wrap_synthetic_factorial.sh replot union

# Per-axis attribution (3 anchor pairs):
sbatch wrap_attribution.sh
```
