# Reporter-informed zero logic for Cohen U6

## Observation taxonomy

```
                    U6 detected?
                   /            \
                 YES              NO
                 |                |
          MPRA obs?          MPRA obs?
         /        \          /        \
       YES         NO      YES         NO
        |           |       |           |
   confirmed    true zero  orphan    not transfected
   nonzero obs  (phantom   nonzero    (no obs, no
   [in fit +    weight)    obs        phantom -- CRE
    zero count             [in fit,   not seen in
    correct]               zero count  this cell]
                           understated
                           -> bug]
```

## The -U6 +MPRA (orphan) case

The bottom-left cell is ambiguous. Two failure modes:

- **(A) Spurious MPRA** -- MPRA barcode is a false positive (index hopping,
  ambient), U6 absence is correct.
- **(B) False-negative U6** -- CRE was genuinely transfected, U6 library was
  too shallow to detect it, MPRA signal is real.

Empirical evidence from Cohen data (2026-04-01):
- 64.5% of U6-confirmed (cell, CRE) pairs have ZERO MPRA obs -- even confirmed
  transfections are frequently missed by MPRA.
- Median MPRA barcode capture rate per confirmed (cell, CRE): 0.01% of library
  barcodes. The MPRA library has ~230K barcodes, severely undersampled per cell.
- U6 is undersequenced

Conclusion: failure mode B dominates. Orphan MPRA obs are predominantly real
transfections where U6 failed to detect.

## Compromise: treat orphan cells as confirmed

`_reporter_zero_counts` (core.py) augments the reporter-confirmed set with
orphan cells before computing phantom zero weights. Orphan cells contribute:
- Their nonzero MPRA obs to the fit (already present in the data).
- Phantom zeros for their unobserved barcodes (same as confirmed cells).

This is a **conservative compromise**: it assumes B dominates and treats all
orphan obs as real. If A were dominant (spurious MPRA), the correct treatment
would be to drop orphan obs entirely. The phantom zero weight for a group cannot
go negative under this treatment (n_total >= n_nonzero by construction).

The scientifically ideal fix would be a probabilistic model of reporter
sensitivity, but that is out of scope for this analysis.
