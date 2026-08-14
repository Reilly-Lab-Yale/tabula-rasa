#!/usr/bin/env python3
"""Compare our per-cell-type activity calls against Lalanne et al.'s published set.

Lalanne et al. call a CRE active in mouse embryoid bodies if its reporter
expression is in significant excess of the basal controls (no promoter and
minimal promoter) by bootstrap resampling, in all three biological
replicates, aggregating over every cell. That yields 58 of 204 endogenous
CREs, of which 10 are further called cell-type specific by a permutation
test. Their per-CRE results are Supplementary Data 5.

We test each (CRE, cell type) pair against the minimal-promoter reference
within that cell type by MWU, with Benjamini-Hochberg applied once across the
whole hypothesis set, so the ten-cell-type multiplicity is already paid for.
A CRE counts as active here if any cell type is significant and up.

The asymmetry is the point rather than a nuisance: a whole-population test
dilutes an element that fires in one lineage, which is the case the stratified
test is meant to catch. So the script reports not just how many extra
elements we call but whether they look like that -- how many cell types each
is significant in, and what aggregate activity the published table assigns
them.

The supplementary file is downloaded on demand and checksummed rather than
vendored, since it is Springer's to distribute.

    python lalanne_active_cres.py
"""
import hashlib
import pathlib
import sys
import urllib.request

import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
CACHE = HERE / "cache"
OUT = HERE / "output"

# Supplementary Data 5 of doi:10.1038/s41592-024-02260-3.
SUPP_URL = (
    "https://static-content.springer.com/esm/"
    "art%3A10.1038%2Fs41592-024-02260-3/MediaObjects/"
    "41592_2024_2260_MOESM7_ESM.xlsx"
)
SUPP_SHA256 = "8cd9c7eb19c9ccad417b2035f1da4f1f03f0fc06f33f1ac11fe368b09863b9ec"
# Row index of the real header in the 'aggregated' sheet; the rows above it
# are a title and a per-column legend.
SUPP_HEADER_ROW = 11

OURS = (HERE.parents[1] / "empirical_testing" / "activity_volcano" / "output"
        / "shendure" / "shendure_activity_mwu.tsv")

ALPHA = 0.05


def fetch_supplement():
    CACHE.mkdir(exist_ok=True)
    dest = CACHE / "lalanne2024_supp_data_5.xlsx"
    if not dest.exists():
        print(f"downloading {SUPP_URL}", flush=True)
        req = urllib.request.Request(
            SUPP_URL, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req) as r, open(dest, "wb") as f:
            f.write(r.read())
    got = hashlib.sha256(dest.read_bytes()).hexdigest()
    assert got == SUPP_SHA256, (
        f"supplementary file checksum changed: {got} != {SUPP_SHA256}. "
        "The published table may have been revised; re-verify the counts "
        "below before trusting this comparison.")
    return dest


def load_published(path):
    d = pd.read_excel(path, sheet_name="aggregated", header=SUPP_HEADER_ROW)
    d = d.dropna(subset=["CRE_id"])
    for col in ("all_rep_act_hit", "all_rep_spec_hit"):
        d[col] = d[col].astype(bool)
    dev = d[d.CRE_class == "devCRE"].set_index("CRE_id")
    # The paper's own headline numbers; if these break, the parse is wrong.
    assert len(dev) == 204, f"expected 204 devCREs, got {len(dev)}"
    assert dev.all_rep_act_hit.sum() == 58, (
        f"expected 58 active, got {dev.all_rep_act_hit.sum()}")
    assert dev.all_rep_spec_hit.sum() == 10, (
        f"expected 10 specific, got {dev.all_rep_spec_hit.sum()}")
    return dev


def load_ours(path):
    d = pd.read_csv(path, sep="\t")
    assert (d.test_type == "mwu").all(), "expected MWU results only"
    assert d.bh_p.notna().all(), "missing BH-adjusted p-values"
    return d


def main():
    pub = load_published(fetch_supplement())
    ours = load_ours(OURS)

    shared = sorted(set(pub.index) & set(ours.comparison_CRE))
    only_ours = sorted(set(ours.comparison_CRE) - set(pub.index))
    only_pub = sorted(set(pub.index) - set(ours.comparison_CRE))
    print(f"devCREs matched by identifier: {len(shared)} of {len(pub)} "
          f"published")
    if only_pub:
        print(f"  published, untested here: {only_pub}")
    if only_ours:
        # Controls and a couple of elements whose identifiers differ by one in
        # the published table; reported rather than dropped silently.
        print(f"  tested here, unmatched: {only_ours}")
    assert len(shared) >= 200, f"identifier join collapsed: {len(shared)}"

    sig = ours[(ours.bh_p < ALPHA) & (ours.log2_fc > 0)
               & ours.comparison_CRE.isin(shared)]
    n_cts = sig.groupby("comparison_CRE").comparison_cell_type.nunique()

    theirs = set(pub.index[pub.all_rep_act_hit]) & set(shared)
    mine = set(n_cts.index)
    recovered = sorted(theirs & mine)
    missed = sorted(theirs - mine)
    extra = sorted(mine - theirs)

    print()
    print(f"published active   : {len(theirs)}")
    print(f"called active here : {len(mine)}")
    print(f"  recovered        : {len(recovered)}")
    print(f"  missed           : {len(missed)}"
          + (f"  {missed}" if missed else ""))
    print(f"  additional       : {len(extra)}")
    print(f"superset of published set: {not missed}")

    print()
    print("cell types significant per CRE:")
    for label, ids in (("recovered", recovered), ("additional", extra)):
        s = n_cts[ids]
        print(f"  {label:11s} n={len(ids):3d}  median {s.median():.0f}  "
              f"mean {s.mean():.2f}")
    print()
    print("published aggregate activity (Activity_all_cells), median:")
    never = [c for c in shared if c not in mine]
    for label, ids in (("recovered", recovered), ("additional", extra),
                       ("never called", never)):
        print(f"  {label:13s} {pub.loc[ids, 'Activity_all_cells'].median():.4f}")

    # The claim the manuscript makes rests on both halves: we lose nothing,
    # and what we add is the lineage-restricted tail rather than noise.
    assert not missed, "no longer a superset of the published calls"
    assert n_cts[extra].median() < n_cts[recovered].median(), (
        "additional calls are not more cell-type-restricted than recovered "
        "ones; the dilution explanation does not hold")

    OUT.mkdir(exist_ok=True)
    rows = pd.DataFrame({
        "CRE_id": shared,
        "published_active": [c in theirs for c in shared],
        "called_here": [c in mine for c in shared],
        "n_cell_types_significant": [int(n_cts.get(c, 0)) for c in shared],
        "published_activity_all_cells": pub.loc[shared,
                                                "Activity_all_cells"].values,
        "published_specific": pub.loc[shared, "all_rep_spec_hit"].values,
    })
    out = OUT / "lalanne_active_cre_comparison.tsv"
    rows.to_csv(out, sep="\t", index=False)
    print(f"\nwrote {out} ({len(rows)} rows)")


if __name__ == "__main__":
    sys.exit(main())
