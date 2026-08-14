#!/usr/bin/env python3
"""Run missing test arms on simulations that already exist, and aggregate.

Testing is a read-only pass over simulated counts, so a sim that was tested
with one test can be tested with another without regenerating anything. This
script does only that: it never simulates, which is the point -- the
per-dataset power scripts decide how many library replicates to make, and
pointing one of those at a restored sim directory would happily start
manufacturing more.

The four arms are {t-test, MWU} x {+reporter, deflated}. The deflated arms
drop every zero-count observation before testing, which is what a no-reporter
experiment sees on disc. Because all four read the same sims, differences
between them are attributable to the test alone.

Both directory layouts are handled: a flat root of sim_* (Zhao et al.) and a
root of per-cell-type directories each holding sim_* (Lalanne et al.).

    python retest_arms.py <sim_root> <out_prefix> [--arms mwu,mwu_deflated]

Writes <out_prefix>_<arm>.parquet next to this script's output/ directory,
with columns (reject_null, fc, cell_type).
"""
import argparse
import pathlib
import sys

import pandas as pd

_repo_root = pathlib.Path(__file__).resolve().parents[3]
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

import scMPRAforge as scm  # noqa: E402

OUTPUT_DIR = pathlib.Path(__file__).resolve().parent / "output"

# (test directory, method name, has_reporter)
ARMS = {
    "ttest":          ("ttest", True),
    "ttest_deflated": ("ttest", False),
    "mwu":            ("mwu",   True),
    "mwu_deflated":   ("mwu",   False),
}
HYPOTHESIS_SET = "hs_all_ct"


def discover(root):
    """Yield (cell_type_or_None, sim_dir) for every simulation under root."""
    subdirs = sorted(p for p in root.iterdir() if p.is_dir())
    assert subdirs, f"no directories under {root}"
    flat = [p for p in subdirs if p.name.startswith("sim_")]
    if flat:
        return [(None, p) for p in flat]
    out = []
    for ct_dir in subdirs:
        out.extend((ct_dir.name, p) for p in sorted(ct_dir.iterdir())
                   if p.is_dir() and p.name.startswith("sim_"))
    assert out, f"no sim_* directories under {root} or its subdirectories"
    return out


def load_tested(root, client):
    """Load sims that already carry at least one arm's results.

    A sim with no tests/ subtree was never completed and is skipped rather
    than tested, since its presence usually means an interrupted run.
    """
    sims = []
    skipped = 0
    for _ct, d in discover(root):
        hs_dir = d / "tests" / HYPOTHESIS_SET
        if not hs_dir.is_dir() or not any(
            (hs_dir / a).is_dir() and any((hs_dir / a).iterdir()) for a in ARMS
        ):
            skipped += 1
            continue
        sims.append(scm.de_novo_simulation(
            location=d.parent, name=d.name, client=client))
    print(f"loaded {len(sims)} tested sims, skipped {skipped} untested",
          flush=True)
    assert sims, "no sims with existing test results; nothing to retest"
    return sims


def run_arm(sims, arm, client):
    method, has_reporter = ARMS[arm]
    done = 0
    for sim in sims:
        d = sim.location / sim.name / "tests" / HYPOTHESIS_SET / arm
        if d.is_dir() and any(d.iterdir()):
            done += 1
            continue
        getattr(sim, method)(HYPOTHESIS_SET, has_reporter=has_reporter)
        sim.save()
    print(f"  {arm}: {len(sims) - done} run, {done} already present",
          flush=True)


def aggregate(sims, arm):
    """Per-(cell type) reject_null and fold change, one row per hypothesis."""
    rows = []
    for sim in sims:
        for i in range(sim.get_state_field("n_sims")):
            m = sim._merge_in_ground_truth(HYPOTHESIS_SET, arm, i)
            m["fc"] = m["comparison_truth"] / m["reference_truth"]
            rows.append(m[["comparison_cell_type", "reject_null", "fc"]])
    combined = pd.concat(rows, ignore_index=True)
    combined = combined.rename(columns={"comparison_cell_type": "cell_type"})
    assert combined["fc"].notna().all(), "fold change has nulls"
    assert (combined["fc"] > 0).all(), "non-positive fold change"
    return combined


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sim_root", type=pathlib.Path)
    ap.add_argument("out_prefix")
    ap.add_argument("--arms", default=",".join(ARMS))
    args = ap.parse_args()

    arms = args.arms.split(",")
    bad = [a for a in arms if a not in ARMS]
    assert not bad, f"unknown arms {bad}; choose from {list(ARMS)}"

    from dask.distributed import Client, LocalCluster
    cluster = LocalCluster(n_workers=1, threads_per_worker=2,
                           memory_limit="48GB")
    client = Client(cluster)

    sims = load_tested(args.sim_root, client)
    for arm in arms:
        run_arm(sims, arm, client)

    OUTPUT_DIR.mkdir(exist_ok=True)
    for arm in arms:
        df = aggregate(sims, arm)
        out = OUTPUT_DIR / f"{args.out_prefix}_{arm}.parquet"
        df.to_parquet(out)
        print(f"wrote {out} ({len(df):,} rows, "
              f"{df['cell_type'].nunique()} cell types)", flush=True)

    client.close()
    print("done")


if __name__ == "__main__":
    main()
