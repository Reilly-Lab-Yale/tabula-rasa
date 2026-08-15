#!/usr/bin/env python3
"""
Write ortho_meta.json onto every fitted ortho.

The record separates what was modelled -- which unobserved combinations became
zeros, at what reporter resolution -- from how it was computed, i.e. whether
those zeros were rows in the input or design-row weights built at fit time.
`fit_mode` conflates the two: a fit whose zeros arrived in the input reads as
`standard`, and its modelling choice then survives only in the directory name.

Every field is derived from an artifact rather than typed in: the fit script
for the call that produced the ortho, core.py's `*_BOUNDS` aliases for which
fit each dataset treats as canonical, the design dicts and the source count
table as independent cross-checks, and the run_stats filename for the date.
A record that cannot be derived is refused, with the missing evidence named.

Schema and rationale: ORTHO_METADATA.md.

    python backfill_ortho_meta.py            # derive and report, write nothing
    python backfill_ortho_meta.py --write
"""
import argparse
import json
import pathlib
import pickle
import re
import sys

DATA_ROOT = pathlib.Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new")
REPO = pathlib.Path(__file__).resolve().parents[2]
FITS_ROOT = REPO / "analyses" / "model_fitting" / "fits"
CORE_PY = REPO / "scMPRAforge" / "core.py"

COUNT_COL = "umis_mpra_bc"

# Zeros a fit never materialises are built from a rule at fit time; zeros that
# arrived in the input were built by data prep and are simply read.
FIT_TIME = {"per_delivery", "per_barcode", "all_combinations", "all_combinations_moi"}


# ---------------------------------------------------------------- evidence

def parse_fit_script(path):
    """Recover the fitting call from a fit.py without importing it."""
    src = path.read_text()

    # The leading dot keeps prose in docstrings ("consider_missing adds full
    # Cartesian zeros...") from reading as a call.
    ev = {
        "consider_missing": bool(
            re.search(r"\.set_consider_missing\s*\(\s*(?:True|enabled\s*=\s*True)", src)),
        "moi_correct_cm": bool(re.search(r"\bmoi_correct_cm\s*=\s*True", src)),
        "phantom_compress": bool(re.search(r"\bphantom_compress\s*=\s*True", src)),
        "coarse_reporter": bool(re.search(r"\.set_coarse_reporter\s*\(", src)),
        "nb_only": bool(re.search(r"\bnb_only\s*=\s*True", src)),
    }

    seen = set(re.findall(r'\breporter_expansion\s*=\s*"([a-z]+)"', src))
    assert len(seen) <= 1, f"{path}: by_cre and by_cell_type disagree on reporter_expansion: {seen}"
    ev["reporter_expansion"] = seen.pop() if seen else None

    m = re.search(r'DATA_DIR\s*=\s*Path\("([^"]+)"\)', src)
    assert m, f"{path}: no DATA_DIR"
    ev["data_dir"] = pathlib.Path(m.group(1))
    ev["dataset"] = ev["data_dir"].name

    m = re.search(r'\.from_(tsv|parquet)\(\s*str\(DATA_DIR\s*/\s*"([^"]+)"\)', src)
    assert m, f"{path}: no from_tsv/from_parquet call"
    ev["source_table"] = ev["data_dir"] / m.group(2)

    return ev


def canonical_orthos():
    """Read the `<DS>_BOUNDS = <ALIAS>` aliases core.py uses to pick a default."""
    src = CORE_PY.read_text()
    alias_to_ortho = dict(re.findall(
        r"^(\w+)\s*=\s*Bounds\.from_tgz\([^)]*?presets/([\w.]+)\.tgz",
        src, re.M))
    canon = {}
    for ds_alias, alias in re.findall(r"^(\w+)_BOUNDS\s*=\s*(\w+_BOUNDS)\s*$", src, re.M):
        ortho = alias_to_ortho.get(alias)
        assert ortho, f"core.py: {ds_alias}_BOUNDS points at {alias}, which loads no preset"
        canon[ortho] = ds_alias.lower()
    return canon


_ZERO_CACHE = {}


def source_has_zeros(path):
    """True/False if the source count table can be read, None if it cannot."""
    key = str(path)
    if key in _ZERO_CACHE:
        return _ZERO_CACHE[key]

    col = None
    if path.name.endswith(".tsv") or path.name.endswith(".tsv.gz"):
        if path.exists():
            import pandas as pd
            col = pd.read_csv(path, sep="\t", usecols=[COUNT_COL])[COUNT_COL]
    elif path.is_dir():
        # A .scmpra is a directory: Parquet plus a members.json manifest.
        import pandas as pd
        col = pd.read_parquet(path / "data.parquet", columns=[COUNT_COL])[COUNT_COL]

    result = None
    if col is not None:
        assert (col >= 0).all(), f"{path}: negative counts in {COUNT_COL}"
        result = bool((col == 0).any())
    _ZERO_CACHE[key] = result
    return result


def design_tags(ortho_dir):
    """fit_mode / reporter_expansion as recorded in the design dicts, if any."""
    for sub in ("by_cell_type_design", "by_cre_design"):
        pkls = sorted((ortho_dir / sub).glob("*.pkl")) if (ortho_dir / sub).is_dir() else []
        if pkls:
            d = pickle.load(open(pkls[0], "rb"))
            return d.get("fit_mode"), d.get("reporter_expansion")
    return None, None


def fitted_at(fit_dir, ortho_dir):
    """Fit date, from the run_stats report if one was kept, else from the save."""
    dates = sorted(re.findall(r"run_stats_(\d{4})(\d{2})(\d{2})_\d{6}\.txt",
                              " ".join(p.name for p in fit_dir.iterdir())))
    if dates:
        return "-".join(dates[-1]), "run_stats filename"
    import datetime
    day = datetime.date.fromtimestamp(ortho_dir.stat().st_mtime).isoformat()
    return day, "ortho directory mtime (no run_stats report kept)"


# ---------------------------------------------------------------- derivation

def derive(ortho_dir, fit_dir, canon):
    ev = parse_fit_script(fit_dir / "fit.py")

    if ev["consider_missing"]:
        # Consider-missing definitionally disregards reporter information, so
        # it never carries a reporter expansion.
        zero_expansion = "all_combinations_moi" if ev["moi_correct_cm"] else "all_combinations"
        reporter_resolution, reporter_source = "none", "absent"
    elif ev["phantom_compress"]:
        assert ev["coarse_reporter"], (
            f"{fit_dir.name}: phantom_compress=True with no set_coarse_reporter; "
            "core.py raises on this, so the ortho cannot have been fit this way")
        zero_expansion = "per_delivery" if ev["reporter_expansion"] == "single" else "per_barcode"
        reporter_resolution, reporter_source = "element", "separate_table"
    else:
        # Neither path taken: no zeros are built, so every zero the fit saw was
        # already a row in the input, at whatever resolution the reporter has.
        zero_expansion = "preexisting"
        reporter_resolution, reporter_source = "barcode", "in_table"

    fit_time = zero_expansion in FIT_TIME
    when, when_from = fitted_at(fit_dir, ortho_dir)
    rec = {
        "zero_expansion": zero_expansion,
        "expansion_stage": "fit_time" if fit_time else "data_prep",
        "model_family": "nb" if ev["nb_only"] else "zinb",
        "reporter_resolution": reporter_resolution,
        "reporter_source": reporter_source,
        "zero_storage": "phantom_compressed" if fit_time else "materialized",
        "dataset": ev["dataset"],
        "source_table": ev["source_table"].name,
        # Recorded rather than inferred, so the `preexisting` classification
        # can be audited without re-reading the count table.
        "source_has_zero_rows": source_has_zeros(ev["source_table"]),
        "canonical": canon.get(ortho_dir.name) == ev["dataset"],
        "fit_script": str(fit_dir.relative_to(REPO) / "fit.py"),
        "fitted_at": when,
        "fitted_at_basis": when_from,
        # No fit recorded the repo state it ran against; a sha inferred from
        # the date would be a guess dressed as provenance.
        "code_version": None,
        "backfilled": True,
        "backfill_basis": "derived from fit.py, core.py *_BOUNDS, design dicts, source table",
    }
    if zero_expansion == "per_barcode" and reporter_resolution == "element":
        # Blanketing an element's whole barcode set from one element-level
        # observation: the comparison Results 2.1 argues against, kept visible.
        rec["counterfactual"] = True

    cross_check(rec, ev, ortho_dir, fit_dir)
    validate(rec, ortho_dir)
    return rec


def cross_check(rec, ev, ortho_dir, fit_dir):
    """Agreement between independent artifacts: name, design dicts, source data."""
    name = ortho_dir.name
    assert (rec["model_family"] == "nb") == bool(re.search(r"_nb_", name)), (
        f"{name}: directory name and nb_only in {fit_dir.name}/fit.py disagree on count family")

    fit_mode, rep_exp = design_tags(ortho_dir)
    if fit_mode is not None:
        expect = {"all_combinations": "cm_phantom",
                  "all_combinations_moi": "cm_phantom_moib",
                  "per_delivery": "obs_phantom",
                  "per_barcode": "obs_phantom",
                  "preexisting": "standard"}[rec["zero_expansion"]]
        assert fit_mode == expect, f"{name}: design dict says fit_mode={fit_mode!r}, derived {expect!r}"
    if rep_exp is not None:
        expect = {"per_delivery": "single", "per_barcode": "coarse"}.get(rec["zero_expansion"])
        assert rep_exp == expect, f"{name}: design dict says reporter_expansion={rep_exp!r}, derived {expect!r}"

    # One-directional: zeros the fit never builds must already be in the input.
    # The converse does not hold -- a table carrying zeros can still be refit
    # under consider-missing, which is how the reporter-free counterfactuals
    # for Lalanne et al. are constructed.
    if rec["zero_expansion"] == "preexisting":
        assert rec["source_has_zero_rows"] is not False, (
            f"{name}: fit builds no zeros and {ev['source_table'].name} has none either, "
            "so the fit saw no zeros at all")


def validate(rec, ortho_dir):
    name = ortho_dir.name
    if rec["reporter_resolution"] == "none":
        assert rec["zero_expansion"] in ("all_combinations", "all_combinations_moi"), \
            f"{name}: no reporter, but zero_expansion={rec['zero_expansion']!r}"
    if rec["zero_expansion"] == "preexisting":
        assert rec["expansion_stage"] == "data_prep" and rec["zero_storage"] == "materialized", \
            f"{name}: preexisting zeros cannot be built or compressed at fit time"
    if rec["zero_storage"] == "phantom_compressed":
        assert rec["expansion_stage"] == "fit_time", \
            f"{name}: nothing compresses zeros that were never built at fit time"


def load(ortho_dir):
    """Read an ortho's metadata, revalidating it. The counterpart to --write."""
    ortho_dir = pathlib.Path(ortho_dir)
    path = ortho_dir / "ortho_meta.json"
    assert path.exists(), f"{ortho_dir.name}: no ortho_meta.json; run backfill_ortho_meta.py"
    rec = json.loads(path.read_text())
    validate(rec, ortho_dir)
    return rec


# ---------------------------------------------------------------- driver

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true", help="write ortho_meta.json (default: report only)")
    ap.add_argument("--data-root", type=pathlib.Path, default=DATA_ROOT)
    args = ap.parse_args()

    canon = canonical_orthos()

    records, skipped = {}, []
    for ds in sorted(args.data_root.iterdir()):
        if not ds.is_dir():
            continue
        for o in sorted(ds.iterdir()):
            # cohen carries a date suffix on the orthos refit on 2026-04-01;
            # the fit script keeps the undated name.
            fit_dir = FITS_ROOT / re.sub(r"_\d{8}$", "", o.name)
            if not (fit_dir / "fit.py").exists():
                continue          # a count table or some other non-ortho entry
            if not ((o / "by_cell_type_design").is_dir() or (o / "by_cre_design").is_dir()):
                # A fit script names it, but the ortho is not readable here.
                # This tree is a partial copy; the fit lives on the cluster.
                why = ("dangling symlink to " + str(o.readlink())) if o.is_symlink() \
                    else "no design directory"
                skipped.append((o.name, why))
                continue
            records[o] = derive(o, fit_dir, canon)
    assert records, f"no fitted orthos found under {args.data_root}"

    by_dataset = {}
    for rec in records.values():
        by_dataset.setdefault(rec["dataset"], []).append(rec["canonical"])
    for ds, flags in sorted(by_dataset.items()):
        n = sum(flags)
        assert n <= 1, f"{ds}: {n} orthos flagged canonical, at most one is allowed"
        if n == 0:
            print(f"note: {ds} has no canonical fit (no {ds.upper()}_BOUNDS alias in core.py)")

    for o, rec in records.items():
        flag = "*" if rec["canonical"] else ("!" if rec.get("counterfactual") else " ")
        print(f"{flag} {o.name:38s} {rec['zero_expansion']:22s} "
              f"{rec['expansion_stage']:10s} {rec['zero_storage']:20s} "
              f"{rec['reporter_resolution']:8s} {rec['model_family']:5s} {rec['fitted_at']}")
        if args.write:
            (o / "ortho_meta.json").write_text(json.dumps(rec, indent=2) + "\n")

    for name, why in skipped:
        print(f"skipped {name}: {why}", file=sys.stderr)

    print(f"\n{len(records)} orthos {'written' if args.write else 'derived'}"
          f"{f', {len(skipped)} skipped' if skipped else ''}"
          f" (* canonical, ! counterfactual)")


if __name__ == "__main__":
    main()
