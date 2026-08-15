"""
The ortho_meta.json sidecar: what a fit modelled, and how it computed it.

`fit_mode` conflates two different things - which unobserved combinations
became zeros (a modelling choice) and whether those zeros were rows in the
input or design-row weights built at fit time (an implementation detail).
A fit whose zeros arrived in the input reads as `standard`, so its modelling
choice survives only in the directory name. This module separates them.

Fits write their own record through `ortho.save()`. Orthos fit before the
sidecar existed are reconstructed from artifacts by
`analyses/model_fitting/backfill_ortho_meta.py`, which is the only thing that
sets `backfilled`. Both go through `classify()` and `validate()` here so the
two paths cannot drift.

Stdlib only, deliberately: the back-fill pass loads this module by path rather
than importing the package, so that reading metadata does not require the
fitting stack.

Schema and rationale: analyses/model_fitting/ORTHO_METADATA.md
"""
import json
import pathlib

FILENAME = "ortho_meta.json"

# Zeros a fit builds from a rule are built at fit time; zeros that arrived in
# the input were built by data prep and are simply read.
FIT_TIME_EXPANSIONS = {"per_delivery", "per_barcode",
                       "all_combinations", "all_combinations_moi"}

# Fields no fit can know about itself. The back-fill pass fills them in; until
# then they are null, which reads as "not yet determined" rather than "no".
DEFERRED = ("canonical", "source_has_zero_rows")


def classify(consider_missing, moi_correct_cm, phantom_compress,
             coarse_reporter, reporter_expansion, nb_only):
    """The six modelling and implementation fields, from the fitting call."""
    if consider_missing:
        # Consider-missing definitionally disregards reporter information, so
        # it never carries a reporter expansion.
        zero_expansion = "all_combinations_moi" if moi_correct_cm else "all_combinations"
        reporter_resolution, reporter_source = "none", "absent"
    elif phantom_compress:
        assert coarse_reporter, (
            "phantom_compress=True with no coarse reporter attached; "
            "standard_fit raises on this")
        zero_expansion = "per_delivery" if reporter_expansion == "single" else "per_barcode"
        reporter_resolution, reporter_source = "element", "separate_table"
    else:
        # Neither path taken: the fit builds no zeros, so every zero it saw was
        # already a row in the input, at whatever resolution the reporter has.
        zero_expansion = "preexisting"
        reporter_resolution, reporter_source = "barcode", "in_table"

    fit_time = zero_expansion in FIT_TIME_EXPANSIONS
    rec = {
        "zero_expansion": zero_expansion,
        "expansion_stage": "fit_time" if fit_time else "data_prep",
        "model_family": "nb" if nb_only else "zinb",
        "reporter_resolution": reporter_resolution,
        "reporter_source": reporter_source,
        "zero_storage": "phantom_compressed" if fit_time else "materialized",
    }
    if zero_expansion == "per_barcode" and reporter_resolution == "element":
        # Blanketing an element's whole barcode set from one element-level
        # observation. Legal, and the comparison Results 2.1 argues against;
        # flagged so it is visible rather than silent.
        rec["counterfactual"] = True
    return rec


def validate(rec, name=""):
    where = f"{name}: " if name else ""
    if rec["reporter_resolution"] == "none":
        assert rec["zero_expansion"] in ("all_combinations", "all_combinations_moi"), \
            f"{where}no reporter, but zero_expansion={rec['zero_expansion']!r}"
    if rec["zero_expansion"] == "preexisting":
        assert rec["expansion_stage"] == "data_prep" and rec["zero_storage"] == "materialized", \
            f"{where}preexisting zeros cannot be built or compressed at fit time"
    if rec["zero_storage"] == "phantom_compressed":
        assert rec["expansion_stage"] == "fit_time", \
            f"{where}nothing compresses zeros that were never built at fit time"
    return rec


def write(ortho_dir, rec):
    ortho_dir = pathlib.Path(ortho_dir)
    validate(rec, ortho_dir.name)
    (ortho_dir / FILENAME).write_text(json.dumps(rec, indent=2) + "\n")


def load(ortho_dir, missing_ok=False):
    """Read an ortho's record, revalidating it."""
    ortho_dir = pathlib.Path(ortho_dir)
    path = ortho_dir / FILENAME
    if not path.exists():
        if missing_ok:
            return None
        raise FileNotFoundError(
            f"{ortho_dir.name} has no {FILENAME}; "
            "run analyses/model_fitting/backfill_ortho_meta.py")
    return validate(json.loads(path.read_text()), ortho_dir.name)
