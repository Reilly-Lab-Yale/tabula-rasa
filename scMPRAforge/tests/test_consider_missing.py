import json
import numpy as np
import pandas as pd
import pytest

pytest.importorskip("dask")
dd = pytest.importorskip("dask.dataframe")
pytest.importorskip("tensorzinb.tensorzinb")

from scMPRAforge.core import (
    _bootstrap_build_bundle,
    _mwu_make_bundle,
    _optimize_mpra_ddf,
    HypothesisSet,
    make_by_celltype_hypotheses,
    make_by_cre_hypotheses,
    scMPRA_data,
)


def _make_obj(pdf: pd.DataFrame, table_type: str = "mpra_umiwise") -> scMPRA_data:
    obj = scMPRA_data()
    obj.data = _optimize_mpra_ddf(dd.from_pandas(pdf, npartitions=2))
    obj.table_type = table_type
    return obj


def _toy_missing_pdf() -> pd.DataFrame:
    # One replicate, two cell types, two CREs/barcodes; each observed in one cell type only.
    return pd.DataFrame(
        {
            "rep_id": ["r1", "r1"],
            "cell_bc": ["c1", "c2"],
            "cell_type": ["ctA", "ctB"],
            "mpra_bc": ["m1", "m2"],
            "cre_id": ["cre1", "cre2"],
            "umis_mpra_bc": [5, 7],
        }
    )


def test_consider_missing_policy_defaults_and_setter_validation():
    obj = _make_obj(_toy_missing_pdf())
    assert obj.consider_missing_enabled is False
    assert obj.consider_missing_max_memory_gb == 100.0
    assert obj.consider_missing_peak_overhead_factor == 4.0
    assert obj.consider_missing_subset_semantics == "full_replicate"

    obj.set_consider_missing(True, max_memory_gb=120.0, peak_overhead_factor=3.0)
    assert obj.consider_missing_enabled is True
    assert obj.consider_missing_max_memory_gb == 120.0
    assert obj.consider_missing_peak_overhead_factor == 3.0

    with pytest.raises(ValueError, match="max_memory_gb"):
        obj.set_consider_missing(True, max_memory_gb=0.0)
    with pytest.raises(ValueError, match="peak_overhead_factor"):
        obj.set_consider_missing(True, peak_overhead_factor=0.0)


def test_get_data_requires_split_context_when_enabled():
    obj = _make_obj(_toy_missing_pdf())
    obj.set_consider_missing(True)

    with pytest.raises(ValueError, match="requires a context"):
        obj.get_data(include_missing=True)
    with pytest.raises(ValueError, match="Global consider_missing inflation is disallowed"):
        obj.get_data(include_missing=True, context={"kind": "global"})

    raw = obj.get_data(include_missing=False).compute()
    assert len(raw) == 2


def test_split_level_inflation_cell_type_and_cre_id():
    obj = _make_obj(_toy_missing_pdf())
    obj.set_consider_missing(True, max_memory_gb=10.0)

    ct = obj.get_data(
        include_missing=True,
        context={"kind": "split_level", "split": "cell_type", "level": "ctA"},
    ).compute().sort_values(["cell_bc", "mpra_bc"]).reset_index(drop=True)
    assert list(ct.columns) == ["cell_bc", "rep_id", "cre_id", "cell_type", "mpra_bc", "umis_mpra_bc"]
    assert len(ct) == 2
    assert isinstance(ct["umis_mpra_bc"].dtype, pd.SparseDtype)
    assert int(ct.loc[ct["mpra_bc"] == "m1", "umis_mpra_bc"].iloc[0]) == 5
    assert int(ct.loc[ct["mpra_bc"] == "m2", "umis_mpra_bc"].iloc[0]) == 0

    cr = obj.get_data(
        include_missing=True,
        context={"kind": "split_level", "split": "cre_id", "level": "cre1"},
    ).compute().sort_values(["cell_bc", "mpra_bc"]).reset_index(drop=True)
    assert len(cr) == 2
    assert int(cr.loc[cr["cell_bc"] == "c1", "umis_mpra_bc"].iloc[0]) == 5
    assert int(cr.loc[cr["cell_bc"] == "c2", "umis_mpra_bc"].iloc[0]) == 0


def test_split_level_inflation_hard_memory_cap():
    obj = _make_obj(_toy_missing_pdf())
    obj.set_consider_missing(True, max_memory_gb=1e-9)
    with pytest.raises(ValueError, match="estimated peak memory"):
        obj.get_data(
            include_missing=True,
            context={"kind": "split_level", "split": "cell_type", "level": "ctA"},
        ).compute()


def test_consider_missing_policy_persistence_and_legacy_defaults(tmp_path):
    obj = _make_obj(_toy_missing_pdf())
    obj.set_consider_missing(True, max_memory_gb=123.0, peak_overhead_factor=2.5)
    out = tmp_path / "policy.scmpra"
    obj.to_parquet(str(out))

    loaded = scMPRA_data.from_parquet(str(out))
    assert loaded.consider_missing_enabled is True
    assert loaded.consider_missing_max_memory_gb == 123.0
    assert loaded.consider_missing_peak_overhead_factor == 2.5
    assert loaded.consider_missing_subset_semantics == "full_replicate"

    members_path = out / "members.json"
    members = json.loads(members_path.read_text())
    for k in [
        "consider_missing_enabled",
        "consider_missing_max_memory_gb",
        "consider_missing_peak_overhead_factor",
        "consider_missing_subset_semantics",
    ]:
        members.pop(k, None)
    members_path.write_text(json.dumps(members))

    legacy = scMPRA_data.from_parquet(str(out))
    assert legacy.consider_missing_enabled is False
    assert legacy.consider_missing_max_memory_gb == 100.0
    assert legacy.consider_missing_peak_overhead_factor == 4.0
    assert legacy.consider_missing_subset_semantics == "full_replicate"


def test_hypothesis_builders_and_mwu_respect_policy():
    obj = _make_obj(_toy_missing_pdf())

    hs_raw = make_by_celltype_hypotheses(
        comparison_cell_type="ctA",
        counts=obj,
        comparison_cres="all",
        reference_cre="reference",
    ).to_dataframe()
    assert set(hs_raw["comparison_CRE"].astype(str)) == {"cre1"}

    hs_cre_raw = make_by_cre_hypotheses(
        comparison_cre="cre1",
        counts=obj,
        comparison_cell_types="all",
        reference_cell_type="ctA",
    ).to_dataframe()
    assert len(hs_cre_raw) == 0

    obj.set_consider_missing(True, max_memory_gb=10.0)
    hs_miss = make_by_celltype_hypotheses(
        comparison_cell_type="ctA",
        counts=obj,
        comparison_cres="all",
        reference_cre="reference",
    ).to_dataframe()
    assert set(hs_miss["comparison_CRE"].astype(str)) == {"cre1", "cre2"}

    hs_cre_miss = make_by_cre_hypotheses(
        comparison_cre="cre1",
        counts=obj,
        comparison_cell_types="all",
        reference_cell_type="ctA",
    ).to_dataframe()
    assert set(hs_cre_miss["comparison_cell_type"].astype(str)) == {"ctB"}

    bundle = _mwu_make_bundle(hypotheses=HypothesisSet.from_dataframe(hs_miss), models_or_counts=obj)
    assert ("ctA", "cre2") in bundle["counts"]
    assert np.all(bundle["counts"][("ctA", "cre2")] == 0.0)


def test_bootstrap_opt_out_warns_and_uses_raw_semantics():
    pdf = pd.DataFrame(
        {
            "rep_id": ["r1", "r1", "r1"],
            "cell_bc": ["c1", "c2", "c2"],
            "transfection_bc": ["t1", "t2", "t2"],
            "cell_type": ["ctA", "ctA", "ctA"],
            "cre_id": ["reference", "creX", "creX"],
            "umis_mpra_bc": [0, 3, 1],
        }
    )
    obj = _make_obj(pdf)
    obj.set_consider_missing(True, max_memory_gb=10.0)

    hs = HypothesisSet.from_dataframe(
        pd.DataFrame(
            [
                {
                    "comparison_CRE": "creX",
                    "comparison_cell_type": "MAX",
                    "reference_CRE": "reference",
                    "reference_cell_type": "MAX",
                    "meta": "bootstrap_activity",
                }
            ]
        )
    )

    with pytest.warns(UserWarning, match="bootstrap uses raw data semantics"):
        bundle = _bootstrap_build_bundle(hypotheses=hs, models_or_counts=obj, n_bootstraps=10)

    integrations = bundle["integrations"]
    # Raw semantics: one integration row for reference (c1,t1) and one for creX (c2,t2)
    assert len(integrations) == 2
