import numpy as np
import pandas as pd
import pytest

pytest.importorskip("dask")
dd = pytest.importorskip("dask.dataframe")
pytest.importorskip("tensorzinb.tensorzinb")

from scMPRAforge.core import _optimize_mpra_ddf, scMPRA_data


def _make_obj(pdf: pd.DataFrame, table_type: str = "mpra_umiwise") -> scMPRA_data:
    obj = scMPRA_data()
    obj.data = _optimize_mpra_ddf(dd.from_pandas(pdf, npartitions=2))
    obj.table_type = table_type
    return obj


def test_consider_missing_expands_replicate_cartesian_and_zero_fills():
    # r1: 2 cells x 2 mpra = 4 rows expected
    # r2: 1 cell x 2 mpra = 2 rows expected
    # total expected = 6
    pdf = pd.DataFrame(
        {
            "rep_id": ["r1", "r1", "r1", "r2", "r2"],
            "cell_bc": ["c1", "c1", "c2", "c3", "c3"],
            "cell_type": ["ct1", "ct1", "ct2", "ct3", "ct3"],
            "mpra_bc": ["m1", "m1", "m2", "m3", "m4"],
            "cre_id": ["creA", "creA", "creB", "creC", "creD"],
            # duplicate (r1,c1,m1) to test summation before filling
            "umis_mpra_bc": [2, 3, 7, 1, 4],
            "reads_mpra_bc": [20, 30, 70, 10, 40],
        }
    )

    obj = _make_obj(pdf)
    obj.consider_missing(max_memory_gb=1.0)
    out = obj.data.compute().sort_values(["rep_id", "cell_bc", "mpra_bc"]).reset_index(drop=True)

    assert list(out.columns) == ["cell_bc", "rep_id", "cre_id", "cell_type", "mpra_bc", "umis_mpra_bc"]
    assert len(out) == 6
    assert isinstance(out["umis_mpra_bc"].dtype, pd.SparseDtype)
    assert out["umis_mpra_bc"].dtype.subtype == np.dtype("int64")

    # Existing summed combo retained
    v = out.loc[(out["rep_id"] == "r1") & (out["cell_bc"] == "c1") & (out["mpra_bc"] == "m1"), "umis_mpra_bc"]
    assert int(v.iloc[0]) == 5

    # Missing combos added with zeros
    z1 = out.loc[(out["rep_id"] == "r1") & (out["cell_bc"] == "c1") & (out["mpra_bc"] == "m2"), "umis_mpra_bc"]
    z2 = out.loc[(out["rep_id"] == "r1") & (out["cell_bc"] == "c2") & (out["mpra_bc"] == "m1"), "umis_mpra_bc"]
    assert int(z1.iloc[0]) == 0
    assert int(z2.iloc[0]) == 0

    # Operation tracking
    op_name, payload = obj.operations[-1]
    assert op_name == "consider_missing"
    assert payload["observed_rows"] == 4
    assert payload["expanded_rows"] == 6
    assert payload["rows_added"] == 2
    assert payload["estimated_peak_memory_gb"] > 0


def test_consider_missing_rejects_non_umiwise_tables():
    pdf = pd.DataFrame(
        {
            "cell_bc": ["c1"],
            "rep_id": ["r1"],
            "cre_id": ["creA"],
            "cell_type": ["ct1"],
            "mpra_bc": ["m1"],
            "umi": ["u1"],
            "reads": [1],
        }
    )
    obj = _make_obj(pdf, table_type="mpra_readwise")
    with pytest.raises(ValueError, match="UMI-wise"):
        obj.consider_missing()


def test_consider_missing_rejects_missing_required_columns():
    pdf = pd.DataFrame(
        {
            "rep_id": ["r1"],
            "cre_id": ["creA"],
            "cell_type": ["ct1"],
            "umis_mpra_bc": [1],
        }
    )
    obj = _make_obj(pdf)
    with pytest.raises(ValueError, match="missing"):
        obj.consider_missing()


def test_consider_missing_rejects_ambiguous_cell_mapping():
    pdf = pd.DataFrame(
        {
            "rep_id": ["r1", "r1"],
            "cell_bc": ["c1", "c1"],
            "cell_type": ["ct1", "ct2"],  # ambiguous within replicate
            "mpra_bc": ["m1", "m2"],
            "cre_id": ["creA", "creB"],
            "umis_mpra_bc": [1, 2],
        }
    )
    obj = _make_obj(pdf)
    with pytest.raises(ValueError, match="cell_bc"):
        obj.consider_missing()


def test_consider_missing_rejects_ambiguous_mpra_mapping():
    pdf = pd.DataFrame(
        {
            "rep_id": ["r1", "r1"],
            "cell_bc": ["c1", "c2"],
            "cell_type": ["ct1", "ct2"],
            "mpra_bc": ["m1", "m1"],  # ambiguous within replicate
            "cre_id": ["creA", "creB"],
            "umis_mpra_bc": [1, 2],
        }
    )
    obj = _make_obj(pdf)
    with pytest.raises(ValueError, match="mpra_bc"):
        obj.consider_missing()


def test_consider_missing_respects_memory_guard():
    pdf = pd.DataFrame(
        {
            "rep_id": ["r1", "r1", "r1", "r1"],
            "cell_bc": ["c1", "c2", "c1", "c2"],
            "cell_type": ["ct1", "ct2", "ct1", "ct2"],
            "mpra_bc": ["m1", "m1", "m2", "m2"],
            "cre_id": ["creA", "creA", "creB", "creB"],
            "umis_mpra_bc": [1, 1, 1, 1],
        }
    )
    obj = _make_obj(pdf)
    with pytest.raises(ValueError, match="max_memory_gb"):
        obj.consider_missing(max_memory_gb=1e-9)
