import numpy as np
import pandas as pd
import pytest

pytest.importorskip("dask")
dd = pytest.importorskip("dask.dataframe")
distributed = pytest.importorskip("dask.distributed")
pytest.importorskip("tensorzinb.tensorzinb")

from scMPRAforge.core import (
    _optimize_mpra_ddf,
    _prepare_subset_for_modeling,
    _smart_matrix,
    get_cell_counts,
    scMPRA_data,
)


def _is_sparse_int64(series: pd.Series) -> bool:
    dtype = series.dtype
    return (
        pd.api.types.is_sparse(dtype)
        and dtype.subtype == np.dtype("int64")
        and dtype.fill_value == 0
    )


def test_optimize_mpra_ddf_makes_umis_sparse_int64():
    pdf = pd.DataFrame(
        {
            "rep_id": ["r1", "r1", "r2", "r2"],
            "cre_id": ["c1", "c1", "c2", "c2"],
            "cell_type": ["t1", "t1", "t2", "t2"],
            "umis_mpra_bc": [1, None, "3", "bad"],
        }
    )
    ddf = dd.from_pandas(pdf, npartitions=2)

    out = _optimize_mpra_ddf(ddf).compute()

    assert _is_sparse_int64(out["umis_mpra_bc"])
    assert out["umis_mpra_bc"].astype("int64").tolist() == [1, 0, 3, 0]


def test_read_wise_to_umi_wise_keeps_sparse_umis():
    readwise = pd.DataFrame(
        {
            "cell_bc": ["cell1", "cell1", "cell1", "cell2"],
            "rep_id": ["r1", "r1", "r1", "r1"],
            "cre_id": ["creA", "creA", "creA", "creB"],
            "cell_type": ["ct1", "ct1", "ct1", "ct2"],
            "mpra_bc": ["m1", "m1", "m1", "m2"],
            "umi": ["u1", "u1", "u2", "u3"],
            "reads": [2, 1, 4, 5],
        }
    )

    obj = scMPRA_data()
    obj.data = dd.from_pandas(readwise, npartitions=1)
    obj.table_type = "mpra_readwise"

    obj.read_wise_to_umi_wise(keep_reads=False)
    out = obj.data.compute().sort_values(["cell_bc", "cre_id"]).reset_index(drop=True)

    assert _is_sparse_int64(out["umis_mpra_bc"])
    assert out["umis_mpra_bc"].astype("int64").tolist() == [2, 1]


def test_parquet_roundtrip_restores_sparse_umis(tmp_path):
    umiwise = pd.DataFrame(
        {
            "cell_bc": ["cell1", "cell2", "cell3"],
            "rep_id": ["r1", "r1", "r2"],
            "cre_id": ["creA", "creB", "creA"],
            "cell_type": ["ct1", "ct2", "ct1"],
            "mpra_bc": ["m1", "m2", "m1"],
            "umis_mpra_bc": [0, 4, 7],
        }
    )

    obj = scMPRA_data()
    obj.data = _optimize_mpra_ddf(dd.from_pandas(umiwise, npartitions=1))
    obj.table_type = "mpra_umiwise"

    out_dir = tmp_path / "toy.scmpra"
    obj.to_parquet(str(out_dir))
    loaded = scMPRA_data.from_parquet(str(out_dir))

    loaded_df = loaded.data.compute().sort_values(["cell_bc"]).reset_index(drop=True)
    assert _is_sparse_int64(loaded_df["umis_mpra_bc"])
    assert loaded_df["umis_mpra_bc"].astype("int64").tolist() == [0, 4, 7]


def test_get_cell_counts_accepts_sparse_umis_without_formula_failure():
    df = pd.DataFrame(
        {
            "rep_id": ["r1", "r1", "r2", "r2"],
            "cell_type": ["ctA", "ctA", "ctB", "ctB"],
            "cre_id": ["cre1", "cre2", "cre1", "cre2"],
            "umis_mpra_bc": [0, 2, 1, 3],
        }
    )
    sparse_df = df.copy()
    sparse_df["umis_mpra_bc"] = sparse_df["umis_mpra_bc"].astype(pd.SparseDtype("int64", 0))

    with distributed.Client(processes=False, n_workers=1, threads_per_worker=1) as client:
        got_sparse = get_cell_counts(client, sparse_df, split="cell_type").sort_index()
        got_dense = get_cell_counts(client, df, split="cell_type").sort_index()

    pd.testing.assert_frame_equal(got_sparse, got_dense)


def test_prepare_subset_for_modeling_densifies_umis_for_standard_fit_boundary():
    raw = pd.DataFrame(
        {
            "rep_id": ["r1", "r1", "r2", "r2"],
            "cell_type": ["ctA", "ctA", "ctB", "ctB"],
            "cre_id": ["reference", "creX", "reference", "creX"],
            "umis_mpra_bc": pd.Series([0, 1, 2, 3], dtype=pd.SparseDtype("int64", 0)),
        }
    )

    prepared = _prepare_subset_for_modeling(raw)
    assert not pd.api.types.is_sparse(prepared["umis_mpra_bc"].dtype)
    assert prepared["umis_mpra_bc"].dtype == np.dtype("int64")

    mats = _smart_matrix(prepared, split="cell_type")
    regressand_dtype = mats["regressand"]["umis_mpra_bc"].dtype
    assert not pd.api.types.is_sparse(regressand_dtype)
    assert pd.api.types.is_numeric_dtype(regressand_dtype)
