"""Tests for io._read_cells_table - the polars fast-load with column projection.

Covers: column projection (only _CELL_COLS materialized, junk dropped, file order kept),
cell_id indexing, parity with a full pandas read, the 'first column becomes the index'
fallback when cell_id is absent, and the pandas fallback on an unreadable/non-parquet path.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import polars as pl
import pytest

from spatialscribe.analysis.io import _read_cells_table


def _write_parquet(path: Path, data: dict) -> Path:
    pl.DataFrame(data).write_parquet(path)
    return path


def test_projection_drops_junk_and_indexes_cell_id(tmp_path):
    """cell_id becomes the index; _CELL_COLS are kept; a non-panel junk column is dropped."""
    path = _write_parquet(
        tmp_path / "cells.parquet",
        {
            "cell_id": ["c0", "c1", "c2"],
            "x_centroid": [1.0, 2.0, 3.0],
            "y_centroid": [4.0, 5.0, 6.0],
            "transcript_counts": [10, 20, 30],
            "junk": ["a", "b", "c"],   # not in _CELL_COLS -> must be projected away
        },
    )

    out = _read_cells_table(path)

    assert isinstance(out, pd.DataFrame)
    assert out.index.name == "cell_id"
    # rows preserved in file order
    assert list(out.index) == ["c0", "c1", "c2"]
    # projected columns present, junk dropped by column projection
    assert "x_centroid" in out.columns
    assert "y_centroid" in out.columns
    assert "transcript_counts" in out.columns
    assert "junk" not in out.columns
    assert "cell_id" not in out.columns   # became the index


def test_parity_with_pandas_full_read(tmp_path):
    """The projected columns' values match a straight pandas.read_parquet of the same file."""
    path = _write_parquet(
        tmp_path / "cells.parquet",
        {
            "cell_id": ["c0", "c1", "c2", "c3"],
            "x_centroid": [1.5, 2.5, 3.5, 4.5],
            "y_centroid": [10.0, 11.0, 12.0, 13.0],
            "transcript_counts": [7, 8, 9, 10],
            "junk": [0, 1, 2, 3],
        },
    )

    out = _read_cells_table(path)
    ref = pd.read_parquet(path)  # full, unprojected read (file order)

    for col in ("x_centroid", "y_centroid", "transcript_counts"):
        assert out[col].tolist() == ref[col].tolist()
    # index parity: the fast path's index equals the cell_id column of the full read
    assert list(out.index) == ref["cell_id"].tolist()


def test_first_column_becomes_index_when_no_cell_id(tmp_path):
    """With no 'cell_id', the first column ('id' here) is used as the index."""
    path = _write_parquet(
        tmp_path / "cells.parquet",
        {
            "id": ["a0", "a1"],           # first column -> index
            "x_centroid": [1.0, 2.0],
            "y_centroid": [3.0, 4.0],
        },
    )

    out = _read_cells_table(path)

    assert out.index.name == "id"
    assert list(out.index) == ["a0", "a1"]
    assert "x_centroid" in out.columns
    assert "id" not in out.columns


def test_pandas_fallback_returns_when_polars_path_fails(tmp_path, monkeypatch):
    """If the polars fast path raises, the except-branch pandas read still returns a valid frame."""
    path = _write_parquet(
        tmp_path / "cells.parquet",
        {
            "cell_id": ["c0", "c1"],
            "x_centroid": [1.0, 2.0],
            "y_centroid": [3.0, 4.0],
        },
    )

    def _boom(*args, **kwargs):
        raise RuntimeError("forced polars failure")

    # _read_cells_table does `import polars as pl; pl.scan_parquet(...)`; break that hop.
    monkeypatch.setattr(pl, "scan_parquet", _boom)

    out = _read_cells_table(path)

    # Came back via the pandas fallback, still indexed by cell_id.
    assert isinstance(out, pd.DataFrame)
    assert out.index.name == "cell_id"
    assert list(out.index) == ["c0", "c1"]


def test_unreadable_file_raises_pandas_error_not_polars(tmp_path):
    """A non-parquet/corrupt path exercises the pandas fallback; the surfaced error is a
    pandas/pyarrow error, never a polars-specific one (the polars error is swallowed)."""
    bad = tmp_path / "not_a_parquet.parquet"
    bad.write_bytes(b"this is plainly not a parquet file")

    with pytest.raises(Exception) as excinfo:
        _read_cells_table(bad)

    # The polars ComputeError was caught; whatever escapes comes from the pandas read.
    assert not type(excinfo.value).__module__.startswith("polars")
