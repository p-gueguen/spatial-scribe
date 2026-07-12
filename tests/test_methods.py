from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")


def _tiny(n=20):
    a = anndata.AnnData(X=np.zeros((n, 3), dtype="float32"))
    a.obs_names = [f"cell{i}" for i in range(n)]
    return a


def test_join_parquet_reindexes_on_cell_id(tmp_path):
    from spatialscribe.analysis import methods

    a = _tiny(20)
    # Parquet covers only the first 10 cells and is out of order.
    df = pd.DataFrame({"cell_id": [f"cell{i}" for i in [3, 1, 0, 2, 4, 5, 6, 7, 8, 9]],
                       "label": ["A"] * 10, "conf": np.linspace(0, 1, 10)})
    p = tmp_path / "m.parquet"
    df.to_parquet(p)
    info = methods.join_parquet(a, p, columns=["label", "conf"], prefix="scanvi")
    assert "scanvi_label" in a.obs and "scanvi_conf" in a.obs
    assert a.obs["scanvi_label"].iloc[0] == "A"          # cell0 present
    assert a.obs["scanvi_label"].iloc[15] is None or (a.obs["scanvi_label"].isna().iloc[15])  # cell15 uncovered
    assert 0.4 < info["coverage"] < 0.6                  # 10/20 covered


def test_join_parquet_missing_column_degrades_to_zero_coverage(tmp_path):
    import anndata, numpy as np, pandas as pd
    from spatialscribe.analysis import methods
    a = anndata.AnnData(X=np.zeros((3, 2), dtype="float32"))
    a.obs_names = ["c0", "c1", "c2"]
    pd.DataFrame({"cell_id": ["c0", "c1", "c2"], "other": [1, 2, 3]}).to_parquet(tmp_path / "m.parquet")
    info = methods.join_parquet(a, tmp_path / "m.parquet", columns=["label", "conf"], prefix="x")
    assert info["coverage"] == 0.0            # neither requested column present -> no crash, 0 coverage
    assert a.obs["x_label"].isna().all()


def test_join_panhumanpy_schema(tmp_path):
    import anndata, numpy as np, pandas as pd
    from spatialscribe.analysis import methods

    a = anndata.AnnData(X=np.zeros((5, 2), dtype="float32"))
    a.obs_names = [f"c{i}" for i in range(5)]
    pd.DataFrame({"cell_id": a.obs_names.tolist(),
                  "broad": ["T"] * 5, "medium": ["CD8 T"] * 5, "fine": ["CD8 Tem"] * 5,
                  "confidence": np.linspace(0, 1, 5)}).to_parquet(tmp_path / "ph.parquet")
    info = methods.join_panhumanpy(a, tmp_path / "ph.parquet")
    assert {"ph_broad", "ph_medium", "ph_fine", "ph_confidence"} <= set(a.obs.columns)
    assert info["coverage"] == 1.0


def test_join_rctd_schema(tmp_path):
    import anndata, numpy as np, pandas as pd
    from spatialscribe.analysis import methods

    a = anndata.AnnData(X=np.zeros((4, 2), dtype="float32"))
    a.obs_names = [f"c{i}" for i in range(4)]
    pd.DataFrame({"cell_id": ["c0", "c1", "c2", "c3"],
                  "first_type": ["T", "Mel", "T", "B"], "second_type": [None, "T", None, None],
                  "spot_class": ["singlet", "doublet_certain", "singlet", "singlet"],
                  "weight": [0.9, 0.6, 0.95, 0.8],
                  "singlet_score": [12.0, 8.0, 15.0, 10.0]}).to_parquet(tmp_path / "r.parquet")
    info = methods.join_rctd(a, tmp_path / "r.parquet")
    assert {"rctd_first_type", "rctd_second_type", "rctd_spot_class", "rctd_weight", "rctd_singlet_score"} <= set(a.obs.columns)
    assert info["coverage"] == 1.0


def test_join_singler_schema(tmp_path):
    import anndata, numpy as np, pandas as pd
    from spatialscribe.analysis import methods

    a = anndata.AnnData(X=np.zeros((4, 2), dtype="float32"))
    a.obs_names = [f"c{i}" for i in range(4)]
    pd.DataFrame({"cell_id": ["c0", "c1", "c2", "c3"],
                  "label": ["T", "Mel", "B", "T"],
                  "delta": [0.4, 0.1, 0.35, 0.2]}).to_parquet(tmp_path / "s.parquet")
    info = methods.join_singler(a, tmp_path / "s.parquet")
    assert {"singler_label", "singler_delta"} <= set(a.obs.columns)
    assert info["coverage"] == 1.0


def test_join_scanvi_schema(tmp_path):
    import anndata, numpy as np, pandas as pd
    from spatialscribe.analysis import methods

    a = anndata.AnnData(X=np.zeros((4, 2), dtype="float32"))
    a.obs_names = [f"c{i}" for i in range(4)]
    pd.DataFrame({"cell_id": ["c0", "c1", "c2", "c3"],
                  "label": ["T", "Mel", "B", "T"],
                  "confidence": [0.9, 0.6, 0.8, 0.7],
                  "entropy": [0.2, 0.7, 0.4, 0.5]}).to_parquet(tmp_path / "sc.parquet")
    info = methods.join_scanvi(a, tmp_path / "sc.parquet")
    assert {"scanvi_label", "scanvi_confidence", "scanvi_entropy"} <= set(a.obs.columns)
    assert info["coverage"] == 1.0
