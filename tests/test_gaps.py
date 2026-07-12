"""Tests for the biology-gap modules: de-novo programs and tumor calling."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")
pytest.importorskip("scanpy")


def _synthetic(n=400, seed=0):
    rng = np.random.default_rng(seed)
    t = ["CD3D", "CD3E", "CD8A", "TRAC"]
    mel = ["MLANA", "PMEL", "TYR", "DCT", "MITF", "SOX10"]
    noise = [f"N{i}" for i in range(30)]
    genes = t + mel + noise
    grp = np.array([0] * (n // 2) + [1] * (n - n // 2))       # 0=T, 1=melanocyte
    X = rng.poisson(0.3, (n, len(genes))).astype("float32")
    gi = {g: i for i, g in enumerate(genes)}
    for g in t:
        X[grp == 0, gi[g]] += rng.poisson(8, (grp == 0).sum())
    for g in mel:
        X[grp == 1, gi[g]] += rng.poisson(8, (grp == 1).sum())
    a = anndata.AnnData(X=X, var=pd.DataFrame(index=genes))
    a.obs["cell_type"] = pd.Categorical(np.where(grp == 0, "T cell", "Malignant/Melanocyte"))
    a.obsm["spatial"] = rng.normal(0, 50, (n, 2))
    return a, grp


def test_discover_programs(monkeypatch):
    monkeypatch.setenv("SPATIALSCRIBE_FORCE_CPU", "1")
    import scanpy as sc

    from spatialscribe.analysis import programs

    a, _ = _synthetic()
    sc.pp.normalize_total(a); sc.pp.log1p(a)
    tab = programs.discover_programs(a, n_programs=4)
    assert "programs" in a.obsm and a.obsm["programs"].shape == (a.n_obs, 4)
    assert "program" in a.obs
    assert len(tab) == 4 and all(tab["top_genes"].map(len) > 0)


def test_malignant_score_higher_in_tumor(monkeypatch):
    monkeypatch.setenv("SPATIALSCRIBE_FORCE_CPU", "1")
    from spatialscribe.analysis import cnv

    a, grp = _synthetic()
    info = cnv.malignant_score(a, tissue="melanoma")
    assert info["status"] == "ok"
    s = np.asarray(a.obs["malignant_score"])
    assert s.min() >= 0 and s.max() <= 1
    # Melanocyte (grp==1) cells should score higher than T cells on average.
    assert s[grp == 1].mean() > s[grp == 0].mean()


def test_cnv_degrades_gracefully_when_env_missing():
    """CNV calling never raises: a missing cnv_env env returns a 'skipped' status.

    (infercnvpy runs in an isolated subprocess env, not in-process, so unavailability is a
    graceful skip - the app falls back to the marker-based malignant_score - not an exception.)
    """
    from spatialscribe.analysis import cnv

    a, _ = _synthetic()
    r = cnv.call_malignant_cnv(a, env_python="/nonexistent/cnv_env/python")
    assert r["status"].startswith("skipped")
    assert r["pct_malignant"] == 0.0


def test_cnv_join_parquet_plumbing(tmp_path):
    """_join_cnv reindexes a {cell_id, cnv_score, is_malignant} parquet onto obs (subset covered)."""
    import numpy as np
    import pandas as pd
    import pyarrow as pa
    import pyarrow.parquet as pq

    from spatialscribe.analysis import cnv

    a, _ = _synthetic()
    ids = a.obs_names.astype(str).tolist()
    sub = ids[: len(ids) // 2][::-1]                     # cover the first half, in reversed order
    df = pd.DataFrame({"cell_id": sub,
                       "cnv_score": np.linspace(0.1, 0.5, len(sub)),
                       "is_malignant": [i % 3 == 0 for i in range(len(sub))]})
    p = tmp_path / "cnv.parquet"
    pq.write_table(pa.Table.from_pandas(df, preserve_index=False), str(p))

    n = cnv._join_cnv(a, p)
    assert n == len(sub)
    assert "cnv_score" in a.obs and "is_malignant" in a.obs
    assert np.isnan(np.asarray(a.obs["cnv_score"], dtype=float)[len(ids) - 1])   # uncovered -> NaN
    assert int(a.obs["is_malignant"].astype(bool).sum()) == int(df["is_malignant"].sum())  # uncovered -> False


def test_cnv_reference_marker_pure_filter():
    """The marker-pure filter drops impure reference cells; tumor cells are never reference."""
    import anndata
    import numpy as np
    import pandas as pd

    from spatialscribe.analysis import cnv

    genes = ["CD3D", "CD3E", "EPCAM", "KRT8"]
    rng = np.random.default_rng(0)
    X = rng.poisson(0.3, size=(200, 4)).astype("float32")
    X[:50, [0, 1]] += 20        # 50 pure T cells (CD3 high)
    X[50:100, [2, 3]] += 20     # 50 'T cell'-labeled cells that are actually tumor (epithelial) -> impure
    X[100:, [2, 3]] += 20       # 100 tumor cells
    a = anndata.AnnData(X=X)
    a.var_names = genes
    a.obs["cell_type"] = pd.Categorical(["T cell"] * 100 + ["Tumor"] * 100)
    a.layers["counts"] = a.X.copy()
    markers = {"T cell": ["CD3D", "CD3E"], "Tumor": ["EPCAM", "KRT8"]}

    lbl0, _ = cnv._build_reference(a, "cell_type", None, None, 0.5, 20)
    assert int((lbl0 == "reference").sum()) == 100                # unfiltered: all 100 'T cell'
    lbl1, _ = cnv._build_reference(a, "cell_type", None, markers, 0.5, min_reference=20)
    n_ref = int((lbl1 == "reference").sum())
    assert 40 <= n_ref < 100                                      # impure 'T cells' dropped
    assert int((lbl1[100:] == "reference").sum()) == 0           # tumor cells never in the reference


def test_cancerfinder_degrades_gracefully_when_unconfigured():
    """Cancer-Finder never raises: a missing env / repo / checkpoint returns a 'skipped' status."""
    from spatialscribe.analysis import cancerfinder

    a, _ = _synthetic()
    r = cancerfinder.call_cancerfinder(a, env_python="/nonexistent/py", repo="/nope", ckpt="/nope")
    assert r["status"].startswith("skipped")
    assert r["pct_malignant"] == 0.0


def test_unconfigured_paths_skip_with_envvar_hint(monkeypatch):
    """With no site paths configured (the public-checkout default), both learned callers skip
    with a clear 'set <ENV_VAR>' message - never a hardcoded /data path and never a raise.

    Guards the self-contained contract: the committed defaults are empty env-var overrides, so a
    fresh clone with SPATIALSCRIBE_CNV_* / CANCERFINDER_* unset degrades to marker-based calling.
    """
    from spatialscribe.analysis import cancerfinder, cnv

    for var in ("SPATIALSCRIBE_CNV_PYTHON", "SPATIALSCRIBE_CNV_GTF", "SPATIALSCRIBE_CNV_LIB"):
        monkeypatch.setattr(cnv, {"SPATIALSCRIBE_CNV_PYTHON": "_GI_INSITUCNV_PY",
                                  "SPATIALSCRIBE_CNV_GTF": "_GTF_DEFAULT",
                                  "SPATIALSCRIBE_CNV_LIB": "_GI_INSITUCNV_LIB"}[var], "")
    for attr in ("_CF_PY", "_CF_REPO", "_CF_CKPT"):
        monkeypatch.setattr(cancerfinder, attr, "")

    a, _ = _synthetic()
    r_cnv = cnv.call_malignant_cnv(a)
    assert r_cnv["status"].startswith("skipped") and "SPATIALSCRIBE_CNV" in r_cnv["status"]
    r_cf = cancerfinder.call_cancerfinder(a)
    assert r_cf["status"].startswith("skipped") and "CANCERFINDER_" in r_cf["status"]


def test_cancerfinder_join_plumbing(tmp_path):
    """_join_cf reindexes {cell_id, cancerfinder_prob} onto obs + thresholds is_malignant."""
    import numpy as np
    import pandas as pd
    import pyarrow as pa
    import pyarrow.parquet as pq

    from spatialscribe.analysis import cancerfinder

    a, _ = _synthetic()
    ids = a.obs_names.astype(str).tolist()
    sub = ids[: len(ids) // 2]
    df = pd.DataFrame({"cell_id": sub, "cancerfinder_prob": np.linspace(0.1, 0.9, len(sub))})
    p = tmp_path / "cf.parquet"
    pq.write_table(pa.Table.from_pandas(df, preserve_index=False), str(p))

    n = cancerfinder._join_cf(a, p, threshold=0.5)
    assert n == len(sub)
    assert "cancerfinder_prob" in a.obs and "cancerfinder_malignant" in a.obs
    assert np.isnan(np.asarray(a.obs["cancerfinder_prob"], dtype=float)[len(ids) - 1])   # uncovered -> NaN
    assert int(a.obs["cancerfinder_malignant"].astype(bool).sum()) == int((df["cancerfinder_prob"] > 0.5).sum())
