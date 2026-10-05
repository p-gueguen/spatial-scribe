"""The 2026-10 speed-ups must not change results: each fast path is checked against the slow one."""
from __future__ import annotations

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp


def test_lut_rgb_matches_per_cell_lookup():
    from backend.app import _lut_rgb

    rng = np.random.default_rng(0)
    vals = rng.choice(["T cell", "B cell", "nan", "Tumour"], 500).astype(str)
    pal = {"T cell": [1, 2, 3], "B cell": [4, 5, 6], "Tumour": [7, 8, 9]}
    slow = np.array([pal.get(v, [0, 0, 0]) for v in vals], dtype=np.uint8).reshape(-1).tolist()
    assert _lut_rgb(vals, pal, [0, 0, 0]) == slow


def test_gzip_and_cache_headers():
    from fastapi.testclient import TestClient

    import backend.app as app_mod

    c = TestClient(app_mod.app)
    r = c.get("/api/demos", headers={"Accept-Encoding": "gzip"})
    assert r.status_code == 200
    if not (app_mod._DIST / "index.html").exists():
        pytest.skip("no built SPA to check cache headers on")
    big = c.get("/", headers={"Accept-Encoding": "gzip"})
    assert big.headers.get("cache-control") == "no-cache"
    asset = next((app_mod._DIST / "assets").glob("*.js")).name
    a = c.get(f"/assets/{asset}", headers={"Accept-Encoding": "gzip"})
    assert "immutable" in a.headers.get("cache-control", "")
    assert a.headers.get("content-encoding") == "gzip"


def test_marker_fidelity_on_marker_slice_is_identical():
    from spatialscribe.analysis import eval_metrics as em

    rng = np.random.default_rng(1)
    a = ad.AnnData(X=sp.csr_matrix(rng.poisson(1.0, (300, 40)).astype(np.float32)))
    a.var_names = [f"g{i}" for i in range(40)]
    a.obs["lab"] = rng.choice(["A", "B", "C"], 300)
    sets = {"A": ["g1", "g5"], "B": ["g7", "g30"], "C": ["g2"]}
    mk = sorted({g for gs in sets.values() for g in gs})
    assert em.marker_program_fidelity(a[:, mk], "lab", sets) == em.marker_program_fidelity(a, "lab", sets)


def test_reference_match_memo_is_a_copy_and_invalidates_on_relabel(monkeypatch):
    from spatialscribe.analysis import reference as ref

    calls = []
    monkeypatch.setattr(ref, "_reference_panel_match",
                        lambda r, g, k, d, t: calls.append(1) or {"global": {"n": len(calls)}})
    monkeypatch.setattr(ref, "_MATCH_MEMO", {})
    r = ad.AnnData(X=np.zeros((4, 2)))
    r.obs["ct"] = ["a", "a", "b", "b"]
    m1 = ref.reference_panel_match(r, ["x"], "ct")
    m1["mutated"] = True                                  # callers annotate the result in place
    m2 = ref.reference_panel_match(r, ["x"], "ct")
    assert len(calls) == 1 and "mutated" not in m2
    r.obs["ct"] = ["a", "b", "b", "b"]                    # in-place relabel -> recompute
    ref.reference_panel_match(r, ["x"], "ct")
    assert len(calls) == 2


def test_celltypist_model_reused_for_same_reference(monkeypatch):
    celltypist = pytest.importorskip("celltypist")
    from spatialscribe.analysis import annotate

    rng = np.random.default_rng(2)
    genes = [f"g{i}" for i in range(20)]
    r = ad.AnnData(X=rng.poisson(2.0, (60, 20)).astype(np.float32))
    r.var_names = genes
    r.obs["ct"] = np.repeat(["A", "B", "C"], 20)
    a = ad.AnnData(X=rng.poisson(2.0, (30, 20)).astype(np.float32))
    a.var_names = genes
    trains = []
    real = celltypist.train
    monkeypatch.setattr(celltypist, "train", lambda *x, **k: trains.append(k["n_jobs"]) or real(*x, **k))
    monkeypatch.setattr(annotate, "_CELLTYPIST_MODEL", {})
    assert annotate.celltypist_transfer(a, r, "ct")["status"] == "ok"
    assert annotate.celltypist_transfer(a, r, "ct")["status"] == "ok"
    assert trains == [16]                                 # trained once, multi-core
    r2 = r.copy()
    annotate.celltypist_transfer(a, r2, "ct")             # a different reference object retrains
    assert len(trains) == 2


def test_rds_cache_never_shares_an_entry_between_paths(tmp_path, monkeypatch):
    import os

    from spatialscribe.analysis import reference as ref

    monkeypatch.setenv("SPATIALSCRIBE_REF_CACHE", str(tmp_path / "cache"))
    monkeypatch.setattr(ref, "_rctd_rds_to_adata",
                        lambda p: ad.AnnData(X=np.zeros((1, 1)), obs={"src": [str(p)]}))
    a, b = tmp_path / "a" / "ref.rds", tmp_path / "b" / "ref.rds"   # same name, size, mtime
    for f in (a, b):
        f.parent.mkdir()
        f.write_bytes(b"x")
        os.utime(f, ns=(1, 1))
    assert ref._rds_cached(a).obs["src"][0] == str(a)
    assert ref._rds_cached(b).obs["src"][0] == str(b)      # not a's cached entry
    assert ref._rds_cached(a).obs["src"][0] == str(a)      # served from the cache, still a's
    assert len(os.listdir(tmp_path / "cache")) == 2
