"""Backend (FastAPI) endpoint tests for the custom-reference + report surface.

Injects the shared processed fixture straight into the session store (skips the slow live pipeline)
and drives the new HTTP endpoints via TestClient. Skips cleanly if fastapi is not installed.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

# backend/ lives at the repo root (a sibling of tests/), not under src - put it on the path.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def _labeled_reference_file(section, tmp_path, n_per=80, k=3, seed=0):
    import anndata as ad
    import pandas as pd
    genes = [str(g) for g in section.var_names]
    ng = len(genes); block = max(1, ng // (k + 1))
    rng = np.random.default_rng(seed)
    X, lab = [], []
    for ci in range(k):
        base = rng.poisson(0.5, (n_per, ng)).astype("float32")
        base[:, ci * block:(ci + 1) * block] += rng.poisson(9, (n_per, block))
        X.append(base); lab += [f"RefType{ci}"] * n_per
    a = ad.AnnData(X=np.vstack(X)); a.var_names = genes
    a.obs_names = [f"r{i}" for i in range(a.n_obs)]; a.layers["counts"] = a.X.copy()
    a.obs["cell_type"] = pd.Categorical(lab)
    p = tmp_path / "ref.h5ad"; a.write_h5ad(p)
    return str(p)


@pytest.fixture
def client_sid(processed_adata):
    pytest.importorskip("fastapi")
    from fastapi.testclient import TestClient
    import backend.app as app
    sid = "pytestsid"
    app._SESSIONS[sid] = {"adata": processed_adata, "tissue": "melanoma", "log": []}
    yield TestClient(app.app), sid, app
    app._SESSIONS.pop(sid, None)


def test_annotation_sig_detects_in_place_mutation(processed_adata):
    # The copilot mutation signal: a merge/relabel (obs cell-type content change) must flip the
    # signature, so the SPA refreshes the map; a no-op (materialising a gene column for colouring)
    # must NOT - a recolour is not a mutation. Fixes the "optimize merging changed nothing on the map".
    import backend.app as app

    key = "cell_type_final" if "cell_type_final" in processed_adata.obs else "cell_type"
    sig0 = app._annotation_sig(processed_adata)
    # colouring the map by a gene writes obs[<gene>] - must be invisible to the mutation signal
    processed_adata.obs["_some_gene"] = 1.0
    assert app._annotation_sig(processed_adata) == sig0
    # a real relabel of the annotation column flips the signature
    import pandas as pd
    s = processed_adata.obs[key].astype(str)
    first = s.iloc[0]
    s.iloc[: (s == first).sum()] = str(first) + " / merged"
    processed_adata.obs[key] = pd.Categorical(s)
    assert app._annotation_sig(processed_adata) != sig0


def test_reference_by_path(client_sid, tmp_path):
    c, sid, _ = client_sid
    # build the reference against the session's adata genes
    import backend.app as app
    refp = _labeled_reference_file(app._SESSIONS[sid]["adata"], tmp_path)
    r = c.post(f"/api/{sid}/reference", json={"path": refp})
    assert r.status_code == 200, r.text
    j = r.json()
    assert j["ok"] and j["label_key"] == "cell_type" and j["n_ref_cells"] == 240
    assert j["match"]["status"] == "ok"
    # the session now carries the reference
    assert app._SESSIONS[sid].get("reference") is not None


def test_reference_bad_path_400(client_sid):
    c, sid, _ = client_sid
    r = c.post(f"/api/{sid}/reference", json={"path": "/does/not/exist.h5ad"})
    assert r.status_code == 400


def test_verify_report(client_sid):
    c, sid, _ = client_sid
    r = c.get(f"/api/{sid}/verify_report")
    assert r.status_code == 200, r.text
    j = r.json()
    assert j["status"] == "ok" and "per_type" in j and "section_agreement" in j


def test_qc_verdict_endpoint_graceful(client_sid, monkeypatch):
    """No API key -> the endpoint still returns 200 with the funnel present and verdict None."""
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    c, sid, _ = client_sid
    r = c.get(f"/api/{sid}/qc_verdict")
    assert r.status_code == 200, r.text
    j = r.json()
    assert j["verdict"] is None and "note" in j                # LLM disabled, graceful


def test_malignant_concordance_endpoint(client_sid):
    c, sid, _ = client_sid
    r = c.post(f"/api/{sid}/run_cap/malignant_concordance", json={"params": {}})
    assert r.status_code == 200, r.text
    j = r.json()
    assert j["ok"] and j["value"]["status"] == "ok"
    names = {cc["name"]: cc["status"] for cc in j["value"]["callers"]}
    assert names.get("marker score") == "ok"


def test_reference_transfer_endpoint(client_sid, tmp_path):
    c, sid, _ = client_sid
    import backend.app as app
    refp = _labeled_reference_file(app._SESSIONS[sid]["adata"], tmp_path)
    c.post(f"/api/{sid}/reference", json={"path": refp})
    r = c.post(f"/api/{sid}/run_cap/reference_transfer", json={"params": {}})
    assert r.status_code == 200, r.text
    j = r.json()
    assert j["ok"] and j["value"]["status"] == "ok"
    # Assert each arm RAN, not merely that it is listed. `capabilities.py` records
    # {"status": "skipped: <exc>"} for a failed arm, so asserting membership alone stays green even
    # when every engine silently no-opped - exactly how the TACCO arm went unrun (see 7d3d20a).
    methods = {a["method"]: a.get("status") for a in j["value"]["arms"]}
    assert methods.get("celltypist") == "ok", f"celltypist arm did not run: {methods!r}"
    assert methods.get("tacco") == "ok", f"tacco arm did not run: {methods!r}"


# --------------------------- bundled demo registry --------------------------- #
# The app ships one processed section: `breast` (Xenium 5K, shallow -> 1/10 typable). Guards the
# registry lookup, not the science; the happy path runs against a temp h5ad injected into DEMO_CACHES.


def test_demos_endpoint_lists_the_section():
    pytest.importorskip("fastapi")
    from fastapi.testclient import TestClient
    import backend.app as app

    demos = TestClient(app.app).get("/api/demos").json()["demos"]
    assert [d["name"] for d in demos] == ["breast"]
    assert all("available" in d for d in demos)


def test_load_demo_unknown_name_is_404_not_500():
    pytest.importorskip("fastapi")
    from fastapi.testclient import TestClient
    import backend.app as app

    r = TestClient(app.app).post("/api/load_demo?name=does_not_exist")
    assert r.status_code == 404
    assert "does_not_exist" in r.json()["detail"]


def test_load_demo_resolves_the_named_cache(processed_adata, tmp_path, monkeypatch):
    """`name` must select its OWN cache path - a regression here would silently serve the breast
    section when the user asked for the positive control (and vice versa)."""
    pytest.importorskip("fastapi")
    from fastapi.testclient import TestClient
    import backend.app as app
    from spatialscribe.analysis import export

    a = processed_adata.copy()
    a.uns["spatialscribe_demo"] = {"tissue": "breast", "source": "fake atera5k",
                                   "role": "positive-control"}
    path = export.export_h5ad(a, tmp_path / "atera5k.h5ad")
    monkeypatch.setitem(app.DEMO_CACHES, "atera5k", str(path))
    monkeypatch.setitem(app.DEMO_CACHES, "breast", str(tmp_path / "absent.h5ad"))

    c = TestClient(app.app)
    j = c.post("/api/load_demo?name=atera5k").json()
    assert j["n_obs"] == a.n_obs
    assert j["demo"] == {"source": "fake atera5k", "role": "positive-control"}
    app._SESSIONS.pop(j["session_id"], None)

    assert c.post("/api/load_demo?name=breast").status_code == 500   # cache missing, not a 404


# ------------------- panel-check marker resolution: CellGuide augmentation ------------------- #
# Panel-check UNIONS CellGuide's canonical markers onto the curated baseline for every category, so
# under-covered curated lineages (Mast 4, NK 4) gain the literature markers the curated dicts never
# named and more land on the panel. CellGuide is a static DB -> deterministic; the LLM is a last
# resort only for labels CellGuide cannot resolve. These pin: curated markers are never discarded,
# curated leads (no dup), abstention labels are excluded, and every real category is grounded.


def _ctx_for(tissue="breast"):
    from spatialscribe.analysis.capabilities import RunContext
    return RunContext(tissue=tissue)


def _adata_with_types(cats):
    import anndata as ad
    import numpy as np
    import pandas as pd
    a = ad.AnnData(X=np.ones((len(cats) * 2, 3), dtype="float32"))
    a.var_names = ["EPCAM", "CD3D", "PECAM1"]
    a.obs["cell_type"] = pd.Categorical(list(cats) * 2)
    return a


def test_curated_categories_are_augmented_with_cellguide(monkeypatch):
    """All-curated section: CellGuide markers are UNIONED onto the curated baseline (curated leads,
    no duplicates), and every non-abstention category is grounded. CellGuide is a static DB, so
    running it for known lineages stays deterministic - the old worry was an LLM's per-run choices."""
    import backend.app as app
    from spatialscribe.analysis import markers as _markers

    curated = _markers.for_tissue(_ctx_for().tissue)
    tcell0 = curated["T cell"][0]
    seen = {}

    def _fake(cell_types, *a, **k):
        seen["cell_types"] = list(cell_types)
        return {"T cell": [tcell0, "IL32", "CD3G"]}   # one curated dup + two CellGuide extras
    monkeypatch.setattr(_markers, "markers_for_types", _fake)

    a = _adata_with_types(["T cell", "Myeloid", "Endothelial", "Unassigned: low quality"])
    ctx = _ctx_for()
    app._apply_llm_markers_for_categories(a, ctx)

    ms = ctx.marker_sets
    assert set(seen["cell_types"]) == {"T cell", "Myeloid", "Endothelial"}   # abstention excluded
    assert ms["T cell"][:len(curated["T cell"])] == curated["T cell"]        # curated FIRST, preserved
    assert "IL32" in ms["T cell"] and "CD3G" in ms["T cell"]                 # CellGuide extras appended
    assert ms["T cell"].count(tcell0) == 1                                   # curated marker not duplicated
    assert "Unassigned: low quality" not in ms                              # abstention never enters


def test_novel_category_grounds_and_merges_onto_curated(monkeypatch):
    """An unanticipated label is grounded (CellGuide/LLM) and MERGED onto the curated baseline.

    Guards the review-caught regression: grounding a novel type must NOT throw away the benchmarked
    curated markers for the known lineages (that regressed panel adequacy 8/10 -> 6/10).
    """
    import backend.app as app
    from spatialscribe.analysis import markers as _markers

    seen = {}

    def _fake(cell_types, *a, **k):
        seen["cell_types"] = list(cell_types)
        return {"Tumor subclone A": ["EPCAM"]}
    monkeypatch.setattr(_markers, "markers_for_types", _fake)

    a = _adata_with_types(["T cell", "Tumor subclone A"])
    ctx = _ctx_for()
    app._apply_llm_markers_for_categories(a, ctx)

    assert set(seen["cell_types"]) == {"T cell", "Tumor subclone A"}   # all non-abstention cats grounded
    assert ctx.marker_sets["Tumor subclone A"] == ["EPCAM"]           # novel grounded
    curated = _markers.for_tissue(ctx.tissue)
    assert ctx.marker_sets["T cell"] == curated["T cell"]            # known lineage keeps curated markers


def test_apply_llm_markers_does_not_raise_on_index_var_names(monkeypatch):
    """Regression (live 500): with the augmentation grounding EVERY category, panel_report passes
    a.var_names (a pandas Index) into markers_for_types -> infer_organism, whose `panel_genes or []`
    raised 'The truth value of a Index is ambiguous'. The REAL path (no mock) must return gracefully
    offline, not 500 - the earlier augmentation tests mocked markers_for_types and so missed it."""
    monkeypatch.setenv("CELLGUIDE_OFFLINE", "1")
    import backend.app as app
    from spatialscribe.analysis import markers as _markers

    a = _adata_with_types(["T cell", "Myeloid"])   # a.var_names is a pandas Index (the crashing input)
    ctx = _ctx_for()
    app._apply_llm_markers_for_categories(a, ctx)   # must not raise
    # curated markers are still applied for a known lineage even with CellGuide offline
    curated = _markers.for_tissue(ctx.tissue)
    assert ctx.marker_sets is not None and ctx.marker_sets["T cell"][:len(curated["T cell"])] == curated["T cell"]


def test_run_pipeline_endpoint_runs_in_background_and_reports(client_sid, monkeypatch):
    # Test the endpoint GLUE (thread spawn, checkbox apply, session update, status merge) with a fast
    # fake - run_pipeline itself is covered by test_pipeline.py, no need to re-run it here.
    import time

    c, sid, app = client_sid
    from spatialscribe.analysis import pipeline as pl

    def _fake(adata, ctx, opts):
        return {"route": "annotate", "stages": [{"name": "compute_qc", "status": "ok"}],
                "action_log": [{"name": "compute_qc", "params": {}}],
                "summary": {"ok": ["compute_qc"], "skipped": [], "failed": []}}

    monkeypatch.setattr(pl, "run_pipeline", _fake)     # late-bound in the thread -> intercepted
    r = c.post(f"/api/{sid}/run_pipeline", json={"tumour": True})
    assert r.status_code == 202 and r.json()["started"] is True

    st = {}
    for _ in range(100):
        st = c.get(f"/api/{sid}/pipeline_status").json()
        if st["done"]:
            break
        time.sleep(0.02)
    assert st["done"] and st["error"] is None, st
    assert st["route"] == "annotate" and st["summary"]["ok"] == ["compute_qc"]
    assert app._SESSIONS[sid]["is_tumour"] is True     # the Data-step checkbox was applied
    assert app._SESSIONS[sid].get("pipeline_running") is False


def test_run_pipeline_guard_409_while_running(client_sid):
    c, sid, app = client_sid
    app._SESSIONS[sid]["pipeline_running"] = True       # simulate an in-flight run (deterministic)
    try:
        assert c.post(f"/api/{sid}/run_pipeline").status_code == 409
        assert c.post(f"/api/{sid}/run_cap/compute_qc", json={"params": {}}).status_code == 409
    finally:
        app._SESSIONS[sid]["pipeline_running"] = False


def test_pipeline_greens_annotate_step_on_any_typing_route():
    # Bug: a full "Run full analysis" that types via reference_transfer / RCTD (e.g. a mouse-brain
    # section labelled from an Allen atlas) left the Annotate rail step - and downstream Report - grey,
    # because is_done("annotate") only fires for the consensus `annotate` cap's full key set
    # (cell_type_final / confidence / verdict / reason). The rail step means "cells are typed", so it
    # must green off cell_type presence, whatever route set it.
    pytest.importorskip("fastapi")
    import anndata as ad
    import pandas as pd
    import backend.app as app
    from spatialscribe.analysis import status

    a = ad.AnnData(X=np.zeros((6, 3), dtype="float32"))
    a.obs["cell_type"] = pd.Categorical(["L6 CT", "VLMC", "L6 CT", "Macrophage", "VLMC", "L6 CT"])
    assert status.is_done(a, "annotate") is False        # reference-transfer route: no confidence keys
    s: dict = {"ran": []}
    app._mark_pipeline_ran(s, a)
    assert "annotate" in s["ran"]                        # greened off the typed outcome

    b = ad.AnnData(X=np.zeros((4, 3), dtype="float32"))  # de-novo cluster route: nothing typed
    s2: dict = {"ran": []}
    app._mark_pipeline_ran(s2, b)
    assert "annotate" not in s2["ran"]                   # stays grey - honest
