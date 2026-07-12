"""Shared, offline, deterministic test substrate for the whole suite.

One golden synthetic section (built once from the demo generator) is driven through the
registry to produce the raw and processed AnnData fixtures every test reuses. Everything
runs CPU-only with no API key, so the suite is fast, hermetic, and reproducible.
"""

from __future__ import annotations

import os

os.environ.setdefault("SPATIALSCRIBE_FORCE_CPU", "1")

import warnings

warnings.filterwarnings("ignore")

import pytest


@pytest.fixture(autouse=True)
def _offline(monkeypatch, tmp_path):
    """Pin the CPU path; guarantee no live LLM OR CellGuide network call leaks into a test."""
    monkeypatch.setenv("SPATIALSCRIBE_FORCE_CPU", "1")
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    # markers_for_types now consults CellGuide; offline it would urlopen the live snapshot. Pin it
    # offline AND point the disk cache at an EMPTY per-test dir: the real cache (populated on the
    # deploy under /data) must NOT leak in, or "offline => no markers" silently becomes "offline =>
    # cached markers" and breaks the curated-fallback contract. _CACHE is read at import, so patch the
    # module attribute, not just the env. test_cellguide.py seeds its own cache to exercise the hit.
    monkeypatch.setenv("CELLGUIDE_OFFLINE", "1")
    from spatialscribe.analysis import cellguide as _cg
    monkeypatch.setattr(_cg, "_CACHE", tmp_path / "_cellguide_empty")
    _cg._reset_for_tests()


@pytest.fixture
def ctx():
    from spatialscribe.analysis import capabilities as cap
    return cap.RunContext(tissue="melanoma", use_llm=False)


@pytest.fixture
def raw_adata():
    """A small synthetic Xenium-like section: raw counts, control mask, obsm['spatial']."""
    from spatialscribe.analysis import demo
    from spatialscribe.analysis.io import build_control_mask

    a = demo.make_demo_adata(n_cells=600, seed=0)
    a.var["control"] = build_control_mask(a)
    return a


@pytest.fixture
def processed_adata(raw_adata, ctx):
    """``raw_adata`` driven through compute_qc -> panel_check -> cluster -> annotate.

    Produced entirely via ``capabilities.run`` (the same path both frontends use), so the
    fixture also exercises the happy-path pipeline on every test that consumes it.
    """
    from spatialscribe.analysis import capabilities as cap

    for name, params in [("compute_qc", {}), ("panel_check", {}),
                         ("cluster", {"resolution": 1.0}), ("annotate", {})]:
        res = cap.run(raw_adata, name, params, ctx)
        assert res.ok, (name, res.error)
    return raw_adata


_LOAD_SECTION_FIXTURE: str | None = None


def _load_section_fixture_path() -> str:
    """A tiny on-disk .h5ad section for the load_section conformance run (written once, memoized).
    load_section loads a NEW section from a path rather than transforming the golden fixture, so the
    registry conformance needs a real, minimal section to point it at."""
    global _LOAD_SECTION_FIXTURE
    if _LOAD_SECTION_FIXTURE is None:
        import tempfile

        import anndata as ad
        import numpy as np
        p = os.path.join(tempfile.mkdtemp(prefix="ss_confsec_"), "section.h5ad")
        a = ad.AnnData(X=np.zeros((10, 3), dtype="float32"))
        a.obs_names = [f"c{i}" for i in range(10)]
        a.var_names = [f"g{i}" for i in range(3)]
        a.obsm["spatial"] = np.tile(np.arange(10.0)[:, None], (1, 2))
        a.write_h5ad(p)
        _LOAD_SECTION_FIXTURE = p
    return _LOAD_SECTION_FIXTURE


def capability_params(name: str, adata) -> dict:
    """Valid params for capability ``name`` given a processed ``adata`` (test helper)."""
    cats = list(adata.obs["cell_type"].astype("category").cat.categories) if "cell_type" in adata.obs else []
    if name == "load_section":
        # A real section path (auto_reference off - selection is exercised in test_load_section).
        return {"path": _load_section_fixture_path(), "auto_reference": False}
    if name == "immune_exclusion":
        a, b = (cats[0], cats[-1]) if len(cats) >= 2 else (cats[0], cats[0]) if cats else ("x", "y")
        return {"type_a": a, "type_b": b}
    if name == "subcluster":
        if "cell_type" in adata.obs:
            return {"cell_type": str(adata.obs["cell_type"].value_counts().index[0])}
        return {"cell_type": "T cell"}   # stub for the prereq test (run() rejects before use)
    if name == "cluster":
        return {"resolution": 1.0}
    if name == "discover_programs":
        return {"n_programs": 4}
    if name == "highlight_cells":
        return {"criterion": "low quality"}   # always resolvable (QC count-floor / lowest-decile)
    if name == "show_spatial":
        return {"color_by": "cell_type"}      # a real field on the processed fixture
    if name == "marker_dotplot":
        return {"genes": [str(g) for g in adata.var_names[:3]]}   # on-panel genes from the fixture
    if name == "expression_violin":
        return {"gene": str(adata.var_names[0])}
    return {}
