"""Tests for the synthetic melanoma demo fixture and the immune-exclusion hero."""

from __future__ import annotations

import pytest

pytest.importorskip("anndata")
pytest.importorskip("squidpy")


def test_load_demo_shape(monkeypatch):
    monkeypatch.setenv("SPATIALSCRIBE_FORCE_CPU", "1")
    from spatialscribe.analysis import demo

    s = demo.load_demo()
    assert s.platform == "xenium-demo"
    assert s.adata.n_obs == 3000
    assert len(s.panel_genes) > 0
    assert s.control_mask.sum() >= 10          # control probes present
    assert "Malignant/Melanocyte" in set(s.adata.obs["true_type"])
    assert "T cell" in set(s.adata.obs["true_type"])


def test_demo_is_immune_excluded(monkeypatch):
    """The demo is built so T cells sit at the tumor margin -> negative enrichment."""
    monkeypatch.setenv("SPATIALSCRIBE_FORCE_CPU", "1")
    from spatialscribe.analysis import demo, spatial

    a = demo.load_demo().adata
    a.obs["cell_type"] = a.obs["true_type"]
    res = spatial.immune_exclusion(a, "Malignant/Melanocyte", "T cell", cluster_key="cell_type")
    assert res["zscore"] < 0
    assert "exclud" in res["verdict"] or res["zscore"] < 0
