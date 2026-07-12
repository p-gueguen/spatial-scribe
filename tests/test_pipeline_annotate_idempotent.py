"""Regression tests for the annotation-consistency fix (2026-07-12).

1. run_pipeline must NOT re-annotate a section that already carries `cell_type` - re-running the
   de-novo cluster-marker annotator clobbered higher-quality reference labels (the atera RCTD cache
   collapsed its 10 per-cell types to a 4-type per-cluster call, and every demo showed different
   labels before vs after "Run full analysis"). The guard mirrors the reference_transfer branch.
2. qc.run_funnel must LOG (not silently swallow) an annotation_quality failure - a bare `except: pass`
   turned a version-skew into a mystery empty AQI tile behind an HTTP 200.
"""
import logging

import pytest

from spatialscribe.analysis import capabilities as cap
from spatialscribe.analysis import demo, pipeline, qc


def _preannotated(n=2000):
    a = demo.make_demo_adata(n_cells=n, seed=0)
    # give it a rich, pre-existing cell_type (as if loaded from a reference-annotated cache)
    src = a.obs["true_type"] if "true_type" in a.obs else a.obs.index.to_series().map(lambda _: "T cell")
    a.obs["cell_type"] = src.astype("category")
    return a


def test_run_pipeline_preserves_existing_cell_type():
    a = _preannotated()
    before = set(a.obs["cell_type"].astype(str).unique())
    assert len(before) >= 3, "fixture should carry several types so a collapse would be visible"
    ctx = cap.RunContext(tissue="melanoma", use_llm=False)
    pipeline.run_pipeline(a, ctx, pipeline.PipelineOptions(export=False))
    after = set(a.obs["cell_type"].astype(str).unique())
    # the de-novo annotator would have overwritten cell_type with a per-cluster call (fewer types);
    # the guard keeps the provided labels verbatim.
    assert after == before, f"cell_type was re-annotated (clobbered): {before} -> {after}"


def test_run_pipeline_annotates_a_raw_section():
    """The guard must not over-fire: a RAW section (no cell_type) is still annotated."""
    a = demo.make_demo_adata(n_cells=2000, seed=1)
    assert "cell_type" not in a.obs
    ctx = cap.RunContext(tissue="melanoma", use_llm=False)
    pipeline.run_pipeline(a, ctx, pipeline.PipelineOptions(export=False))
    assert "cell_type" in a.obs and a.obs["cell_type"].nunique() >= 1


def test_run_funnel_logs_instead_of_swallowing_aqi_failure(monkeypatch, caplog):
    from spatialscribe.analysis import eval_metrics as _em
    from spatialscribe.analysis import markers

    a = demo.make_demo_adata(n_cells=2000, seed=0)
    ctx = cap.RunContext(tissue="melanoma", use_llm=False)
    for step in ("compute_qc", "cluster", "annotate"):
        cap.run(a, step, {}, ctx)

    def _boom(*args, **kw):
        raise TypeError("annotation_quality() got an unexpected keyword argument (simulated skew)")

    monkeypatch.setattr(_em, "annotation_quality", _boom)
    with caplog.at_level(logging.WARNING):
        headline = qc.run_funnel(a, cluster_key="cell_type", marker_sets=markers.for_tissue("melanoma"))
    # degrade-graceful: the funnel still returns, the AQI is just absent...
    assert "annotation_quality" not in headline
    assert "quality_metrics" not in headline.get("layers_run", [])
    # ...but it is NOT silent anymore.
    assert any("AQI tile will be empty" in r.message for r in caplog.records), \
        "annotation_quality failure was swallowed silently (the bug)"


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
