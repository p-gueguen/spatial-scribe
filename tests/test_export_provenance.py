"""Tests for export.export_script - the re-runnable provenance script.

Covers: structured {name, params} records render to registry-backed cap.run(...) calls and
the emitted source is valid Python; legacy {label, code} records still emit their snippet
(backward compat); and the tissue kwarg is interpolated into the RunContext header.
"""

from __future__ import annotations

import pytest

from spatialscribe.analysis import export


def test_structured_records_render_and_compile(tmp_path):
    """Structured provenance -> compilable source with a registry cap.run() call + header."""
    action_log = [
        {"name": "compute_qc", "params": {}},
        {"name": "cluster", "params": {"resolution": 1.0}},
    ]
    path = export.export_script(action_log, tmp_path / "rerun.py")
    src = path.read_text()

    # Must be syntactically valid Python.
    compile(src, str(path), "exec")

    # Registry-backed call assembled from the record (params repr'd exactly).
    assert "cap.run(adata, 'cluster', {'resolution': 1.0}, ctx)" in src
    assert "cap.run(adata, 'compute_qc', {}, ctx)" in src
    # Header wiring.
    assert 'RunContext(tissue="melanoma")' in src
    assert "capabilities as cap" in src


def test_legacy_record_still_renders_code_snippet(tmp_path):
    """A legacy {label, code} record falls back to emitting its stored code snippet."""
    action_log = [{"label": "load", "code": "adata = load(...)"}]
    path = export.export_script(action_log, tmp_path / "legacy.py")
    src = path.read_text()

    assert "adata = load(...)" in src
    assert "# --- load ---" in src


def test_export_script_threads_source_path(tmp_path):
    """With a source_path the emitted rerun.py loads the REAL section, not a placeholder - so the
    CLI-written script is one-shot re-runnable (the backend endpoint already did this; the CLI now too)."""
    path = export.export_script([{"name": "compute_qc", "params": {}}],
                                tmp_path / "rerun.py", source_path="/data/my_section")
    src = path.read_text()
    assert "io.load('/data/my_section')" in src
    assert "PATH/TO" not in src                       # no placeholder when the path is known


def test_report_discloses_pipeline_stages(tmp_path):
    """The HTML report lists which optional arms ran vs skipped, so a marker-only run is not
    mistaken for a decontaminated / malignant-called one (honest by disclosure, not omission)."""
    ad = pytest.importorskip("anndata")
    import numpy as np

    a = ad.AnnData(X=np.ones((5, 3), dtype="float32"))
    a.uns["pipeline"] = {"route": "annotate", "stages": [
        {"name": "annotate", "status": "ok"},
        {"name": "split_purify", "status": "skipped:no reference"},
    ]}
    out = export.render_analysis_report(a, [], "breast", tmp_path / "report.html")
    html = out.read_text()
    assert "Pipeline stages" in html
    assert "split_purify" in html and "skipped" in html and "no reference" in html


def test_tissue_kwarg_interpolated_into_header(tmp_path):
    """The tissue kwarg flows into the RunContext line in the generated header."""
    path = export.export_script(
        [{"name": "compute_qc", "params": {}}], tmp_path / "breast.py", tissue="breast"
    )
    src = path.read_text()

    assert 'RunContext(tissue="breast")' in src
    assert 'RunContext(tissue="melanoma")' not in src
