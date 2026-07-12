"""run_pipeline: the full analysis spine as one shared function (CLI + app), CPU / no-LLM."""
from __future__ import annotations

import pytest

pytest.importorskip("scanpy")
pytest.importorskip("squidpy")

from spatialscribe.analysis import capabilities as cap
from spatialscribe.analysis import demo, pipeline


def _status(res):
    return {s["name"]: s["status"] for s in res["stages"]}


def test_run_pipeline_runs_the_spine_and_gates_malignant_on_tumour():
    adata = demo.load_demo().adata
    ctx = cap.RunContext(tissue="melanoma", use_llm=False, is_tumour=True)
    res = pipeline.run_pipeline(adata, ctx, pipeline.PipelineOptions(export=False))
    st = _status(res)
    # the always-on spine must actually RUN (status ok), not merely be listed (silent-skip gotcha)
    for stage in ("compute_qc", "panel_check", "annotation_strategy", "cluster", "niches"):
        assert st.get(stage) == "ok", (stage, st)
    assert "cell_type" in adata.obs and "niche" in adata.obs
    # honest record: every status is in the allowed vocabulary, and nothing hard-failed on the demo
    assert all(s["status"].split(":", 1)[0] in {"ok", "skipped", "failed"} for s in res["stages"])
    assert res["summary"]["failed"] == [], res["stages"]
    # tumour context -> malignant was NOT gated out at the pipeline level
    assert not st.get("malignant_concordance", "").startswith("skipped:non-tumour")
    assert adata.uns["pipeline"]["route"] in {"annotate", "reference_transfer", "cluster"}
    # the qc_funnel stage RAN (not merely listed) -> per-cell confidence + the section AQI are
    # produced, so a headless report carries them (both were absent before this stage existed).
    assert st.get("qc_funnel") == "ok", st
    assert {"annotation_confidence", "annotation_verdict", "cell_type_final"} <= set(adata.obs.columns)
    aq = adata.uns.get("annotation_quality")
    assert isinstance(aq, dict) and isinstance(aq.get("aqi"), dict), aq
    assert aq["aqi"].get("aqi") is not None


def test_run_pipeline_progress_is_monotonic():
    # The green-bar bug: the stage-start tick (idx/N) and a capability's own local 0..1 ticks fought
    # each other, so the fraction jumped to ~1.0 inside a stage then snapped back at the next stage.
    # Every reported fraction must now be monotonically non-decreasing (each stage's local progress
    # maps into its global band), starting at 0.0 and ending at 1.0.
    adata = demo.load_demo().adata
    seen: list[float] = []
    ctx = cap.RunContext(tissue="melanoma", use_llm=False, is_tumour=True)
    ctx.progress = lambda frac, label: seen.append(float(frac))
    pipeline.run_pipeline(adata, ctx, pipeline.PipelineOptions(export=False))
    assert seen, "no progress reported"
    assert seen[0] <= 0.02 and seen[-1] == 1.0                      # starts ~0, ends at done=1.0
    drops = [(i, round(seen[i - 1], 4), round(seen[i], 4))
             for i in range(1, len(seen)) if seen[i] < seen[i - 1] - 1e-9]
    assert not drops, f"progress went backwards at {drops[:6]}"


def test_run_pipeline_non_tumour_gate_and_completes_cleanly():
    adata = demo.load_demo().adata
    ctx = cap.RunContext(tissue="melanoma", use_llm=False, is_tumour=False)
    res = pipeline.run_pipeline(adata, ctx, pipeline.PipelineOptions(export=False, split=True))
    st = _status(res)
    # non-tumour -> the malignant stage is honestly gated out at the pipeline level
    assert st.get("malignant_concordance", "").startswith("skipped:non-tumour")
    assert st.get("cluster") == "ok"
    # optional arms (here split=True) never hard-fail the run - it completes clean
    assert res["summary"]["failed"] == [], res["stages"]
    # split was REQUESTED: it appears in the record (ran, or honestly skipped) - never silently absent.
    # (Whether it ran real SPLIT vs its in-env fallback is surfaced by the honesty-messaging track.)
    assert st.get("split_purify", "").split(":", 1)[0] in {"ok", "skipped"}
