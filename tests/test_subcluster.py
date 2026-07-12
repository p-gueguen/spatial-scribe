"""Subcluster: LLM naming + confidence populate (not all-'n/a'), and a lower default resolution.

Regression for three linked bugs: (1) too many subclusters (resolution 0.5 default), (2) the CONF
column stuck at 'n/a', and (3) no LLM subtype names - (2) and (3) shared a root cause: name_subtypes
used a FIXED token budget that truncated the JSON for many subclusters, so every row fell back to
"subcluster N" / "n/a". The LLM boundary is monkeypatched (no network).
"""

from __future__ import annotations


def test_name_subtypes_scales_token_budget(monkeypatch):
    from spatialscribe.analysis import llm

    captured = {}

    def fake_complete(system, user, max_tokens=1024, model=None, json_mode=False):
        captured["max_tokens"] = max_tokens
        return "{" + ",".join(
            f'"{i}": {{"label": "CD8 T", "confidence": "high", "rationale": "x"}}' for i in range(8)
        ) + "}"

    monkeypatch.setattr(llm, "complete", fake_complete)
    out = llm.name_subtypes("T cell", {str(i): ["CD8A", "GZMB"] for i in range(8)})
    # scaled with the number of subclusters (200*8+400=2000), not the old fixed 1200 that truncated.
    assert captured["max_tokens"] >= 200 * 8 + 400
    assert out["0"]["confidence"] == "high" and out["7"]["label"] == "CD8 T"


def test_name_subtypes_small_case_keeps_a_floor(monkeypatch):
    from spatialscribe.analysis import llm

    captured = {}

    def fake_complete(system, user, max_tokens=1024, model=None, json_mode=False):
        captured["mt"] = max_tokens
        return "{}"

    monkeypatch.setattr(llm, "complete", fake_complete)
    llm.name_subtypes("T cell", {"0": ["CD8A"]})
    assert captured["mt"] >= 1200            # floor preserved for the tiny case


def test_subcluster_populates_llm_labels_and_confidence(processed_adata, monkeypatch):
    from spatialscribe.analysis import llm, subcluster as sc

    ct = str(processed_adata.obs["cell_type"].value_counts().index[0])   # most common type (>=20 cells)

    def fake_name(parent, top, context="x"):
        return {g: {"label": f"{parent} / CD8 T", "confidence": "high", "rationale": "r"} for g in top}

    monkeypatch.setattr(llm, "name_subtypes", fake_name)
    _, rows = sc.subcluster(processed_adata, ct, use_llm=True)
    assert rows
    assert all(r["confidence"] == "high" for r in rows)          # not 'n/a'
    assert all("CD8 T" in r["label"] for r in rows)              # LLM name, not "subcluster N"


def test_subcluster_default_resolution_is_low():
    import inspect

    from spatialscribe.analysis import subcluster as sc

    # "Low", not an exact constant: the default was tuned 0.5 -> 0.3 -> 0.1 (commit d2a7da8) to avoid
    # over-splitting subtypes. Assert the property (a coarse first pass), not a brittle literal that
    # re-breaks on every retune - which is exactly what left this test asserting a stale 0.3.
    assert inspect.signature(sc.subcluster).parameters["resolution"].default <= 0.3


def test_subcluster_capability_exposes_resolution():
    from spatialscribe.analysis import capabilities as cap

    assert "resolution" in cap.REGISTRY["subcluster"].params
