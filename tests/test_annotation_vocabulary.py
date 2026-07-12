"""The LLM names clusters; it does not get to invent labels or overrule marker evidence.

Two defects this pins, both measured on the 2026-07 overnight benchmark (5 real sections, the two
arms differing ONLY in which LLM named the same leiden clusters):

1. OPEN VOCABULARY. ``llm._ANNOTATE_SYS`` asked for ``{"label": str}`` with no allowed set, and
   ``consensus_annotate`` wrote that free text straight into ``obs['cell_type']``. 44 of the 59
   cluster-label disagreements between the two arms were NAMING-ONLY (same lineage, different
   string). Constrain the model at GENERATION time and validate on the way back; never string-snap
   an off-vocabulary answer onto a curated key ("pDC"/"RBC" are near nothing and would snap wrong).

2. SILENT OVERRIDE. ``final = claude if use_llm and claude else marker`` meant the LLM won every
   genuine disagreement, and ``agreement`` was a loose lowercase substring test recorded for display
   that never gated anything. The remaining 15 disagreements were real biology, and there is NO
   evidence the LLM is more accurate than the marker argmax. A true lineage conflict must be flagged,
   not silently resolved in the LLM's favour.
"""
from __future__ import annotations

import numpy as np
import pytest

anndata = pytest.importorskip("anndata")


NOVEL = "Novel / unknown"


@pytest.fixture
def clustered(raw_adata, ctx):
    """raw_adata through compute_qc + cluster (no annotation yet)."""
    from spatialscribe.analysis import capabilities as cap

    for name, params in [("compute_qc", {}), ("cluster", {"resolution": 1.0})]:
        res = cap.run(raw_adata, name, params, ctx)
        assert res.ok, (name, res.error)
    return raw_adata


@pytest.fixture
def marker_sets(clustered):
    from spatialscribe.analysis import markers as mk

    present = mk.present(list(clustered.var_names), mk.for_tissue("melanoma"))
    sets = {k: v for k, v in present.items() if v}
    assert len(sets) >= 2, "fixture needs >=2 lineages on the panel"
    return sets


def _fake_llm(monkeypatch, label_for):
    """Patch llm.annotate_clusters; return the captured kwargs dict."""
    from spatialscribe.analysis import llm

    captured: dict = {}

    def _fake(markers_per_cluster, context="", **kwargs):
        captured.update(kwargs)
        captured["clusters"] = list(markers_per_cluster)
        return {cl: {"label": label_for(cl), "confidence": "high", "rationale": "r"}
                for cl in markers_per_cluster}

    monkeypatch.setattr(llm, "annotate_clusters", _fake)
    return captured


def test_allowed_labels_and_the_escape_hatch_are_sent_to_the_model(clustered, marker_sets, monkeypatch):
    """Closing the vocabulary happens at GENERATION time, not by post-hoc snapping."""
    from spatialscribe.analysis import annotate

    captured = _fake_llm(monkeypatch, lambda cl: next(iter(marker_sets)))
    annotate.consensus_annotate(clustered, cluster_key="leiden", use_llm=True,
                                marker_sets=marker_sets)

    allowed = captured.get("allowed_labels")
    assert allowed is not None, "the model must be told the allowed label set"
    assert set(marker_sets).issubset(set(allowed))
    assert NOVEL in allowed, "the 'none of these fit' escape hatch must stay reachable"


def test_off_vocabulary_label_falls_back_to_marker_and_is_never_snapped(clustered, marker_sets, monkeypatch):
    """An invented label ('pDC') must not enter cell_type, and must not be snapped to a curated key."""
    from spatialscribe.analysis import annotate

    _fake_llm(monkeypatch, lambda cl: "pDC")
    df = annotate.consensus_annotate(clustered, cluster_key="leiden", use_llm=True,
                                     marker_sets=marker_sets)

    cats = set(map(str, clustered.obs["cell_type"].astype("category").cat.categories))
    assert "pDC" not in cats
    assert cats.issubset(set(marker_sets) | {NOVEL, "Unassigned", "Unknown"})
    assert (df["final"] == df["marker_label"]).all(), "off-vocabulary must fall back to marker evidence"
    assert df["off_vocabulary"].all()


def test_novel_unknown_escape_hatch_is_preserved(clustered, marker_sets, monkeypatch):
    """'None of these fit' is an honest answer and must survive validation."""
    from spatialscribe.analysis import annotate

    _fake_llm(monkeypatch, lambda cl: NOVEL)
    df = annotate.consensus_annotate(clustered, cluster_key="leiden", use_llm=True,
                                     marker_sets=marker_sets)

    assert (df["final"] == NOVEL).all()
    assert not df["off_vocabulary"].any()


def test_agreeing_label_is_used_and_flags_no_conflict(clustered, marker_sets, monkeypatch):
    """When the LLM confirms the marker argmax, that label is the answer and nothing is flagged."""
    from spatialscribe.analysis import annotate

    m_lab = annotate.marker_labels(clustered, "leiden", marker_sets=marker_sets)
    _fake_llm(monkeypatch, lambda cl: m_lab[str(cl)])
    df = annotate.consensus_annotate(clustered, cluster_key="leiden", use_llm=True,
                                     marker_sets=marker_sets)

    assert (df["final"] == df["marker_label"]).all()
    assert df["agreement"].all()
    assert not df["conflict"].any()
    assert not clustered.obs["label_conflict"].to_numpy().any()


def test_llm_cannot_silently_overrule_marker_evidence(clustered, marker_sets, monkeypatch):
    """A genuine lineage conflict is FLAGGED, not resolved in the LLM's favour."""
    from spatialscribe.analysis import annotate

    m_lab = annotate.marker_labels(clustered, "leiden", marker_sets=marker_sets)
    keys = list(marker_sets)

    def _other(cl):                       # a DIFFERENT, valid, in-vocabulary lineage
        mine = m_lab[str(cl)]
        return next(k for k in keys if k != mine)

    _fake_llm(monkeypatch, _other)
    df = annotate.consensus_annotate(clustered, cluster_key="leiden", use_llm=True,
                                     marker_sets=marker_sets)

    assert (df["final"] == df["marker_label"]).all(), "the LLM must not win a real disagreement"
    assert (df["claude_label"] != df["marker_label"]).all()
    assert not df["agreement"].any()
    assert df["conflict"].all()
    assert clustered.obs["label_conflict"].to_numpy().all()


def _confident_adata(n=40, conflict=False):
    """A section engineered to yield PASS cells, so the conflict cap is testable.

    Each cell expresses ONLY its own lineage's markers (the other lineage stays 0), otherwise
    ``purity.crisp_purity`` marks every cell impure -> ``contam = 0.5`` -> confidence caps at exactly
    0.5, one ulp below the PASS threshold, and the control arm has no PASS cells to lose.
    Fillers keep >= 5 genes detected so the independent low-signal gate does not fire.
    """
    import pandas as pd

    rng = np.random.default_rng(0)
    genes = ["CD3D", "CD3E", "MS4A1", "CD79A", "F1", "F2", "F3", "F4"]
    X = np.zeros((n, len(genes)), dtype="float32")
    X[:, 4:] = rng.poisson(2, size=(n, 4)) + 1
    lab = np.where(np.arange(n) < n // 2, "T", "B")
    for i in range(n):
        own = (0, 1) if lab[i] == "T" else (2, 3)
        X[i, list(own)] = 50
    a = anndata.AnnData(X=X)
    a.var_names = genes
    a.layers["counts"] = a.X.copy()
    a.obs["cell_type"] = pd.Categorical(lab)
    a.obs["pct_counts_control"] = 0.0
    a.obs["label_conflict"] = bool(conflict)
    return a


def test_conflicted_cluster_never_reads_as_a_confident_call():
    """A flagged conflict caps the per-cell verdict below PASS (flag, don't silently pick).

    Non-vacuous by construction: the unconflicted control arm must produce PASS cells, else
    "no PASS under conflict" would hold trivially and prove nothing.
    """
    from spatialscribe.analysis import annotate

    markers = {"T": ["CD3D", "CD3E"], "B": ["MS4A1", "CD79A"]}

    control = _confident_adata(conflict=False)
    res = annotate.apply_confidence(control, marker_sets=markers)
    assert res["pct_pass"] > 0, "control arm must contain PASS cells or this test is vacuous"
    assert res["pct_label_conflict"] == 0.0

    flagged = _confident_adata(conflict=True)
    res2 = annotate.apply_confidence(flagged, marker_sets=markers)
    verdict = flagged.obs["annotation_verdict"].astype(str).to_numpy()
    assert not (verdict == "PASS").any(), "a conflicted cluster must not contain PASS cells"
    assert res2["pct_label_conflict"] == 1.0
    # Capped, not abstained: the marker evidence still stands, it is merely uncorroborated.
    assert (verdict == "WARN").any()


def test_vocabulary_is_every_scored_lineage_not_just_the_argmax_winners(clustered, monkeypatch):
    """A lineage no cluster happened to win must still be nameable, or the LLM cannot fix a miss.

    Subset to a SINGLE leiden cluster so exactly one lineage wins the argmax while many are scored:
    with the old ``vocab = argmax winners`` the model would have been offered a 1-label vocabulary.
    """
    from spatialscribe.analysis import annotate

    first = str(clustered.obs["leiden"].astype(str).iloc[0])
    one = clustered[clustered.obs["leiden"].astype(str) == first].copy()

    captured = _fake_llm(monkeypatch, lambda cl: NOVEL)
    annotate.consensus_annotate(one, cluster_key="leiden", use_llm=True, marker_sets=None)

    scored = set(annotate.score_marker_sets(one, None))
    winners = set(annotate.marker_labels(one, "leiden").values())
    allowed = set(captured["allowed_labels"])
    assert len(winners) == 1, "single-cluster subset must have exactly one argmax winner"
    assert len(scored) > 1, "many lineages must be scored for this test to discriminate"
    assert scored.issubset(allowed), "every scored lineage must be offered to the model"


def test_skin_vocabulary_has_a_bucket_for_mast_cells():
    """Closing the vocabulary makes a MISSING lineage a real cost: mast cells exist in skin.

    Before, dermal mast cells could only be named 'Novel / unknown' (or misfiled into Myeloid) because
    the skin/melanoma set had no Mast key, while the sibling EPITHELIAL_LINEAGES did.
    """
    from spatialscribe.analysis import markers as mk

    skin = mk.resolve_markers("melanoma")
    assert "Mast" in skin
    assert {"TPSAB1", "CPA3"}.issubset(set(skin["Mast"]))
    # the immune lineages that share markers with Mast must stay distinct
    assert "Myeloid" in skin and "TPSAB1" not in skin["Myeloid"]


def test_a_lineage_absent_from_the_panel_never_enters_the_vocabulary(clustered):
    """present() filters marker sets to the panel, so an off-panel lineage cannot be offered."""
    from spatialscribe.analysis import markers as mk

    panel = [g for g in clustered.var_names if g not in ("TPSAB1", "TPSB2", "CPA3", "MS4A2")]
    on_panel = mk.present(panel, mk.resolve_markers("melanoma"))
    assert not on_panel.get("Mast"), "a lineage with no on-panel markers must not be nameable"


def test_apply_confidence_tolerates_a_missing_conflict_flag():
    """Objects from an older build, or from reference_transfer, carry no label_conflict column."""
    import pandas as pd

    from spatialscribe.analysis import annotate

    a = _confident_adata(conflict=False)
    del a.obs["label_conflict"]
    res = annotate.apply_confidence(a, marker_sets={"T": ["CD3D", "CD3E"], "B": ["MS4A1", "CD79A"]})
    assert res["pct_label_conflict"] == 0.0
    assert res["pct_pass"] > 0                       # nothing is capped when nothing is flagged
    assert pd is not None


def test_conflicted_cell_that_also_fails_reports_the_conflict_reason():
    """When a conflicted cell fails on its own merits, 'label_conflict' is the actionable reason."""
    from spatialscribe.analysis import annotate

    markers = {"T": ["CD3D", "CD3E"], "B": ["MS4A1", "CD79A"]}
    a = _confident_adata(conflict=True)
    a.obs["pct_counts_control"] = 12.0            # contaminated -> confidence 0 -> FAIL
    annotate.apply_confidence(a, marker_sets=markers)

    reasons = set(a.obs["annotation_reason"].astype(str))
    assert "label_conflict" in reasons
    finals = set(a.obs["cell_type_final"].astype(str))
    assert annotate._ABSTAIN["label_conflict"] in finals


def test_no_llm_keeps_the_marker_label_and_flags_nothing(clustered, marker_sets):
    """use_llm=False is unchanged: marker argmax wins, no conflict, no off-vocabulary."""
    from spatialscribe.analysis import annotate

    df = annotate.consensus_annotate(clustered, cluster_key="leiden", use_llm=False,
                                     marker_sets=marker_sets)
    assert (df["final"] == df["marker_label"]).all()
    assert not df["conflict"].any()
    assert not df["off_vocabulary"].any()


def test_label_validation_is_case_and_whitespace_tolerant(clustered, marker_sets, monkeypatch):
    """'  t cell ' is the curated 'T cell', not an off-vocabulary invention - canonicalise, don't reject."""
    from spatialscribe.analysis import annotate

    m_lab = annotate.marker_labels(clustered, "leiden", marker_sets=marker_sets)
    _fake_llm(monkeypatch, lambda cl: f"  {m_lab[str(cl)].lower()} ")
    df = annotate.consensus_annotate(clustered, cluster_key="leiden", use_llm=True,
                                     marker_sets=marker_sets)

    assert not df["off_vocabulary"].any()
    assert (df["final"] == df["marker_label"]).all()      # canonical casing restored
    assert not df["conflict"].any()
