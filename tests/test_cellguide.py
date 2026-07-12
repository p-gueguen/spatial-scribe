"""CellGuide marker provider - hermetic (no network).

Seeds a temp cache with real CellGuide payloads (tests/fixtures/cellguide/, trimmed from the live
snapshot) and pins CELLGUIDE_OFFLINE, so the suite never reaches the network. The real live endpoint
is exercised separately by an internal benchmark harness
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

FIX = Path(__file__).parent / "fixtures" / "cellguide"
TCELL = "CL_0000084"


@pytest.fixture
def cg(monkeypatch, tmp_path):
    """The cellguide module wired to a seeded offline cache under `<tmp>/offline/`."""
    from spatialscribe.analysis import cellguide

    snap = tmp_path / "offline"
    (snap / "canonical_marker_genes").mkdir(parents=True)
    (snap / "computational_marker_genes").mkdir(parents=True)
    shutil.copy(FIX / "celltype_metadata.sample.json", snap / "celltype_metadata.json")
    shutil.copy(FIX / f"canonical_{TCELL}.json", snap / "canonical_marker_genes" / f"{TCELL}.json")
    shutil.copy(FIX / f"computational_{TCELL}.json",
                snap / "computational_marker_genes" / f"{TCELL}.json")

    monkeypatch.setattr(cellguide, "_CACHE", tmp_path)
    monkeypatch.setenv("CELLGUIDE_OFFLINE", "1")
    cellguide._reset_for_tests()
    return cellguide


def test_resolve_cl_maps_names_synonyms_and_aliases(cg):
    assert cg.resolve_cl("T cell") == "CL:0000084"
    assert cg.resolve_cl("t  CELL") == "CL:0000084"                 # case/punctuation-insensitive
    assert cg.resolve_cl("Stromal/CAF") == "CL:0000057"            # composite alias -> fibroblast
    assert cg.resolve_cl("Malignant/Melanocyte") == "CL:0000148"   # -> melanocyte
    assert cg.resolve_cl("not a real cell type xyz") is None       # exact only, never fuzzy-snap


def test_canonical_is_preferred_and_human(cg):
    out = cg.markers_for_label("T cell", organism="Homo sapiens")
    assert out["source"] == "canonical"
    assert out["cl_id"] == "CL:0000084"
    assert "CD3D" in out["markers"]                                # literature-curated, human symbols
    assert len(out["markers"]) <= 12


def test_falls_back_to_computational_when_no_canonical(cg, tmp_path):
    # remove the seeded canonical file -> must fall back to the computational aggregate
    (tmp_path / "offline" / "canonical_marker_genes" / f"{TCELL}.json").unlink()
    cg._reset_for_tests()
    out = cg.markers_for_label("T cell", organism="Homo sapiens")
    assert out["source"] == "computational"
    # tissue-agnostic human aggregate, ranked by marker_score ("Effect Size") = real T-cell markers
    assert out["markers"][:3] == ["CD3E", "CD3D", "CD52"]
    assert "SKAP1" not in out["markers"]                           # a tissue-specific census non-marker


def test_mouse_uses_computational_with_mouse_case_symbols(cg):
    """Canonical is human-only (no organism, human symbols); mouse resolves to computational."""
    out = cg.markers_for_label("T cell", organism="mouse")         # 'mouse' folds to Mus musculus
    assert out["source"] == "computational"
    assert out["cl_id"] == "CL:0000084"
    assert "Cd3g" in out["markers"] or "Cd3d" in out["markers"]    # mouse-case
    assert "CD3D" not in out["markers"]                            # never the human-canonical symbols


def test_distilled_cache_is_written_small_and_capped(cg, tmp_path):
    (tmp_path / "offline" / "canonical_marker_genes" / f"{TCELL}.json").unlink()
    cg._reset_for_tests()
    cg.markers_for_label("T cell", organism="Homo sapiens", n=20)

    dist = tmp_path / "offline" / "computational_top" / f"{TCELL}.json"
    assert dist.exists(), "the distilled (not the raw) file must be cached"
    d = json.loads(dist.read_text())
    assert set(d) <= {"Homo sapiens", "Mus musculus"}
    assert all(len(v) <= cg._CACHE_TOP_N for v in d.values())     # capped at 20 to save space
    assert dist.stat().st_size < 4000                             # a few hundred bytes, not ~1.3 MB


def test_cache_survives_an_online_to_offline_transition(monkeypatch, tmp_path):
    """Files cached under a real snapshot id must still be served offline (via <cache>/.latest)."""
    from spatialscribe.analysis import cellguide

    snap_id = "1764612212"
    snap = tmp_path / snap_id
    (snap / "canonical_marker_genes").mkdir(parents=True)
    shutil.copy(FIX / "celltype_metadata.sample.json", snap / "celltype_metadata.json")
    shutil.copy(FIX / f"canonical_{TCELL}.json", snap / "canonical_marker_genes" / f"{TCELL}.json")
    (tmp_path / ".latest").write_text(snap_id)                  # what an online run leaves behind

    monkeypatch.setattr(cellguide, "_CACHE", tmp_path)
    monkeypatch.setenv("CELLGUIDE_OFFLINE", "1")                # now the network is gone
    cellguide._reset_for_tests()

    assert cellguide.snapshot() == snap_id                      # reused the last-known id, not None
    out = cellguide.markers_for_label("T cell")
    assert out["source"] == "canonical" and "CD3D" in out["markers"]


def test_offline_with_empty_cache_degrades_to_none(monkeypatch, tmp_path):
    """No network, no cache -> the caller must be able to keep its curated markers."""
    from spatialscribe.analysis import cellguide

    monkeypatch.setattr(cellguide, "_CACHE", tmp_path / "empty")
    monkeypatch.setenv("CELLGUIDE_OFFLINE", "1")
    cellguide._reset_for_tests()
    out = cellguide.markers_for_label("T cell")
    assert out == {"cl_id": None, "source": "none", "markers": []}
    assert cellguide.resolve_cl("T cell") is None                 # no metadata -> no crash


def test_markers_for_labels_batches_and_dedupes(cg):
    out = cg.markers_for_labels(["T cell", "T cell", "not real xyz"])
    assert set(out) == {"T cell", "not real xyz"}                  # deduped
    assert out["T cell"]["source"] == "canonical"
    assert out["not real xyz"] == {"cl_id": None, "source": "none", "markers": []}


def test_markers_for_types_prefers_cellguide_over_the_llm(cg, monkeypatch):
    """markers.markers_for_types grounds a type via CellGuide without ever calling the LLM."""
    from spatialscribe.analysis import markers as mk

    called = {"llm": False}

    def _boom(*a, **k):
        called["llm"] = True
        return {}
    monkeypatch.setattr(mk, "_generate_llm_panel", _boom)

    panel = ["CD3D", "CD3E", "MS4A1", "FOXP3", "GZMB", "CCL5"]   # some canonical T-cell markers on-panel
    out = mk.markers_for_types(["T cell"], panel_genes=panel, tissue="blood", use_llm=True)
    assert "T cell" in out and "CD3D" in out["T cell"]
    assert called["llm"] is False, "CellGuide grounded it, so the LLM must not be consulted"


def test_canonical_is_ranked_by_effect_size_not_alphabetical(cg):
    """Flagship high-effect markers must lead; the alphabetical head must not bury them (review fix)."""
    out = cg.markers_for_label("T cell", organism="Homo sapiens", n=6)
    assert out["source"] == "canonical"
    # CD3D (a top-effect T-cell marker in the computational aggregate) beats ANXA2/AQP3, which head
    # the canonical file alphabetically but are not effect-size leaders.
    assert "CD3D" in out["markers"]
    assert out["markers"].index("CD3D") < 3


def test_thin_canonical_falls_back_to_computational(cg, tmp_path, monkeypatch):
    """< _MIN_CANONICAL canonical symbols -> use the richer computational set (review fix)."""
    canon = tmp_path / "offline" / "canonical_marker_genes" / f"{TCELL}.json"
    canon.write_text(json.dumps([{"tissue": "blood", "symbol": "CD3D"}]))   # only 1 symbol
    cg._reset_for_tests()
    out = cg.markers_for_label("T cell", organism="Homo sapiens")
    assert out["source"] == "computational"
    assert out["markers"][:2] == ["CD3E", "CD3D"]


def test_non_dict_record_never_raises(cg, tmp_path):
    """A malformed (non-dict) computational record must not break the 'never raises' contract."""
    comp = tmp_path / "offline" / "computational_marker_genes" / f"{TCELL}.json"
    recs = json.loads(comp.read_text())
    comp.write_text(json.dumps(["oops-not-a-dict", *recs]))
    (tmp_path / "offline" / "canonical_marker_genes" / f"{TCELL}.json").unlink()
    cg._reset_for_tests()
    out = cg.markers_for_label("T cell", organism="Homo sapiens")   # must not raise
    assert out["source"] == "computational" and out["markers"]


def test_resolve_cl_names_win_over_earlier_synonyms(cg, tmp_path):
    """A later type's canonical NAME must not be shadowed by an earlier type's synonym (review fix)."""
    meta = tmp_path / "offline" / "celltype_metadata.json"
    m = json.loads(meta.read_text())
    # earlier record carries a synonym that normalises to a LATER record's canonical name
    m["CL:0000001"] = {"name": "widget cell", "id": "CL:0000001", "synonyms": ["ghost type"]}
    m["CL:0000002"] = {"name": "ghost type", "id": "CL:0000002", "synonyms": []}
    meta.write_text(json.dumps(m))
    cg._reset_for_tests()
    assert cg.resolve_cl("ghost type") == "CL:0000002"          # the NAME wins, not the earlier synonym


def test_markers_for_types_grounds_with_no_llm(cg, monkeypatch):
    """The whole point: panel-check grounds novel types via CellGuide even with the LLM OFF."""
    from spatialscribe.analysis import markers as mk

    monkeypatch.setattr(mk, "_generate_llm_panel", lambda *a, **k: {})   # LLM unavailable
    out = mk.markers_for_types(["T cell"], panel_genes=["CD3D", "CD3E", "MS4A1"],
                               tissue="blood", use_llm=False)
    assert out.get("T cell") and "CD3D" in out["T cell"]


_MINI_OBO = """format-version: 1.2
[Term]
id: CL:0000988
name: hematopoietic cell
[Term]
id: CL:0000738
name: leukocyte
is_a: CL:0000988 ! hematopoietic cell
[Term]
id: CL:0000542
name: lymphocyte
is_a: CL:0000738 ! leukocyte
[Term]
id: CL:0000084
name: T cell
is_a: CL:0000542 ! lymphocyte {is_inferred="true"}
[Term]
id: CL:0000624
name: CD4 T cell
is_a: CL:0000084 ! T cell
[Term]
id: CL:0000625
name: CD8 T cell
is_a: CL:0000084 ! T cell
[Term]
id: CL:0000097
name: mast cell
is_a: CL:0000738 ! leukocyte
[Term]
id: CL:0000499
name: stromal cell
[Term]
id: CL:0000057
name: fibroblast
is_a: CL:0000499 ! stromal cell
[Term]
id: CL:0000066
name: epithelial cell
"""


def test_cell_ontology_mergeability(cg, tmp_path):
    """The Cell Ontology decides which types may be MERGED: same major lineage -> yes, different -> no.
    CD4/CD8 T (both leukocytes) merge; a mast cell (immune) never merges into epithelium or stroma."""
    meta = tmp_path / "offline" / "celltype_metadata.json"
    m = json.loads(meta.read_text())
    m.update({
        "CL:0000624": {"name": "CD4 T cell", "id": "CL:0000624", "synonyms": []},
        "CL:0000625": {"name": "CD8 T cell", "id": "CL:0000625", "synonyms": []},
        "CL:0000097": {"name": "mast cell", "id": "CL:0000097", "synonyms": []},
        "CL:0000066": {"name": "epithelial cell", "id": "CL:0000066", "synonyms": []},
        "CL:0000057": {"name": "fibroblast", "id": "CL:0000057", "synonyms": []},
    })
    meta.write_text(json.dumps(m))
    (tmp_path / "cl-basic.obo").write_text(_MINI_OBO)   # default _CL_OBO = _CACHE/'cl-basic.obo'
    cg._reset_for_tests()

    # is_a closure ignores the {is_inferred=...} qualifier; anchors resolve to compartments
    assert cg.cl_compartments("CL:0000625") == {"immune", "hematopoietic"}   # CD8 T -> leukocyte
    assert cg.cl_compartments("CL:0000066") == {"epithelial"}
    assert cg.cl_compartments("CL:0000057") == {"stromal"}

    assert cg.labels_mergeable("CD4 T cell", "CD8 T cell")["mergeable"] is True     # same lineage
    assert cg.labels_mergeable("mast cell", "epithelial cell")["mergeable"] is False  # immune vs epithelial
    assert cg.labels_mergeable("mast cell", "fibroblast")["mergeable"] is False       # immune vs stromal
    # an unresolved label -> None (undecidable), so the caller falls back to its heuristic
    assert cg.labels_mergeable("mast cell", "not a real type xyz")["mergeable"] is None


def test_keyword_compartment_fallback_for_composite_labels(cg):
    # Study-specific composite labels (RCTD / CosMx custom vocab) that do NOT resolve to a CL id still
    # classify by lineage keyword, so the merge guard works beyond canonical names. No obo needed.
    assert cg._label_compartments("MFAP5_IGFBP6_fibroblast") == {"stromal"}
    assert cg._label_compartments("LYVE1_macrophage") == {"immune"}
    assert cg._label_compartments("Capillary_EC") == {"endothelial"}
    assert cg._label_compartments("PIP_mammary_luminal_cell") == {"epithelial"}
    assert cg._label_compartments("Vascular_smooth_muscle") == {"muscle"}
    # cross-lineage composite pair is blocked; same-lineage composite pair merges
    assert cg.labels_mergeable("T_cell", "PIP_mammary_luminal_cell")["mergeable"] is False
    assert cg.labels_mergeable("MFAP5_IGFBP6_fibroblast", "GPC3_fibroblast")["mergeable"] is True


def test_cell_ontology_degrades_without_obo(cg):
    """No obo cached -> the is-a graph is empty, but the keyword fallback still decides for named
    lineages; only a label with NO lineage keyword is undecidable (None). Never crashes."""
    cg._reset_for_tests()
    assert cg._cl_graph() == {}                                    # no obo -> empty is-a graph
    # keyword fallback still classifies canonical lineages (both immune) -> decidable, not a crash
    assert cg.labels_mergeable("T cell", "mast cell")["mergeable"] is True
    # a label with no CL id and no lineage keyword is undecidable -> caller falls back to its heuristic
    assert cg.labels_mergeable("widget", "gadget")["mergeable"] is None


def test_markers_for_types_infers_mouse_from_symbol_case(cg):
    from spatialscribe.analysis import markers as mk

    assert mk.infer_organism(["Cd3d", "Cd3e", "Ms4a1", "Gzmb"]) == "mouse"
    assert mk.infer_organism(["CD3D", "CD3E", "MS4A1"]) == "human"
    out = mk.markers_for_types(["T cell"], panel_genes=["Cd3d", "Cd3e", "Cd3g"], use_llm=False)
    # mouse -> computational, mouse-case symbols (never the human canonical set)
    assert out.get("T cell")
    assert any(g in out["T cell"] for g in ("Cd3g", "Cd3d"))
    assert "CD3D" not in out["T cell"]
