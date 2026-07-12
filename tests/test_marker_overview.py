"""Grouped marker dot-plot overview (CONTRACT section 4).

Covers: the offline fallback grouping is deterministic and total (every input type lands in exactly
one category); ``llm.group_cell_types`` is grounded (drops types the model invented, keeps every input
once, unassigned -> "Other"); ``category_overview`` ranks columns by expression (not the LLM); and
``dotplot_view`` honours an explicit ``group_order``.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")


def _section(types, per=8, seed=0, genes=("g0", "g1", "g2", "g3", "g4", "g5")):
    rng = np.random.default_rng(seed)
    n = len(types) * per
    a = anndata.AnnData(X=rng.poisson(3, size=(n, len(genes))).astype("float32"))
    a.obs_names = [f"c{i}" for i in range(n)]
    a.var_names = list(genes)
    a.obsm["spatial"] = rng.normal(size=(n, 2)) * 100
    a.obs["cell_type"] = pd.Categorical([t for t in types for _ in range(per)])
    return a


def test_fallback_grouping_total_and_deterministic():
    from spatialscribe.analysis import views

    types = ["CD8 T", "CD4 T", "Treg", "B cell", "Plasma", "Macrophage",
             "Tumor epithelial", "Endothelial", "CAF", "Weirdtype 7"]
    g1 = views._fallback_group(types)
    assert g1 == views._fallback_group(types)                    # deterministic

    placed = [t for members in g1.values() for t in members]
    assert sorted(placed) == sorted(types)                       # total: every input assigned
    assert len(placed) == len(set(placed))                       # exactly once each
    assert "Weirdtype 7" in g1.get("Other", [])                  # unmatched -> Other
    assert set(g1["T/NK lymphoid"]) == {"CD8 T", "CD4 T", "Treg"}  # related types adjacent


def test_group_cell_types_drops_invented_and_is_total(monkeypatch):
    from spatialscribe.analysis import llm

    monkeypatch.setattr(llm, "available", lambda: True)
    # the model invents "Ghost cell" and forgets "Treg" - grounding must drop the stray and
    # sweep the unassigned real type into "Other".
    monkeypatch.setattr(
        llm, "complete",
        lambda *a, **k: '{"Lymphoid": ["CD8 T", "Ghost cell"], "Myeloid": ["Macrophage"]}')
    llm._GROUP_CACHE.clear()

    out = llm.group_cell_types(["CD8 T", "Treg", "Macrophage"], tissue="breast")
    placed = [t for m in out.values() for t in m]
    assert "Ghost cell" not in placed                            # invented type dropped
    assert sorted(placed) == ["CD8 T", "Macrophage", "Treg"]     # total over the inputs
    assert placed.count("CD8 T") == 1                            # each input exactly once
    assert "Treg" in out.get("Other", [])                        # unassigned -> Other


def test_group_cell_types_empty_without_llm(monkeypatch):
    from spatialscribe.analysis import llm

    monkeypatch.setattr(llm, "available", lambda: False)
    assert llm.group_cell_types(["A", "B"], tissue="x") == {}    # callers fall back


def test_category_overview_ranks_columns_by_expression():
    """The per-category top-2 marker columns come from a computed one-vs-rest expression
    difference, never the LLM: the gene that is high in a category's cells is what gets picked."""
    from spatialscribe.analysis import capabilities as cap, views

    rng = np.random.default_rng(0)
    genes = ["CD3D", "PECAM1", "NOISE1", "NOISE2"]
    n = 40
    X = rng.poisson(1, (n, 4)).astype("float32")
    X[:20, 0] += 25                       # CD3D high in the T cells
    X[20:, 1] += 25                       # PECAM1 high in the Endothelial cells
    a = anndata.AnnData(X=X)
    a.var_names = genes
    a.obs_names = [f"c{i}" for i in range(n)]
    a.obsm["spatial"] = rng.normal(size=(n, 2))
    a.obs["cell_type"] = pd.Categorical(["T cell"] * 20 + ["Endothelial"] * 20)

    ctx = cap.RunContext(tissue="melanoma", use_llm=False)        # curated LINEAGE_MARKERS, no LLM
    genes_out, group_order, note = views.category_overview(a, ctx)
    assert "CD3D" in genes_out and "PECAM1" in genes_out
    assert set(group_order) == {"T cell", "Endothelial"}
    # Offline (use_llm=False): grouping is the deterministic keyword fallback, recorded structurally
    # in uns (the prose "rows grouped by ..." caption was dropped as UI noise, so `note` stays quiet).
    assert a.uns.get("marker_overview_source") == "keyword fallback"   # offline honesty, structural
    assert "LLM" not in note                                          # never falsely claims an LLM grouped


def test_dotplot_view_honours_group_order():
    from spatialscribe.analysis import views

    a = _section(["A", "B", "C"], per=8)
    fig = views.dotplot_view(a, ["g0", "g1"], group_order=["C", "A", "ZZ_absent"])
    yt = list(fig.layout.yaxis.ticktext)
    assert yt[0] == "C" and yt[1] == "A"                          # requested order first
    assert "ZZ_absent" not in yt                                  # absent-from-section skipped
    assert set(yt) == {"A", "B", "C"} and yt[2] == "B"            # present-but-unlisted appended


def test_dotplot_colorbar_labels_relative_not_absolute():
    """The colour is mean expression scaled 0-1 PER GENE (relative across types), so the colorbar must
    NOT say 'mean expr' - which reads as absolute and made a spillover-driven low gene (T-cell MS4A1)
    look like a strong marker just because it was the per-gene argmax."""
    from spatialscribe.analysis import views

    a = _section(["A", "B", "C"], per=8)
    fig = views.dotplot_view(a, ["g0", "g1"])
    label = fig.data[0].marker.colorbar.title.text
    assert "relative" in label.lower() and "mean expr" not in label.lower()


def test_marker_dotplot_auto_shows_split_pair_when_layer_present():
    """Once a split_corrected layer exists, ANY dot-plot (incl. the auto overview, called with corrected
    unset) draws the raw + purified PAIR - so 'after SPLIT' stops silently rendering raw X."""
    from spatialscribe.analysis import capabilities as cap

    a = _section(["T cell", "B/Plasma", "Myeloid"], per=8)
    ctx = cap.RunContext(tissue="melanoma", use_llm=False)

    # no split layer -> a single figure
    r1 = cap.run(a, "marker_dotplot", {"genes": ["g0", "g1"]}, ctx)
    figs1 = [x for x in ctx.artifacts if x.get("kind") == "figure"]
    assert r1.ok and len(figs1) == 1 and r1.value.get("split_available") is False

    # add a purified layer -> the SAME call (corrected unset) now emits raw + purified
    a.layers["split_corrected"] = a.X.copy()
    ctx.artifacts.clear()
    r2 = cap.run(a, "marker_dotplot", {"genes": ["g0", "g1"]}, ctx)
    figs2 = [x for x in ctx.artifacts if x.get("kind") == "figure"]
    assert r2.ok and r2.value.get("compared") is True and len(figs2) == 2
    titles = {x["title"] for x in figs2}
    assert titles == {"marker dot-plot (raw)", "marker dot-plot (SPLIT-purified)"}


def test_fallback_group_never_files_an_abstention_label_as_a_lineage():
    """Regression: the keyword "nk" is a bare substring of "u-nk-nown", so an untyped cell used to be
    grouped under T/NK lymphoid in the offline dot-plot. Abstention labels name the ABSENCE of a cell
    type, so they belong in "Other" - and no keyword may match inside a word."""
    from spatialscribe.analysis import views

    types = ["NK cell", "CD8 T cell", "Unknown", "Unassigned: low quality", "Ambiguous: mixed",
             "Uncertain", "Unresolvable", "Novel state 3", "Endothelial", "conventional DC"]
    g = views._fallback_group(types)

    assert set(g.get("Other", [])) == {"Unknown", "Unassigned: low quality", "Ambiguous: mixed",
                                       "Uncertain", "Unresolvable", "Novel state 3"}
    assert g["T/NK lymphoid"] == ["NK cell", "CD8 T cell"]        # real lymphoid types still match
    assert g["Vascular"] == ["Endothelial"]                       # "epithel" must not match "endothelial"
    assert g["Myeloid"] == ["conventional DC"]                    # space-prefixed keyword " dc" still works
    flat = [t for v in g.values() for t in v]
    assert sorted(flat) == sorted(types)                          # still total: nothing dropped
