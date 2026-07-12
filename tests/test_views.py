"""Copilot-driven spatial views: select_cells criteria + the highlight_cells capability."""
from __future__ import annotations

import numpy as np
import pytest

anndata = pytest.importorskip("anndata")


def _section(n=200, seed=0):
    rng = np.random.default_rng(seed)
    a = anndata.AnnData(X=rng.poisson(3, size=(n, 6)).astype("float32"))
    a.obs_names = [f"c{i}" for i in range(n)]
    a.var_names = [f"g{i}" for i in range(6)]
    a.obsm["spatial"] = rng.normal(size=(n, 2)) * 100
    import pandas as pd
    a.obs["cell_type"] = pd.Categorical(rng.choice(["T cell", "Myeloid", "Tumor"], n))
    return a


def test_dotplot_bubbles_fit_within_cell_pitch():
    """A 100%-expressing bubble must fit inside its own grid cell so it never bleeds onto the
    neighbouring category. Marker `size` is in pixels, so the max dot is derived from the row/column
    pitch - regression for a fixed 26 px max that overlapped rows on a tall (22-type) dot-plot."""
    import pandas as pd
    from spatialscribe.analysis import views

    rng = np.random.default_rng(0)
    n = 22 * 20
    a = anndata.AnnData(X=rng.poisson(2.0, size=(n, 40)).astype("float32"))
    a.var_names = [f"g{i}" for i in range(40)]
    a.obs_names = [f"c{i}" for i in range(n)]
    a.obs["cell_type"] = pd.Categorical([f"T{i % 22}" for i in range(n)])
    genes = [f"g{i}" for i in range(36)]
    fig = views.dotplot_view(a, genes, height=460)
    sizes = list(fig.data[0].marker.size)
    row_pitch = (460 - 40.0) / 22
    assert max(sizes) < row_pitch, (max(sizes), row_pitch)     # no vertical overlap onto neighbour rows
    assert min(sizes) > 0


def test_select_cells_celltype_and_lowquality_fallback():
    from spatialscribe.analysis import views

    a = _section()
    # cell-type match (case-insensitive)
    m, resolved, note = views.select_cells(a, "t cell")
    assert resolved == "T cell" and m.sum() == (a.obs["cell_type"].astype(str) == "T cell").sum()

    # low-quality with a rich-count section (nothing under the floor) -> honest lowest-decile fallback
    a.obs["total_counts"] = np.full(a.n_obs, 500.0)   # all well above any floor
    m2, resolved2, note2 = views.select_cells(a, "low quality")
    assert "lowest-count decile" in resolved2.lower()
    assert "no cells were flagged" in note2.lower()      # states the fallback honestly
    assert 0 < m2.sum() <= a.n_obs                        # highlights *something* (the tail), not everything

    # unresolvable -> empty mask + explanatory note, never raises
    m3, _, note3 = views.select_cells(a, "qwerty")
    assert m3.sum() == 0 and "could not resolve" in note3.lower()


def test_highlight_cell_type_drives_the_main_map_not_an_embedded_figure():
    """A "where are the T cells" query lights up the MAIN specimen canvas (map_view + emphasise),
    rather than dumping a thumbnail into the chat."""
    from spatialscribe.analysis import capabilities as cap

    a = _section()
    ctx = cap.RunContext(tissue="melanoma")
    res = cap.run(a, "highlight_cells", {"criterion": "Myeloid"}, ctx)
    assert res.ok and res.value["rendered_view"] is True and res.value["drove_map"] is True
    assert res.value["n_matched"] == int((a.obs["cell_type"].astype(str) == "Myeloid").sum())
    # a map_view directive that recolours by cell_type and spotlights the matched category - no figure
    assert len(ctx.artifacts) == 1 and ctx.artifacts[0]["kind"] == "map_view"
    assert ctx.artifacts[0]["color_by"] == "cell_type" and ctx.artifacts[0]["highlight"] == "Myeloid"


def test_highlight_non_celltype_criterion_drives_the_main_map():
    """A criterion that is NOT a cell-type category (e.g. low quality) now DRIVES the MAIN map -
    writes a transient obs highlight + a map_view directive - rather than dumping an inline figure
    (same 'drive the real plot' UX as a cell-type highlight; see 'drive the MAIN map for subset
    highlights')."""
    from spatialscribe.analysis import capabilities as cap

    a = _section()
    ctx = cap.RunContext(tissue="melanoma")
    res = cap.run(a, "highlight_cells", {"criterion": "low quality"}, ctx)
    assert res.ok and res.value["rendered_view"] is True and res.value.get("drove_map") is True
    assert len(ctx.artifacts) == 1 and ctx.artifacts[0]["kind"] == "map_view"
    assert ctx.artifacts[0]["color_by"] == "_copilot_highlight"


def test_highlight_cells_is_copilot_exposed():
    from spatialscribe.analysis import capabilities as cap

    assert "highlight_cells" in {t["name"] for t in cap.copilot_tools()}


def test_show_spatial_gene_field_and_unknown_are_grounded():
    from spatialscribe.analysis import capabilities as cap

    a = _section()
    # by a real (on-panel) gene -> drives the MAIN canvas (map_view), not an inline figure: its
    # expression is written to a persistent obs column so /api/points can colour by it, and the
    # column is recorded in uns so the app lists it in the colour dropdown.
    ctx = cap.RunContext()
    r = cap.run(a, "show_spatial", {"color_by": "g0"}, ctx)
    assert r.ok and r.value["rendered_view"] and len(ctx.artifacts) == 1
    assert ctx.artifacts[0]["kind"] == "map_view" and ctx.artifacts[0]["color_by"] == "g0"
    assert "g0" in a.obs and "g0" in list(a.uns.get("gene_color_fields", []))
    # by a field -> map_view
    ctx2 = cap.RunContext()
    r2 = cap.run(a, "show_spatial", {"color_by": "cell type"}, ctx2)
    assert r2.ok and r2.value["rendered_view"] and len(ctx2.artifacts) == 1
    assert ctx2.artifacts[0]["kind"] == "map_view"
    # unknown key -> grounded: no figure, honest note, no error
    ctx3 = cap.RunContext()
    r3 = cap.run(a, "show_spatial", {"color_by": "not_a_gene_zzz"}, ctx3)
    assert r3.ok and r3.value["rendered_view"] is False and len(ctx3.artifacts) == 0
    assert "neither an obs field nor an on-panel gene" in r3.value["note"]


def test_dotplot_violin_composition_emit_figures():
    from spatialscribe.analysis import capabilities as cap

    a = _section()
    for name, params in [("marker_dotplot", {"genes": ["g0", "g1", "g2"]}),
                         ("expression_violin", {"gene": "g0"}),
                         ("composition_chart", {})]:
        ctx = cap.RunContext()
        r = cap.run(a, name, params, ctx)
        assert r.ok, (name, r.error)
        assert len(ctx.artifacts) == 1 and ctx.artifacts[0]["kind"] == "figure"

    # off-panel gene in a violin is grounded: no figure, a note, still ok
    ctx = cap.RunContext()
    r = cap.run(a, "expression_violin", {"gene": "NOPE"}, ctx)
    assert r.ok and r.value["rendered_view"] is False and len(ctx.artifacts) == 0


def test_segmentation_load_and_capability(tmp_path):
    import anndata
    import pandas as pd
    from spatialscribe.analysis import capabilities as cap, views

    rows = []
    for cid in ["c0", "c1", "c2"]:                 # 3 unit-square cells
        for dx, dy in [(0, 0), (1, 0), (1, 1), (0, 1)]:
            rows.append({"cell_id": cid, "vertex_x": float(dx), "vertex_y": float(dy)})
    p = tmp_path / "cell_boundaries.parquet"
    pd.DataFrame(rows).to_parquet(p)
    assert views.load_boundaries(str(tmp_path)) is not None    # found via a directory
    assert views.load_boundaries(str(p)) is not None           # or a direct file

    a = anndata.AnnData(X=np.zeros((3, 2), dtype="float32"))
    a.obs_names = ["c0", "c1", "c2"]; a.obsm["spatial"] = np.zeros((3, 2))
    a.obs["cell_type"] = pd.Categorical(["A", "B", "A"])
    # no boundaries recorded -> honest note, no figure, still ok
    ctx = cap.RunContext()
    r = cap.run(a, "show_segmentation", {}, ctx)
    assert r.ok and r.value["rendered_view"] is False and "cell_boundaries" in r.value["note"]
    # boundaries recorded -> renders polygons
    a.uns["boundaries_path"] = str(p)
    ctx2 = cap.RunContext()
    r2 = cap.run(a, "show_segmentation", {}, ctx2)
    assert r2.ok and r2.value["rendered_view"] and len(ctx2.artifacts) == 1
