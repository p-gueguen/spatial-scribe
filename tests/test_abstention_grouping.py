"""Abstention labels are not cell types, and every consumer must group on the SAME column.

Measured defect (2026-07): ``export.report_figures`` grouped the DELIVERED composition + spatial
figures on ``cell_type_final`` - which carries the five abstention pseudo-labels as ordinary
categories of the same pd.Categorical - while ``niches`` / ``spatial.nhood_enrichment`` /
``states`` grouped on ``cell_type`` (carrying plain "Unassigned"/"Unknown"). So a delivered legend
showed "Ambiguous: mixed" as if it were a lineage, and its composition percentages had a different
denominator from the enrichment z-scores printed beside them.

Fix: one source of truth for what an abstention label IS (``annotate.ABSTENTION_LABELS``), one
``annotate.annotation_key`` that every consumer groups on, abstained cells excluded from
composition (with the abstained fraction reported, not hidden), and collapsed into a single
"Not assigned" category wherever they must stay in the spatial graph.
"""
from __future__ import annotations

import numpy as np
import pytest

anndata = pytest.importorskip("anndata")


def _labelled(n=60, seed=0):
    """Synthetic section: 2 real lineages + every abstention label, with spatial coords."""
    import pandas as pd

    from spatialscribe.analysis import annotate

    rng = np.random.default_rng(seed)
    a = anndata.AnnData(X=rng.poisson(5, size=(n, 6)).astype("float32"))
    a.var_names = [f"G{i}" for i in range(6)]
    a.obsm["spatial"] = rng.uniform(0, 100, size=(n, 2))

    abst = sorted(annotate.ABSTENTION_LABELS)
    labels = ["T cell"] * 25 + ["Myeloid"] * 20
    labels += [abst[i % len(abst)] for i in range(n - len(labels))]
    a.obs["cell_type"] = pd.Categorical(labels)
    a.obs["cell_type_final"] = pd.Categorical(labels)
    return a


def test_marker_dotplot_derives_discriminative_genes_for_any_vocabulary():
    """The report marker dot-plot must pick each type's own discriminative gene from the section's
    DEGs (not a curated dict), so it works on arbitrary reference-transfer labels."""
    import pandas as pd

    from spatialscribe.analysis import export

    rng = np.random.default_rng(0)
    n_per, n_genes = 40, 20
    blocks, labels = [], []
    for ti, name in enumerate(("osteoblast", "myotube", "erythroid progenitor cell")):
        M = rng.poisson(1.0, size=(n_per, n_genes)).astype("float32")
        M[:, ti] += rng.poisson(15.0, size=n_per)          # a private high gene per type
        blocks.append(M)
        labels += [name] * n_per
    a = anndata.AnnData(X=np.vstack(blocks))
    a.var_names = [f"g{i}" for i in range(n_genes)]
    a.obs["cell_type"] = pd.Categorical(labels)

    fig = export.fig_marker_dotplot(a, "cell_type")
    assert fig is not None
    ax = fig.axes[0]
    genes = {t.get_text() for t in ax.get_xticklabels()}
    types = {t.get_text() for t in ax.get_yticklabels()}
    assert {"osteoblast", "myotube", "erythroid progenitor cell"} <= types
    assert {"g0", "g1", "g2"} <= genes                     # each type's private gene was selected


def test_marker_dotplot_absent_without_labels():
    from spatialscribe.analysis import export

    a = anndata.AnnData(X=np.random.default_rng(0).poisson(3, size=(30, 6)).astype("float32"))
    assert export.fig_marker_dotplot(a, "cell_type") is None   # no cell_type -> honest skip, no raise


def test_annotation_quality_section_renders_aqi_and_caveat():
    from spatialscribe.analysis import export

    a = anndata.AnnData(X=np.zeros((10, 3), dtype="float32"))
    # index_only regime (no >=3-method ensemble, G is None) -> the caveat must appear
    a.uns["annotation_quality"] = {"aqi": {"aqi": 0.25, "argmin": "contamination",
                                           "components": {"A": 0.72, "C": 0.23, "M": 0.85, "G": None}}}
    frag = export._section_annotation_quality(a)
    assert frag and "0.25" in frag and "0.72" in frag and "0.23" in frag
    assert "ensemble" in frag                                   # the honest index_only caveat
    # no AQI stamped -> section is omitted, never raises
    assert export._section_annotation_quality(anndata.AnnData(X=np.zeros((4, 2), dtype="float32"))) is None


def test_abstention_labels_are_a_single_source_of_truth():
    """The five _ABSTAIN strings plus the two plain ones, defined once and reused by the backend."""
    from spatialscribe.analysis import annotate

    assert set(annotate._ABSTAIN.values()).issubset(annotate.ABSTENTION_LABELS)
    assert {"Unassigned", "Unknown"}.issubset(annotate.ABSTENTION_LABELS)
    assert annotate.is_abstention("Ambiguous: mixed")
    assert annotate.is_abstention("Unassigned")
    assert not annotate.is_abstention("T cell")
    # A real lineage that merely starts with a scary word must not be swallowed by prefix matching.
    assert not annotate.is_abstention("Novel epithelial subtype X")


def test_backend_reuses_the_same_constant():
    """backend/app.py must not keep its own private copy of the abstention vocabulary."""
    pytest.importorskip("fastapi")
    from backend import app

    from spatialscribe.analysis import annotate

    assert app._ABSTENTION_PREFIXES is annotate.ABSTENTION_PREFIXES


def test_annotation_key_prefers_cell_type_final():
    from spatialscribe.analysis import annotate

    a = _labelled()
    assert annotate.annotation_key(a) == "cell_type_final"
    del a.obs["cell_type_final"]
    assert annotate.annotation_key(a) == "cell_type"


def test_composition_excludes_abstained_cells_and_reports_the_fraction():
    """Composition percentages must be over TYPED cells, and the abstained share must be surfaced."""
    from spatialscribe.analysis import annotate, export

    a = _labelled()
    vc, abstained_frac = export.composition_table(a)

    assert not any(annotate.is_abstention(str(i)) for i in vc.index)
    assert set(vc.index) == {"T cell", "Myeloid"}
    assert vc.sum() == 45                                  # 25 + 20 typed cells
    expected = (a.n_obs - 45) / a.n_obs
    assert abs(abstained_frac - expected) < 1e-9
    assert abstained_frac > 0


def test_delivered_figures_never_show_an_abstention_label_as_a_cell_type():
    """The composition figure's tick labels are lineages only; the abstained share goes in the title."""
    from spatialscribe.analysis import annotate, export

    figs = dict(export.report_figures(_labelled()))
    fig = figs["composition"]
    ax = fig.axes[0]
    ticks = [t.get_text() for t in ax.get_yticklabels()]
    assert ticks, "composition figure must have labelled bars"
    assert not any(annotate.is_abstention(t) for t in ticks)
    assert "abstain" in ax.get_title().lower()


def test_spatial_figure_collapses_abstention_into_one_not_assigned_category():
    """Five pseudo-labels must not become five legend entries; abstained cells stay visible as one."""
    from spatialscribe.analysis import annotate, export

    figs = dict(export.report_figures(_labelled()))
    ax = figs["spatial_celltypes"].axes[0]
    legend = ax.get_legend()
    assert legend is not None
    texts = [t.get_text() for t in legend.get_texts()]
    assert annotate.NOT_ASSIGNED in texts
    assert sum(t == annotate.NOT_ASSIGNED for t in texts) == 1
    assert not any(annotate.is_abstention(t) and t != annotate.NOT_ASSIGNED for t in texts)


def test_collapse_abstention_merges_every_pseudo_label():
    from spatialscribe.analysis import annotate

    a = _labelled()
    out = annotate.collapse_abstention(a.obs["cell_type_final"])
    cats = set(map(str, out.cat.categories))
    assert cats == {"T cell", "Myeloid", annotate.NOT_ASSIGNED}


def test_nhood_enrichment_does_not_report_abstention_as_a_cell_type():
    """Abstained cells stay in the spatial graph (they occupy space) but are not an enrichment class."""
    pytest.importorskip("squidpy")
    from spatialscribe.analysis import annotate, spatial

    a = _labelled(n=120)
    res = spatial.nhood_enrichment(a, cluster_key=annotate.annotation_key(a))
    assert not any(annotate.is_abstention(c) for c in res["categories"])
    z = np.asarray(res["zscore"])
    assert z.shape == (len(res["categories"]), len(res["categories"]))


def test_nhood_enrichment_leaves_no_temp_column_or_stray_uns_key():
    """The collapse is an implementation detail: it must not leak into obs or uns (h5ad round-trip)."""
    pytest.importorskip("squidpy")
    from spatialscribe.analysis import annotate, keys, spatial

    a = _labelled(n=120)
    key = annotate.annotation_key(a)
    spatial.nhood_enrichment(a, cluster_key=key)

    assert not [c for c in a.obs.columns if c.startswith("_nhood")], "temp obs column leaked"
    assert "_nhood_collapsed_nhood_enrichment" not in a.uns, "temp uns key leaked"
    # The capability's declared produces-key must exist, and its matrix must match its categories.
    published = keys.nhood_enrichment_key(key).name
    assert published in a.uns
    z = np.asarray(a.uns[published]["zscore"])
    cats = a.uns[published]["categories"]
    assert z.shape == (len(cats), len(cats))
    assert not any(annotate.is_abstention(c) for c in cats)


def test_neighborhood_composition_collapses_abstention():
    pytest.importorskip("squidpy")
    from spatialscribe.analysis import annotate, niches

    a = _labelled(n=120)
    _, cats = niches.neighborhood_composition(a, cluster_key=annotate.annotation_key(a))
    assert not any(annotate.is_abstention(c) and c != annotate.NOT_ASSIGNED for c in cats)


# --------------------------------------------------------------------------- #
# Consumers found by the adversarial review: fixing composition/nhood/niches while leaving states,
# the interactive chart, immune_exclusion and the app's spatial colouring on the raw column would
# have re-created the very inconsistency this change set exists to remove.
# --------------------------------------------------------------------------- #
def test_state_by_celltype_has_no_abstention_rows():
    """'mean IFN score of the Ambiguous: mixed cells' is not a biological statement."""
    from spatialscribe.analysis import annotate, states

    a = _labelled(n=80)
    a.var_names = ["MKI67", "ISG15", "IFI6", "HSPA1A", "VIM", "CD3D"]
    mat = states.state_by_celltype(a, cluster_key=annotate.annotation_key(a))
    assert not mat.empty
    assert not any(annotate.is_abstention(str(i)) for i in mat.index)


def test_interactive_composition_view_shares_the_report_denominator():
    """The plotly artifact and the delivered PNG must not disagree about what fraction is 'T cell'."""
    pytest.importorskip("plotly")
    from spatialscribe.analysis import annotate, export, views

    a = _labelled()
    counts, abstained = export.composition_table(a)
    fig = views.composition_view(a)
    bar = fig.data[0]
    assert not any(annotate.is_abstention(str(x)) for x in bar.x)
    assert int(sum(bar.y)) == int(counts.sum())
    assert "abstained" in fig.layout.title.text


def test_nhood_enrichment_on_a_fully_abstained_section_returns_empty_not_a_traceback():
    """A barren panel stamps every cell 'Unresolvable: panel'; squidpy would raise on 1 category."""
    pytest.importorskip("squidpy")
    import pandas as pd

    from spatialscribe.analysis import annotate, spatial

    a = _labelled(n=60)
    a.obs["cell_type_final"] = pd.Categorical([annotate._ABSTAIN["unresolvable_panel"]] * a.n_obs)
    res = spatial.nhood_enrichment(a, cluster_key="cell_type_final")
    assert res["categories"] == [] and res["zscore"] == []
    assert res["pct_abstained_in_graph"] == 1.0
    assert "note" in res
    assert not [c for c in a.obs.columns if c.startswith("_nhood")]


def test_nhood_enrichment_discloses_abstained_cells_left_in_the_permutation_null():
    pytest.importorskip("squidpy")
    from spatialscribe.analysis import annotate, spatial

    a = _labelled(n=120)
    res = spatial.nhood_enrichment(a, cluster_key=annotate.annotation_key(a))
    assert 0 < res["pct_abstained_in_graph"] < 1


def test_niches_are_never_named_after_the_abstention_class():
    """'Niche 3: Not assigned + T cell' describes annotation failure, not biology."""
    pytest.importorskip("squidpy")
    from spatialscribe.analysis import annotate, niches

    a = _labelled(n=150)
    df = niches.call_niches(a, cluster_key=annotate.annotation_key(a), n_niches=3)
    for row in df.to_dict("records"):
        assert annotate.NOT_ASSIGNED not in row["niche"]
        assert not any(annotate.is_abstention(t) for t in row["top_types"])
