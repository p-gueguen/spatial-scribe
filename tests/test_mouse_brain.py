"""Mouse-brain support: species-adaptive marker matching + a brain (CNS) marker context.

Regression guard for the two things that blocked a mouse brain Xenium section:
  1. gene symbols are mouse title-case (``Slc17a7``), but every built-in marker set is human
     uppercase (``SLC17A7``) and the matcher was exact/case-sensitive -> 0 markers resolved;
  2. the app had no CNS cell-type vocabulary (neurons / glia), only tumour-immune-stromal.

All synthetic - no on-disk data, no network.
"""

from __future__ import annotations


def _mouse_brain_adata(n: int = 400, seed: int = 0):
    """Tiny title-case mouse-brain-like section: a few real CNS markers + spanning background.

    The background rates span a wide range on purpose - ``sc.tl.score_genes`` bins genes by mean
    expression to pick control sets, and a uniformly-low background yields "no control genes".
    """
    import anndata as ad
    import numpy as np

    rng = np.random.default_rng(seed)
    markers = ["Slc17a7", "Satb2", "Gad1", "Gad2", "Aqp4", "Gfap",
               "Sox10", "Opalin", "Pdgfra", "Cspg4", "Cldn5", "Pecam1"]
    n_bg = 100
    bg = [f"Bg{i}" for i in range(n_bg)]
    M = np.zeros((n, len(markers) + n_bg), dtype="float32")
    M[:, :len(markers)] = rng.poisson(3.0, size=(n, len(markers)))
    bg_rates = rng.uniform(0.5, 20.0, size=n_bg)                 # spanning -> control genes exist
    M[:, len(markers):] = rng.poisson(np.broadcast_to(bg_rates, (n, n_bg)))
    a = ad.AnnData(X=M)
    a.var_names = markers + bg
    a.obs_names = [f"c{i}" for i in range(n)]
    return a


# --------------------------------------------------------------------------- #
# species-adaptive matcher
# --------------------------------------------------------------------------- #
def test_present_maps_human_markers_onto_mouse_panel():
    from spatialscribe.analysis import markers as m

    panel = ["Acta2", "Aqp4", "Gad1", "Foobar"]                 # mouse title-case panel
    got = m.present(panel, {"x": ["ACTA2", "AQP4", "GAD1", "NOPE"]})
    # resolved to the PANEL's actual symbols (so adata[:, genes] indexing is valid), NOPE dropped.
    assert got["x"] == ["Acta2", "Aqp4", "Gad1"]


def test_on_panel_flat_helper():
    from spatialscribe.analysis import markers as m

    assert m.on_panel(["Aqp4", "Gfap"], ["AQP4", "GFAP", "ZZZ"]) == ["Aqp4", "Gfap"]


def test_present_exact_first_and_ambiguity_guard():
    from spatialscribe.analysis import markers as m

    # A human panel is unchanged: exact case-sensitive match hits first.
    assert m.present(["ACTA2", "CD3D", "GFAP"], {"x": ["ACTA2", "CD3D"]})["x"] == ["ACTA2", "CD3D"]
    # If two panel genes share a lowercase form the fold refuses to guess (returns nothing), so it
    # can never silently conflate distinct symbols.
    assert m.present(["Abc", "ABC"], {"x": ["abc"]})["x"] == []


# --------------------------------------------------------------------------- #
# brain marker context + routing
# --------------------------------------------------------------------------- #
def test_for_tissue_routes_brain():
    from spatialscribe.analysis import markers as m

    for t in ["mouse brain", "adult mouse cortex", "hippocampus CA1", "cerebellum",
              "CNS", "neural tissue", "cortical neurons"]:
        assert m.for_tissue(t) is m.BRAIN_MARKERS, t
    assert m.for_tissue("breast") is m.EPITHELIAL_LINEAGES
    assert m.for_tissue("melanoma") is m.LINEAGE_MARKERS


def test_brain_markers_resolve_and_score_on_mouse_panel():
    from spatialscribe.analysis import annotate, markers as m

    a = _mouse_brain_adata()
    pres = m.present(a.var_names, m.BRAIN_MARKERS)
    resolved = {k: v for k, v in pres.items() if v}
    assert len(resolved) >= 4                                   # neurons + glia resolve, not zero
    assert "Slc17a7" in pres["Excitatory neuron"]               # returns the panel's title-case symbol
    assert all(g in set(a.var_names) for gs in pres.values() for g in gs)
    # end-to-end: the annotation scorer produces one column per resolvable lineage, no error.
    colmap = annotate.score_marker_sets(a, m.BRAIN_MARKERS)
    assert len(colmap) >= 4


# --------------------------------------------------------------------------- #
# reference registry
# --------------------------------------------------------------------------- #
def test_mouse_brain_reference_registered_and_ranked():
    from spatialscribe.analysis import reference as r

    assert "mouse_brain" in r.REFERENCE_REGISTRY
    # a mouse-specific query surfaces it above the human brain-tumour entry (no path needed).
    ranked = r.choose_reference("mouse cortex")
    assert ranked[0]["tissue_key"] == "mouse_brain"


# --------------------------------------------------------------------------- #
# annotate threads the tissue's markers (regression: mouse brain must NOT get melanoma labels)
# --------------------------------------------------------------------------- #
def test_label_context_is_species_aware():
    from spatialscribe.analysis import capabilities as cap

    assert cap.RunContext(tissue="mouse brain").label_context() == "mouse brain"
    assert cap.RunContext(tissue="melanoma").label_context() == "human melanoma"   # demo default unchanged


def test_consensus_annotate_honours_supplied_markers():
    import numpy as np

    from spatialscribe.analysis import annotate, markers as m

    a = _mouse_brain_adata()
    a.obs["leiden"] = np.where(np.arange(a.n_obs) % 2 == 0, "0", "1").astype(object)
    a.obs["leiden"] = a.obs["leiden"].astype("category")
    annotate.consensus_annotate(a, cluster_key="leiden", use_llm=False, marker_sets=m.BRAIN_MARKERS)
    labels = set(a.obs["cell_type"].astype(str))
    assert labels <= set(m.BRAIN_MARKERS) | {"Unknown"}      # brain vocabulary
    assert not (labels & set(m.LINEAGE_MARKERS))             # never the melanoma default
