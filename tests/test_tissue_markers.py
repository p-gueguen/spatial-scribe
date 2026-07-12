"""Free-text tissue -> marker panel resolution (any tissue), with the LLM boundary mocked.

Curated tissues use their offline set; an unrecognised free-text tissue with a panel + API key gets
a Claude-generated panel grounded to the panel's genes (cached). No network in these tests - the
``llm`` boundary is monkeypatched.
"""

from __future__ import annotations


def _marked(labels, counts, n_genes=30, seed=0):
    """Synthetic COUNTS section/reference: label i is marked by its private high gene ``g{i}``."""
    import anndata
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(seed)
    blocks, labs = [], []
    for i, (name, n) in enumerate(zip(labels, counts)):
        M = rng.poisson(1.0, size=(n, n_genes)).astype("float32")
        M[:, i] += rng.poisson(30.0, size=n)
        blocks.append(M)
        labs += [name] * n
    a = anndata.AnnData(X=np.vstack(blocks).astype("float32"))
    a.var_names = [f"g{i}" for i in range(n_genes)]
    a.obs["cell_type"] = pd.Categorical(labs)
    a.obsm["spatial"] = rng.uniform(0, 100, size=(a.n_obs, 2))
    return a


# The fine reference-transfer vocabulary (MOCA/Allen atlas labels) that no curated human-tissue
# marker dict anticipates - the case where the AQI purity term collapsed to a coverage artifact.
_ATLAS_LABELS = ("osteoblast", "myotube", "erythroid progenitor cell", "fibroblast")


def test_ctx_markers_grounds_on_reference_when_curated_is_low():
    # A section carrying an atlas vocabulary the curated dicts miss, WITH an in-memory reference:
    # ctx.markers must ground the markers on the reference's own per-label DEGs (not the generic
    # EPITHELIAL_LINEAGES fallback), so downstream purity/fidelity can actually be computed.
    from spatialscribe.analysis import capabilities as cap

    ref = _marked(_ATLAS_LABELS, [60] * 4, seed=0)
    sec = _marked(_ATLAS_LABELS, [40] * 4, seed=1)
    ctx = cap.RunContext(tissue="mouse embryo", use_llm=False, reference=ref, ref_label_key="cell_type")
    ms = ctx.markers(sec)

    covered = [lab for lab in _ATLAS_LABELS if ms.get(lab)]
    assert len(covered) >= 3                                # grounded on the reference, not left uncovered
    assert "g0" in ms["osteoblast"]                         # the reference-derived private marker


def test_reference_grounded_markers_measure_purity_not_coverage():
    # The payoff + negative control: with reference-grounded markers the PMP is DEFINED for ~all cells
    # (retention -> 1, killing the C=coverage artifact) AND it still drops when the labels are wrong -
    # so it measures real purity, not just coverage.
    import numpy as np

    from spatialscribe.analysis import capabilities as cap
    from spatialscribe.analysis import purity as p

    ref = _marked(_ATLAS_LABELS, [60] * 4, seed=0)
    sec = _marked(_ATLAS_LABELS, [50] * 4, seed=2)
    ctx = cap.RunContext(tissue="mouse embryo", use_llm=False, reference=ref, ref_label_key="cell_type")
    ms = ctx.markers(sec)

    p.pmp(sec, assigned_label_key="cell_type", lineage_markers=ms)
    pmp_clean = np.asarray(sec.obs["pmp"], dtype=float)
    assert np.mean(~np.isnan(pmp_clean)) > 0.9              # retention ~1 (was ~0.15 on generic fallback)

    scr = sec.copy()                                        # same labels+coverage, shuffled onto cells
    scr.obs["cell_type"] = (scr.obs["cell_type"].sample(frac=1.0, random_state=0).to_numpy())
    p.pmp(scr, assigned_label_key="cell_type", lineage_markers=ms)
    pmp_scr = np.asarray(scr.obs["pmp"], dtype=float)
    assert np.mean(~np.isnan(pmp_scr)) > 0.9               # coverage identical...
    assert np.nanmedian(pmp_clean) > np.nanmedian(pmp_scr) + 0.1   # ...but real purity is higher when correct


def test_curated_tissue_wins_without_calling_llm(monkeypatch):
    from spatialscribe.analysis import llm, markers as m

    def boom(*a, **k):
        raise AssertionError("the LLM must not be called for a curated tissue")

    monkeypatch.setattr(llm, "generate_marker_panel", boom)
    assert m.resolve_markers("mouse brain", panel_genes=["Aqp4"], use_llm=True) is m.BRAIN_MARKERS
    assert m.resolve_markers("human breast carcinoma", panel_genes=["EPCAM"], use_llm=True) is m.EPITHELIAL_LINEAGES
    assert m.resolve_markers("cutaneous melanoma", panel_genes=["MLANA"], use_llm=True) is m.LINEAGE_MARKERS


def test_generic_marker_fallback_stamps_a_wrong_tissue_caveat():
    # An unrecognised organ (kidney) falls back to the generic breast/carcinoma EPITHELIAL_LINEAGES
    # set; that must be surfaced (uns['marker_set_warning']) so a kidney is not silently typed
    # 'Epithelial/Tumor' with no tell. A curated tissue (breast) must NOT warn.
    from spatialscribe.analysis import capabilities as cap
    from spatialscribe.analysis import demo
    from spatialscribe.analysis import markers as m

    assert not m.tissue_has_curated_set("kidney")
    assert m.tissue_has_curated_set("mouse brain") and m.tissue_has_curated_set("human breast")

    a = demo.make_demo_adata(n_cells=200, seed=0)
    res = cap.RunContext(tissue="kidney", use_llm=False).markers(a)
    assert res is m.EPITHELIAL_LINEAGES                      # the generic fallback
    assert "marker_set_warning" in a.uns
    assert "kidney" in a.uns["marker_set_warning"]["message"]

    b = demo.make_demo_adata(n_cells=200, seed=1)
    cap.RunContext(tissue="human breast", use_llm=False).markers(b)   # curated -> no caveat
    assert "marker_set_warning" not in b.uns


def test_unknown_tissue_generates_grounded_panel_and_caches(monkeypatch):
    from spatialscribe.analysis import llm, markers as m

    m._TISSUE_PANEL_CACHE.clear()
    panel = ["Alb", "Apoa1", "Ttr", "Clec4f", "Csf1r", "Pecam1", "Dcn", "Col1a1", "Bg1", "Bg2"]
    calls = {"n": 0}

    def fake(tissue, panel_genes, max_types=15, cell_types=None):
        calls["n"] += 1
        return {"Hepatocyte": ["ALB", "APOA1", "TTR", "NOTONPANEL"],   # 3 resolve (case-insensitive)
                "Kupffer cell": ["CLEC4F", "CSF1R"],
                "Ghost": ["ZZZ1", "ZZZ2"]}                             # 0 on panel -> dropped

    monkeypatch.setattr(llm, "generate_marker_panel", fake)
    r1 = m.resolve_markers("mouse liver", panel_genes=panel, use_llm=True)
    assert set(r1) == {"Hepatocyte", "Kupffer cell"}                  # zero-coverage lineage dropped
    assert r1["Hepatocyte"] == ["ALB", "APOA1", "TTR", "NOTONPANEL"]  # raw list kept for panel_check
    r2 = m.resolve_markers("mouse liver", panel_genes=panel, use_llm=True)
    assert r2 == r1 and calls["n"] == 1                               # cached: LLM invoked once


def test_offline_and_empty_reply_fall_back_to_epithelial(monkeypatch):
    from spatialscribe.analysis import llm, markers as m

    m._TISSUE_PANEL_CACHE.clear()
    # no key / use_llm False -> deterministic epithelial fallback (== for_tissue behaviour)
    assert m.resolve_markers("mouse liver", panel_genes=["Alb"], use_llm=False) is m.EPITHELIAL_LINEAGES
    # LLM returns nothing usable -> same fallback, nothing cached
    monkeypatch.setattr(llm, "generate_marker_panel", lambda *a, **k: {})
    assert m.resolve_markers("axolotl gill", panel_genes=["Alb"], use_llm=True) is m.EPITHELIAL_LINEAGES


def test_generate_marker_panel_parses_and_normalizes(monkeypatch):
    from spatialscribe.analysis import llm

    reply = '{"Hepatocyte": ["ALB","APOA1"], "Kupffer": ["CLEC4F"], "bad": "notalist", "x": [1,2]}'
    monkeypatch.setattr(llm, "complete", lambda *a, **k: reply)
    out = llm.generate_marker_panel("liver", ["ALB", "APOA1", "CLEC4F"])
    assert out["Hepatocyte"] == ["ALB", "APOA1"]
    assert out["Kupffer"] == ["CLEC4F"]
    assert "bad" not in out           # non-list value skipped
    assert out["x"] == []             # non-str genes filtered


def test_generate_marker_panel_bad_reply_is_empty(monkeypatch):
    from spatialscribe.analysis import llm

    monkeypatch.setattr(llm, "complete", lambda *a, **k: "sorry, no JSON here")
    assert llm.generate_marker_panel("liver", ["ALB"]) == {}


def test_generate_marker_panel_honours_requested_cell_types(monkeypatch):
    from spatialscribe.analysis import llm

    seen = {}

    def fake_complete(system, user, **k):
        seen["user"] = user
        return '{"My weird type": ["EPCAM"], "T cell": ["CD3D"]}'

    monkeypatch.setattr(llm, "complete", fake_complete)
    out = llm.generate_marker_panel("breast", ["EPCAM", "CD3D"], cell_types=["My weird type", "T cell"])
    assert "My weird type" in seen["user"] and "EXACTLY" in seen["user"]   # types went into the prompt
    assert set(out) == {"My weird type", "T cell"}


def test_markers_for_types_rekeys_to_current_categories(monkeypatch):
    from spatialscribe.analysis import llm, markers as m

    m._TYPES_PANEL_CACHE.clear()
    panel = ["EPCAM", "CD3D", "PECAM1"]
    calls = {"n": 0}

    def fake(tissue, panel_genes, max_types=15, cell_types=None):
        calls["n"] += 1
        # model renames one type and invents a stray - both must be dropped/normalised
        return {"epithelial/tumor": ["EPCAM"], "T cell": ["CD3D"], "Ghost": ["ZZZ"]}

    monkeypatch.setattr(llm, "generate_marker_panel", fake)
    cats = ["Epithelial/Tumor", "T cell", "Endothelial"]
    r = m.markers_for_types(cats, panel_genes=panel, tissue="breast", use_llm=True)
    assert set(r) == {"Epithelial/Tumor", "T cell"}       # case-normalised key kept; Ghost dropped
    assert "Endothelial" not in r                          # model omitted it -> absent, not fabricated
    m.markers_for_types(cats, panel_genes=panel, tissue="breast", use_llm=True)
    assert calls["n"] == 1                                  # cached per (tissue, panel, types)
    # off / unannotated -> empty, so the caller falls back to the curated resolver
    assert m.markers_for_types(cats, panel_genes=panel, tissue="breast", use_llm=False) == {}
    assert m.markers_for_types([], panel_genes=panel, tissue="breast", use_llm=True) == {}
