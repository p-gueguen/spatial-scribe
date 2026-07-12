"""Reference<->panel matching (analysis/reference.py) - synthetic only, no /data, no network.

Covers label detection, Ensembl->symbol var-name repair, keyword + panel-overlap reference
ranking, the global match score + per-type resolvability + clustering nudge, and the
reference_match capability (no-reference skip + with-path path).
"""
from __future__ import annotations

import json

import anndata
import numpy as np
import pandas as pd


def _panel_ref(n=80, n_genes=40, seed=0):
    """A synthetic reference: 3 well-separated types (A/B/C), each with a private high gene."""
    rng = np.random.default_rng(seed)
    blocks, labels = [], []
    for t, hi in [("A", 0), ("B", 1), ("C", 2)]:
        M = rng.poisson(1.0, size=(n, n_genes)).astype("float32")
        M[:, hi] += rng.poisson(20.0, size=n)
        blocks.append(M)
        labels += [t] * n
    a = anndata.AnnData(X=np.vstack(blocks))
    a.var_names = [f"g{i}" for i in range(n_genes)]
    a.obs["cell_type"] = pd.Categorical(labels)
    return a


# --------------------------------------------------------------------------- #
# detect_label_key + load_reference
# --------------------------------------------------------------------------- #
def test_detect_label_key_priority_and_rejections():
    from spatialscribe.analysis import reference as ref

    a = anndata.AnnData(X=np.zeros((60, 3), dtype="float32"))
    a.obs["cell_type"] = pd.Categorical((["T"] * 30) + (["B"] * 30))
    a.obs["author_cell_type"] = pd.Categorical((["x"] * 20) + (["y"] * 40))
    a.obs["cell_type_ontology_term_id"] = pd.Categorical((["CL:1"] * 30) + (["CL:2"] * 30))
    a.obs["barcode"] = [f"c{i}" for i in range(60)]      # near-unique -> not a label
    assert ref.detect_label_key(a) == "cell_type"        # curated priority beats author_cell_type

    b = anndata.AnnData(X=np.zeros((40, 2), dtype="float32"))
    b.obs["annotation"] = pd.Categorical((["p"] * 20) + (["q"] * 20))
    assert ref.detect_label_key(b) == "annotation"       # only sane label column

    c = anndata.AnnData(X=np.zeros((10, 2), dtype="float32"))
    c.obs["id"] = [f"c{i}" for i in range(10)]            # every value unique
    assert ref.detect_label_key(c) is None


def test_load_reference_symbolizes_ensembl(tmp_path):
    from spatialscribe.analysis import reference as ref

    a = _panel_ref()
    a.var_names = [f"ENSG00000{i:05d}" for i in range(a.n_vars)]     # Ensembl ids
    a.var["feature_name"] = [f"g{i}" for i in range(a.n_vars)]       # symbols live here
    p = tmp_path / "ref.h5ad"
    a.write_h5ad(p)
    ad, key = ref.load_reference(p)
    assert key == "cell_type"
    assert list(ad.var_names)[:3] == ["g0", "g1", "g2"]              # var_names became the symbols


def test_load_reference_autodetects_annotation(tmp_path):
    from spatialscribe.analysis import reference as ref

    a = _panel_ref()
    a.obs = a.obs.rename(columns={"cell_type": "annotation"})
    p = tmp_path / "ref.h5ad"
    a.write_h5ad(p)
    _, key = ref.load_reference(p)
    assert key == "annotation"


def test_load_reference_missing_and_unlabelled(tmp_path):
    import pytest

    from spatialscribe.analysis import reference as ref

    with pytest.raises(FileNotFoundError):
        ref.load_reference(tmp_path / "nope.h5ad")

    a = anndata.AnnData(X=np.zeros((10, 2), dtype="float32"))       # no label column
    a.obs["id"] = [f"c{i}" for i in range(10)]
    p = tmp_path / "nolabel.h5ad"
    a.write_h5ad(p)
    with pytest.raises(ValueError):
        ref.load_reference(p)


# --------------------------------------------------------------------------- #
# choose_reference
# --------------------------------------------------------------------------- #
def test_choose_reference_ranks_by_keyword():
    from spatialscribe.analysis import reference as ref

    fake = {
        "breast": {"path": None, "label_key": "cell_type", "gene_name_col": None,
                   "description": "breast atlas", "keywords": ["breast", "mammary", "carcinoma"]},
        "kidney": {"path": None, "label_key": "cell_type", "gene_name_col": None,
                   "description": "kidney atlas", "keywords": ["kidney", "renal"]},
    }
    out = ref.choose_reference("breast carcinoma", registry=fake)
    assert out[0]["tissue_key"] == "breast"
    assert out[0]["available"] is False and out[0]["gene_overlap_frac"] is None
    assert out[0]["score"] >= out[-1]["score"]


def test_choose_reference_reranks_by_panel_overlap(tmp_path):
    from spatialscribe.analysis import reference as ref

    panel = [f"g{i}" for i in range(30)]
    hi = _panel_ref(n=20, n_genes=40); hi.var_names = [f"g{i}" for i in range(40)]   # covers the panel
    lo = _panel_ref(n=20, n_genes=40); lo.var_names = [f"z{i}" for i in range(40)]   # disjoint
    hp, lp = tmp_path / "hi.h5ad", tmp_path / "lo.h5ad"
    hi.write_h5ad(hp); lo.write_h5ad(lp)
    fake = {
        "hi": {"path": str(hp), "label_key": "cell_type", "gene_name_col": None,
               "description": "tissue", "keywords": ["tissue"]},
        "lo": {"path": str(lp), "label_key": "cell_type", "gene_name_col": None,
               "description": "tissue", "keywords": ["tissue"]},
    }
    out = ref.choose_reference("tissue", panel_genes=panel, registry=fake)
    assert out[0]["tissue_key"] == "hi"
    assert out[0]["gene_overlap_frac"] > out[1]["gene_overlap_frac"]
    # score IS the keyword match (overlap is a positive TIE-BREAK, not blended in / a demotion): the
    # measured-overlap ref must not be scored below an unmeasured one that ties on keyword.
    assert out[0]["score"] == out[0]["keyword_score"]
    assert out[0]["score"] == out[1]["score"]                          # tie on keyword; overlap breaks it


# --------------------------------------------------------------------------- #
# reference_panel_match
# --------------------------------------------------------------------------- #
def test_reference_panel_match_good_fit():
    from spatialscribe.analysis import reference as ref

    a = _panel_ref()
    panel = [f"g{i}" for i in range(30)]                 # >= 25 shared, spans the 3 private markers
    m = ref.reference_panel_match(a, panel, "cell_type")
    assert m["status"] == "ok"
    assert m["global"]["verdict"] in {"good", "fair"}
    assert m["clustering_nudge"] is None or m["global"]["verdict"] == "fair"
    assert m["global"]["mean_f1"] is not None
    assert all(isinstance(d["resolvable"], bool) for d in m["per_type"].values())
    json.dumps(m)                                        # fully JSON-able


def test_detect_and_tissue_consistency():
    from spatialscribe.analysis import reference as ref

    a = _panel_ref()
    assert ref.detect_reference_tissue(a) is None                      # no tissue column -> None
    a.obs["tissue"] = pd.Categorical(["kidney"] * a.n_obs)
    assert ref.detect_reference_tissue(a) == "kidney"
    assert ref.tissue_consistency("human breast", a)["consistent"] is False   # breast section, kidney ref
    a.obs["tissue"] = pd.Categorical(["breast"] * a.n_obs)
    assert ref.tissue_consistency("human breast carcinoma", a)["consistent"] is True
    # synonyms do not false-alarm: an eye/uveal query is consistent with an eye reference
    a.obs["tissue"] = pd.Categorical(["eye"] * a.n_obs)
    assert ref.tissue_consistency("uveal melanoma", a)["consistent"] is True
    del a.obs["tissue"]
    assert ref.tissue_consistency("human breast", a)["consistent"] is None     # no metadata -> caveat


def test_reference_panel_match_wrong_tissue_overrides_to_clustering():
    # A well-separated reference that DECLARES a different tissue is a wrong-tissue reference no matter
    # how cleanly its types resolve internally: the verdict must flag the mismatch and route away from
    # transfer (the silent-nonsense guard - a kidney atlas on a breast section).
    from spatialscribe.analysis import reference as ref

    a = _panel_ref()
    a.obs["tissue"] = pd.Categorical(["kidney"] * a.n_obs)
    panel = [f"g{i}" for i in range(30)]
    m = ref.reference_panel_match(a, panel, "cell_type", tissue="human breast")
    assert m["global"]["tissue_mismatch"] is True
    assert m["tissue_check"]["reference_tissue"] == "kidney"
    assert m["recommendation"]["action"] == "unsupervised_clustering"
    assert "MISMATCH" in m["clustering_nudge"]
    # same reference, matching tissue -> no override
    m2 = ref.reference_panel_match(a, panel, "cell_type", tissue="human kidney")
    assert m2["global"]["tissue_mismatch"] is False


def test_plan_annotation_strategy_routes_to_cluster_on_wrong_tissue():
    from spatialscribe.analysis import reference as ref

    a = _panel_ref()
    a.obs["tissue"] = pd.Categorical(["kidney"] * a.n_obs)
    plan = ref.plan_annotation_strategy(a, [f"g{i}" for i in range(30)], "cell_type",
                                        tissue="human breast")
    assert plan["tissue_mismatch"] is True
    assert plan["recommended_mode"] == "cluster"


def test_reference_panel_match_poor_fit_fires_nudge():
    from spatialscribe.analysis import reference as ref

    a = _panel_ref()
    m = ref.reference_panel_match(a, ["zzz1", "zzz2"], "cell_type")     # < 25 shared genes
    assert m["status"] == "insufficient_overlap"
    assert m["global"]["verdict"] == "poor"
    assert m["clustering_nudge"] and "unsupervised clustering" in m["clustering_nudge"]


def test_reference_panel_match_confusable_pair_not_resolvable():
    from spatialscribe.analysis import reference as ref

    # C and D differ only OFF-panel -> at least one is not resolvable on the panel.
    rng = np.random.default_rng(0)
    blocks, labels = [], []
    for t, hi in [("A", 0), ("B", 1), ("C", 2), ("D", 2)]:
        M = rng.poisson(1.0, size=(100, 40)).astype("float32")
        M[:, hi] += rng.poisson(20.0, size=100)
        if t == "D":
            M[:, 39] += rng.poisson(20.0, size=100)      # D's only distinct gene is g39 (off-panel)
        blocks.append(M); labels += [t] * 100
    a = anndata.AnnData(X=np.vstack(blocks))
    a.var_names = [f"g{i}" for i in range(40)]
    a.obs["cell_type"] = pd.Categorical(labels)
    m = ref.reference_panel_match(a, [f"g{i}" for i in range(39)], "cell_type")   # panel omits g39
    assert m["status"] == "ok"
    assert any(not d["resolvable"] for d in m["per_type"].values())


# --------------------------------------------------------------------------- #
# reference_match capability
# --------------------------------------------------------------------------- #
def test_reference_match_capability_no_reference(processed_adata, ctx):
    from spatialscribe.analysis import capabilities as cap

    res = cap.run(processed_adata, "reference_match", {}, ctx)       # ctx.reference is None
    assert res.ok
    assert processed_adata.uns["reference_match"]["status"] == "no_reference"
    json.dumps(res.value)


def test_reference_match_capability_with_path(processed_adata, ctx, tmp_path):
    from spatialscribe.analysis import capabilities as cap

    panel = processed_adata.var_names[~processed_adata.var["control"].to_numpy(bool)].tolist()
    genes = [str(g) for g in panel[:40]]                             # guaranteed overlap with the panel
    rng = np.random.default_rng(0)
    blocks, labels = [], []
    for t, hi in [("A", 0), ("B", 1), ("C", 2)]:
        M = rng.poisson(1.0, size=(60, len(genes))).astype("float32")
        M[:, hi] += rng.poisson(20.0, size=60)
        blocks.append(M); labels += [t] * 60
    r = anndata.AnnData(X=np.vstack(blocks))
    r.var_names = genes
    r.obs["cell_type"] = pd.Categorical(labels)
    p = tmp_path / "ref.h5ad"
    r.write_h5ad(p)

    res = cap.run(processed_adata, "reference_match", {"reference_path": str(p)}, ctx)
    assert res.ok, res.error
    rm = processed_adata.uns["reference_match"]
    assert rm["status"] == "ok"
    assert "global" in res.value and "clustering_nudge" in res.value
    assert rm["global"]["target_depth"] > 0             # derived from obs['total_counts']


def test_choose_reference_short_token_does_not_hijack():
    """Regression: a skin/kidney query must not be won by the breast atlas via a stray 1-2 char
    token (the pre-fix "k"-from-"K=29" substring bug). Uses the built-in registry (env paths unset
    -> all unavailable -> keyword-only ranking, the exact path where nothing corrected the tie)."""
    from spatialscribe.analysis import reference as ref

    assert ref.choose_reference("skin")[0]["tissue_key"] == "skin"
    assert ref.choose_reference("kidney")[0]["tissue_key"] == "kidney"


def test_reference_panel_match_nan_depth_falls_back():
    """Regression: a NaN target_depth (empty/corrupt total_counts) must fall back to 50, not crash."""
    from spatialscribe.analysis import reference as ref

    a = _panel_ref()
    m = ref.reference_panel_match(a, [f"g{i}" for i in range(30)], "cell_type", target_depth=float("nan"))
    assert m["status"] == "ok"
    assert m["global"]["target_depth"] == 50.0


def test_detect_label_key_accepts_granular_reference():
    """Regression: a small/granular reference (more types than half its cells) is not rejected."""
    from spatialscribe.analysis import reference as ref

    a = anndata.AnnData(X=np.zeros((100, 3), dtype="float32"))
    a.obs["cell_type"] = pd.Categorical([f"t{i % 60}" for i in range(100)])   # 60 types / 100 cells
    assert ref.detect_label_key(a) == "cell_type"


def test_reference_key_and_capability_registered():
    from spatialscribe.analysis import capabilities as cap
    from spatialscribe.analysis import keys

    assert keys.Uns.REFERENCE_MATCH in keys.all_keys()
    assert "reference_match" in cap.copilot_names()


# --------------------------------------------------------------------------- #
# CELLxGENE routing: select_reference_cells + resolve_census_tissue + _assay_rank
# --------------------------------------------------------------------------- #
def _census_obs(seed=0):
    """A synthetic census obs frame: 3 common types across 10x v3/v2, a spatial (Visium) block,
    and a rare type - the shape select_reference_cells must route well."""
    rows = []
    jid = 0
    for ct in ("T cell", "B cell", "macrophage"):
        for assay, n in (("10x 3' v3", 100), ("10x 3' v2", 100)):
            for _ in range(n):
                rows.append((jid, ct, assay)); jid += 1
    for _ in range(80):                                   # Visium spatial block - must be excluded
        rows.append((jid, "epithelial cell", "Visium Spatial Gene Expression")); jid += 1
    for _ in range(10):                                   # rare type - must survive stratification
        rows.append((jid, "mast cell", "10x 3' v2")); jid += 1
    df = pd.DataFrame(rows, columns=["soma_joinid", "cell_type", "assay"])
    return df.sample(frac=1.0, random_state=seed).reset_index(drop=True)   # shuffle input order


def test_assay_rank_orders_recent_chemistry_first():
    from spatialscribe.analysis import reference as ref

    assert ref._assay_rank("10x 3' v3") > ref._assay_rank("10x 3' v2") > ref._assay_rank("Smart-seq2")
    assert ref._assay_rank("Smart-seq2") > ref._assay_rank("microwell-seq")
    assert ref._assay_rank("10x Flex") > ref._assay_rank("10x 3' v3")     # newest wins once present
    assert ref._is_spatial_assay("Visium Spatial Gene Expression")
    assert ref._is_spatial_assay("Slide-seqV2")
    assert not ref._is_spatial_assay("10x 3' v3")


def test_select_reference_cells_stratifies_and_prefers_recent_chemistry():
    from spatialscribe.analysis import reference as ref

    obs = _census_obs()
    picked = ref.select_reference_cells(obs, target_cells=300, min_cells_per_type=25)

    # spatial excluded entirely
    assert not picked["assay"].map(ref._is_spatial_assay).any()
    assert "epithelial cell" not in set(picked["cell_type"])
    # every non-spatial type represented, including the rare one
    assert {"T cell", "B cell", "macrophage", "mast cell"} <= set(picked["cell_type"])
    # rare type fully kept (10 cells, under any cap)
    assert (picked["cell_type"] == "mast cell").sum() == 10
    # within a common type the recent chemistry (v3) is preferred over v2
    tcell = picked[picked["cell_type"] == "T cell"]
    assert (tcell["assay"] == "10x 3' v3").sum() >= (tcell["assay"] == "10x 3' v2").sum()
    # target respected (never blows past it)
    assert len(picked) <= 300


def test_select_reference_cells_keeps_spatial_only_tissue():
    from spatialscribe.analysis import reference as ref

    obs = pd.DataFrame({"soma_joinid": range(30), "cell_type": ["x"] * 30,
                        "assay": ["Visium Spatial Gene Expression"] * 30})
    picked = ref.select_reference_cells(obs, target_cells=100)
    assert len(picked) == 30                              # never strip the tissue down to nothing


def test_resolve_census_tissue_maps_free_text():
    from spatialscribe.analysis import reference as ref

    vocab = ["kidney", "lung", "brain", "embryo", "skin of body"]
    assert ref.resolve_census_tissue("kidney", vocab) == "kidney"
    assert ref.resolve_census_tissue("KIDNEY", vocab) == "kidney"           # case-insensitive
    assert ref.resolve_census_tissue("mouse embryo", vocab) == "embryo"     # two-way substring
    assert ref.resolve_census_tissue("skin", vocab) == "skin of body"       # token overlap
    assert ref.resolve_census_tissue("pancreas", vocab) is None             # honest miss


# --------------------------------------------------------------------------- #
# reference_marker_sets: per-label markers from a REFERENCE's own DEGs (fixes the AQI
# purity-coverage artifact on fine reference-transfer vocabularies).
# --------------------------------------------------------------------------- #
def _marked_reference(labels=("osteoblast", "myotube", "erythroid progenitor cell", "fibroblast"),
                      counts=None, n_genes=30, seed=0):
    """A synthetic COUNTS reference: label i is marked by its own private high gene ``g{i}``."""
    rng = np.random.default_rng(seed)
    counts = counts or [60] * len(labels)
    blocks, labs = [], []
    for i, (name, n) in enumerate(zip(labels, counts)):
        M = rng.poisson(1.0, size=(n, n_genes)).astype("float32")
        M[:, i] += rng.poisson(30.0, size=n)               # private marker gene g{i}
        blocks.append(M)
        labs += [name] * n
    a = anndata.AnnData(X=np.vstack(blocks).astype("float32"))
    a.var_names = [f"g{i}" for i in range(n_genes)]
    a.obs["cell_type"] = pd.Categorical(labs)
    return a


def test_reference_marker_sets_derives_panel_markers_per_label():
    from spatialscribe.analysis import reference as ref

    r = _marked_reference()
    panel = [f"g{i}" for i in range(20)]                    # panel covers g0..g3 (the markers) + others
    sets = ref.reference_marker_sets(r, "cell_type", panel, top_n=5, min_cells=20)

    assert set(sets) == {"osteoblast", "myotube", "erythroid progenitor cell", "fibroblast"}
    assert sets["osteoblast"][0] == "g0"                   # each label's private gene ranks first
    assert sets["myotube"][0] == "g1"
    assert all(all(g in panel for g in gs) for gs in sets.values())        # only panel genes


def test_reference_marker_sets_drops_tiny_labels_and_caches():
    from spatialscribe.analysis import reference as ref

    r = _marked_reference(labels=("A", "B", "C", "rare"), counts=[60, 60, 60, 8])
    panel = [f"g{i}" for i in range(20)]
    sets = ref.reference_marker_sets(r, "cell_type", panel, top_n=4, min_cells=20)
    assert "rare" not in sets                              # < min_cells -> dropped (unresolvable, not faked)
    assert {"A", "B", "C"} <= set(sets)
    # second call is served from the reference.uns cache (no recompute)
    sets2 = ref.reference_marker_sets(r, "cell_type", panel, top_n=4, min_cells=20)
    assert sets2 == sets


def test_reference_marker_sets_empty_without_labels_or_overlap():
    from spatialscribe.analysis import reference as ref

    r = _marked_reference()
    assert ref.reference_marker_sets(r, "cell_type", ["zzz1", "zzz2"]) == {}    # no panel overlap
    assert ref.reference_marker_sets(None, "cell_type", ["g0"]) == {}           # no reference
    assert ref.reference_marker_sets(r, "missing_key", ["g0"]) == {}            # no such label col
