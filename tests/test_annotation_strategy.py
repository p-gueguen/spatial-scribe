"""Reference<->panel self-verify + auto-rerun ladder and the supervised-vs-clustering route.

Covers reference.plan_annotation_strategy (the ladder: score -> coarsen confusable labels ->
reselect a better reference -> route) and the annotation_strategy capability. Synthetic only,
no network, no /data (registry references are unset in tests, so the reselect rung is a no-op).
"""

from __future__ import annotations


def _ref_adata(n_per: int = 80, n_genes: int = 40, *, confusable: bool = False,
               disjoint: bool = False, n_types: int = 3, seed: int = 0):
    """A synthetic single-cell reference: each type gets one private high-expressed gene.

    ``confusable`` makes types A and B share the SAME private gene (so the panel confuses them ->
    the coarsen rung should merge them). ``disjoint`` renames genes so they do NOT overlap a
    ``g*`` panel (-> insufficient_overlap -> cluster).
    """
    import anndata as ad
    import numpy as np

    rng = np.random.default_rng(seed)
    types = [chr(ord("A") + i) for i in range(n_types)]
    n = n_per * len(types)
    X = rng.poisson(1.0, size=(n, n_genes)).astype("float32")
    labels: list[str] = []
    for ti, t in enumerate(types):
        priv = 0 if (confusable and t in ("A", "B")) else ti
        X[ti * n_per:(ti + 1) * n_per, priv] += rng.poisson(25.0, size=(n_per,))
        labels += [t] * n_per
    a = ad.AnnData(X=X)
    a.var_names = [(f"REF{i}" if disjoint else f"g{i}") for i in range(n_genes)]
    a.obs["cell_type"] = labels
    a.obs["cell_type"] = a.obs["cell_type"].astype("category")
    return a


_PANEL = [f"g{i}" for i in range(40)]     # matches the non-disjoint reference's gene space


def test_merge_groups_from_confusion_unions_transitively():
    from spatialscribe.analysis import reference as R

    per = {"A": {"confused_with": "B", "confused_frac": 0.6},
           "B": {"confused_with": "A", "confused_frac": 0.5},
           "C": {"confused_with": None, "confused_frac": 0.0},
           "D": {"confused_with": "E", "confused_frac": 0.1}}   # below min_frac -> not merged
    m = R.merge_groups_from_confusion(per, min_frac=0.25)
    assert m["A"] == m["B"] == "A + B"        # merged
    assert m["C"] == "C" and m["D"] == "D"    # unmerged (identity)


def test_good_reference_recommends_supervised_transfer():
    from spatialscribe.analysis import reference as R

    r = R.plan_annotation_strategy(_ref_adata(), _PANEL, "cell_type", target_depth=50.0)
    assert r["recommended_mode"] == "reference_transfer"
    assert r["final_verdict"] == "good"
    assert r["coarsen_map"] is None                     # a good fit needs no coarsening
    assert r["ladder"][0]["step"] == "initial"


def test_confusable_reference_coarsens_then_transfers():
    from spatialscribe.analysis import reference as R

    r = R.plan_annotation_strategy(_ref_adata(confusable=True), _PANEL, "cell_type", target_depth=50.0)
    steps = [s["step"] for s in r["ladder"]]
    assert "coarsen" in steps                           # the auto-rerun fired
    assert r["coarsen_map"] and r["coarsen_map"].get("A") == "A + B"
    # coarsening must not make the verdict worse, and it should end supervised (panel CAN resolve the
    # coarser labels), not clustering.
    assert r["recommended_mode"] == "reference_transfer"


def test_no_gene_overlap_routes_to_clustering():
    from spatialscribe.analysis import reference as R

    r = R.plan_annotation_strategy(_ref_adata(disjoint=True), _PANEL, "cell_type", target_depth=50.0)
    assert r["status"] != "ok"                          # insufficient_overlap
    assert r["recommended_mode"] == "cluster"
    assert "clustering" in (r["reason"] or "").lower()


def test_plan_never_raises_on_degenerate_reference():
    from spatialscribe.analysis import reference as R

    r = R.plan_annotation_strategy(_ref_adata(n_types=1), _PANEL, "cell_type", target_depth=50.0)
    assert r["recommended_mode"] in {"reference_transfer", "cluster"}   # too_few_types -> cluster


def test_capability_no_reference_recommends_marker_annotation(processed_adata, ctx):
    import json

    from spatialscribe.analysis import capabilities as cap, keys

    res = cap.run(processed_adata, "annotation_strategy", {}, ctx)
    assert res.ok
    assert res.value["recommended_mode"] == "annotate"          # no reference -> marker supervised
    assert processed_adata.uns["annotation_route"]["status"] == "no_reference"
    assert keys.Uns.ANNOTATION_ROUTE in keys.all_keys()
    assert "annotation_strategy" in cap.copilot_names()
    json.dumps(res.value, default=str)                          # JSON-able


def test_capability_with_reference_file_routes(processed_adata, ctx, tmp_path):
    from spatialscribe.analysis import capabilities as cap

    processed_adata.var_names = [f"g{i}" for i in range(processed_adata.n_vars)]   # align to the ref
    rp = tmp_path / "ref.h5ad"
    _ref_adata().write_h5ad(rp)
    res = cap.run(processed_adata, "annotation_strategy", {"reference_path": str(rp)}, ctx)
    assert res.ok
    assert res.value["recommended_mode"] in {"reference_transfer", "cluster"}
    assert processed_adata.uns["annotation_route"]["ladder"][0]["step"] == "initial"


def test_plan_bad_label_key_never_raises():
    # Contract: plan_annotation_strategy is pure and never raises - a label_key that is not a
    # reference column (a common LLM mistake) must degrade, not KeyError out of the ladder.
    from spatialscribe.analysis import reference as R

    r = R.plan_annotation_strategy(_ref_adata(), _PANEL, "NOT_A_COLUMN", target_depth=50.0)
    assert r["recommended_mode"] in {"reference_transfer", "annotate", "cluster"}
    assert r["ladder"][0]["status"] in {"error", "insufficient_overlap", "too_few_types", "ok"}


def test_capability_bad_label_key_auto_detects_and_writes_route(processed_adata):
    # A wrong ctx.ref_label_key / label_key param must NOT crash the capability or skip its produces
    # key: resolve_reference validates the key against the in-memory reference and auto-detects.
    from spatialscribe.analysis import capabilities as cap, keys

    ref = _ref_adata()
    ref.obs = ref.obs.rename(columns={"cell_type": "my_label"})   # non-standard label column name
    processed_adata.var_names = [f"g{i}" for i in range(processed_adata.n_vars)]
    ctx = cap.RunContext(tissue="melanoma", use_llm=False, reference=ref, ref_label_key="WRONG_COL")
    res = cap.run(processed_adata, "annotation_strategy", {"label_key": "ALSO_WRONG"}, ctx)
    assert res.ok                                                 # auto-detected 'my_label'
    assert keys.Uns.ANNOTATION_ROUTE.present(processed_adata)     # produces contract held
