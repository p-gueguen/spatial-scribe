"""merge_types: APPLY the panel-driven merges (collapse types the panel cannot separate into one
coarser label) and report the annotation-quality delta - the executor half of the merge advice that
self_verify / panel_check only surface. Also its fold into self_heal (merge_confusable)."""
from __future__ import annotations

import numpy as np
import pytest

anndata = pytest.importorskip("anndata")
pytest.importorskip("scanpy")

# A and B are indistinguishable (identical markers, no private gene -> a panel gap); C is private.
_MS = {"A": ["g0", "g1", "g2"], "B": ["g0", "g1", "g2"], "C": ["g10", "g11", "g12"]}


def _confusable_section(n_per=60, seed=0):
    import pandas as pd

    rng = np.random.default_rng(seed)
    n = n_per * 3
    # continuous per-gene baseline means so genes spread across expression bins (else scanpy
    # score_genes finds no control genes when a boosted marker is a lone outlier in its bin).
    means = rng.uniform(0.5, 12.0, size=60)
    X = rng.poisson(np.tile(means, (n, 1))).astype("float32")
    types = np.array(["A"] * n_per + ["B"] * n_per + ["C"] * n_per)
    ab = np.isin(types, ["A", "B"])
    c = types == "C"
    X[ab, 0:3] += rng.poisson(6.0, size=(int(ab.sum()), 3)).astype("float32")    # shared A/B program
    X[c, 10:13] += rng.poisson(6.0, size=(int(c.sum()), 3)).astype("float32")    # C's private program
    a = anndata.AnnData(X=X)
    a.obs_names = [f"c{i}" for i in range(n)]
    a.var_names = [f"g{i}" for i in range(60)]
    a.obsm["spatial"] = rng.normal(size=(n, 2)) * 100
    a.obs["cell_type"] = pd.Categorical(types)
    return a


def _multi_lineage_section(seed=0):
    """T cell + NK cell (both immune) + Endothelial (vascular), each with its own on-panel markers."""
    import pandas as pd
    import scanpy as sc

    from spatialscribe.analysis import markers as mk

    rng = np.random.default_rng(seed)
    lin = {k: mk.LINEAGE_MARKERS[k] for k in ("T cell", "NK cell", "Endothelial")}
    genes = list(dict.fromkeys(sum(lin.values(), []))) + [f"BG{i}" for i in range(70)]
    gi = {g: i for i, g in enumerate(genes)}
    nbg = 70

    def block(n, name):
        M = rng.poisson(0.3, size=(n, len(genes))).astype("float32")
        M[:, len(genes) - nbg:] += rng.poisson(rng.uniform(0.5, 15, nbg), size=(n, nbg)).astype("float32")
        for g in lin[name]:
            if g in gi:
                M[:, gi[g]] += rng.poisson(25, size=n)
        return M

    blocks, labs = [], []
    for name in lin:
        blocks.append(block(60, name)); labs += [name] * 60
    a = anndata.AnnData(X=np.vstack(blocks))
    a.var_names = genes
    a.obs["cell_type"] = pd.Categorical(labs)
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a); sc.pp.log1p(a)
    return a, lin


def test_merge_never_crosses_major_lineage():
    # The live "optimize merging" bug: panel_check blobs any two markerless types into one group, so
    # an immune type got merged into neural/epithelial/stromal (Oligodendrocyte + B cell). The lineage
    # guard must split a mixed group into lineage-coherent subgroups: the two immune types coarsen
    # together, but the vascular type is NEVER folded in. Holds offline via the compartment fallback.
    from spatialscribe.analysis import verify

    a, lin = _multi_lineage_section()
    res = verify.merge_confusable_types(a, [["T cell", "NK cell", "Endothelial"]],
                                        cluster_key="cell_type", marker_sets=lin)
    assert res["status"] == "ok"
    cats = set(a.obs["cell_type"].astype(str))
    assert "Endothelial" in cats                                    # the vascular type stays separate
    assert not any("Endothelial" in c and ("T cell" in c or "NK cell" in c) for c in cats)
    assert any("T cell" in c and "NK cell" in c for c in cats)      # the two immune types coarsen


def test_merge_confusable_types_applies_labels_and_improves_quality():
    from spatialscribe.analysis import verify

    a = _confusable_section()
    res = verify.merge_confusable_types(a, [["A", "B"]], cluster_key="cell_type", marker_sets=_MS)

    assert res["status"] == "ok" and res["n_groups_merged"] == 1
    assert res["merged"][0]["label"] == "A / B" and res["merged"][0]["n_cells"] == 120
    # the label actually changed: A and B collapsed into one category, C untouched
    cats = set(a.obs["cell_type"].astype(str))
    assert "A / B" in cats and "A" not in cats and "B" not in cats and "C" in cats
    # before, one of the indistinguishable pair fails marker-agreement; merging removes that failure
    assert res["quality_before"]["n_failed"] >= 1
    assert res["n_failed_delta"] <= 0
    # and the honest quality rate does not drop (the merged union program is scored, not dropped)
    assert res["agreement_delta"] is not None and res["agreement_delta"] >= 0


def test_merge_dry_run_previews_without_mutating():
    # The Curate "Suggest merges" button: preview the ontology-partitioned merges WITHOUT applying.
    from spatialscribe.analysis import capabilities as cap
    from spatialscribe.analysis import verify

    a = _confusable_section()
    before = set(a.obs["cell_type"].astype(str))
    prev = verify.merge_confusable_types(a, [["A", "B"]], cluster_key="cell_type",
                                         marker_sets=_MS, dry_run=True)
    assert prev["status"] == "preview"
    assert prev["would_merge"][0]["label"] == "A / B" and prev["would_merge"][0]["n_cells"] == 120
    assert set(a.obs["cell_type"].astype(str)) == before          # obs NOT mutated
    assert "annotation_verification" not in a.uns                 # the before/after verify was skipped

    # same through the capability (what the button calls) - preview, then apply
    ctx = cap.RunContext(tissue="gut", marker_sets=_MS)
    rp = cap.run(a, "merge_types", {"dry_run": True}, ctx)
    assert rp.ok and rp.value["status"] == "preview" and "A" in set(a.obs["cell_type"].astype(str))
    ra = cap.run(a, "merge_types", {}, ctx)
    assert ra.ok and ra.value["n_groups_merged"] == 1 and "A / B" in set(a.obs["cell_type"].astype(str))


def test_merge_types_capability_defaults_to_panel_check_groups():
    from spatialscribe.analysis import capabilities as cap

    a = _confusable_section()
    ctx = cap.RunContext(tissue="gut", marker_sets=_MS)
    # no explicit groups -> derive them from panel_check (A+B have no private marker)
    r = cap.run(a, "merge_types", {}, ctx)
    assert r.ok, r.error
    assert r.value["n_groups_merged"] == 1 and "A / B" in set(a.obs["cell_type"].astype(str))

    # a panel with no confusable groups -> honest no-op, never raises
    a2 = _confusable_section()
    ctx2 = cap.RunContext(tissue="gut", marker_sets={"C": ["g10", "g11", "g12"]})
    r2 = cap.run(a2, "merge_types", {}, ctx2)
    assert r2.ok and r2.value["status"] == "no_merges" and r2.value["n_groups_merged"] == 0


def test_self_heal_merge_confusable_applies_the_merge_and_reports_it():
    from spatialscribe.analysis import capabilities as cap

    # default self_heal leaves the merge advisory (merge is None); opting in applies it. Fresh
    # section each time - self_heal mutates labels (abstain/subcluster), so never share state.
    r0 = cap.run(_confusable_section(), "self_heal", {}, cap.RunContext(tissue="gut", marker_sets=_MS))
    assert r0.ok and r0.value.get("merge") is None

    a = _confusable_section()
    r = cap.run(a, "self_heal", {"merge_confusable": True}, cap.RunContext(tissue="gut", marker_sets=_MS))
    assert r.ok, r.error
    assert r.value["merge"] is not None and r.value["merge"]["n_groups_merged"] >= 1
    assert "A / B" in set(a.obs["cell_type"].astype(str))
