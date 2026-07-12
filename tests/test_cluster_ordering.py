"""Leiden cluster ordering (numeric, not lexical) + the lowered default resolution."""
from __future__ import annotations

import inspect

import numpy as np
import pytest

anndata = pytest.importorskip("anndata")


def test_order_clusters_numeric_sorts_by_value():
    """Integer-string labels order 0,1,2,...,10 - not the lexical 0,1,10,11,2."""
    import pandas as pd

    from spatialscribe.analysis import cluster

    a = anndata.AnnData(X=np.zeros((6, 2), dtype="float32"))
    a.obs["leiden"] = pd.Categorical(["10", "2", "1", "0", "11", "2"])
    # sanity: pandas' default category order is lexical - the bug we are fixing
    assert list(a.obs["leiden"].cat.categories) == ["0", "1", "10", "11", "2"]

    cluster.order_clusters_numeric(a, "leiden")

    assert list(a.obs["leiden"].cat.categories) == ["0", "1", "2", "10", "11"]
    # only the category ORDER changes; the per-cell labels are untouched
    assert a.obs["leiden"].astype(str).tolist() == ["10", "2", "1", "0", "11", "2"]


def test_order_clusters_numeric_non_integer_falls_back_to_natsort():
    """Non-integer labels natural-sort instead of raising ValueError."""
    import pandas as pd

    from spatialscribe.analysis import cluster

    a = anndata.AnnData(X=np.zeros((3, 2), dtype="float32"))
    a.obs["leiden"] = pd.Categorical(["c10", "c2", "c1"])
    cluster.order_clusters_numeric(a, "leiden")
    assert list(a.obs["leiden"].cat.categories) == ["c1", "c2", "c10"]


def test_order_clusters_numeric_noop_without_key():
    """No cluster column -> no error, no change (defensive no-op)."""
    from spatialscribe.analysis import cluster

    a = anndata.AnnData(X=np.zeros((2, 2), dtype="float32"))
    cluster.order_clusters_numeric(a, "leiden")  # must not raise
    assert "leiden" not in a.obs


def test_ordered_group_labels_honors_categorical_order():
    """The plot legend/trace helper follows a Categorical's numeric order, not lexical."""
    import pandas as pd

    from spatialscribe.analysis.plots import _ordered_group_labels

    col = pd.Series(pd.Categorical(["2", "10", "0"], categories=["0", "1", "2", "10"]))
    # '1' has no cells present -> skipped, remaining order stays numeric (not 0,10,2)
    assert _ordered_group_labels(col, {"0", "2", "10"}) == ["0", "2", "10"]
    # a present label outside the declared categories is appended, never dropped
    assert _ordered_group_labels(col, {"0", "10", "cell"}) == ["0", "10", "cell"]
    # a plain (non-categorical) column falls back to a lexical sort (prior behavior)
    assert _ordered_group_labels(pd.Series(["b", "a"], dtype=object), {"b", "a"}) == ["a", "b"]


def test_default_resolution_lowered_across_surfaces():
    """Every default-resolution surface is the lowered 0.5 - no lingering 1.0."""
    from spatialscribe.analysis import backend, capabilities, cluster

    assert inspect.signature(cluster.cluster).parameters["resolution"].default == 0.5
    # the wizard-button / copilot adapter default
    assert inspect.signature(capabilities._cap_cluster).parameters["resolution"].default == 0.5
    # whichever backend is active (CPU in the test env)
    assert inspect.signature(backend.get_backend().leiden).parameters["resolution"].default == 0.5


def test_relabel_clusters_by_size_makes_zero_the_largest():
    """Seurat convention: after relabel, cluster 0 is the biggest, and NaN low-signal cells stay NaN."""
    import pandas as pd

    from spatialscribe.analysis import cluster

    a = anndata.AnnData(X=np.zeros((10, 2), dtype="float32"))
    # sizes: "7"->1, "3"->4, "9"->3, plus 2 low-signal NaN cells
    a.obs["leiden"] = pd.Categorical(
        ["3", "3", "3", "3", "9", "9", "9", "7", pd.NA, pd.NA])
    cluster.relabel_clusters_by_size(a, "leiden")

    vc = a.obs["leiden"].value_counts()
    assert int(vc["0"]) == 4 and int(vc["1"]) == 3 and int(vc["2"]) == 1   # 3->0, 9->1, 7->2
    assert list(a.obs["leiden"].cat.categories) == ["0", "1", "2"]
    assert int(a.obs["leiden"].isna().sum()) == 2                          # NaN preserved


def test_relabel_rekeys_rank_genes_without_recompute():
    """A pre-existing rank_genes_groups is re-keyed through the same mapping, so cluster 0's genes are
    the OLD largest cluster's genes - no recompute, and never a mislabelled ranking."""
    import numpy as np
    import pandas as pd

    from spatialscribe.analysis import cluster

    a = anndata.AnnData(X=np.zeros((7, 3), dtype="float32"))
    a.obs["leiden"] = pd.Categorical(["5", "5", "5", "5", "2", "2", "8"])   # sizes 5->4, 2->2, 8->1
    # ranking keyed by the OLD labels; each cluster's marker is a distinct gene
    names = np.array([("g5", "g2", "g8")], dtype=[("5", "O"), ("2", "O"), ("8", "O")])
    a.uns["rank_genes_groups"] = {"names": names}

    cluster.relabel_clusters_by_size(a, "leiden")

    rg = a.uns["rank_genes_groups"]["names"]
    assert list(rg.dtype.names) == ["0", "1", "2"]        # numeric-ordered
    assert rg["0"][0] == "g5"                             # new 0 = old largest (5) -> its gene g5
    assert rg["1"][0] == "g2" and rg["2"][0] == "g8"


def test_cluster_markers_reuses_the_ranking(monkeypatch):
    """cluster_markers is a READ when rank_genes_groups already matches the clusters - it must NOT call
    the backend ranking again (that is the whole 'make it fast' point)."""
    import numpy as np
    import pandas as pd

    from spatialscribe.analysis import backend, capabilities as cap

    a = anndata.AnnData(X=np.random.default_rng(0).poisson(1.0, (12, 4)).astype("float32"))
    a.var_names = [f"g{i}" for i in range(4)]
    a.obs["leiden"] = pd.Categorical(["0"] * 6 + ["1"] * 6)
    a.uns["rank_genes_groups"] = {
        "names": np.array([("g0", "g1"), ("g2", "g3")], dtype=[("0", "O"), ("1", "O")])}

    called = {"n": 0}
    orig = backend.get_backend().rank_genes_groups
    def spy(self, *a_, **k_):
        called["n"] += 1
        return orig(*a_, **k_)
    monkeypatch.setattr(type(backend.get_backend()), "rank_genes_groups", spy)

    out = cap.run(a, "cluster_markers", {"n_genes": 2}, cap.RunContext(tissue="x", use_llm=False))
    assert out.ok and out.value["source"].startswith("reused")
    assert called["n"] == 0                               # the fast path did not recompute
