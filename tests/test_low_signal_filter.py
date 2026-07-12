"""Low-signal cells are KEPT in the object but held OUT of the embedding (the UMAP-disc fix).

Regression cover for the cluster.preprocess/cluster rework: near-empty cells must not be dropped
(that corrupted the annotatability/QC reporting the app exists to produce), but must be excluded
from PCA/neighbors/Leiden/UMAP so they don't form a phantom UMAP disc or spawn junk clusters.
"""
from __future__ import annotations

import numpy as np
import pytest

anndata = pytest.importorskip("anndata")
sparse = pytest.importorskip("scipy.sparse")


def _demo_with_low(n: int = 260, seed: int = 0):
    """A synthetic section with 20 truly-empty (0-gene) + 10 three-gene cells at the tail."""
    import scipy.sparse as sp

    from spatialscribe.analysis import demo
    from spatialscribe.analysis.io import build_control_mask

    a = demo.make_demo_adata(n_cells=n, seed=seed)
    a.var["control"] = build_control_mask(a)
    a.uns["n_panel_genes"] = 2000  # rich-panel regime -> low-signal threshold is qc.RICH_PANEL_MIN_GENES
    a.layers.pop("counts", None)   # let cluster() re-derive counts from the (edited) X
    X = (a.X.toarray() if sp.issparse(a.X) else np.asarray(a.X)).astype("float32")
    X[-30:, :] = 0.0        # last 30 cells wiped...
    X[-30:-20, :3] = 1.0    # ...but 10 of them keep 3 genes (still below the rich-panel floor -> low)
    a.X = sp.csr_matrix(X)
    return a, 30            # 30 low-signal cells (20 empty + 10 three-gene)


def test_low_signal_cells_kept_but_excluded_from_embedding():
    from spatialscribe.analysis import cluster

    a, n_low = _demo_with_low()
    n0 = a.n_obs
    cluster.cluster(a, resolution=1.0)

    assert a.n_obs == n0                                    # nothing is dropped
    low = a.obs["low_signal"].to_numpy().astype(bool)
    assert int(low.sum()) == n_low
    # leiden is NaN exactly for the low-signal cells (so groupby/nunique/annotate skip them)
    leiden_nan = a.obs["leiden"].isna().to_numpy()
    assert np.array_equal(leiden_nan, low)
    # the cluster count counts only real clusters (NaN excluded)
    assert a.obs["leiden"].nunique() == len(a.obs["leiden"].cat.categories) > 0
    # X_umap is NaN exactly for low-signal cells; X_pca is finite everywhere (0 for low)
    um_nan = ~np.isfinite(np.asarray(a.obsm["X_umap"])).all(1)
    assert np.array_equal(um_nan, low)
    assert np.isfinite(np.asarray(a.obsm["X_pca"])).all()
    assert a.uns["low_signal_filter"]["n_excluded"] == n_low


def test_panel_aware_threshold_keeps_low_gene_cells_on_small_panels():
    """On a small targeted panel (<1000 genes) only truly-empty cells are excluded, not 3-gene ones."""
    from spatialscribe.analysis import cluster

    a, _ = _demo_with_low()
    a.uns["n_panel_genes"] = 200
    cluster.cluster(a, resolution=1.0)
    assert int(a.obs["low_signal"].sum()) == 20   # only the 20 zero-gene cells


def test_rich_panel_floor_is_absolute_not_a_fraction_of_median():
    """The rich-panel floor is an ABSOLUTE gene count, so a DEEP section keeps its many-gene cells
    (a fraction-of-median floor wrongly dropped 9% of the atera5k positive control), while a 10-gene
    cell - kept by the old flat <5 floor - is now held out (it seeded 1-gene artefact micro-clusters)."""
    import scipy.sparse as sp

    from spatialscribe.analysis import cluster, qc
    from spatialscribe.analysis.demo import make_demo_adata
    from spatialscribe.analysis.io import build_control_mask

    assert qc.RICH_PANEL_MIN_GENES == 15

    a = make_demo_adata(n_cells=200, seed=1)
    a.var["control"] = build_control_mask(a)
    a.uns["n_panel_genes"] = 5101                       # rich panel
    a.layers.pop("counts", None)
    X = (a.X.toarray() if sp.issparse(a.X) else np.asarray(a.X)).astype("float32")
    X[-10:, :] = 0.0
    X[-10:, :10] = 1.0                                  # 10 cells with exactly 10 genes (old floor kept them)
    a.X = sp.csr_matrix(X)
    cluster.cluster(a, resolution=1.0)
    low = a.obs["low_signal"].to_numpy().astype(bool)
    assert low[-10:].all()                              # 10-gene cells now excluded (10 < 15)
    # the demo's normal cells are well above 15 genes -> an absolute floor does NOT sweep them in
    assert not low[:-10].any()


def test_recluster_is_idempotent_on_n_obs():
    from spatialscribe.analysis import cluster

    a, _ = _demo_with_low()
    cluster.cluster(a, resolution=0.5)
    n1 = a.n_obs
    cluster.cluster(a, resolution=1.5)            # re-cluster at a new resolution must not drop cells
    assert a.n_obs == n1


def test_subcluster_without_counts_layer_labels_every_masked_cell():
    """Reference-transfer path: cell_type present, no counts layer, a low-gene cell in the target."""
    from spatialscribe.analysis import subcluster

    a, _ = _demo_with_low(n=140)
    a.layers.pop("counts", None)
    # give every cell a single cell_type so the whole section is the subcluster target
    a.obs["cell_type"] = "T cell"
    sub, rows = subcluster.subcluster(a, "T cell", resolution=0.5, use_llm=False)
    assert rows and "subtype" in a.obs
    # every originally-selected cell (all of them here) got a subtype - none left null
    assert a.obs["subtype"].notna().all()
    assert a.n_obs == 140
