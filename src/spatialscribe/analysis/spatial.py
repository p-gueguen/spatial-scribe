"""Spatial statistics (squidpy) - neighborhoods, immune exclusion, spatial genes.

What it does
------------
`spatial_neighbors` builds the spatial graph (Delaunay, KNN fallback).
`nhood_enrichment` computes cell-type neighborhood enrichment with min-cell guards.
`immune_exclusion` reads that enrichment to answer "are T cells excluded from the
tumor?" - the copilot hero moment - returning the z-score and a plain verdict.
`spatially_variable_genes` ranks Moran's I.

Depends on
----------
squidpy, numpy. Operates on ``adata`` with ``obsm['spatial']`` and a categorical
cell-type column.

Guards (from the repo's fS_nhood_enrichment_2_squidpy.py): skip cell types with
< ``MIN_CELLS_PER_CELLTYPE`` cells to avoid unstable enrichment on sparse groups.
"""

from __future__ import annotations

MIN_CELLS_PER_CELLTYPE = 20


def spatial_neighbors(adata, method: str = "delaunay", n_neighs: int = 6) -> None:
    """Build the spatial neighbor graph in place (Delaunay by default). Routed through the backend
    (squidpy on both legs - rsc 0.14 has no GPU spatial_neighbors; it runs on obsm['spatial'])."""
    from .backend import get_backend

    be = get_backend()
    if method == "delaunay":
        be.spatial_neighbors(adata, coord_type="generic", delaunay=True)
    else:
        be.spatial_neighbors(adata, coord_type="generic", n_neighs=n_neighs)


def _valid_groups(adata, cluster_key: str):
    vc = adata.obs[cluster_key].value_counts()
    keep = vc[vc >= MIN_CELLS_PER_CELLTYPE].index.tolist()
    return keep, [g for g in vc.index if g not in keep]


def nhood_enrichment(adata, cluster_key: str = "cell_type", seed: int = 0) -> dict:
    """Cell-type neighborhood enrichment (z-scores). Returns a JSON-able summary.

    Cell types below ``MIN_CELLS_PER_CELLTYPE`` are FLAGGED in ``dropped_sparse_types`` (they are
    still scored - squidpy runs over every category), so a z-score resting on a handful of cells can
    be discounted rather than silently trusted.

    Abstained cells KEEP their place in the neighbour graph (they occupy real tissue, and deleting
    them would silently rewire every neighbourhood) but collapse into one ``Not assigned`` class that
    is not reported: an abstention is not a lineage to be enriched near. They therefore remain in the
    permutation null, which ``pct_abstained_in_graph`` discloses.
    """
    import numpy as np

    from . import annotate as _an
    from .backend import get_backend

    if "spatial_connectivities" not in adata.obsp:
        spatial_neighbors(adata)

    abstained = _an.abstention_mask(adata.obs[cluster_key].astype(str).to_numpy())
    key = cluster_key
    if abstained.any():
        key = "_nhood_collapsed"
        adata.obs[key] = _an.collapse_abstention(adata.obs[cluster_key]).to_numpy()
        adata.obs[key] = adata.obs[key].astype("category")

    # squidpy raises `Expected at least 2 clusters, found 1` when fewer than two REAL cell types
    # survive - reachable on a barren panel, where apply_confidence stamps every cell
    # 'Unresolvable: panel'. Return an empty, honest result instead of a traceback.
    real = [c for c in map(str, adata.obs[key].astype("category").cat.categories)
            if not _an.is_abstention(c)]
    if len(real) < 2:
        if key != cluster_key:
            del adata.obs[key]
        return {"cluster_key": cluster_key, "categories": real, "zscore": [],
                "dropped_sparse_types": [], "pct_abstained_in_graph": float(abstained.mean()),
                "note": "fewer than 2 non-abstained cell types; neighborhood enrichment is undefined"}
    try:
        _, dropped = _valid_groups(adata, key)
        # 1000 permutations (squidpy default) is ~6 s on 100k cells; 200 gives stable z-scores much
        # faster (the enrichment/immune-exclusion buttons reuse this).
        n_perms = 200 if adata.n_obs > 20_000 else 1000
        get_backend().nhood_enrichment(adata, cluster_key=key, seed=seed, n_perms=n_perms,
                                       show_progress_bar=False)
        raw = adata.uns.pop(f"{key}_nhood_enrichment")
        z = np.asarray(raw["zscore"])
        cats = list(map(str, adata.obs[key].cat.categories))
    finally:
        if key != cluster_key and key in adata.obs:
            del adata.obs[key]

    keep = [i for i, c in enumerate(cats) if not _an.is_abstention(c)]
    z = z[np.ix_(keep, keep)] if keep else z[:0, :0]
    cats = [cats[i] for i in keep]
    # Republish under the CALLER's key (keys.nhood_enrichment_key), never the temp one, and carry the
    # category names so the matrix can never be re-indexed against a differently-ordered obs column.
    out = {"zscore": z, "categories": cats}
    if "count" in raw:
        cnt = np.asarray(raw["count"])
        out["count"] = cnt[np.ix_(keep, keep)] if keep else cnt[:0, :0]
    adata.uns[f"{cluster_key}_nhood_enrichment"] = out

    return {
        "cluster_key": cluster_key,
        "categories": cats,
        "zscore": z.tolist(),
        # NOTE: these are FLAGGED as low-n, not removed - squidpy scores every category. Reported so a
        # z-score backed by a handful of cells can be discounted.
        "dropped_sparse_types": [d for d in dropped if not _an.is_abstention(str(d))],
        # Abstained cells remain in the graph AND in the permutation null (as one inert class) even
        # though their row/col is not reported. Disclose it: a z-score reads differently when a fifth
        # of the neighbourhood was uncallable.
        "pct_abstained_in_graph": float(abstained.mean()),
    }


def immune_exclusion(adata, tumor_label: str, tcell_label: str,
                     cluster_key: str = "cell_type") -> dict:
    """Is ``tcell_label`` excluded from ``tumor_label``? (the hero readout).

    Uses the neighborhood-enrichment z-score between the two labels: strongly negative
    => spatially segregated (immune-excluded); positive => infiltrated/co-localized.
    """
    res = nhood_enrichment(adata, cluster_key=cluster_key)
    cats = res["categories"]
    # Resolve the labels tolerantly (the copilot passes "T cells"/"tumor", not the exact category):
    # an exact-match here would error on the near-miss and the copilot loop would silently recover
    # into a different tool. See analysis.labels.match_cell_type.
    from .labels import match_cell_type
    a, b = match_cell_type(tumor_label, cats), match_cell_type(tcell_label, cats)
    if a is None or b is None:
        missing = [lbl for lbl, m in ((tumor_label, a), (tcell_label, b)) if m is None]
        return {"error": f"cell type(s) {missing} not found in this section. Available: {cats}"}
    tumor_label, tcell_label = a, b
    i, j = cats.index(tumor_label), cats.index(tcell_label)
    zval = res["zscore"][i][j]
    verdict = (
        "immune-excluded (T cells segregated from the tumor)" if zval < -2
        else "infiltrated (T cells co-localize with the tumor)" if zval > 2
        else "no strong spatial preference"
    )
    return {"tumor": tumor_label, "tcell": tcell_label, "zscore": round(float(zval), 2),
            "verdict": verdict}


def co_occurrence(adata, cluster_key: str = "cell_type") -> dict:
    """Cell-type co-occurrence vs distance (squidpy). Returns a JSON-able summary."""
    import numpy as np

    from .backend import get_backend

    get_backend().co_occurrence(adata, cluster_key=cluster_key, show_progress_bar=False)
    occ = adata.uns[f"{cluster_key}_co_occurrence"]
    return {
        "cluster_key": cluster_key,
        "categories": list(adata.obs[cluster_key].cat.categories),
        "occ": np.asarray(occ["occ"]).tolist(),
        "interval": np.asarray(occ["interval"]).tolist(),
    }


def spatially_variable_genes(adata, n_top: int = 20) -> list[str]:
    """Top spatially variable genes by Moran's I."""
    from .backend import get_backend

    if "spatial_connectivities" not in adata.obsp:
        spatial_neighbors(adata)
    get_backend().spatial_autocorr(adata, mode="moran", show_progress_bar=False)
    moran = adata.uns["moranI"].sort_values("I", ascending=False)
    return moran.head(n_top).index.tolist()


def spatial_coherence(adata, label_key: str = "cell_type", k: int = 15,
                      pas_threshold: float = 0.2) -> dict:
    """Layer 6: per-cell same-label spatial-neighbor fraction + dataset PAS (writes
    ``obs['spatial_coherence']``). Thin shim over ``spatial_anno_metrics.spatial.spatial_coherence``
    (the single source of truth); builds the graph with SpatialScribe's ``spatial_neighbors`` first to
    preserve the pipeline's kNN choice. DOWN-WEIGHT ONLY - rare infiltrating cells are legitimately
    incoherent, so this never drops cells."""
    from spatial_anno_metrics import spatial as _sp
    if "spatial_connectivities" not in adata.obsp:
        spatial_neighbors(adata, method="knn", n_neighs=k)
    return _sp.spatial_coherence(adata, label_key=label_key, k=k, pas_threshold=pas_threshold)


def neighborhood_sanity(adata, cluster_key: str = "cell_type") -> dict:
    """Layer 6 dataset-level sanity: flags any cell type whose strongest neighborhood-enrichment partner
    is a DIFFERENT type (a systematic-mislabeling smell). Thin shim over
    ``spatial_anno_metrics.spatial.neighborhood_sanity``; returns the self-enrichment diagonal + flags."""
    from spatial_anno_metrics import spatial as _sp
    return _sp.neighborhood_sanity(adata, label_key=cluster_key)
