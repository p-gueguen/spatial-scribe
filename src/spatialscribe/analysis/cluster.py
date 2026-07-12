"""Normalize + dimensionality reduction + Leiden clustering (GPU/CPU via backend).

What it does
------------
`preprocess` normalizes/log-transforms and runs HVG + PCA + neighbors on whatever object it
is given (no cell dropping). `cluster` flags near-empty cells, embeds + Leiden-clusters the
good cells only, and scatters the results back onto the full object so low-signal cells stay
present (for QC and hard-abstention) yet never pollute the embedding.

Depends on
----------
scanpy (+ rapids_singlecell on GPU), via ``backend``. Reads/writes ``adata`` in place.
"""

from __future__ import annotations


def preprocess(adata, n_neighbors: int = 15, n_pcs: int = 50, progress=None) -> None:
    """Store counts, normalize+log, HVG, PCA, neighbor graph (in place, idempotent).

    Does NOT drop cells - excluding low-signal cells from the embedding is :func:`cluster`'s job.
    ``progress(frac, label)`` (optional) reports coarse checkpoints for the app's progress bar.
    """
    from .backend import get_backend

    be = get_backend()
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    if "log1p" not in adata.uns:  # normalize only once
        if progress:
            progress(0.05, "normalizing")
        be.normalize_log(adata)
        adata.uns["log1p"] = {"base": None}
    # Drop genes with zero expression before PCA. A targeted panel skips HVG (n_vars <= 2000), so PCA
    # sees EVERY panel gene, and the GPU (rapids) PCA REJECTS all-zero genes ("There are genes with zero
    # expression. Please remove them before running PCA") - the observed "clustering does nothing / empty
    # leiden" on a Xenium panel carrying genes not expressed in this tissue. A 0-everywhere gene is
    # neither a principal component nor a marker, so nothing is lost. On the usual path preprocess runs
    # on a `good` COPY (adata[~low]), so the caller's full object keeps every panel gene.
    be.filter_genes(adata, min_cells=1)
    if progress:
        progress(0.20, "selecting variable genes")
    be.hvg(adata)
    if progress:
        progress(0.35, "PCA")
    be.pca(adata, n_comps=n_pcs)
    if progress:
        progress(0.50, "building the neighbour graph")
    be.neighbors(adata, n_neighbors=n_neighbors)


def _low_signal_mask(adata):
    """Boolean mask of near-empty cells to EXCLUDE FROM THE EMBEDDING (never from the object).

    Near-empty cells have no expression structure: they collapse to near-identical vectors, so
    UMAP scatters them into a spurious uniform "disc" and Leiden spawns junk clusters from them
    (the breast 5K section carries ~2k zero-count + ~17k <10-count cells whose distributional
    count-floor is 0, so nothing else removes them). We KEEP them in the object - annotate's
    low-signal gate hard-abstains them as ``Unassigned:low_quality`` and QC still counts them -
    and only hold them out of PCA/neighbors/Leiden/UMAP. The threshold is panel-aware:
    ``< qc.RICH_PANEL_MIN_GENES`` detected genes on a rich panel (>= 1000 genes) - an ABSOLUTE floor
    (a fraction-of-median one wrongly cut deep, fine sections), shared with annotate's abstain gate so
    a cell held out of the embedding is also typed as low_quality - and only truly-empty (0-gene) cells
    on a small targeted panel where a flat floor would gut legitimate low-plex cells.
    """
    import numpy as np

    from . import qc as _qc

    raw = adata.layers["counts"] if "counts" in adata.layers else adata.X
    genes = np.asarray((raw > 0).sum(1)).ravel()
    n_panel = _qc.resolve_panel_size(adata)
    thr = _qc.RICH_PANEL_MIN_GENES if n_panel >= 1000 else 1
    return genes < thr


def _leiden_umap_markers(adata, resolution: float, key: str, progress=None) -> None:
    """Leiden + size relabelling + UMAP + numeric ordering + t-test marker genes (in place)."""
    from .backend import get_backend

    be = get_backend()
    if progress:
        progress(0.70, "Leiden communities")
    be.leiden(adata, resolution=resolution, key_added=key)
    # Seurat convention: cluster 0 is the largest. Relabel BEFORE rank_genes_groups, or the ranking
    # table would be keyed by the pre-relabel names.
    relabel_clusters_by_size(adata, key)
    if progress:
        progress(0.85, "UMAP embedding")
    be.umap(adata)
    order_clusters_numeric(adata, key)  # numeric (not lexical) cluster order, before rank_genes
    if progress:
        progress(0.95, "ranking cluster markers")
    # Marker genes per cluster. t-test (not wilcoxon): wilcoxon is a per-cell rank test that runs for
    # minutes on a big CPU section and blocks the single-process backend (niche/programs queue behind
    # it); t-test is ~10-50x faster with equivalent top-marker ordering for naming/subtyping. Routed
    # through the backend so it is the GPU t-test kernel on rsc 0.14, else scanpy on CPU.
    be.rank_genes_groups(adata, groupby=key, method="t-test")


def _rekey_rank_genes(adata, mapping: dict, key: str = "rank_genes_groups") -> None:
    """Rename the per-group columns of ``uns['rank_genes_groups']`` through ``mapping`` (in place).

    ``rank_genes_groups`` is a dict of structured arrays whose FIELD names are the cluster labels
    (``names``, ``scores``, ``pvals``, ...). After a size relabel the labels changed, so the ranking
    must be re-keyed or it points at the wrong clusters. Renaming fields keeps the computed ranking -
    no recompute - because the field that WAS old-label ``c`` still holds ``c``'s genes, and ``c`` is
    exactly the cluster now called ``mapping[c]``. No-op when the table is absent."""
    import numpy as np

    rg = adata.uns.get(key)
    if not isinstance(rg, dict):
        return
    for field, arr in list(rg.items()):
        names = getattr(getattr(arr, "dtype", None), "names", None)
        if not names:
            continue
        # Rebuild column-by-column: an in-place `arr.dtype = ...` refuses on object columns (the gene
        # names are Python-string references). Ordered numerically so downstream iteration is tidy.
        new_dtype = sorted(((mapping.get(n, n), arr.dtype[n]) for n in names),
                           key=lambda t: (not t[0].isdigit(), int(t[0]) if t[0].isdigit() else t[0]))
        out = np.empty(arr.shape, dtype=new_dtype)
        for n in names:
            out[mapping.get(n, n)] = arr[n]
        rg[field] = out


def relabel_clusters_by_size(adata, key: str = "leiden") -> None:
    """Rename the clusters so ``0`` is the LARGEST, ``1`` the next, ... (in place, Seurat convention).

    Leiden's own label numbers carry no meaning - they fall out of the community-detection order and
    differ between the igraph (CPU) and cuGraph (GPU) backends for the same section. Sorting by cell
    count makes the number informative and stable, so "cluster 0" always names the dominant
    population and a legend truncated at the top still shows the clusters that matter.

    Ties are broken by the original label so the mapping is deterministic. Cells with a NaN label
    (``cluster.cluster`` holds low-signal cells out of the embedding) stay NaN. If a
    ``rank_genes_groups`` table is already present (e.g. a pre-clustered cache), it is re-keyed through
    the same mapping so the ranking still matches the labels - so this is safe to call BEFORE
    ``rank_genes_groups`` (in the cluster step) OR after it (on a cache load)."""
    import pandas as pd

    if key not in adata.obs:
        return
    s = adata.obs[key].astype(str)
    counts = s[s != "nan"].value_counts()
    ordered = sorted(counts.index, key=lambda c: (-int(counts[c]), str(c)))
    mapping = {old: str(i) for i, old in enumerate(ordered)}
    _rekey_rank_genes(adata, mapping)          # keep a pre-existing ranking consistent, no recompute
    new = s.map(mapping)                       # unmapped ("nan") -> NaN, which is what we want
    adata.obs[key] = pd.Categorical(new, categories=[str(i) for i in range(len(ordered))])


def order_clusters_numeric(adata, key: str = "leiden") -> None:
    """Re-order the cluster categorical numerically (0, 1, 2, ..., 10), not lexically (in place).

    Leiden labels are integer strings, so pandas' default lexical category order sorts
    ``'10'`` before ``'2'``. This makes ``'2'`` precede ``'10'`` everywhere the *category
    order* is read - the ``rank_genes_groups`` group order, plot legends, the
    neighborhood-enrichment axes, and the annotation group iteration - all of which take
    ``adata.obs[key].cat.categories`` at face value. Non-integer labels fall back to a
    natural sort. No-op if ``key`` is absent.
    """
    import pandas as pd

    if key not in adata.obs:
        return
    s = adata.obs[key].astype(str)
    uniq = s.unique().tolist()
    try:
        ordered = sorted(uniq, key=int)
    except ValueError:  # not all-integer labels -> natural sort
        from natsort import natsorted
        ordered = natsorted(uniq)
    adata.obs[key] = pd.Categorical(s, categories=ordered)


def cluster(adata, resolution: float = 0.5, key: str = "leiden", progress=None) -> None:
    """Leiden clustering + per-cluster marker genes (in place), excluding low-signal cells.

    The default ``resolution`` is deliberately low (0.5) so the first pass yields a few coarse,
    interpretable clusters rather than an over-split map; raise it to subdivide. Near-empty cells
    are flagged (``obs['low_signal']``) and held OUT of the embedding so they don't form a UMAP
    disc or spawn junk clusters, but they remain in the object: their ``leiden`` is left NaN
    (excluded from every ``cat.categories`` / ``groupby`` / ``nunique`` consumer, so annotate and
    the cluster count skip them cleanly), their ``X_umap`` is NaN (the UMAP view drops them), and
    their ``X_pca`` is 0 (finite, so cluster-confidence distance math never sees NaN). Marker
    genes land in ``adata.uns['rank_genes_groups']`` for annotation/subtyping downstream.
    """
    import numpy as np
    import pandas as pd

    from .backend import get_backend

    be = get_backend()
    if progress:
        progress(0.05, "normalizing")
    # Normalize the FULL object once (annotate marker scoring + NMF read adata.X for every cell).
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    if "log1p" not in adata.uns:
        be.normalize_log(adata)
        adata.uns["log1p"] = {"base": None}

    low = _low_signal_mask(adata)
    if low.all():  # degenerate section: nothing would remain - keep everything rather than blank out
        low = np.zeros(adata.n_obs, dtype=bool)
    adata.obs["low_signal"] = low
    adata.uns["low_signal_filter"] = {"n_excluded": int(low.sum()), "n_kept": int((~low).sum())}

    if not low.any():  # no low-signal cells: embed the whole object in place
        if "neighbors" not in adata.uns:
            preprocess(adata, progress=progress)
        _leiden_umap_markers(adata, resolution, key, progress=progress)
        if progress:
            progress(1.0, "done")
        return

    # Embed + cluster the good cells only, then scatter results back to the full object. `good`
    # inherits counts + the log1p flag from the parent, so preprocess() only (re)builds HVG/PCA/
    # neighbors on the subset - it does not re-normalize an already-normalized X.
    good = adata[~low].copy()
    preprocess(good, progress=progress)
    _leiden_umap_markers(good, resolution, key, progress=progress)

    gi = np.where(~low)[0]
    # leiden: real labels for good cells; low-signal cells stay NaN (outside cat.categories).
    lab = pd.Series(pd.NA, index=adata.obs_names, dtype="object")
    lab.iloc[gi] = good.obs[key].astype(str).to_numpy()
    adata.obs[key] = pd.Categorical(lab, categories=list(good.obs[key].cat.categories))
    # X_umap: real coords for good, NaN for low (points() drops non-finite rows on the UMAP basis).
    um = np.full((adata.n_obs, good.obsm["X_umap"].shape[1]), np.nan, dtype=float)
    um[gi] = np.asarray(good.obsm["X_umap"], dtype=float)
    adata.obsm["X_umap"] = um
    # X_pca: real coords for good, 0 for low (finite - cluster-confidence distances never NaN).
    pc = np.zeros((adata.n_obs, good.obsm["X_pca"].shape[1]), dtype=float)
    pc[gi] = np.asarray(good.obsm["X_pca"], dtype=float)
    adata.obsm["X_pca"] = pc
    adata.uns["rank_genes_groups"] = good.uns["rank_genes_groups"]
    if progress:
        progress(1.0, "done")


def top_markers(adata, group: str, n: int = 15, key: str = "rank_genes_groups") -> list[str]:
    """Return the top-``n`` marker gene symbols for one cluster label."""
    import numpy as np

    names = adata.uns[key]["names"]
    col = group if group in names.dtype.names else str(group)
    return list(np.asarray(names[col])[:n])
