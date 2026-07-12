"""Click-to-subcluster + marker+LLM subtyping (feature H2).

What it does
------------
Take one annotated cell type (e.g. "T cell"), re-cluster just those cells, compute
per-subcluster markers, and let Claude name the subtypes (CD4/CD8/Treg/exhausted...).
Writes the refined labels back to ``adata.obs['subtype']`` for the selected cells and
returns a per-subcluster table for the UI.

Depends on: scanpy (+ backend), :mod:`cluster`, :mod:`llm`.
"""

from __future__ import annotations

from . import cluster as _cluster
# Tolerant label resolution now lives in .labels (shared with immune_exclusion et al.). Re-exported
# under the old names so callers and tests that use subcluster._match_cell_type keep working.
from .labels import match_cell_type as _match_cell_type  # noqa: F401
from .labels import norm_label as _norm_label  # noqa: F401


def subcluster(adata, cell_type: str, cluster_key: str = "cell_type",
               resolution: float = 0.1, use_llm: bool = True,
               context: str = "human skin / melanoma"):
    """Re-cluster cells of ``cell_type`` and name subtypes.

    Returns ``(sub_adata, table)`` where ``table`` is a list of
    ``{subcluster, n_cells, top_markers, label, confidence, rationale}`` and the parent
    ``adata.obs['subtype']`` is updated in place for the selected cells.
    """
    cats = list(dict.fromkeys(adata.obs[cluster_key].astype(str)))
    resolved = _match_cell_type(cell_type, cats)
    if resolved is None:
        raise ValueError(f"'{cell_type}' is not a cell type in this section. "
                         f"Available cell types: {sorted(cats)}.")
    cell_type = resolved
    mask = adata.obs[cluster_key].astype(str) == cell_type
    if mask.sum() < 20:
        raise ValueError(f"Too few '{cell_type}' cells to subcluster ({int(mask.sum())}).")

    sub = adata[mask].copy()
    # Re-embed + Leiden on the subset only, re-deriving normalization from RAW counts: the parent's
    # X is already log-normalized, so reset X from the inherited counts layer and clear the derived
    # flags before cluster() (which would otherwise double-normalize, or reuse the parent's median
    # scale-factor + full-panel HVGs that are stale for a single cell type).
    if "counts" in sub.layers:
        sub.X = sub.layers["counts"].copy()
    for k in ("neighbors", "log1p", "pca", "hvg"):
        sub.uns.pop(k, None)
    sub.obsm.pop("X_pca", None)
    _cluster.cluster(sub, resolution=resolution, key="subleiden")

    # cluster() leaves low-signal cells' subleiden NaN (outside cat.categories) - they are not a
    # subtype, so skip them here and let them keep the parent cell_type label below.
    groups = list(sub.obs["subleiden"].cat.categories)
    top = {g: _cluster.top_markers(sub, g, n=12, key="rank_genes_groups") for g in groups}

    labels = {}
    if use_llm:
        from . import llm

        named = llm.name_subtypes(cell_type, top, context=context)
        labels = {g: named.get(g, {}) for g in groups}

    rows = []
    mapping = {}
    for g in groups:
        info = labels.get(g, {})
        label = info.get("label") or f"{cell_type} subcluster {g}"
        mapping[g] = label
        rows.append({
            "subcluster": g,
            "n_cells": int((sub.obs["subleiden"] == g).sum()),
            "top_markers": top[g],
            "label": label,
            "confidence": info.get("confidence", "n/a"),
            "rationale": info.get("rationale", ""),
        })

    # Low-signal cells (subleiden NaN) map to NaN -> keep the parent cell_type label, not a subtype.
    sub.obs["subtype"] = (
        sub.obs["subleiden"].astype(str).map(mapping).fillna(cell_type).astype("category")
    )
    import pandas as pd
    # Init `subtype` EMPTY (NaN) so cells NOT subclustered stay uncoloured CONTEXT on the map instead
    # of inheriting cell_type -- otherwise colouring by `subtype` recolours the whole tissue. Init only
    # if absent so successive subcluster runs on different types accumulate.
    if "subtype" not in adata.obs:
        adata.obs["subtype"] = pd.Series(pd.NA, index=adata.obs_names, dtype="object")
    adata.obs["subtype"] = adata.obs["subtype"].astype("object")
    # Align by obs_names (robust to any row filtering inside cluster()) rather than positional mask.
    adata.obs.loc[sub.obs_names, "subtype"] = sub.obs["subtype"].astype(str).values
    adata.obs["subtype"] = adata.obs["subtype"].astype("category")
    return sub, rows
