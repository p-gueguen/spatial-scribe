"""Spatial niches / TME domains - the tissue neighborhoods cells live in.

Fills the "what are the TME niches?" gap (tumor core / immune-excluded margin / TLS-like).
Two backends:

- ``method="composition"`` (default, no extra deps): for each cell, build the cell-type
  composition of its spatial neighborhood, then cluster those composition vectors into
  niches and name each by its dominant types. CPU-light, uses squidpy + scikit-learn.
- ``method="novae"`` (optional): the scverse graph foundation model (pretrained on ~30M
  iST cells) for zero-shot niches; used if ``novae`` is installed.

Writes ``adata.obs['niche']`` and returns a per-niche composition table.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas


def neighborhood_composition(adata, cluster_key: str = "cell_type", n_neighs: int = 15):
    """Return an (n_cells x n_types) matrix of each cell's neighborhood composition.

    The abstention pseudo-labels collapse into ONE ``Not assigned`` feature. They are kept (a cell
    surrounded by uncallable tissue is a real, informative niche) but they are not five distinct
    lineages, which is what the raw ``cell_type_final`` categories would assert.
    """
    import numpy as np

    from . import annotate as _an
    from .backend import get_backend

    if "spatial_connectivities" not in adata.obsp:
        get_backend().spatial_neighbors(adata, coord_type="generic", n_neighs=n_neighs)
    A = adata.obsp["spatial_connectivities"]
    labels = _an.collapse_abstention(adata.obs[cluster_key]).astype("category")
    cats = list(map(str, labels.cat.categories))
    codes = labels.cat.codes.to_numpy()
    # one-hot of cell types, then average over neighbors (row-normalized adjacency).
    onehot = np.zeros((adata.n_obs, len(cats)), dtype=float)
    onehot[np.arange(adata.n_obs), codes] = 1.0
    deg = np.asarray(A.sum(1)).ravel()
    deg[deg == 0] = 1
    comp = (A @ onehot) / deg[:, None]
    return comp, cats


def call_niches(adata, cluster_key: str = "cell_type", n_niches: int = 8,
                n_neighs: int = 15, method: str = "composition", progress=None) -> "pandas.DataFrame":  # noqa: F821
    """Assign each cell a spatial niche. Writes ``adata.obs['niche']``; returns a table.

    ``progress(frac, label)`` (optional) reports coarse checkpoints for the app's progress bar.
    """
    import numpy as np
    import pandas as pd

    if method == "novae":
        return _novae_niches(adata, n_niches)

    if progress:
        progress(0.20, "spatial neighbour graph")
    comp, cats = neighborhood_composition(adata, cluster_key, n_neighs)
    if progress:
        progress(0.55, "neighbourhood profiles")
    from sklearn.cluster import KMeans, MiniBatchKMeans

    # The old KMeans(n_init=10) default clusters the whole section 10x -> ~10 s on 100k cells (this
    # was the slow "Spatial+niches" load). Scale n_init down with size, and switch to MiniBatchKMeans
    # past 100k so niche calling stays responsive up to ~10^6 cells.
    n_k = min(n_niches, adata.n_obs)
    if adata.n_obs > 300_000:
        km = MiniBatchKMeans(n_clusters=n_k, n_init=3, batch_size=4096, random_state=0)
    else:
        # sklearn's old n_init=10 default clusters the whole section 10x (~7 s on 100k cells). On
        # these low-dimensional composition vectors 3 restarts match the 10-restart inertia to
        # within 0.02% at ~1/3 the cost; keep 10 only on small sections where it is already cheap.
        km = KMeans(n_clusters=n_k, n_init=(3 if adata.n_obs > 20_000 else 10), random_state=0)
    if progress:
        progress(0.85, "clustering niches")
    lab = km.fit_predict(comp)

    # Name each niche by its top-2 enriched cell types (vs the global mean).
    glob = comp.mean(0)
    rows, names = [], {}
    from . import annotate as _an

    # "Not assigned" stays a composition FEATURE (a niche of uncallable tissue is real), but never
    # NAMES a niche: "Niche 3: Not assigned + T cell" describes annotation failure, not biology.
    nameable = [i for i, c in enumerate(cats) if not _an.is_abstention(c)] or list(range(len(cats)))
    for k in range(km.n_clusters):
        m = lab == k
        prof = comp[m].mean(0)
        enr = prof - glob
        order = sorted(nameable, key=lambda i: -enr[i])
        top = [cats[i] for i in order[:2]]
        name = f"Niche {k}: {' + '.join(top)}"
        names[k] = name
        rows.append({"niche": name, "n_cells": int(m.sum()),
                     "top_types": top, "dominant_frac": float(prof.max())})
    adata.obs["niche"] = pd.Categorical([names[k] for k in lab])
    if progress:
        progress(1.0, "done")
    return pd.DataFrame(rows)


def _novae_niches(adata, n_niches: int = 8):
    """Optional zero-shot niches via the novae graph foundation model.

    Verified against novae 1.1.0 on a 100k-cell Xenium demo (8 niches, transcriptomics-only, no H&E).
    """
    try:
        import novae
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError(
            "method='novae' needs the optional `novae` package (pip install novae). "
            "Use method='composition' for the dependency-free caller."
        ) from exc
    # novae needs a spatial graph and RAW counts in X (a processed section may ship log1p X; the raw
    # counts are kept in obsm/layers). obsm['spatial'] must be present.
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    novae.spatial_neighbors(adata)                                  # required before representations
    model = novae.Novae.from_pretrained("prism-oncology/novae-human-0")  # org renamed from MICS-Lab
    model.compute_representations(adata, zero_shot=True)            # -> obsm['novae_latent']
    # assign_domains RETURNS the obs column name it wrote; `level` is novae's hierarchy level
    # (n_domains=/resolution= are alternatives - resolution is recommended for zero-shot).
    dom_key = model.assign_domains(adata, level=n_niches) or next(
        c for c in adata.obs.columns if c.startswith("novae_domains"))  # fallback if a version returns None
    # Cast to str (not left categorical): a categorical niche column mapped to RGBA tuples downstream
    # builds a tuple-Index -> MultiIndex, and pandas isna() is not implemented for MultiIndex.
    adata.obs["niche"] = adata.obs[dom_key].astype(str).astype("category")
    return adata.obs["niche"].value_counts().rename_axis("niche").reset_index(name="n_cells")
