"""Cluster confidence - data-driven over/under-clustering merge/split nudges.

What it does
------------
The top rung of the app's uncertainty ladder (cluster -> cell -> panel). Instead of forcing
the user to guess a Leiden resolution, it tells them, grounded in real numbers:

  * ``merge_test``  - which cluster PAIRS are statistically indistinguishable. For each adjacent
    pair it trains a RandomForest to separate the two clusters on the PCA embedding and compares
    its cross-validated balanced accuracy to a permuted-label null (``permutation_test_score``).
    A pair the classifier cannot beat chance on (``p > alpha`` or accuracy below a floor) is an
    over-split -> **merge**. This is a lightweight reimplementation of the sc-SHC / CHOIR idea
    (Grabski & Purdom, Genome Biology 2023; CHOIR, Nature Genetics 2025) in scikit-learn - it does
    NOT wrap their R packages.
  * ``split_test``  - which single clusters hide substructure. It projects each cluster's cells
    onto their main within-cluster PC and computes the bimodality coefficient
    ``BC = (skew^2 + 1) / kurtosis`` (uniform -> 1, normal -> 1/3); ``BC`` above ~5/9 means the
    dominant axis is bimodal -> **split** candidate. When de-novo NMF programs are present
    (``obs['program']``) it additionally requires the cluster to span >= 2 programs, so a single
    program's spread is not mistaken for two cell states.
  * ``cluster_confidence`` - the orchestrator that runs both, ties each cluster to its dominant
    annotation (a merge pair sharing one label is a stronger nudge), and returns the unified panel
    the wizard/report/copilot consume.

Everything here is ADVISORY: it never merges or splits on its own; it hands the user a ranked,
plain-language nudge with the underlying statistics. Reference-free, CPU-only, no API key.

How to use it
-------------
>>> from spatialscribe.analysis import cluster_confidence as ccf
>>> ccf.cluster_confidence(adata, cluster_key="leiden", annotation_key="cell_type")

Depends on: numpy, pandas, scipy.stats (skew/kurtosis), scikit-learn (RandomForest,
permutation_test_score, PCA). Reads :mod:`spatialscribe.analysis.config` for tunable thresholds.
Heavy imports live inside the functions so importing the module stays cheap.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _embedding(adata, key: str):
    """Return the cell x dim embedding matrix (``obsm[key]``, else a lazy PCA of ``X``)."""
    import numpy as np

    if key in adata.obsm:
        return np.asarray(adata.obsm[key], dtype=float)
    from sklearn.decomposition import PCA

    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X, dtype=float)
    n_comp = int(min(30, X.shape[1], max(2, X.shape[0] - 1)))
    return PCA(n_components=n_comp, random_state=0).fit_transform(X)


def _clusters(adata, cluster_key: str):
    """Ordered list of cluster labels (respects an existing categorical order)."""
    s = adata.obs[cluster_key]
    if hasattr(s, "cat"):
        return [str(c) for c in s.cat.categories if (s.astype(str) == str(c)).any()]
    return sorted(str(c) for c in s.unique())


def _adjacent_pairs(E, labels, clusters, k: int):
    """Pairs to test: all pairs when few clusters, else each cluster with its ``k`` nearest
    centroids (keeps the O(n_clusters^2) RF-training cost bounded on finely-clustered sections)."""
    from itertools import combinations

    import numpy as np

    if len(clusters) <= k + 1:
        return list(combinations(clusters, 2))
    cents = np.vstack([E[labels == c].mean(0) for c in clusters])
    d = np.linalg.norm(cents[:, None, :] - cents[None, :, :], axis=2)
    pairs = set()
    for i in range(len(clusters)):
        for j in np.argsort(d[i])[1 : k + 1]:
            pairs.add(tuple(sorted((clusters[i], clusters[int(j)]))))
    return sorted(pairs)


def _bimodality_coefficient(x) -> float:
    """Sample bimodality coefficient ``(g1^2 + 1) / (g2 + 3(n-1)^2/((n-2)(n-3)))`` (SAS form).

    ``g1`` = sample skewness, ``g2`` = sample EXCESS kurtosis. Uniform -> 1, normal -> ~1/3,
    a clean two-mode mixture -> > 5/9. Returns 0 for degenerate input."""
    import numpy as np
    from scipy.stats import kurtosis, skew

    x = np.asarray(x, dtype=float)
    n = x.size
    if n < 4 or np.allclose(x, x[0]):
        return 0.0
    g1 = float(skew(x, bias=False))
    g2 = float(kurtosis(x, fisher=True, bias=False))
    denom = g2 + 3.0 * (n - 1) ** 2 / ((n - 2) * (n - 3))
    if denom <= 0:
        return 0.0
    return float((g1 ** 2 + 1.0) / denom)


# --------------------------------------------------------------------------- #
# merge
# --------------------------------------------------------------------------- #
def merge_test(adata, cluster_key: str = "leiden", embedding: str = "X_pca", *,
               accuracy_floor: float | None = None, max_cells: int = 120,
               n_permutations: int | None = None, cv: int = 3, n_estimators: int = 50,
               adjacency_k: int = 3, max_dims: int = 20, random_state: int = 0) -> "pandas.DataFrame":  # noqa: F821
    """Pairwise "are these two clusters the same population?" test (sc-SHC / CHOIR-style).

    Returns a DataFrame ``[cluster_a, cluster_b, rf_accuracy, null_mean, p_value, verdict]`` sorted
    most-mergeable first. The verdict is ACCURACY-primary (CHOIR's core): ``'merge'`` when a
    RandomForest cannot separate the pair above ``accuracy_floor`` balanced accuracy, else
    ``'distinct'``. ``null_mean`` (permuted-label accuracy, ~0.5) and ``p_value`` are reported for
    interpretability but do NOT drive the verdict, so a coarse (small-``n_permutations``) p never
    misclassifies a clearly-separable pair. Empty DataFrame when there are < 2 clusters. Compute is
    bounded (top ``max_dims`` PCs, ``max_cells``/cluster, ``adjacency_k`` nearest neighbours) so a
    synchronous wizard/copilot call stays fast. Reference-free, CPU.
    """
    import numpy as np
    import pandas as pd

    from . import config

    if accuracy_floor is None:
        accuracy_floor = float(config.get("cluster_confidence", "merge_accuracy_floor", default=0.65))
    if n_permutations is None:
        n_permutations = int(config.get("cluster_confidence", "merge_n_permutations", default=10))

    cols = ["cluster_a", "cluster_b", "rf_accuracy", "null_mean", "p_value", "verdict"]
    clusters = _clusters(adata, cluster_key)
    if len(clusters) < 2:
        return pd.DataFrame(columns=cols)

    E = _embedding(adata, embedding)[:, :max_dims]
    labels = adata.obs[cluster_key].astype(str).to_numpy()
    rng = np.random.default_rng(random_state)

    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import permutation_test_score

    rows = []
    for a, b in _adjacent_pairs(E, labels, clusters, adjacency_k):
        ia, ib = np.where(labels == a)[0], np.where(labels == b)[0]
        n_sub = int(min(len(ia), len(ib), max_cells))
        eff_cv = int(min(cv, n_sub))
        if n_sub < 2 or eff_cv < 2:
            continue
        sa = rng.choice(ia, n_sub, replace=False)
        sb = rng.choice(ib, n_sub, replace=False)
        X = np.vstack([E[sa], E[sb]])
        y = np.r_[np.zeros(n_sub), np.ones(n_sub)]
        clf = RandomForestClassifier(n_estimators=n_estimators, random_state=random_state,
                                     class_weight="balanced", n_jobs=1)
        acc, perm, p = permutation_test_score(
            clf, X, y, cv=eff_cv, n_permutations=n_permutations,
            scoring="balanced_accuracy", random_state=random_state, n_jobs=1)
        verdict = "merge" if acc < accuracy_floor else "distinct"
        rows.append({"cluster_a": a, "cluster_b": b, "rf_accuracy": round(float(acc), 4),
                     "null_mean": round(float(np.mean(perm)), 4), "p_value": round(float(p), 4),
                     "verdict": verdict})

    df = pd.DataFrame(rows, columns=cols)
    return df.sort_values("rf_accuracy", ascending=True).reset_index(drop=True) if not df.empty else df


# --------------------------------------------------------------------------- #
# split
# --------------------------------------------------------------------------- #
def split_test(adata, cluster_key: str = "leiden", embedding: str = "X_pca", *,
               bimodality_threshold: float | None = None, min_cells: int | None = None,
               program_key: str = "program", program_min_frac: float = 0.15) -> "pandas.DataFrame":  # noqa: F821
    """Per-cluster "does this cluster hide substructure?" flag.

    Returns a DataFrame ``[cluster, n_cells, bimodality, n_programs, verdict]``. ``verdict ==
    'split'`` when the cluster's main within-cluster axis is bimodal (``bimodality >
    threshold``); when NMF programs are present (``obs[program_key]``) the cluster must ALSO span
    >= 2 programs (each held by >= ``program_min_frac`` of its cells) to avoid over-splitting a
    single program's spread. Clusters below ``min_cells`` are reported as ``'ok'`` (too small).
    """
    import numpy as np
    import pandas as pd
    from sklearn.decomposition import PCA

    from . import config

    if bimodality_threshold is None:
        bimodality_threshold = float(config.get("cluster_confidence", "split_bimodality", default=0.555))
    if min_cells is None:
        min_cells = int(config.get("cluster_confidence", "min_cluster_cells", default=20))

    E = _embedding(adata, embedding)
    labels = adata.obs[cluster_key].astype(str).to_numpy()
    has_prog = program_key in adata.obs
    prog = adata.obs[program_key].astype(str).to_numpy() if has_prog else None

    rows = []
    for c in _clusters(adata, cluster_key):
        idx = np.where(labels == c)[0]
        n = int(idx.size)
        n_prog = 0
        if has_prog and n:
            vc = pd.Series(prog[idx]).value_counts(normalize=True)
            n_prog = int((vc >= program_min_frac).sum())
        if n < min_cells:
            rows.append({"cluster": c, "n_cells": n, "bimodality": 0.0,
                         "n_programs": n_prog, "verdict": "ok"})
            continue
        proj = PCA(n_components=1, random_state=0).fit_transform(E[idx]).ravel()
        bc = _bimodality_coefficient(proj)
        prog_ok = (n_prog >= 2) if has_prog else True
        verdict = "split" if (bc > bimodality_threshold and prog_ok) else "ok"
        rows.append({"cluster": c, "n_cells": n, "bimodality": round(bc, 4),
                     "n_programs": n_prog, "verdict": verdict})
    return pd.DataFrame(rows, columns=["cluster", "n_cells", "bimodality", "n_programs", "verdict"])


# --------------------------------------------------------------------------- #
# orchestrator
# --------------------------------------------------------------------------- #
def cluster_confidence(adata, cluster_key: str = "leiden", annotation_key: str = "cell_type",
                       embedding: str = "X_pca", **kw) -> dict:
    """Unified cluster-confidence panel: merge suggestions + split suggestions + per-cluster verdict.

    Ties each cluster to its dominant ``annotation_key`` label (when present) so a merge pair that
    shares one label reads as a stronger nudge (``same_label=True``). ``**kw`` (e.g.
    ``n_permutations``) is forwarded to :func:`merge_test`. Returns a JSON-able dict the wizard,
    report and copilot render. Advisory only - it applies nothing.
    """
    merges = merge_test(adata, cluster_key, embedding, **kw)
    splits = split_test(adata, cluster_key, embedding)
    clusters = _clusters(adata, cluster_key)

    label_of: dict[str, str] = {}
    if annotation_key and annotation_key in adata.obs:
        sub = adata.obs[[cluster_key, annotation_key]].astype(str)
        for c, g in sub.groupby(cluster_key, observed=True):
            m = g[annotation_key].mode()
            label_of[str(c)] = str(m.iloc[0]) if len(m) else ""

    merge_suggestions = []
    for r in merges[merges["verdict"] == "merge"].itertuples(index=False):
        la, lb = label_of.get(r.cluster_a, ""), label_of.get(r.cluster_b, "")
        merge_suggestions.append({
            "cluster_a": r.cluster_a, "cluster_b": r.cluster_b,
            "label_a": la, "label_b": lb,
            "same_label": bool(la and lb and la == lb),
            "rf_accuracy": r.rf_accuracy, "p_value": r.p_value,
        })

    split_suggestions = [
        {"cluster": r.cluster, "label": label_of.get(r.cluster, ""),
         "bimodality": r.bimodality, "n_programs": int(r.n_programs)}
        for r in splits[splits["verdict"] == "split"].itertuples(index=False)
    ]

    merged_clusters = {c for m in merge_suggestions for c in (m["cluster_a"], m["cluster_b"])}
    split_clusters = {s["cluster"] for s in split_suggestions}
    per_cluster = [
        {"cluster": c, "label": label_of.get(c, ""),
         "verdict": ("merge" if c in merged_clusters else "split" if c in split_clusters else "ok")}
        for c in clusters
    ]

    return {
        "n_clusters": len(clusters),
        "merge_suggestions": merge_suggestions,
        "split_suggestions": split_suggestions,
        "cluster_verdicts": per_cluster,
        "nudge": _nudge(merge_suggestions, split_suggestions),
    }


def _nudge(merge_suggestions: list, split_suggestions: list) -> str:
    """One plain-language sentence summarizing the panel (grounded in the counts)."""
    bits = []
    if merge_suggestions:
        pairs = ", ".join(f"{m['cluster_a']}+{m['cluster_b']}" for m in merge_suggestions[:4])
        same = sum(m["same_label"] for m in merge_suggestions)
        extra = f" ({same} share a cell-type label)" if same else ""
        bits.append(f"{len(merge_suggestions)} cluster pair(s) look statistically indistinguishable "
                    f"- consider merging: {pairs}{extra}.")
    if split_suggestions:
        cl = ", ".join(s["cluster"] for s in split_suggestions[:4])
        bits.append(f"{len(split_suggestions)} cluster(s) show internal substructure "
                    f"- consider splitting: {cl}.")
    return " ".join(bits) if bits else "Clustering looks well resolved - no merge/split suggested."
