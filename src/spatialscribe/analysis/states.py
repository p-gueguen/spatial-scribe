"""Cell-type x cell-state view (H6) - programs orthogonal to lineage identity.

Scores cell-state signatures (cycling, IFN/ISG, hypoxia, stress/HSP, EMT, T-exhaustion,
T-cytotoxicity, ECM remodeling, and antigen presentation - the default `markers.CELL_STATES`)
per cell and summarizes them as a cell-type x state matrix - so a scientist sees, e.g., which
lineage is proliferating or interferon-activated, on top of who each cell is. Panel-restricted;
small control set (the Xenium maxRank=150 spirit).

Depends on: scanpy, pandas, numpy; :mod:`markers`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from . import markers as _m

if TYPE_CHECKING:
    import pandas


def score_states(adata, state_sets: dict[str, list[str]] | None = None) -> list[str]:
    """Add one ``state_<name>`` score column per state signature. Returns the names."""
    from .markers import present

    from .backend import get_backend

    be = get_backend()
    # Default = the generic tumour-agnostic programs (CELL_STATES) ONLY. The Myeloid-only TAM states
    # were dropped from the default (they still live in `markers.TAM_STATES` for callers that want them).
    sets = state_sets or _m.CELL_STATES
    panel = set(adata.var_names)
    added = []
    for name, genes in present(panel, sets).items():
        if len(genes) < 2:                       # need >=2 on-panel genes to score a program
            continue
        col = "state_" + name.replace("/", "_").replace(" ", "_")
        be.score_genes(adata, gene_list=genes, score_name=col, ctrl_size=min(50, adata.n_vars))
        added.append((name, col))
    adata.uns["state_columns"] = dict(added)
    return [c for _, c in added]


def state_by_celltype(adata, cluster_key: str = "cell_type",
                      state_sets: dict[str, list[str]] | None = None) -> "pandas.DataFrame":  # noqa: F821
    """Return a cell-type (rows) x state (cols) mean-score matrix for the H6 heatmap.

    The score->obs-column mapping (``{state_name: 'state_<name>'}``, the same one in
    ``uns['state_columns']``) is also attached to ``df.attrs['score_fields']`` so a caller can
    recolour the canvas by a program without re-deriving the column names.
    """
    import numpy as np
    import pandas as pd

    from . import annotate as _an

    score_states(adata, state_sets)
    colmap = adata.uns["state_columns"]           # {state_name: column}
    if not colmap:
        return pd.DataFrame()
    rows = {}
    for ct, sub in adata.obs.groupby(cluster_key, observed=True):
        # An abstention is not a cell type: "mean IFN score of the Ambiguous: mixed cells" is not a
        # biological statement, and it would read as a lineage row in the delivered heatmap.
        if _an.is_abstention(str(ct)):
            continue
        rows[str(ct)] = {name: float(np.mean(sub[col])) for name, col in colmap.items()}
    df = pd.DataFrame(rows).T                       # cell types x states
    df.attrs["score_fields"] = dict(colmap)
    return df


def assign_cell_states(adata, cluster_key: str = "cell_type",
                       state_sets: dict[str, list[str]] | None = None,
                       min_z: float = 1.0) -> dict:
    """Type each cell with a DOMINANT cell-state program (CyteType-style state typing).

    Scores the default state signatures (cycling / IFN-ISG / hypoxia / stress-HSP / EMT / T-exhaustion /
    T-cytotoxicity / ECM-remodeling / antigen-presentation - `markers.CELL_STATES`, minus any program
    with <2 on-panel genes), z-scores each program across cells so they are comparable, and labels each cell with the single
    highest-scoring program IF it clears ``min_z`` (else ``"None"`` - no dominant program). Writes
    ``obs['cell_state']`` (categorical, colourable on the map). Returns the state distribution and,
    per cell type, its dominant (non-baseline) state and how prevalent it is.

    States are not mutually exclusive biologically; this reports the single strongest program per cell
    for an interpretable, colourable label. The per-program z-scores stay in ``obs['state_*']`` (used
    by the heatmap and the LLM naming), so nothing is lost.
    """
    import numpy as np
    import pandas as pd

    score_states(adata, state_sets)
    colmap = adata.uns.get("state_columns", {})    # {state_name: column}
    names = list(colmap)
    n = adata.n_obs
    if not names:
        adata.obs["cell_state"] = pd.Categorical(["None"] * n, categories=["None"])
        return {"states": [], "distribution": {"None": n}, "per_celltype": {}, "min_z": float(min_z),
                "note": "no state signatures had >=2 on-panel genes"}

    Z = np.zeros((n, len(names)), dtype=float)
    for j, name in enumerate(names):
        v = np.asarray(adata.obs[colmap[name]], dtype=float)
        Z[:, j] = (v - v.mean()) / (v.std() + 1e-9)
    top = Z.argmax(1)
    topz = Z[np.arange(n), top]
    labels = np.array([names[i] for i in top], dtype=object)
    labels[topz < float(min_z)] = "None"           # below threshold on every program -> baseline

    present_cats = [c for c in [*names, "None"] if c in set(labels.tolist())] or ["None"]
    adata.obs["cell_state"] = pd.Categorical(labels, categories=present_cats)

    dist = {str(k): int(v) for k, v in pd.Series(labels).value_counts().items()}
    per_ct: dict[str, dict] = {}
    if cluster_key in adata.obs:
        from . import annotate as _an

        df = pd.DataFrame({"ct": adata.obs[cluster_key].astype(str).to_numpy(), "state": labels})
        for ct, g in df.groupby("ct"):
            if _an.is_abstention(str(ct)):     # not a cell type -> no per-cell-type state summary
                continue
            non_none = g["state"][g["state"] != "None"].value_counts()
            top_state = str(non_none.index[0]) if len(non_none) else "None"
            per_ct[str(ct)] = {
                "top_state": top_state,
                "frac_top": round(float((g["state"] == top_state).mean()), 3),
                "pct_stateful": round(float((g["state"] != "None").mean()), 3)}
    return {"states": names, "distribution": dist, "per_celltype": per_ct, "min_z": float(min_z)}


def name_states_llm(adata, cluster_key: str = "cell_type", context: str = "human tumour",
                    state_sets: dict[str, list[str]] | None = None, n_markers: int = 12) -> dict:
    """CyteType-style LLM naming of each cell type's dominant functional STATE, grounded in its mean
    program z-scores + its top markers. Returns ``{cell_type: {state, confidence, rationale}}`` (or
    ``{}`` if there is nothing to score / the LLM is unavailable). Never raises."""
    from . import llm as _llm

    mat = state_by_celltype(adata, cluster_key, state_sets)     # cell-type x state mean scores
    if mat.empty:
        return {}
    z = (mat - mat.mean(0)) / (mat.std(0) + 1e-9)               # z per program across cell types
    per: dict[str, dict] = {}
    for ct in mat.index:
        try:
            from . import cluster as _cl
            markers = _cl.top_markers(adata, str(ct), n=n_markers) if cluster_key in adata.obs else []
        except Exception:
            markers = []
        per[str(ct)] = {"markers": markers,
                        "state_scores": {str(s): round(float(z.loc[ct, s]), 2) for s in mat.columns}}
    try:
        return _llm.name_cell_states(per, context=context)
    except Exception:
        return {}
