"""Panel-adequacy / cell-type resolvability check (feature H3).

What it does
------------
Given the panel gene list (controls already stripped by ``io.load``), tell the user, in
plain language, which cell types the panel can resolve and which are indistinguishable:

  (a) marker coverage   - fraction of each cell type's canonical markers on the panel
                          (traffic light: >=3 present = green, 1-2 = amber, 0 = red);
  (b) discriminability  - pairs of cell types with no private on-panel marker either way
                          are flagged "cannot separate"; the confusability graph's
                          connected components suggest merge groups.

`check_panel` returns a JSON-able dict stored in ``adata.uns['panel_check']`` and
resurfaced at annotation time. Claude turns it into a verdict (see ``llm.panel_verdict``)
grounded strictly in these tables - gene *presence* is necessary, not sufficient
(imaging probes drop out), so the numbers are always shown alongside.

Depends on
----------
Optional: ``decoupler`` v2 (``dc.op.resource``) for a marker DB. Falls back to the
curated repo panels in :mod:`spatialscribe.analysis.markers` when a DB is unavailable
(keeps the demo offline-reproducible).
"""

from __future__ import annotations

from itertools import combinations

from . import markers as _m

# One-vs-rest AUC (and gene-presence coverage) are noise for a type with only a handful of cells:
# a 6-49-positive AUC is not evidence of separability, and a 0-cell "green coverage" type is not
# something the section confirms can be typed. Below this floor a typability row is flagged
# weak_evidence and is NOT counted confidently_typable (documented in CLAUDE.md; now enforced in code).
MIN_TYPABLE_CELLS = 50


def _load_marker_db(marker_db: str, organism: str) -> tuple[dict[str, list[str]], str]:
    """Return ({cell_type: [markers]}, source) from decoupler, else curated repo panels."""
    try:
        import decoupler as dc

        net = dc.op.resource(marker_db, organism=organism, license="academic")
        # Column names vary across resources; find the (cell_type, gene) columns.
        cols = {c.lower(): c for c in net.columns}
        src = cols.get("cell_type") or cols.get("source") or cols.get("celltype")
        tgt = cols.get("genesymbol") or cols.get("target") or cols.get("gene")
        if "canonical_marker" in cols:  # keep canonical markers when the flag exists
            net = net[net[cols["canonical_marker"]].astype(str).str.lower().isin({"true", "1"})]
        out: dict[str, list[str]] = {}
        for ct, g in zip(net[src], net[tgt]):
            out.setdefault(str(ct), []).append(str(g))
        if out:
            return out, marker_db
    except Exception:
        pass
    # Fallback: curated lineages + TAM states from the repo.
    return {**_m.LINEAGE_MARKERS, **_m.TAM_STATES}, "curated-fallback"


def _traffic(n_present: int) -> str:
    return "green" if n_present >= 3 else "amber" if n_present >= 1 else "red"


def check_panel(panel_genes, marker_sets: dict[str, list[str]] | None = None,
                marker_db: str = "PanglaoDB", organism: str = "human",
                cell_types: list[str] | None = None) -> dict:
    """Assess whether ``panel_genes`` can resolve/distinguish cell types.

    Parameters
    ----------
    marker_sets: explicit ``{cell_type: [markers]}`` (takes precedence - use the curated
        melanoma panels for the demo); otherwise loaded from ``marker_db`` via decoupler.

    Returns a dict: ``{marker_db, coverage: {ct: {...}}, confusable_pairs: [...],
    merge_groups: [[...]], n_cell_types}``.
    """
    panel = {g.upper() for g in panel_genes}
    if marker_sets is not None:
        db, source = dict(marker_sets), "custom"
    else:
        db, source = _load_marker_db(marker_db, organism)
    if cell_types:
        db = {k: v for k, v in db.items() if k in set(cell_types)}

    coverage: dict[str, dict] = {}
    present_by_ct: dict[str, set[str]] = {}
    for ct, mk in db.items():
        mk_up = [g.upper() for g in mk]
        present = [g for g in mk_up if g in panel]
        missing = [g for g in mk_up if g not in panel]
        present_by_ct[ct] = set(present)
        coverage[ct] = {
            "n_markers": len(mk_up),
            "n_present": len(present),
            "coverage_frac": round(len(present) / len(mk_up), 3) if mk_up else 0.0,
            "present": present,
            "missing": missing,
            "status": _traffic(len(present)),
        }

    # (b) pairwise discriminability on on-panel markers.
    confusable = []
    for a, b in combinations(present_by_ct, 2):
        pa, pb = present_by_ct[a], present_by_ct[b]
        private_a, private_b = pa - pb, pb - pa
        if not private_a and not private_b and (pa or pb):
            confusable.append({
                "pair": [a, b],
                "shared_present": sorted(pa & pb),
                "reason": "no private on-panel marker for either type",
            })

    merge_groups = _connected_components(
        list(present_by_ct), [tuple(c["pair"]) for c in confusable]
    )
    return {
        "marker_db": source,
        "n_cell_types": len(coverage),
        "coverage": coverage,
        "confusable_pairs": confusable,
        "merge_groups": [g for g in merge_groups if len(g) > 1],
    }


def is_valid(pc) -> bool:
    """True if ``pc`` is a usable check_panel result.

    Guards against the h5ad round-trip that splits '/'-containing cell-type keys (e.g.
    'Epithelial/Tumor') into nested groups - recompute (cheap) when this returns False.
    """
    try:
        cov = pc["coverage"]
        cov_ok = len(cov) > 0 and all(isinstance(d, dict) and "status" in d for d in cov.values())
        # An h5ad round-trip also turns the confusable_pairs list into a numpy array; require a list.
        return cov_ok and isinstance(pc.get("confusable_pairs", []), list)
    except Exception:
        return False


def _connected_components(nodes, edges) -> list[list[str]]:
    """Undirected connected components (merge groups) from confusable pairs."""
    parent = {n: n for n in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in edges:
        parent[find(a)] = find(b)
    groups: dict[str, list[str]] = {}
    for n in nodes:
        groups.setdefault(find(n), []).append(n)
    return list(groups.values())


def identifiability_auc(adata, cluster_key: str = "cell_type",
                        marker_sets: dict[str, list[str]] | None = None) -> dict:
    """Per-cell-type identifiability GIVEN THE RESTRICTED PANEL.

    For each cell type, one-vs-rest ROC-AUC of that type's on-panel marker score at separating it
    from all other cells (using the section's own labels). High AUC (~1) => the panel can
    confidently pick the type out; ~0.5 => the panel's markers for that type don't distinguish it
    here (a probe-set limitation, complementing the marker-coverage traffic light). Returns
    ``{cell_type: {auc, n_markers, n_cells}}``.
    """
    import numpy as np
    from sklearn.metrics import roc_auc_score

    marker_sets = marker_sets or _m.LINEAGE_MARKERS
    labels = adata.obs[cluster_key].astype(str).to_numpy()
    out: dict[str, dict] = {}
    for ct, mk in marker_sets.items():
        genes = _m.on_panel(adata.var_names, mk)   # species/case-adaptive (human markers -> mouse panel)
        y = (labels == ct).astype(int)
        n_pos = int(y.sum())
        if not genes or n_pos == 0 or n_pos == len(y):
            out[ct] = {"auc": None, "n_markers": len(genes), "n_cells": n_pos}
            continue
        from .backend import get_backend
        get_backend().score_genes(adata, gene_list=genes, score_name="_idf", ctrl_size=min(50, adata.n_vars))
        score = np.asarray(adata.obs["_idf"], dtype=float)
        try:
            auc = float(roc_auc_score(y, score))
        except Exception:
            auc = None
        out[ct] = {"auc": auc, "n_markers": len(genes), "n_cells": n_pos}
    if "_idf" in adata.obs.columns:
        del adata.obs["_idf"]
    return out


def typability_table(adata, cluster_key: str = "cell_type",
                     marker_sets: dict[str, list[str]] | None = None,
                     reference_match: dict | None = None, auc_min: float = 0.7,
                     panel_check_result: dict | None = None,
                     min_typable_cells: int | None = None) -> list[dict]:
    """A COMPUTABLE per-cell-type "can this type be confidently typed on this panel?" decision.

    Deliberately does NOT decide on marker PRESENCE alone (a marker being on-panel is necessary but
    not sufficient). It prefers the strongest computed discriminability available, per cell type:

      1. **depth-matched per-class F1** from a reference (``reference_match['per_type']``) - the gold
         standard (panel-size aware; supersedes one-vs-rest AUC). ``confidently_typable`` = tier
         ``resolvable``.
      2. else **one-vs-rest identifiability AUC** on the section's own labels
         (:func:`identifiability_auc`). ``confidently_typable`` = ``auc >= auc_min``.
      3. else **marker coverage** (gene presence) - returned with ``basis='coverage_only'`` and
         ``weak_evidence=True`` so the caller (and the LLM) knows this row is a hand-wave, not a
         computed separability.

    Returns a list of ``{cell_type, confidently_typable, basis, score, threshold, n_present,
    n_markers, ...}``, hardest-first (not-typable, then lowest score).
    """
    marker_sets = marker_sets or _m.LINEAGE_MARKERS
    min_n = int(min_typable_cells if min_typable_cells is not None else MIN_TYPABLE_CELLS)
    pc = panel_check_result or adata.uns.get("panel_check") or {}
    cov = pc.get("coverage", {}) if is_valid(pc) else {}
    ref_per = (reference_match or {}).get("per_type", {}) or {}
    auc: dict = {}
    if not ref_per and cluster_key in getattr(adata, "obs", {}):
        try:
            auc = identifiability_auc(adata, cluster_key, marker_sets)
        except Exception:
            auc = {}

    # Per-type cell counts (only when the section carries labels). Used to discount a tiny-n AUC (noise)
    # and a 0-cell "green coverage" type (nothing present to type) - the documented <~50-cell guard. When
    # no labels exist (n_cells is None) the coverage hypothetical is preserved (the panel COULD resolve it).
    counts: dict[str, int] = {}
    if cluster_key in getattr(adata, "obs", {}):
        vc = adata.obs[cluster_key].astype(str).value_counts()
        counts = {str(k): int(v) for k, v in vc.items()}

    types = set(cov) | set(ref_per) | set(auc)
    rows: list[dict] = []
    for ct in types:
        d = cov.get(ct, {}) if isinstance(cov.get(ct), dict) else {}
        n_present, n_markers = d.get("n_present"), d.get("n_markers")
        n_cells = counts.get(str(ct)) if counts else None
        enough_cells = (n_cells is None) or (n_cells >= min_n)   # None = no labels -> do not gate
        if ct in ref_per:
            f1 = float(ref_per[ct].get("f1", 0.0)); tier = ref_per[ct].get("tier")
            rows.append({"cell_type": str(ct), "confidently_typable": tier == "resolvable",
                         "basis": "depth_matched_f1", "score": round(f1, 3), "threshold": None,
                         "tier": tier, "confused_with": ref_per[ct].get("confused_with"),
                         "n_present": n_present, "n_markers": n_markers, "n_cells": n_cells,
                         "weak_evidence": False})
        elif auc.get(ct, {}).get("auc") is not None:
            a = float(auc[ct]["auc"])
            rows.append({"cell_type": str(ct), "confidently_typable": (a >= auc_min) and enough_cells,
                         "basis": "identifiability_auc", "score": round(a, 3), "threshold": auc_min,
                         "n_present": n_present, "n_markers": (n_markers or auc[ct].get("n_markers")),
                         "n_cells": (n_cells if n_cells is not None else auc[ct].get("n_cells")),
                         "weak_evidence": not enough_cells})
        else:
            status = d.get("status")
            # gene-presence only: a 0-cell (or tiny-n) type is NOT confidently typable, only a hypothetical.
            rows.append({"cell_type": str(ct),
                         "confidently_typable": (status == "green") and enough_cells,
                         "basis": "coverage_only",
                         "score": (round(float(d.get("coverage_frac", 0)), 3) if d else None),
                         "threshold": None, "n_present": n_present, "n_markers": n_markers,
                         "n_cells": n_cells, "status": status, "weak_evidence": True})
    # Panel-reality floor: a type with ZERO of its canonical markers ON THIS PANEL cannot be confidently
    # typed on it, whatever a reference's depth-matched F1 says - that F1 lives in the reference's gene
    # space, not the restricted panel's. Without this, a foreign type dragged in by a wrong-tissue
    # reference (e.g. Schwann cell / melanocyte / OPC on a gut panel) read "0/6 markers, typable: yes".
    for r in rows:
        if r.get("n_present") == 0 and r.get("confidently_typable"):
            r["confidently_typable"] = False
            r["zero_on_panel"] = True   # the panel measures none of this type's markers
    rows.sort(key=lambda r: (r["confidently_typable"], r["score"] if r["score"] is not None else 1.0))
    return rows
