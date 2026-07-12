"""Trust ledger - three INDEPENDENT per-cell-type verdicts, and the disagreements between them.

The AQI is one honest *section* number, but it cannot see a coherent whole-cluster MISLABEL (the colon
case: internally tidy + marker-consistent, yet wrong). No single reference-free number can. The fix is
not a better number - it is to line up three signals that fail in different ways and flag where they
DISAGREE:

  * **resolvable?**  - can the panel even separate this type? (``panel_check.typability_table``)
  * **coherent?**    - do this type's cells actually express its markers? (``verify.verify_annotation``)
  * **agreed?**      - do the reference methods back the label? (per-cell ``consensus_agreement``, >=3 voters)

Agreement holds when all three agree; the informative rows are the contradictions - *resolvable but the
markers dispute it*, or (the mislabel signature) *coherent but the reference methods dispute it*. That
contradiction is the one thing the index alone can't give, and it is what a bench scientist should look at.
"""
from __future__ import annotations


def trust_ledger(adata, marker_sets: dict[str, list[str]] | None = None, cluster_key: str = "cell_type",
                 reference_match: dict | None = None, *, argmax_min: float = 0.5, agree_min: float = 0.5,
                 min_cells: int = 50) -> dict:
    """Per cell type: resolvable / coherent / agreed + a list of disagreement flags.

    Reuses the three existing computations (typability, verify, consensus_agreement) - no new science.
    Types with ``< min_cells`` cells are skipped (their per-type stats are noise). ``agreed`` is ``None``
    when no reference ensemble voted (``consensus_agreement`` absent); the disputed-mislabel flag then
    cannot fire and the ``note`` says to run the reference methods. Returns
    ``{per_type, n_types, n_flagged, has_ensemble, note}``; guarded, never raises the funnel.
    """
    import numpy as np

    from . import panel_check as _pc
    from . import verify as _vf

    key = cluster_key
    if key not in adata.obs:
        return {"per_type": [], "n_types": 0, "n_flagged": 0, "has_ensemble": False,
                "note": "no cell_type - annotate first"}

    # 1) resolvable (panel adequacy per type) - depth-matched F1 with a reference, else AUC, else coverage.
    try:
        tyt = {str(r.get("cell_type")): r for r in _pc.typability_table(
            adata, cluster_key=key, reference_match=reference_match)}
    except Exception:
        tyt = {}
    # 2) coherent (do the cells express the type's markers?)
    try:
        vf = _vf.verify_annotation(adata, marker_sets=marker_sets, cluster_key=key)
        per = vf.get("per_type", {}) or {}
    except Exception:
        per = {}
    # 3) agreed (reference-ensemble agreement per cell), only when >=3 methods voted upstream.
    has_ens = "consensus_agreement" in adata.obs
    labels = adata.obs[key].astype(str)
    agree_col = (np.asarray(adata.obs["consensus_agreement"], dtype=float) if has_ens else None)

    rows: list[dict] = []
    for ct, n in labels.value_counts().items():
        ct = str(ct)
        n = int(n)
        if n < min_cells:
            continue
        resolvable = tyt.get(ct, {}).get("confidently_typable")            # True | False | None
        v = per.get(ct, {})
        coh = v.get("argmax_agreement")                                    # 0-1 | None
        auc = v.get("auc")
        cause = v.get("cause")
        agreed = None
        if agree_col is not None:
            a = agree_col[(labels == ct).to_numpy()]
            a = a[~np.isnan(a)]
            agreed = float(a.mean()) if a.size else None

        coherent_ok = coh is not None and coh >= argmax_min
        flags: list[str] = []
        if resolvable is False:
            flags.append("panel cannot resolve this type")
        if coh is not None and coh < argmax_min:
            flags.append(f"markers dispute the label ({cause or 'low marker agreement'})")
        if agreed is not None and agreed < agree_min:
            flags.append("reference methods dispute the label")
        # THE mislabel signature: internally coherent, but the ensemble disagrees -> a whole-cluster mislabel
        # that AQI's C/M would score high and only cross-method agreement catches.
        if coherent_ok and agreed is not None and agreed < agree_min:
            flags.append("coherent but DISPUTED - possible whole-cluster mislabel; check the markers + reference")

        rows.append({
            "cell_type": ct, "n_cells": n,
            "resolvable": resolvable,
            "coherent": (None if coh is None else round(float(coh), 2)),
            "auc": (None if auc is None else round(float(auc), 2)),
            "agreed": (None if agreed is None else round(agreed, 2)),
            "cause": cause,
            "flags": flags,
        })
    # Contradictions first (most flags), then the biggest types.
    rows.sort(key=lambda r: (-len(r["flags"]), -r["n_cells"]))
    return {
        "per_type": rows, "n_types": len(rows),
        "n_flagged": sum(1 for r in rows if r["flags"]),
        "has_ensemble": has_ens,
        "note": ("reference ensemble voted - the 'agreed' column + the disputed-mislabel check are live"
                 if has_ens else
                 "no reference ensemble - 'agreed' is blank; run the reference methods (>=3) so a "
                 "coherent-but-wrong cluster can be caught, not just a marker-inconsistent one"),
    }


def _demo() -> None:
    """Self-check: a type that is marker-coherent but the ensemble disputes is flagged as a mislabel."""
    import anndata as ad
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(0)
    genes = ["CD3D", "MS4A1", "PECAM1", "N1", "N2"]
    n = 300
    X = rng.poisson(1, (n, 5)).astype("float32")
    ct = np.array(["T cell"] * 100 + ["B cell"] * 100 + ["Endothelial"] * 100)
    X[ct == "T cell", 0] += 30            # T cells express CD3D -> coherent
    X[ct == "B cell", 1] += 30            # B cells express MS4A1 -> coherent
    X[ct == "Endothelial", 2] += 30       # "Endothelial" cells express PECAM1 -> coherent
    a = ad.AnnData(X=X)
    a.var_names = genes
    a.obs["cell_type"] = pd.Categorical(ct)
    # ensemble disputes ONLY the Endothelial cluster (methods say it is something else) -> mislabel signature
    agree = np.where(ct == "Endothelial", 0.2, 0.9)
    a.obs["consensus_agreement"] = agree
    markers = {"T cell": ["CD3D"], "B cell": ["MS4A1"], "Endothelial": ["PECAM1"]}
    out = trust_ledger(a, marker_sets=markers, min_cells=10)
    assert out["has_ensemble"], out
    byct = {r["cell_type"]: r for r in out["per_type"]}
    endo = byct["Endothelial"]
    assert any("DISPUTED" in f for f in endo["flags"]), endo          # coherent + disputed -> flagged
    assert not byct["T cell"]["flags"], byct["T cell"]                # coherent + agreed -> clean
    assert out["n_flagged"] >= 1
    print("trust._demo OK:", out["n_flagged"], "flagged;", endo["flags"])


if __name__ == "__main__":
    _demo()
