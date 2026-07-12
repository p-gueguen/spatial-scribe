"""Trust ledger - the coherent-but-disputed mislabel flag, and clean rows when the signals agree."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")


def _section(agree_endo: float = 0.2):
    """Three marker-coherent types; the ensemble backs T/B but DISPUTES Endothelial (a mislabel)."""
    rng = np.random.default_rng(0)
    markers = ["CD3D", "MS4A1", "PECAM1"]
    filler = [f"G{i}" for i in range(40)]                 # spanning-range background so score_genes has controls
    genes = markers + filler
    ct = np.array(["T cell"] * 100 + ["B cell"] * 100 + ["Endothelial"] * 100)
    lam = rng.uniform(0.5, 15, len(genes))
    X = rng.poisson(lam, (len(ct), len(genes))).astype("float32")
    X[ct == "T cell", 0] += 30
    X[ct == "B cell", 1] += 30
    X[ct == "Endothelial", 2] += 30                       # each type expresses its marker -> coherent
    a = anndata.AnnData(X=X)
    a.var_names = genes
    a.obs["cell_type"] = pd.Categorical(ct)
    a.obs["consensus_agreement"] = np.where(ct == "Endothelial", agree_endo, 0.9)
    return a


def test_coherent_but_disputed_is_flagged_as_mislabel():
    from spatialscribe.analysis import trust
    m = {"T cell": ["CD3D"], "B cell": ["MS4A1"], "Endothelial": ["PECAM1"]}
    out = trust.trust_ledger(_section(agree_endo=0.2), marker_sets=m, min_cells=10)
    assert out["has_ensemble"]
    by = {r["cell_type"]: r for r in out["per_type"]}
    # coherent (expresses PECAM1) but the ensemble disputes -> the whole-cluster-mislabel flag AQI can't give
    assert any("DISPUTED" in f for f in by["Endothelial"]["flags"]), by["Endothelial"]
    # coherent AND agreed -> clean (no flags)
    assert by["T cell"]["flags"] == [] and by["B cell"]["flags"] == []
    assert out["n_flagged"] == 1


def test_no_ensemble_disables_agreed_and_says_so():
    from spatialscribe.analysis import trust
    a = _section()
    del a.obs["consensus_agreement"]                      # no reference ensemble voted
    out = trust.trust_ledger(
        a, marker_sets={"T cell": ["CD3D"], "B cell": ["MS4A1"], "Endothelial": ["PECAM1"]}, min_cells=10)
    assert not out["has_ensemble"]
    assert all(r["agreed"] is None for r in out["per_type"])          # 'agreed' column unavailable
    assert "run the reference methods" in out["note"]                 # honest: a coherent mislabel can't be caught
