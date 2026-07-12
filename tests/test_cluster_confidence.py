# tests/test_cluster_confidence.py
"""Cluster-confidence: data-driven over/under-clustering merge/split nudges.

The top rung of the uncertainty ladder (cluster -> cell -> panel): tell the user which
Leiden cluster PAIRS are statistically indistinguishable (merge) and which single clusters
hide substructure (split), grounded in real numbers (RF-vs-permutation p-value; bimodality
coefficient), never a forced call. Reference-free, CPU, no API key.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")


def _ad(embedding, labels):
    """AnnData carrying a PCA embedding + a Leiden labelling (X is a dummy non-negative matrix)."""
    a = anndata.AnnData(X=np.abs(embedding).astype("float32"))
    a.var_names = [f"G{i}" for i in range(embedding.shape[1])]
    a.obs_names = [f"c{i}" for i in range(embedding.shape[0])]
    a.obsm["X_pca"] = np.asarray(embedding, dtype="float32")
    a.obs["leiden"] = pd.Categorical([str(x) for x in labels])
    return a


# --------------------------------------------------------------------------- #
# merge_test
# --------------------------------------------------------------------------- #
def test_merge_test_flags_two_identical_clusters():
    """Two clusters drawn from the SAME distribution can't be separated above chance -> merge."""
    from spatialscribe.analysis import cluster_confidence as ccf

    rng = np.random.default_rng(0)
    E = rng.normal(0, 1, size=(120, 10))
    a = _ad(E, ["0"] * 60 + ["1"] * 60)          # arbitrary split of one blob
    df = ccf.merge_test(a, n_permutations=25)
    row = df.iloc[0]
    assert row["verdict"] == "merge"
    assert row["p_value"] > 0.05
    assert row["rf_accuracy"] < 0.65             # near chance


def test_merge_test_keeps_well_separated_clusters():
    """Two well-separated blobs are perfectly separable -> distinct, never merged."""
    from spatialscribe.analysis import cluster_confidence as ccf

    rng = np.random.default_rng(1)
    E = np.vstack([rng.normal(-6, 0.5, size=(60, 10)), rng.normal(6, 0.5, size=(60, 10))])
    a = _ad(E, ["0"] * 60 + ["1"] * 60)
    df = ccf.merge_test(a, n_permutations=25)
    row = df.iloc[0]
    assert row["verdict"] == "distinct"
    assert row["rf_accuracy"] > 0.9


def test_merge_test_empty_for_single_cluster():
    from spatialscribe.analysis import cluster_confidence as ccf

    rng = np.random.default_rng(2)
    a = _ad(rng.normal(0, 1, size=(40, 10)), ["0"] * 40)
    assert ccf.merge_test(a, n_permutations=10).empty


# --------------------------------------------------------------------------- #
# split_test
# --------------------------------------------------------------------------- #
def test_split_test_flags_bimodal_cluster():
    """A cluster whose cells form two sub-blobs on their main axis -> split candidate."""
    from spatialscribe.analysis import cluster_confidence as ccf

    rng = np.random.default_rng(3)
    sub = np.vstack([rng.normal(-6, 0.5, size=(60, 10)), rng.normal(6, 0.5, size=(60, 10))])
    a = _ad(sub, ["0"] * 120)
    df = ccf.split_test(a)
    row = df[df["cluster"] == "0"].iloc[0]
    assert row["verdict"] == "split"
    assert row["bimodality"] > 0.555


def test_split_test_keeps_unimodal_cluster():
    """A single Gaussian blob is unimodal (BC ~ 0.33) -> ok, not split."""
    from spatialscribe.analysis import cluster_confidence as ccf

    rng = np.random.default_rng(4)
    a = _ad(rng.normal(0, 1, size=(120, 10)), ["0"] * 120)
    df = ccf.split_test(a)
    assert df[df["cluster"] == "0"].iloc[0]["verdict"] == "ok"


def test_split_test_requires_two_programs_when_programs_present():
    """When NMF programs exist, a bimodal cluster carried by ONE program is NOT split (guards
    against over-splitting a single program's spread); two programs confirms the split."""
    from spatialscribe.analysis import cluster_confidence as ccf

    rng = np.random.default_rng(5)
    sub = np.vstack([rng.normal(-6, 0.5, size=(60, 10)), rng.normal(6, 0.5, size=(60, 10))])
    a = _ad(sub, ["0"] * 120)
    a.obs["program"] = pd.Categorical(["Program 0"] * 120)                 # single program
    df = ccf.split_test(a)
    assert df[df["cluster"] == "0"].iloc[0]["verdict"] == "ok"
    a.obs["program"] = pd.Categorical(["Program 0"] * 60 + ["Program 1"] * 60)  # two programs
    df = ccf.split_test(a)
    row = df[df["cluster"] == "0"].iloc[0]
    assert row["verdict"] == "split"
    assert row["n_programs"] == 2


# --------------------------------------------------------------------------- #
# cluster_confidence orchestrator
# --------------------------------------------------------------------------- #
def test_cluster_confidence_panel_merge_split_and_annotation_tie():
    from spatialscribe.analysis import cluster_confidence as ccf

    rng = np.random.default_rng(6)
    c0 = rng.normal(0, 1, size=(60, 10))                                    # cluster 0
    c1 = rng.normal(0, 1, size=(60, 10))                                    # cluster 1 == 0 -> merge
    c2 = np.vstack([rng.normal([-20] + [0] * 9, 0.5, size=(30, 10)),
                    rng.normal([20] + [0] * 9, 0.5, size=(30, 10))])        # cluster 2 bimodal + far -> split
    E = np.vstack([c0, c1, c2])
    a = _ad(E, ["0"] * 60 + ["1"] * 60 + ["2"] * 60)
    # Annotation: clusters 0 and 1 got the SAME label (strengthens the merge nudge).
    a.obs["cell_type"] = pd.Categorical(["T cell"] * 120 + ["Other"] * 60)

    out = ccf.cluster_confidence(a, n_permutations=25)
    merges = {tuple(sorted((m["cluster_a"], m["cluster_b"]))) for m in out["merge_suggestions"]}
    assert ("0", "1") in merges
    m01 = next(m for m in out["merge_suggestions"]
               if sorted((m["cluster_a"], m["cluster_b"])) == ["0", "1"])
    assert m01["same_label"] is True                     # both annotated 'T cell'
    assert "2" in {s["cluster"] for s in out["split_suggestions"]}
    assert out["n_clusters"] == 3
    assert isinstance(out["nudge"], str) and out["nudge"]


# --------------------------------------------------------------------------- #
# capability wiring
# --------------------------------------------------------------------------- #
def test_cluster_confidence_capability_runs_and_is_copilot_exposed():
    from spatialscribe.analysis import capabilities as cap

    rng = np.random.default_rng(7)
    E = rng.normal(0, 1, size=(120, 10))
    a = _ad(E, ["0"] * 60 + ["1"] * 60)
    res = cap.run(a, "cluster_confidence", {}, cap.RunContext())
    assert res.ok, res.error
    assert "merge_suggestions" in res.value and "split_suggestions" in res.value
    assert "cluster_confidence" in cap.copilot_names()


def test_cluster_confidence_capability_needs_leiden():
    """Prerequisite gate: without a Leiden clustering the capability returns a structured error."""
    from spatialscribe.analysis import capabilities as cap

    a = anndata.AnnData(X=np.abs(np.random.default_rng(8).normal(size=(20, 5))).astype("float32"))
    res = cap.run(a, "cluster_confidence", {}, cap.RunContext())
    assert not res.ok
    assert res.error["error_type"] == "prerequisite_missing"
