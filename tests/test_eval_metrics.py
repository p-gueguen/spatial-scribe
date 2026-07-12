"""Tests for analysis/eval_metrics.py - the annotation-quality metric battery.

Covers the three regimes from docs/research/cell-annotation-quality-metrics.md:
  * internal (reference-free): scTypeEval cluster-validity family + inter-sample consistency;
  * marker-program fidelity (reference-free, needs marker sets): AUC-ROC / Cohen's d;
  * external (needs ground truth): F1/ARI/kappa, ECS, hierarchical & composition accuracy.

Fixtures are deterministic well-separated blobs (clean) vs shuffled labels (noise); every metric
must score the clean labeling above the shuffled one, which is the whole point of the module.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")


def _blobs(n_per=80, k=4, dim=10, sep=8.0, seed=0):
    """k well-separated Gaussian blobs in a 10-d embedding, correctly labeled T0..T{k-1}."""
    rng = np.random.default_rng(seed)
    centers = rng.normal(0, sep, size=(k, dim))
    X, lab = [], []
    for t in range(k):
        X.append(rng.normal(centers[t], 1.0, size=(n_per, dim)))
        lab += [f"T{t}"] * n_per
    X = np.vstack(X).astype("float32")
    a = anndata.AnnData(X=X)
    a.obs_names = [f"c{i}" for i in range(a.n_obs)]
    a.obsm["X_pca"] = X
    a.obs["cell_type"] = pd.Categorical(lab)
    return a


def _shuffle_labels(a, seed=1):
    a2 = a.copy()
    rng = np.random.default_rng(seed)
    a2.obs["cell_type"] = pd.Categorical(rng.permutation(a2.obs["cell_type"].to_numpy()))
    return a2


# --------------------------------------------------------------------------- #
# Internal cluster-validity (scTypeEval family)
# --------------------------------------------------------------------------- #
def test_internal_validity_high_for_clean_low_for_shuffled():
    from spatialscribe.analysis import eval_metrics as em

    a = _blobs()
    clean = em.internal_validity(a, label_key="cell_type", embedding="X_pca", k=15)
    shuf = em.internal_validity(_shuffle_labels(a), label_key="cell_type", embedding="X_pca", k=15)

    for m in ("silhouette", "silhouette_2label", "neighborhood_purity",
              "orbital_medoid", "ward_propmatch", "integrated"):
        assert clean[m] > shuf[m], (m, clean[m], shuf[m])
    assert clean["neighborhood_purity"] > 0.9      # clean blobs: neighbors share label
    assert 0.0 <= clean["integrated"] <= 1.0


def test_internal_validity_computes_pca_when_no_embedding():
    from spatialscribe.analysis import eval_metrics as em

    a = _blobs()
    del a.obsm["X_pca"]                              # force the PCA fallback
    out = em.internal_validity(a, label_key="cell_type", embedding="X_pca", k=15)
    assert out["neighborhood_purity"] > 0.9


# --------------------------------------------------------------------------- #
# Inter-sample consistency (ISC)
# --------------------------------------------------------------------------- #
def _isc_adata(consistent: bool, n_per=40, n_types=3, n_samples=4, n_genes=30, seed=0):
    """Each (sample, type) is a pseudobulk. Consistent: a type's signature block is the SAME
    across samples. Inconsistent: each sample uses a DIFFERENT block for a given type."""
    rng = np.random.default_rng(seed)
    block = n_genes // n_types
    X, types, samples = [], [], []
    for s in range(n_samples):
        for t in range(n_types):
            sig = t if consistent else rng.integers(0, n_types)  # which block is high
            base = rng.normal(0.2, 0.1, size=(n_per, n_genes))
            base[:, sig * block:(sig + 1) * block] += 5.0
            X.append(base)
            types += [f"T{t}"] * n_per
            samples += [f"S{s}"] * n_per
    a = anndata.AnnData(X=np.vstack(X).astype("float32"))
    a.obs_names = [f"c{i}" for i in range(a.n_obs)]
    a.obs["cell_type"] = pd.Categorical(types)
    a.obs["sample"] = pd.Categorical(samples)
    return a


def test_inter_sample_consistency_high_for_reproducible_types():
    from spatialscribe.analysis import eval_metrics as em

    good = em.inter_sample_consistency(_isc_adata(consistent=True), "cell_type", "sample")
    bad = em.inter_sample_consistency(_isc_adata(consistent=False), "cell_type", "sample")
    assert good["consistency"] > bad["consistency"]


# --------------------------------------------------------------------------- #
# Marker-program fidelity (reference-free, needs marker sets)
# --------------------------------------------------------------------------- #
def _marker_adata(correct=True, seed=0):
    genes = ["CD3D", "CD3E", "TRAC", "MLANA", "SOX10", "DCT"] + [f"N{i}" for i in range(10)]
    gi = {g: i for i, g in enumerate(genes)}
    n = 120
    rng = np.random.default_rng(seed)
    X = rng.poisson(0.2, size=(n, len(genes))).astype("float32")
    for g in ("CD3D", "CD3E", "TRAC"):
        X[:60, gi[g]] += 10
    for g in ("MLANA", "SOX10", "DCT"):
        X[60:, gi[g]] += 10
    a = anndata.AnnData(X=X)
    a.var_names = genes
    a.obs_names = [f"c{i}" for i in range(n)]
    labels = ["T cell"] * 60 + ["Mel"] * 60
    if not correct:
        labels = list(rng.permutation(labels))
    a.obs["cell_type"] = pd.Categorical(labels)
    return a


def test_marker_program_fidelity_high_auc_for_correct_labels():
    from spatialscribe.analysis import eval_metrics as em

    markers = {"T cell": ["CD3D", "CD3E", "TRAC"], "Mel": ["MLANA", "SOX10", "DCT"]}
    good = em.marker_program_fidelity(_marker_adata(correct=True), "cell_type", markers)
    bad = em.marker_program_fidelity(_marker_adata(correct=False), "cell_type", markers)
    assert good["mean_auc"] > 0.9
    assert good["mean_auc"] > bad["mean_auc"]
    assert good["per_type"]["T cell"]["auc"] > 0.9
    assert good["per_type"]["T cell"]["cohens_d"] > 0.5           # large, positive separation


# --------------------------------------------------------------------------- #
# External harness (needs ground truth)
# --------------------------------------------------------------------------- #
def test_external_scores_perfect():
    from spatialscribe.analysis import eval_metrics as em

    truth = np.array(["A"] * 90 + ["B"] * 10)
    out = em.external_scores(truth.copy(), truth)
    for m in ("accuracy", "balanced_accuracy", "macro_f1", "ari", "ami", "kappa", "ecs"):
        assert out[m] == pytest.approx(1.0), (m, out[m])


def test_external_scores_macro_f1_penalizes_rare_miss():
    from spatialscribe.analysis import eval_metrics as em

    truth = np.array(["A"] * 90 + ["B"] * 10)
    pred = np.array(["A"] * 100)                       # miss every rare-B cell
    out = em.external_scores(pred, truth)
    assert out["accuracy"] == pytest.approx(0.9)       # inflated by the majority class
    assert out["macro_f1"] < out["accuracy"]           # macro-F1 exposes the dropped class
    assert out["balanced_accuracy"] < 0.9
    assert out["per_class_f1"]["B"] == pytest.approx(0.0)


# --------------------------------------------------------------------------- #
# Element-Centric Similarity (ECS)
# --------------------------------------------------------------------------- #
def test_ecs_identical_is_one_and_resolution_aware():
    from spatialscribe.analysis import eval_metrics as em

    a = np.array([0, 0, 1, 1, 2, 2])
    assert em.element_centric_similarity(a, a) == pytest.approx(1.0)
    allsame = np.zeros(6, dtype=int)
    refine = np.array([0, 3, 1, 1, 2, 2])              # split cluster 0 into two singletons
    assert em.element_centric_similarity(a, allsame) < 1.0
    # a near-refinement is more similar to `a` than collapsing everything to one cluster
    assert em.element_centric_similarity(a, refine) > em.element_centric_similarity(a, allsame)


# --------------------------------------------------------------------------- #
# Hierarchical & composition accuracy
# --------------------------------------------------------------------------- #
def test_hierarchical_accuracy_partial_credit():
    from spatialscribe.analysis import eval_metrics as em

    hier = {"CD4 T": "T", "CD8 T": "T", "B": "B"}
    truth = np.array(["CD4 T", "CD8 T", "B"])
    assert em.hierarchical_accuracy(truth.copy(), truth, hier)["hierarchical_accuracy"] == pytest.approx(1.0)

    partial = em.hierarchical_accuracy(np.array(["CD8 T", "CD8 T", "B"]), truth, hier, partial=0.5)
    assert partial["hierarchical_accuracy"] == pytest.approx((0.5 + 1 + 1) / 3)   # right lineage, wrong subtype
    assert partial["subtype_accuracy"] == pytest.approx(2 / 3)

    wrong = em.hierarchical_accuracy(np.array(["B", "CD8 T", "B"]), truth, hier, partial=0.5)
    assert wrong["hierarchical_accuracy"] == pytest.approx((0 + 1 + 1) / 3)       # wrong lineage
    assert partial["hierarchical_accuracy"] > wrong["hierarchical_accuracy"]


def test_composition_accuracy():
    from spatialscribe.analysis import eval_metrics as em

    truth = np.array(["A"] * 60 + ["B"] * 40)
    same = em.composition_accuracy(truth.copy(), truth)
    assert same["l1"] == pytest.approx(0.0, abs=1e-9)
    assert same["jsd"] == pytest.approx(0.0, abs=1e-9)
    diff = em.composition_accuracy(np.array(["A"] * 90 + ["B"] * 10), truth)
    assert diff["l1"] > 0


def test_deconvolution_metrics_perfect_and_degraded():
    """OpenProblems spatial-decomposition metrics (R2 uniform-average + JSD axis=0) + RMSE on
    predicted vs true cell-type PROPORTION matrices (spots x types)."""
    from spatialscribe.analysis import eval_metrics as em

    rng = np.random.default_rng(0)
    true = rng.dirichlet(np.ones(4), size=60)                 # 60 spots x 4 types, rows sum to 1
    perfect = em.deconvolution_metrics(true, true.copy())
    assert perfect["r2"] == pytest.approx(1.0)
    assert perfect["jsd"] == pytest.approx(0.0, abs=1e-9)
    assert perfect["jsd_per_spot"] == pytest.approx(0.0, abs=1e-9)
    assert perfect["rmse"] == pytest.approx(0.0, abs=1e-9)
    assert perfect["pearson"] == pytest.approx(1.0)

    noisy = true + rng.uniform(0, 0.3, size=true.shape)
    noisy /= noisy.sum(1, keepdims=True)                      # still valid proportions
    deg = em.deconvolution_metrics(true, noisy)
    assert deg["r2"] < 1.0
    assert deg["jsd"] > 0 and deg["rmse"] > 0


# --------------------------------------------------------------------------- #
# Orchestrator + QC-funnel wiring
# --------------------------------------------------------------------------- #
def test_annotation_quality_orchestrator():
    from spatialscribe.analysis import eval_metrics as em

    a = _blobs()
    a.var["control"] = False
    out = em.annotation_quality(a, label_key="cell_type", embedding="X_pca")
    assert "internal_validity" in out
    assert out["internal_validity"]["neighborhood_purity"] > 0.9


def test_run_funnel_includes_annotation_quality():
    """The reference-free battery is wired into the QC-funnel headline when labels are present."""
    from spatialscribe.analysis import markers, qc

    a = _blobs(n_per=40, k=3)
    a.var_names = [f"g{i}" for i in range(a.n_vars)]
    a.var["control"] = False
    a.obsm["spatial"] = np.random.default_rng(0).uniform(0, 100, size=(a.n_obs, 2))
    head = qc.run_funnel(a, cluster_key="cell_type", marker_sets=markers.LINEAGE_MARKERS)
    assert "annotation_quality" in head
    assert "internal_validity" in head["annotation_quality"]


def test_run_funnel_aqi_matches_mouse_titlecase_panel():
    """Regression: curated marker sets are human-UPPERCASE, but a mouse panel ships title-case symbols
    (Gfap vs GFAP). run_funnel must re-key markers to the panel's casing (markers.on_panel) so the
    AQI's purity (C) + marker-fidelity (M) still compute - otherwise every marker reads off-panel,
    no_dict_frac hits 1.0, and the AQI silently degrades to the internal-validity-only fallback (the
    observed failure on a real Xenium Mouse Brain section)."""
    from spatialscribe.analysis import qc

    genes = ["Cd3d", "Cd3e", "Trac", "Mlana", "Sox10", "Dct"] + [f"N{i}" for i in range(14)]  # mouse title-case panel
    gi = {g: i for i, g in enumerate(genes)}
    n = 160
    rng = np.random.default_rng(0)
    X = rng.poisson(0.3, size=(n, len(genes))).astype("float32")
    for g in ("Cd3d", "Cd3e", "Trac"):
        X[:80, gi[g]] += 12
    for g in ("Mlana", "Sox10", "Dct"):
        X[80:, gi[g]] += 12
    X[:, 6:] += rng.poisson(rng.uniform(0.5, 15.0, len(genes) - 6), size=(n, len(genes) - 6)).astype("float32")
    a = anndata.AnnData(X=X)
    a.var_names = genes
    a.var["control"] = False
    a.obs_names = [f"c{i}" for i in range(n)]
    a.obs["cell_type"] = pd.Categorical(["T cell"] * 80 + ["Mel"] * 80)
    a.obsm["spatial"] = rng.uniform(0, 100, size=(n, 2))
    markers_up = {"T cell": ["CD3D", "CD3E", "TRAC"], "Mel": ["MLANA", "SOX10", "DCT"]}   # human UPPERCASE (as ctx.markers returns)
    head = qc.run_funnel(a, cluster_key="cell_type", marker_sets=markers_up)
    aqi = head["annotation_quality"]["aqi"]
    # markers matched via case-fold -> the marker-based terms populate (both were None on a mouse panel before)
    assert aqi["components"]["M"] is not None, aqi
    assert aqi["components"]["C"] is not None, aqi
    assert (aqi.get("no_dict_frac") or 0) < 0.5, aqi


def _panel_ref(n=120, n_genes=40, seed=0):
    """Reference where A/B have distinct on-panel markers but C/D differ ONLY off-panel (last gene)."""
    rng = np.random.default_rng(seed)
    off = n_genes - 1
    blocks, labels = [], []
    for t, hi in [("A", 0), ("B", 1), ("C", 2), ("D", 2)]:      # C,D share on-panel marker g2
        M = rng.poisson(1.0, size=(n, n_genes)).astype("float32")
        M[:, hi] += rng.poisson(20.0, size=n)
        if t == "D":
            M[:, off] += rng.poisson(20.0, size=n)              # D's only distinct gene is OFF-panel
        blocks.append(M); labels += [t] * n
    a = anndata.AnnData(X=np.vstack(blocks))
    a.var_names = [f"g{i}" for i in range(n_genes)]
    a.obs["cell_type"] = labels
    return a


def test_panel_resolvability_flags_off_panel_confusable_pair():
    """A panel that omits a type-pair's only distinguishing gene cannot resolve that pair."""
    from spatialscribe.analysis import eval_metrics as em

    a = _panel_ref()
    panel = {f"g{i}" for i in range(39)}                         # 39 shared genes; excludes g39 (D's marker)
    res = em.panel_resolvability(a, "cell_type", panel, target_depth=40)
    assert res["status"] == "ok" and res["n_types"] == 4
    assert res["per_type"]["A"]["tier"] == "resolvable"         # private on-panel markers -> resolvable
    assert res["per_type"]["B"]["tier"] == "resolvable"
    cd = {res["per_type"]["C"]["tier"], res["per_type"]["D"]["tier"]}
    assert "not_resolvable" in cd                                # identical on the panel -> confusable
    assert (res["per_type"]["C"]["confused_with"] == "D"
            or res["per_type"]["D"]["confused_with"] == "C")
    assert res["frac_resolvable"] < 1.0


def test_panel_resolvability_insufficient_overlap():
    """A reference sharing < 25 genes with the panel returns an honest skip, not a bogus score."""
    from spatialscribe.analysis import eval_metrics as em

    a = _panel_ref()
    res = em.panel_resolvability(a, "cell_type", {"zzz1", "zzz2"}, target_depth=40)
    assert res["status"] == "insufficient_overlap"
