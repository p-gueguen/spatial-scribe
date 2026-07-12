# tests/test_qc_layers.py
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")


def _tiny(n=200, n_genes=40, seed=0):
    """Small AnnData with spatial coords + segmentation columns for the Layer tests."""
    rng = np.random.default_rng(seed)
    X = rng.poisson(0.5, size=(n, n_genes)).astype("float32")
    a = anndata.AnnData(X=X)
    a.var_names = [f"G{i}" for i in range(n_genes)]
    a.obs_names = [f"c{i}" for i in range(n)]
    a.obs["cell_area"] = rng.uniform(20, 120, n)
    a.obs["nucleus_area"] = a.obs["cell_area"].to_numpy() * rng.uniform(0.2, 0.6, n)
    a.obs.loc[a.obs_names[0], "nucleus_area"] = 0.0  # one anucleate cell
    a.obsm["spatial"] = rng.uniform(0, 100, size=(n, 2))
    a.var["control"] = False
    return a


def test_qc_summary_control_pct_survives_empty_cells():
    """A single empty cell must not NaN out the section's neg-control metric.

    ``pct_counts_control`` is 0/0 -> NaN for a zero-count cell, and ``np.median`` of any array
    containing NaN is NaN. Every downstream ``metric > threshold`` comparison against NaN is then
    False, so the neg-control flag SILENTLY passes. Measured on the 2026-07 overnight benchmark:
    section ``pct_counts_control`` was NaN on 4 of 5 real sections while the per-cell median was
    ~15% - "input QC ok" on a section whose annotation was unusable.
    """
    from spatialscribe.analysis import qc

    a = _tiny(n=50)
    a.obs["total_counts"] = np.r_[0.0, np.full(49, 50.0)]        # one empty cell
    a.obs["n_genes_by_counts"] = np.r_[0.0, np.full(49, 20.0)]
    a.obs["pct_counts_control"] = np.r_[np.nan, np.full(49, 1.5)]

    out = qc.qc_summary(a)
    assert not np.isnan(out["pct_counts_control"]), "NaN silently disables the neg-control flag"
    assert abs(out["pct_counts_control"] - 1.5) < 1e-9


def test_annotation_completeness_flags_resolvable_but_absent_lineage():
    # A lineage the PANEL can resolve (green coverage) that annotates to 0 cells is surfaced - so a
    # "usability ok" headline cannot hide an annotation that recovered only some resolvable lineages
    # (the CosMx case: green T-cell coverage, 0 T cells annotated). Present + non-green types excluded.
    import anndata as ad
    import numpy as np
    import pandas as pd

    from spatialscribe.analysis import annotate

    a = ad.AnnData(X=np.zeros((6, 2), dtype="float32"))
    a.obs["cell_type"] = pd.Categorical(["T cell", "T cell", "Myeloid", "Myeloid", "Myeloid", "T cell"])
    a.uns["panel_check"] = {
        "coverage": {"T cell": {"status": "green", "n_present": 4, "n_markers": 7},
                     "Myeloid": {"status": "green", "n_present": 4, "n_markers": 7},
                     "Endothelial": {"status": "green", "n_present": 4, "n_markers": 4},   # resolvable, 0 cells
                     "Mast": {"status": "amber", "n_present": 1, "n_markers": 4}},          # not green -> ignored
        "confusable_pairs": [], "merge_groups": []}
    comp = annotate.annotation_completeness(a, "cell_type")
    assert comp["resolvable_absent"] == ["Endothelial"]   # green + absent only; present/amber excluded


def test_segmentation_qc_writes_columns():
    from spatialscribe.analysis import qc

    a = _tiny()
    out = qc.segmentation_qc(a)
    assert out["status"] == "ok"
    assert {"nucleus_present", "nucleus_to_cell_ratio", "seg_area_flag"} <= set(a.obs.columns)
    assert bool(a.obs["nucleus_present"].iloc[0]) is False  # the anucleate cell
    assert 0.0 <= out["pct_area_outlier"] <= 1.0


def test_segmentation_qc_skips_without_columns():
    from spatialscribe.analysis import qc

    a = _tiny()
    del a.obs["cell_area"]
    out = qc.segmentation_qc(a)
    assert out["status"] == "skipped"


def test_count_floor_is_fixed_on_small_panel():
    from spatialscribe.analysis import qc

    a = _tiny(n=300, n_genes=40)
    qc.compute_qc(a)
    out = qc.suggest_count_floor(a)
    assert out["mode"] == "fixed"
    assert out["floor"] == 10


def test_count_floor_is_distributional_on_large_panel():
    from spatialscribe.analysis import qc

    a = _tiny(n=300, n_genes=40)
    qc.compute_qc(a)
    # Pretend this is a 5K-panel section: index by panel size, not gene count in X.
    out = qc.suggest_count_floor(a, n_panel_genes=5000)
    assert out["mode"] == "distributional"
    assert 0.0 <= out["pct_removed"] <= 0.10  # a percentile floor removes a bounded fraction


def test_pmp_high_for_marker_pure_cells():
    from spatialscribe.analysis import purity

    rng = np.random.default_rng(1)
    genes = ["CD3D", "CD3E", "TRAC", "MLANA", "SOX10", "COL1A1"]
    n = 60
    X = rng.poisson(0.2, size=(n, len(genes))).astype("float32")
    # First 30 cells are T cells with strong T-marker counts and clean background.
    gi = {g: i for i, g in enumerate(genes)}
    for g in ("CD3D", "CD3E", "TRAC"):
        X[:30, gi[g]] += 10
    a = anndata.AnnData(X=X)
    a.var_names = genes
    a.obs_names = [f"c{i}" for i in range(n)]
    a.obs["cell_type"] = pd.Categorical(["T cell"] * 30 + ["Stromal/CAF"] * 30)
    lineage = {"T cell": ["CD3D", "CD3E", "TRAC"], "Stromal/CAF": ["COL1A1"]}
    purity.pmp(a, assigned_label_key="cell_type", lineage_markers=lineage)
    assert "pmp" in a.obs
    assert a.obs["pmp"].iloc[:30].mean() > 0.5  # T cells are marker-pure


def test_spatial_coherence_high_for_segregated_labels():
    from spatialscribe.analysis import spatial

    rng = np.random.default_rng(2)
    n = 200
    X = rng.poisson(0.5, size=(n, 10)).astype("float32")
    a = anndata.AnnData(X=X)
    a.var_names = [f"G{i}" for i in range(10)]
    a.obs_names = [f"c{i}" for i in range(n)]
    grp = np.array([0] * 100 + [1] * 100)
    # Two spatially separated blobs => neighbors share labels => high coherence.
    a.obsm["spatial"] = np.column_stack([rng.normal(np.where(grp == 0, 0, 60), 3), rng.normal(0, 3, n)])
    a.obs["cell_type"] = pd.Categorical(np.where(grp == 0, "A", "B"))
    out = spatial.spatial_coherence(a, label_key="cell_type", k=10)
    assert "spatial_coherence" in a.obs
    assert out["mean_coherence"] > 0.8
    assert 0.0 <= out["pas"] <= 1.0


def test_neighborhood_sanity_reports_self_enrichment():
    from spatialscribe.analysis import spatial

    rng = np.random.default_rng(3)
    n = 240
    X = rng.poisson(0.5, size=(n, 10)).astype("float32")
    a = anndata.AnnData(X=X)
    a.var_names = [f"G{i}" for i in range(10)]
    a.obs_names = [f"c{i}" for i in range(n)]
    grp = np.array([0] * 120 + [1] * 120)
    a.obsm["spatial"] = np.column_stack([rng.normal(np.where(grp == 0, 0, 60), 3), rng.normal(0, 3, n)])
    a.obs["cell_type"] = pd.Categorical(np.where(grp == 0, "A", "B"))
    out = spatial.neighborhood_sanity(a, cluster_key="cell_type")
    assert set(out["categories"]) == {"A", "B"}
    assert len(out["self_enrichment"]) == 2
    # Well-separated blobs self-associate => nothing suspicious.
    assert out["suspicious"] == []


def test_annotation_stability_low_for_strong_markers():
    from spatialscribe.analysis import annotate

    rng = np.random.default_rng(4)
    genes = ["CD3D", "CD3E", "TRAC", "MLANA", "SOX10", "DCT"] + [f"N{i}" for i in range(10)]
    n = 80
    gi = {g: i for i, g in enumerate(genes)}
    X = rng.poisson(0.2, size=(n, len(genes))).astype("float32")
    for g in ("CD3D", "CD3E", "TRAC"):
        X[:40, gi[g]] += 15
    for g in ("MLANA", "SOX10", "DCT"):
        X[40:, gi[g]] += 15
    a = anndata.AnnData(X=X)
    a.var_names = genes
    a.obs_names = [f"c{i}" for i in range(n)]
    a.layers["counts"] = a.X.copy()
    lineage = {"T cell": ["CD3D", "CD3E", "TRAC"], "Malignant/Melanocyte": ["MLANA", "SOX10", "DCT"]}
    annotate.annotation_stability(a, marker_sets=lineage, reps=3, seed=0)
    assert "annotation_stability" in a.obs
    # Strong, well-separated markers => labels rarely flip under 20% dropout.
    assert a.obs["annotation_stability"].mean() < 0.3


def test_apply_confidence_survives_corrupted_panel_check():
    """Regression: an h5ad round-trip splits '/'-keyed cell types so coverage values lose their
    'status' field; apply_confidence must guard it (not KeyError) and never return NaN confidence."""
    import math

    from spatialscribe.analysis import annotate

    rng = np.random.default_rng(7)
    genes = ["CD3D", "CD3E", "TRAC", "MLANA", "SOX10", "DCT"] + [f"N{i}" for i in range(14)]
    gi = {g: i for i, g in enumerate(genes)}
    n = 60
    X = rng.poisson(0.2, size=(n, len(genes))).astype("float32")
    for g in ("CD3D", "CD3E", "TRAC"):
        X[:30, gi[g]] += 12
    for g in ("MLANA", "SOX10", "DCT"):
        X[30:, gi[g]] += 12
    a = anndata.AnnData(X=X)
    a.var_names = genes
    a.obs_names = [f"c{i}" for i in range(n)]
    a.var["control"] = False
    a.obs["total_counts"] = np.asarray(X.sum(1)).ravel()
    a.obs["n_genes_by_counts"] = (X > 0).sum(1)
    a.obs["cell_type"] = pd.Categorical(["T cell"] * 30 + ["Malignant/Melanocyte"] * 30)
    lineage = {"T cell": ["CD3D", "CD3E", "TRAC"], "Malignant/Melanocyte": ["MLANA", "SOX10", "DCT"]}
    # Corrupted panel_check: coverage dicts WITHOUT the 'status' key (the round-trip failure mode).
    bad_pc = {"coverage": {"T cell": {"n_present": 3}, "Malignant/Melanocyte": {"nested": {"x": 1}}}}
    out = annotate.apply_confidence(a, cluster_key="cell_type", marker_sets=lineage, panel_check_result=bad_pc)
    assert abs(out["pct_pass"] + out["pct_warn"] + out["pct_abstain"] - 1.0) < 1e-6
    assert not math.isnan(out["mean_confidence"])


def test_apply_confidence_uses_optional_penalties():
    from spatialscribe.analysis import annotate, markers

    rng = np.random.default_rng(5)
    markers_ = ["CD3D", "CD3E", "TRAC", "MLANA", "SOX10", "DCT"]
    noise = [f"NOISE{i}" for i in range(60)]
    genes = markers_ + noise
    gi = {g: i for i, g in enumerate(genes)}
    n = 60
    X = rng.poisson(0.2, size=(n, len(genes))).astype("float32")
    for g in ("CD3D", "CD3E", "TRAC"):
        X[:30, gi[g]] += 12
    for g in ("MLANA", "SOX10", "DCT"):
        X[30:, gi[g]] += 12
    a = anndata.AnnData(X=X)
    a.var_names = genes
    a.obs_names = [f"c{i}" for i in range(n)]
    a.var["control"] = False
    a.obs["total_counts"] = np.asarray(X.sum(1)).ravel()
    a.obs["n_genes_by_counts"] = (X > 0).sum(1)
    a.obs["cell_type"] = pd.Categorical(["T cell"] * 30 + ["Malignant/Melanocyte"] * 30)
    lineage = {"T cell": ["CD3D", "CD3E", "TRAC"], "Malignant/Melanocyte": ["MLANA", "SOX10", "DCT"]}

    # Baseline (no optional columns) still returns a valid summary.
    base = annotate.apply_confidence(a, cluster_key="cell_type", marker_sets=lineage)
    assert abs(base["pct_pass"] + base["pct_warn"] + base["pct_abstain"] - 1.0) < 1e-6

    # Inject a low-coherence column for half the cells, corroborated by a VSI doublet flag, and
    # confirm mean confidence drops. Spatial incoherence only penalizes when an independent
    # mixed-signal cue is present (see test_incoherence_alone_does_not_penalize_a_clean_cell);
    # without the VSI corroboration these clean cells would keep full confidence.
    a.obs["spatial_coherence"] = np.r_[np.full(30, 0.05), np.full(30, 0.95)]
    a.obs["vsi"] = np.r_[np.full(30, 0.3), np.full(30, 0.95)]
    penalized = annotate.apply_confidence(a, cluster_key="cell_type", marker_sets=lineage)
    assert penalized["mean_confidence"] < base["mean_confidence"]


def test_incoherence_alone_does_not_penalize_a_clean_cell():
    """A clean, high-signal cell (strong markers, low contamination, no doublet) must KEEP its
    confidence even when spatially isolated: a lone infiltrating T cell in a tumor nest is real
    biology, not an annotation error. Spatial incoherence is a CORROBORATOR - it lowers
    confidence only when an independent mixed-signal cue (CRISP impurity or an ovrlpy VSI
    vertical-overlap doublet) is ALSO present. Regression for the immune-infiltration false
    penalty (location is a prior, not evidence about a cell's transcriptomic identity)."""
    from spatialscribe.analysis import annotate

    rng = np.random.default_rng(11)
    markers_ = ["CD3D", "CD3E", "TRAC", "MLANA", "SOX10", "DCT"]
    noise = [f"NOISE{i}" for i in range(60)]
    genes = markers_ + noise
    gi = {g: i for i, g in enumerate(genes)}
    n = 60
    X = rng.poisson(0.2, size=(n, len(genes))).astype("float32")
    # Make the two populations genuinely CRISP-pure: each expresses ONLY its own lineage's
    # markers (zero the other lineage's), so "clean" really means clean and no cell is flagged
    # impure by Poisson noise - the test then isolates the coherence-corroborator behaviour.
    for g in ("CD3D", "CD3E", "TRAC"):
        X[:30, gi[g]] += 12
        X[30:, gi[g]] = 0
    for g in ("MLANA", "SOX10", "DCT"):
        X[30:, gi[g]] += 12
        X[:30, gi[g]] = 0

    lineage = {"T cell": ["CD3D", "CD3E", "TRAC"], "Malignant/Melanocyte": ["MLANA", "SOX10", "DCT"]}

    def _mk():
        a = anndata.AnnData(X=X.copy())
        a.var_names = genes
        a.obs_names = [f"c{i}" for i in range(n)]
        a.var["control"] = False
        a.obs["total_counts"] = np.asarray(X.sum(1)).ravel()
        a.obs["n_genes_by_counts"] = (X > 0).sum(1)
        a.obs["cell_type"] = pd.Categorical(["T cell"] * 30 + ["Malignant/Melanocyte"] * 30)
        return a

    base = annotate.apply_confidence(_mk(), cluster_key="cell_type", marker_sets=lineage)

    # (A) The 30 T cells are spatially isolated (coherence ~0) but clean -> NOT penalized.
    a_clean = _mk()
    a_clean.obs["spatial_coherence"] = np.r_[np.full(30, 0.03), np.full(30, 0.97)]
    clean = annotate.apply_confidence(a_clean, cluster_key="cell_type", marker_sets=lineage)
    assert clean["mean_confidence"] == pytest.approx(base["mean_confidence"])

    # (B) Same low coherence, but now corroborated by an ovrlpy VSI doublet -> confidence drops.
    a_corrob = _mk()
    a_corrob.obs["spatial_coherence"] = np.r_[np.full(30, 0.03), np.full(30, 0.97)]
    a_corrob.obs["vsi"] = np.r_[np.full(30, 0.3), np.full(30, 0.95)]
    corrob = annotate.apply_confidence(a_corrob, cluster_key="cell_type", marker_sets=lineage)
    assert corrob["mean_confidence"] < base["mean_confidence"]


def test_run_funnel_returns_headline():
    from spatialscribe.analysis import qc

    a = _tiny(n=150, n_genes=40)
    qc.compute_qc(a)
    out = qc.run_funnel(a)  # no cell_type yet -> pre-annotation layers only
    assert out["n_cells"] == a.n_obs
    assert "section" in out and "count_floor" in out and "segmentation" in out
    assert "layers_run" in out


def test_run_funnel_needs_tissue_matched_markers():
    """Regression for the qc_deepdive tool: run_funnel must be called with the section's own
    tissue marker set, not silently fall back to the melanoma default. On an epithelial
    (breast/carcinoma) section, the melanoma dict has no 'Epithelial/Tumor' key, so
    purity.pmp leaves those cells at pmp=0 and annotate.apply_confidence halves their
    confidence for no biological reason - exactly the bug qc_deepdive had before it threaded
    ``markers.for_tissue(tissue)`` into this call.
    """
    from spatialscribe.analysis import markers, qc

    rng = np.random.default_rng(6)
    epi_markers = markers.EPITHELIAL_LINEAGES["Epithelial/Tumor"]
    tcell_markers = markers.EPITHELIAL_LINEAGES["T cell"]
    noise = [f"NOISE{i}" for i in range(60)]
    genes = epi_markers + tcell_markers + noise
    gi = {g: i for i, g in enumerate(genes)}
    n = 60
    X = rng.poisson(0.2, size=(n, len(genes))).astype("float32")
    for g in epi_markers:
        X[:30, gi[g]] += 12
    for g in tcell_markers:
        X[30:, gi[g]] += 12

    a = anndata.AnnData(X=X)
    a.var_names = genes
    a.obs_names = [f"c{i}" for i in range(n)]
    a.var["control"] = False
    a.obs["cell_type"] = pd.Categorical(["Epithelial/Tumor"] * 30 + ["T cell"] * 30)
    a.obsm["spatial"] = rng.uniform(0, 100, size=(n, 2))

    a_epi = a.copy()
    out_epi = qc.run_funnel(a_epi, cluster_key="cell_type", marker_sets=markers.EPITHELIAL_LINEAGES)

    a_mel = a.copy()
    out_mel = qc.run_funnel(a_mel, cluster_key="cell_type", marker_sets=markers.LINEAGE_MARKERS)

    # Tissue-correct markers resolve the epithelial compartment (finite, panel-invariant PMP);
    # the melanoma default cannot (no 'Epithelial/Tumor' key), so PMP is undefined -> NaN for
    # those cells - still the smoking gun. (NaN is now treated as neutral, so mel confidence
    # rides on the weak melanoma margin, not a spurious pmp=0 penalty; the gap below persists.)
    assert (a_epi.obs.loc[a_epi.obs["cell_type"] == "Epithelial/Tumor", "pmp"] > 0).any()
    assert a_mel.obs.loc[a_mel.obs["cell_type"] == "Epithelial/Tumor", "pmp"].isna().all()
    assert out_epi["annotatability"]["mean_confidence"] > out_mel["annotatability"]["mean_confidence"]

    # Now exercise the actual buggy code path: the qc_deepdive agent tool must thread
    # markers.for_tissue(tissue) into run_funnel, not call it bare (which silently falls
    # back to the melanoma LINEAGE_MARKERS default regardless of the section's tissue).
    # A tissue="breast" section run through the real tool should match the epithelial-marker
    # result above, not the melanoma-default one.
    # The copilot's qc_deepdive is now the registry capability 'qc_funnel'; dispatched via
    # capabilities.run with a RunContext(tissue=...), which threads markers.for_tissue(tissue)
    # into run_funnel (replacing the old private agent.tools._run_tool path).
    from spatialscribe.analysis import capabilities as cap

    a_tool = a.copy()
    out_tool = cap.run(a_tool, "qc_funnel", {}, cap.RunContext(tissue="breast")).value
    assert out_tool["annotatability"]["mean_confidence"] == pytest.approx(
        out_epi["annotatability"]["mean_confidence"])
    assert out_tool["annotatability"]["mean_confidence"] > out_mel["annotatability"]["mean_confidence"]


def test_identifiability_auc_separates_clean_types():
    """Per-type identifiability AUC is high for cleanly-separated types on their on-panel markers."""
    from spatialscribe.analysis import panel_check

    rng = np.random.default_rng(9)
    genes = ["CD3D", "CD3E", "TRAC", "MLANA", "SOX10", "DCT"]
    gi = {g: i for i, g in enumerate(genes)}
    n = 80
    X = rng.poisson(0.2, size=(n, len(genes))).astype("float32")
    for g in ("CD3D", "CD3E", "TRAC"):
        X[:40, gi[g]] += 12
    for g in ("MLANA", "SOX10", "DCT"):
        X[40:, gi[g]] += 12
    a = anndata.AnnData(X=X)
    a.var_names = genes
    a.obs_names = [f"c{i}" for i in range(n)]
    a.obs["cell_type"] = pd.Categorical(["T cell"] * 40 + ["Malignant/Melanocyte"] * 40)
    lineage = {"T cell": ["CD3D", "CD3E", "TRAC"], "Malignant/Melanocyte": ["MLANA", "SOX10", "DCT"]}
    auc = panel_check.identifiability_auc(a, cluster_key="cell_type", marker_sets=lineage)
    assert auc["T cell"]["auc"] > 0.9
    assert auc["Malignant/Melanocyte"]["auc"] > 0.9


def test_consensus_annotate_per_cell_majority_and_agreement(monkeypatch):
    import anndata, numpy as np, pandas as pd
    from spatialscribe.analysis import annotate

    n = 6
    a = anndata.AnnData(X=np.zeros((n, 2), dtype="float32"))
    a.obs_names = [f"c{i}" for i in range(n)]
    a.obs["leiden"] = pd.Categorical(["0"] * n)
    a.obs["rctd_first_type"] = ["T", "T", "Mel", "T", "Mel", "Mel"]
    a.obs["singler_label"] = ["T", "Mel", "Mel", "T", "Mel", None]
    # This synthetic AnnData has no real marker signal (2 dummy zero genes), so the
    # existing marker-scoring path would resolve every cluster to "Unknown" - not the
    # deterministic label this test needs to check the per-cell reconciliation math in
    # isolation. Pin the cluster-labelling step's output instead of fighting it, per the
    # task brief: adapt the test to read the *post-cluster* cell_type ("T" for cluster "0"),
    # matching the same value the reconciliation votes over.
    monkeypatch.setattr(annotate, "marker_labels", lambda *a_, **k_: {"0": "T"})
    # cell0: T(cluster),T(rctd),T(singler) -> T (3/3); cell2: T(cluster),Mel(rctd),Mel(singler) -> Mel (2/3)
    df = annotate.consensus_annotate(
        a, cluster_key="leiden", use_llm=False,
        method_label_cols=["rctd_first_type", "singler_label"])
    assert a.obs["cell_type"].iloc[0] == "T"
    assert a.obs["cell_type"].iloc[2] == "Mel"
    assert "consensus_agreement" in a.obs
    assert abs(float(a.obs["consensus_agreement"].iloc[0]) - 1.0) < 1e-6      # 3/3 agree
    assert abs(float(a.obs["consensus_agreement"].iloc[2]) - (2 / 3)) < 1e-6  # 2/3 agree


def test_reference_qc_skips_without_reference_and_scanvi_conf_lowers_confidence():
    """(a) purity.nmp/ncp are reference-guarded no-ops without a reference. (b) a per-cell
    scanvi_confidence column (reference-posterior signal) is an ADDITIONAL guarded soft
    penalty in apply_confidence, following the pmp/spatial_coherence/annotation_stability
    pattern - present only when the obs column exists, so the 34 pre-existing tests
    (which never set it) are unaffected.
    """
    from spatialscribe.analysis import annotate, purity

    # (a) reference-based QC no-ops without a reference.
    a = anndata.AnnData(X=np.abs(np.random.default_rng(0).normal(size=(40, 6)).astype("float32")))
    a.var_names = [f"g{i}" for i in range(6)]
    a.obs_names = [f"c{i}" for i in range(40)]
    assert purity.nmp(a, reference=None)["status"] == "skipped"
    assert purity.ncp(a, reference=None)["status"] == "skipped"

    # (b) a low scanvi_confidence column drags apply_confidence's mean confidence DOWN.
    # Deterministic two-run comparison (simplified per the task brief) instead of comparing
    # against a re-derived obs column, which is awkward to express cleanly.
    rng = np.random.default_rng(8)
    genes = ["CD3D", "CD3E", "TRAC", "MLANA", "SOX10", "DCT"]
    gi = {g: i for i, g in enumerate(genes)}
    n = 60
    X = rng.poisson(0.2, size=(n, len(genes))).astype("float32")
    for g in ("CD3D", "CD3E", "TRAC"):
        X[:30, gi[g]] += 12
    for g in ("MLANA", "SOX10", "DCT"):
        X[30:, gi[g]] += 12
    b = anndata.AnnData(X=X)
    b.var_names = genes
    b.obs_names = [f"c{i}" for i in range(n)]
    b.var["control"] = False
    b.obs["total_counts"] = np.asarray(X.sum(1)).ravel()
    b.obs["n_genes_by_counts"] = (X > 0).sum(1)
    b.obs["cell_type"] = pd.Categorical(["T cell"] * 30 + ["Malignant/Melanocyte"] * 30)
    lineage = {"T cell": ["CD3D", "CD3E", "TRAC"], "Malignant/Melanocyte": ["MLANA", "SOX10", "DCT"]}

    base = annotate.apply_confidence(b.copy(), cluster_key="cell_type", marker_sets=lineage)

    b2 = b.copy()
    b2.obs["scanvi_confidence"] = np.r_[np.full(30, 0.05), np.full(30, 0.95)]
    withconf = annotate.apply_confidence(b2, cluster_key="cell_type", marker_sets=lineage)

    assert withconf["mean_confidence"] < base["mean_confidence"]


def test_ncp_tiny_reference_returns_quickly_without_hanging():
    """purity.ncp must not enumerate the full O(n_genes^2) itertools.combinations space on a
    large panel (e.g. ~13M pairs on a 5K-gene panel hangs). Regression guard: a tiny synthetic
    reference (30 cells x 50 genes) with small max_genes/max_pairs caps returns a status dict
    promptly; no assertion on the exact ncp value (per the task brief).
    """
    from spatialscribe.analysis import purity

    ref = anndata.AnnData(
        X=np.abs(np.random.default_rng(1).normal(size=(30, 50)).astype("float32")))
    ref.var_names = [f"g{i}" for i in range(50)]
    ref.obs_names = [f"r{i}" for i in range(30)]
    ref.obs["cell_type"] = "typeA"

    a = anndata.AnnData(
        X=np.abs(np.random.default_rng(2).normal(size=(40, 50)).astype("float32")))
    a.var_names = ref.var_names
    a.obs_names = [f"c{i}" for i in range(40)]

    out = purity.ncp(a, reference=ref, max_genes=20, max_pairs=10)
    assert isinstance(out, dict) and "status" in out
