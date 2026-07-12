# tests/test_consensus.py
"""popV-style consensus uncertainty: trust cross-method AGREEMENT, not any single method's self-score.

Ports popV's load-bearing metric (Kimmel/Ergen/Yosef, Nat Genet 2024): the number of ensemble methods
that agree on a cell's label tracks accuracy far better than any method's own certainty. Tests the
consensus scoring, its integration into the confidence call (trusted only with a diverse ensemble), and
the plain-language 'methods disagree' rejection reason.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")


# --------------------------------------------------------------------------- #
# consensus_metrics + reliability bins
# --------------------------------------------------------------------------- #
def test_reliability_bin_thresholds():
    from spatialscribe.analysis import consensus

    assert consensus.reliability_bin(1.0, 4) == "very_high"
    assert consensus.reliability_bin(0.66, 3) == "moderate"   # recalibrated: 0.66 < high (0.7)
    assert consensus.reliability_bin(0.80, 5) == "high"       # 0.80 >= high (0.7)
    assert consensus.reliability_bin(0.40, 5) == "low"        # recalibrated: 0.40 < moderate (0.5)
    assert consensus.reliability_bin(0.20, 5) == "low"


def test_consensus_metrics_scores_agreement_with_winner():
    from spatialscribe.analysis import consensus

    a = anndata.AnnData(X=np.zeros((3, 2), dtype="float32"))
    a.obs_names = ["c0", "c1", "c2"]
    a.obs["cell_type"] = pd.Categorical(["T", "T", "T"])          # the consensus/winner label
    a.obs["rctd_first_type"] = ["T", "T", "B"]
    a.obs["singler_label"] = ["T", "Mel", "Mel"]
    a.obs["scanvi_label"] = ["T", "T", "T"]
    out = consensus.consensus_metrics(a, ["rctd_first_type", "singler_label", "scanvi_label"])
    assert out["status"] == "ok"
    assert list(a.obs["consensus_n_methods"]) == [3, 3, 3]
    assert list(a.obs["consensus_score"]) == [3, 2, 1]           # voters agreeing with winner 'T'
    assert a.obs["consensus_agreement"].tolist() == pytest.approx([1.0, 2 / 3, 1 / 3])
    assert a.obs["consensus_reliability"].iloc[0] == "very_high"
    assert a.obs["consensus_reliability"].iloc[1] == "moderate"   # 2/3 = 0.67 < high (0.7) after recalibration
    assert a.obs["consensus_reliability"].iloc[2] == "low"


def test_consensus_metrics_derives_winner_when_no_winner_col():
    from spatialscribe.analysis import consensus

    a = anndata.AnnData(X=np.zeros((2, 2), dtype="float32"))
    a.obs_names = ["c0", "c1"]
    a.obs["rctd_first_type"] = ["T", "T"]
    a.obs["singler_label"] = ["T", "B"]
    a.obs["scanvi_label"] = ["T", "B"]                            # c1 majority = B (2/3)
    consensus.consensus_metrics(a, ["rctd_first_type", "singler_label", "scanvi_label"],
                                winner_col=None)
    assert a.obs["consensus_score"].tolist() == [3, 2]           # c0 all T; c1 winner B, 2 agree
    assert a.obs["consensus_agreement"].tolist() == pytest.approx([1.0, 2 / 3])


def test_consensus_metrics_ignores_none_votes():
    from spatialscribe.analysis import consensus

    a = anndata.AnnData(X=np.zeros((1, 2), dtype="float32"))
    a.obs_names = ["c0"]
    a.obs["cell_type"] = pd.Categorical(["T"])
    a.obs["rctd_first_type"] = ["T"]
    a.obs["singler_label"] = [None]                              # method did not run for this cell
    a.obs["scanvi_label"] = ["T"]
    consensus.consensus_metrics(a, ["rctd_first_type", "singler_label", "scanvi_label"])
    assert int(a.obs["consensus_n_methods"].iloc[0]) == 2        # None is not a voter
    assert int(a.obs["consensus_score"].iloc[0]) == 2


def test_consensus_metrics_skips_without_method_columns():
    from spatialscribe.analysis import consensus

    a = anndata.AnnData(X=np.zeros((3, 2), dtype="float32"))
    a.obs_names = ["c0", "c1", "c2"]
    out = consensus.consensus_metrics(a, ["nonexistent_a", "nonexistent_b"])
    assert out["status"] == "skipped"
    assert "consensus_agreement" not in a.obs


# --------------------------------------------------------------------------- #
# apply_confidence: consensus is the TRUSTED factor, but only with a diverse ensemble
# --------------------------------------------------------------------------- #
def _clean_two_type(seed=0):
    rng = np.random.default_rng(seed)
    genes = ["CD3D", "CD3E", "TRAC", "MLANA", "SOX10", "DCT"] + [f"NOISE{i}" for i in range(40)]
    gi = {g: i for i, g in enumerate(genes)}
    n = 60
    X = rng.poisson(0.2, size=(n, len(genes))).astype("float32")
    for g in ("CD3D", "CD3E", "TRAC"):
        X[:30, gi[g]] += 12
        X[30:, gi[g]] = 0
    for g in ("MLANA", "SOX10", "DCT"):
        X[30:, gi[g]] += 12
        X[:30, gi[g]] = 0
    a = anndata.AnnData(X=X)
    a.var_names = genes
    a.obs_names = [f"c{i}" for i in range(n)]
    a.var["control"] = False
    a.obs["total_counts"] = np.asarray(X.sum(1)).ravel()
    a.obs["n_genes_by_counts"] = (X > 0).sum(1)
    a.obs["cell_type"] = pd.Categorical(["T cell"] * 30 + ["Malignant/Melanocyte"] * 30)
    return a


def test_apply_confidence_trusts_consensus_only_with_enough_methods():
    from spatialscribe.analysis import annotate

    lineage = {"T cell": ["CD3D", "CD3E", "TRAC"], "Malignant/Melanocyte": ["MLANA", "SOX10", "DCT"]}
    base = annotate.apply_confidence(_clean_two_type(), cluster_key="cell_type", marker_sets=lineage)

    # (A) Low consensus on half the cells, but only 2 methods voted -> NOT trusted (popV needs a
    #     diverse ensemble), so confidence is unchanged from baseline.
    a2 = _clean_two_type()
    a2.obs["consensus_agreement"] = np.r_[np.full(30, 0.2), np.full(30, 1.0)]
    a2.obs["consensus_n_methods"] = np.full(60, 2)
    two = annotate.apply_confidence(a2, cluster_key="cell_type", marker_sets=lineage)
    assert two["mean_confidence"] == pytest.approx(base["mean_confidence"])

    # (B) Same low consensus, but now 4 methods voted -> TRUSTED -> confidence drops.
    a4 = _clean_two_type()
    a4.obs["consensus_agreement"] = np.r_[np.full(30, 0.2), np.full(30, 1.0)]
    a4.obs["consensus_n_methods"] = np.full(60, 4)
    four = annotate.apply_confidence(a4, cluster_key="cell_type", marker_sets=lineage)
    assert four["mean_confidence"] < base["mean_confidence"]


def test_abstention_basis_reports_whether_it_is_ensemble_backed_or_heuristic_only():
    """The 2026-07 GT benchmark: the confidence heuristic barely ranks correct cells (within-lineage
    AUC ~0.54), only cross-method agreement does (~0.77). So apply_confidence must SAY whether its
    abstention is ensemble-backed (trustworthy) or heuristic-only (advisory) - honest, not implied."""
    from spatialscribe.analysis import annotate

    lineage = {"T cell": ["CD3D", "CD3E", "TRAC"], "Malignant/Melanocyte": ["MLANA", "SOX10", "DCT"]}

    # No ensemble voted -> heuristic-only, flagged untrustworthy, note tells the user how to upgrade.
    a0 = _clean_two_type()
    h0 = annotate.apply_confidence(a0, cluster_key="cell_type", marker_sets=lineage)
    assert h0["abstention_basis"] == "confidence_heuristic"
    assert a0.uns["abstention_basis"]["trustworthy"] is False
    assert "annotation_methods" in a0.uns["abstention_basis"]["note"]

    # A 2-method consensus is NOT a trusted ensemble (popV needs >=3) -> still heuristic.
    a2 = _clean_two_type()
    a2.obs["consensus_agreement"] = np.full(60, 0.5)
    a2.obs["consensus_n_methods"] = np.full(60, 2)
    assert annotate.apply_confidence(a2, cluster_key="cell_type",
                                     marker_sets=lineage)["abstention_basis"] == "confidence_heuristic"

    # >=3 diverse methods voted -> ensemble-backed, trustworthy.
    a4 = _clean_two_type()
    a4.obs["consensus_agreement"] = np.full(60, 0.8)
    a4.obs["consensus_n_methods"] = np.full(60, 4)
    h4 = annotate.apply_confidence(a4, cluster_key="cell_type", marker_sets=lineage)
    assert h4["abstention_basis"] == "ensemble_agreement"
    assert a4.uns["abstention_basis"]["trustworthy"] is True
    assert a4.uns["abstention_basis"]["frac_ensemble_backed"] == pytest.approx(1.0)


# --------------------------------------------------------------------------- #
# rejection: plain-language 'methods disagree' reason
# --------------------------------------------------------------------------- #
def test_rejection_flags_method_disagreement():
    from spatialscribe.analysis import rejection

    n, g = 10, 10
    a = anndata.AnnData(X=np.zeros((n, g), dtype="float32"))
    a.obs_names = [f"c{i}" for i in range(n)]
    a.obs["cell_type"] = pd.Categorical(["T cell"] * n)
    # All cells clean on counts/genes (no low-signal reason wins), but WARN verdict + low consensus.
    a.obs["total_counts"] = np.full(n, 2000.0)
    a.obs["n_genes_by_counts"] = np.full(n, g)
    a.layers["counts"] = np.full((n, g), 200.0, dtype="float32")     # genes=10 (>=5), counts=2000
    a.obs["annotation_verdict"] = np.array(["WARN"] * n, dtype=object)
    a.obs["annotation_reason"] = np.array([""] * n, dtype=object)
    a.obs["consensus_agreement"] = np.full(n, 0.3)               # methods disagree
    a.obs["consensus_n_methods"] = np.full(n, 5)
    rejection.assign_rejection_reasons(a)
    codes = a.obs["rejection_reason"].astype(str)
    assert (codes == "method_disagreement").any()
    det = a.obs.loc[codes == "method_disagreement", "rejection_detail"].iloc[0]
    assert "methods" in det.lower()


# --------------------------------------------------------------------------- #
# Reliability-weighted consensus (opt-in): majority loses when ONE method dominates
#
# popV/LatchBio showed weighting by a method's own SELF-confidence is futile - those scores are
# calibrated differently and are not comparable across methods. Reliability weighting is a DIFFERENT
# quantity: each method's accuracy measured against labels (or a documented prior). The annotator_bench
# (3k CosMx subsample vs independent GT, lineage axis) is the motivating evidence: naive majority scored
# 0.858, WORSE than RCTD alone (0.923), because SingleR (0.667) and panhumanpy (0.272) dilute it.
# --------------------------------------------------------------------------- #
def test_weighted_vote_flips_a_diluted_majority_to_the_reliable_label():
    from spatialscribe.analysis import consensus

    # RCTD says T; the two weaker methods both say B. Naive majority = B (2 votes vs 1) and is WRONG.
    labels = {"rctd_first_type": "T", "singler_label": "B", "ph_fine": "B"}
    win, agreement = consensus.weighted_vote(labels, consensus.DEFAULT_RELIABILITY)
    assert win == "T"
    # agreement keeps its popV meaning: the UNWEIGHTED fraction of voters backing the winner.
    assert agreement == pytest.approx(1 / 3)


def test_linear_accuracy_weights_would_still_pick_the_wrong_label():
    """Locks in WHY the weight is a log-odds transform, not raw accuracy.

    Linear: RCTD 0.83 < SingleR 0.67 + panhuman 0.45 = 1.12 -> still B (wrong).
    Log-odds: logit(.83)=1.59 > logit(.67)=0.71 + clamp0(logit(.45))=0 -> T (right).
    Guards against a future 'simplification' that drops the logit. (weights = the 2026-07-10 measured prior)
    """
    from spatialscribe.analysis import consensus

    labels = {"rctd_first_type": "T", "singler_label": "B", "ph_fine": "B"}
    w = consensus.DEFAULT_RELIABILITY

    linear = {}
    for voter, lab in labels.items():
        linear[lab] = linear.get(lab, 0.0) + w[voter]
    assert max(linear, key=linear.get) == "B"                       # naive linear weighting fails

    assert consensus.weighted_vote(labels, w)[0] == "T"             # log-odds weighting succeeds


def test_weighted_vote_ignores_a_method_worse_than_a_coin_flip():
    from spatialscribe.analysis import consensus

    # A single reliable voter must beat any number of sub-0.5 voters (they clamp to weight 0).
    labels = {"rctd_first_type": "T", "a": "B", "b": "B", "c": "B", "d": "B"}
    w = {"rctd_first_type": 0.92, "a": 0.3, "b": 0.3, "c": 0.3, "d": 0.3}
    assert consensus.weighted_vote(labels, w)[0] == "T"


def test_weighted_vote_flat_weights_reduce_to_majority():
    from spatialscribe.analysis import consensus

    labels = {"m1": "T", "m2": "B", "m3": "B"}
    flat = {v: 0.7 for v in labels}
    assert consensus.weighted_vote(labels, flat)[0] == "B"          # equal weights -> plain majority


def test_weighted_vote_ignores_missing_votes_and_prefers_the_cluster_label_on_ties():
    from spatialscribe.analysis import consensus

    labels = {"rctd_first_type": "T", "singler_label": None, "ph_fine": float("nan")}
    win, agreement = consensus.weighted_vote(labels, consensus.DEFAULT_RELIABILITY)
    assert win == "T" and agreement == pytest.approx(1.0)           # one voter, one agreement

    tied = {"m1": "T", "m2": "B"}
    w = {"m1": 0.8, "m2": 0.8}
    assert consensus.weighted_vote(tied, w, prefer="B")[0] == "B"   # cluster label breaks the tie
    assert consensus.weighted_vote(tied, w)[0] == "B"               # else alphabetical, deterministic


def test_reliability_from_labels_measures_accuracy_and_degrades_without_gt():
    from spatialscribe.analysis import consensus

    a = anndata.AnnData(X=np.zeros((4, 2), dtype="float32"))
    a.obs_names = [f"c{i}" for i in range(4)]
    a.obs["_gt"] = ["T", "T", "B", "B"]
    a.obs["rctd_first_type"] = ["T", "T", "B", "B"]      # 4/4 correct
    a.obs["singler_label"] = ["T", "B", "B", "T"]        # 2/4 correct
    cols = ["rctd_first_type", "singler_label"]

    assert consensus.reliability_from_labels(a, cols, gt_col="absent") == {}   # honest degradation
    w = consensus.reliability_from_labels(a, cols, gt_col="_gt")
    assert w["rctd_first_type"] == pytest.approx(1.0)
    assert w["singler_label"] == pytest.approx(0.5)


def test_consensus_annotate_default_stays_naive_majority(processed_adata):
    """Backward-compat: reliability_weights=None must reproduce the existing majority vote exactly."""
    from spatialscribe.analysis import annotate, markers

    a = processed_adata
    a.obs["rctd_first_type"] = "T cell"
    a.obs["singler_label"] = "B cell"
    a.obs["scanvi_label"] = "B cell"
    a.obs["ph_fine"] = "B cell"
    cols = ["rctd_first_type", "singler_label", "scanvi_label", "ph_fine"]

    annotate.consensus_annotate(a, cluster_key="leiden", use_llm=False,
                                marker_sets=markers.for_tissue("melanoma"), method_label_cols=cols)
    # Three B votes outnumber RCTD's single T vote even when the cluster label also says T, so the naive
    # majority is B cell everywhere - the dilution the benchmark measured (consensus 0.858 < RCTD 0.923).
    assert set(a.obs["cell_type"].astype(str)) == {"B cell"}


def test_consensus_annotate_reliability_weighted_recovers_the_dominant_method(processed_adata):
    from spatialscribe.analysis import annotate, consensus, markers

    a = processed_adata
    a.obs["rctd_first_type"] = "T cell"
    a.obs["singler_label"] = "B cell"
    a.obs["scanvi_label"] = "B cell"
    a.obs["ph_fine"] = "B cell"
    cols = ["rctd_first_type", "singler_label", "scanvi_label", "ph_fine"]

    annotate.consensus_annotate(a, cluster_key="leiden", use_llm=False,
                                marker_sets=markers.for_tissue("melanoma"), method_label_cols=cols,
                                reliability_weights=consensus.DEFAULT_RELIABILITY)
    # RCTD's log-odds weight (1.59) outweighs SingleR (0.71) + scANVI/panhuman (0, sub-coin-flip clamp);
    # the cluster label is only the tie-break (prefer), and there is no tie. The same 3-vs-1 split the
    # naive majority got backwards. (weights = the 2026-07-10 measured prior)
    assert set(a.obs["cell_type"].astype(str)) == {"T cell"}
    # consensus_agreement still reports the UNWEIGHTED fraction; weighting moves the winner, not the metric.
    assert a.obs["consensus_agreement"].max() <= 0.5
