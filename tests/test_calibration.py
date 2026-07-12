# tests/test_calibration.py
"""Post-hoc calibration of the per-cell annotation confidence.

Motivating evidence (an internal benchmark harness independent 49-type GT on a 12k CosMx breast section coarsened
to a shared lineage axis): ``annotation_confidence`` is mis-calibrated (ECE 0.33; 0.45 on the 100k demo)
AND non-monotonic - lineage accuracy peaks near confidence 0.44 then FALLS through the 0.5-0.75 band.
So the score neither means P(correct) nor ranks cells by correctness.

Two separate defects, and only one of them is fixable post-hoc:

* **Calibration** (does the number mean P(correct)?) - fixed by a monotone isotonic fit against real
  labels. That is what :mod:`spatialscribe.analysis.calibration` does.
* **Discrimination** (does the number RANK correct cells above incorrect ones?) - NOT fixable by any
  monotone map, because a monotone map preserves order. This is why the benchmark's abstention gate did
  not help. :func:`calibration.auc` measures it so the report can say so instead of implying otherwise.

These tests pin both facts.
"""
from __future__ import annotations

import numpy as np
import pytest

anndata = pytest.importorskip("anndata")
pytest.importorskip("sklearn")


def _miscalibrated(n: int = 4000, seed: int = 0):
    """Reproduce the benchmark's shape: accuracy is NON-MONOTONIC in raw confidence.

    Low band [0.35, 0.5) -> 85% correct. High band [0.5, 0.75] -> 30% correct. A confident cell is
    therefore LESS likely to be right, exactly the inversion the benchmark measured.
    """
    rng = np.random.default_rng(seed)
    half = n // 2
    conf = np.concatenate([rng.uniform(0.35, 0.5, half), rng.uniform(0.5, 0.75, half)])
    correct = np.concatenate([(rng.random(half) < 0.85).astype(int),
                              (rng.random(half) < 0.30).astype(int)])
    return conf, correct


def _split(n, frac=0.5, seed=0):
    idx = np.random.default_rng(seed).permutation(n)
    cut = int(round(n * frac))
    return idx[:cut], idx[cut:]


# --------------------------------------------------------------------------- #
# metrics
# --------------------------------------------------------------------------- #
def test_ece_zero_when_perfectly_calibrated():
    from spatialscribe.analysis import calibration

    # 1000 cells at conf 0.2 of which exactly 20% correct; 1000 at 0.8 of which exactly 80%.
    conf = np.concatenate([np.full(1000, 0.2), np.full(1000, 0.8)])
    correct = np.concatenate([np.repeat([1, 0], [200, 800]), np.repeat([1, 0], [800, 200])])
    assert calibration.ece(conf, correct) == pytest.approx(0.0, abs=1e-9)


def test_brier_of_constant_half_on_balanced_labels():
    from spatialscribe.analysis import calibration

    correct = np.repeat([1, 0], [500, 500])
    assert calibration.brier(np.full(1000, 0.5), correct) == pytest.approx(0.25)


def test_reliability_curve_shape_and_drops_empty_bins():
    from spatialscribe.analysis import calibration

    conf, correct = _miscalibrated()
    cur = calibration.reliability_curve(conf, correct, n_bins=10)
    keys = {"bin_mid", "bin_conf", "bin_acc", "bin_count"}
    assert keys <= set(cur)
    lens = {len(cur[k]) for k in keys}
    assert len(lens) == 1                              # all aligned
    assert all(c > 0 for c in cur["bin_count"])        # empty bins dropped
    assert min(cur["bin_count"]) > 0


# --------------------------------------------------------------------------- #
# the two defects: calibration is fixable, discrimination is not
# --------------------------------------------------------------------------- #
def test_isotonic_lowers_ece_on_nonmonotonic_input():
    """THE load-bearing guarantee: fit on a held-out split, ECE must drop on the complement."""
    from spatialscribe.analysis import calibration

    conf, correct = _miscalibrated()
    fit, ev = _split(len(conf))
    cal = calibration.fit_isotonic(conf[fit], correct[fit])
    assert cal is not None

    before = calibration.ece(conf[ev], correct[ev])
    after = calibration.ece(np.clip(cal.predict(conf[ev]), 0, 1), correct[ev])
    assert before > 0.20, f"fixture is not miscalibrated enough (ECE {before:.3f})"
    assert after < 0.10, f"isotonic failed to calibrate (ECE {before:.3f} -> {after:.3f})"
    assert after < before


def test_isotonic_is_monotone():
    from spatialscribe.analysis import calibration

    conf, correct = _miscalibrated()
    cal = calibration.fit_isotonic(conf, correct)
    grid = np.linspace(0.0, 1.0, 101)
    assert np.all(np.diff(cal.predict(grid)) >= -1e-9)


def test_auc_exposes_the_inverted_gate_that_calibration_cannot_fix():
    """Confidence ANTI-ranks correctness here (AUC < 0.5), so abstaining on low confidence
    discards the cells most likely to be right. No monotone recalibration can repair a ranking."""
    from spatialscribe.analysis import calibration

    conf, correct = _miscalibrated()
    assert calibration.auc(conf, correct) < 0.45          # worse than a coin flip: the gate is inverted
    assert calibration.auc(conf, 1 - correct) > 0.55      # sanity: flipping the label flips the AUC


def test_fit_isotonic_returns_none_when_unfittable():
    from spatialscribe.analysis import calibration

    assert calibration.fit_isotonic(np.full(50, 0.5), np.ones(50, dtype=int)) is None   # 1 distinct conf
    assert calibration.fit_isotonic(np.array([0.1]), np.array([1])) is None


# --------------------------------------------------------------------------- #
# apply(): writes a NEW field, never overwrites the raw heuristic
# --------------------------------------------------------------------------- #
def _adata_with_conf(n=200, seed=0):
    conf, correct = _miscalibrated(n=n, seed=seed)
    a = anndata.AnnData(X=np.zeros((n, 2), dtype="float32"))
    a.obs_names = [f"c{i}" for i in range(n)]
    a.obs["annotation_confidence"] = conf
    a.obs["cell_type"] = ["T" if c else "B" for c in correct]
    a.obs["_gt"] = ["T"] * n                       # correct cells match, incorrect ones do not
    return a, conf, correct


def test_apply_writes_new_field_and_never_overwrites():
    from spatialscribe.analysis import calibration

    a, conf, correct = _adata_with_conf()
    cal = calibration.fit_isotonic(conf, correct)
    n = calibration.apply(a, cal)
    assert n == a.n_obs
    assert "annotation_confidence_calibrated" in a.obs
    assert np.array_equal(np.asarray(a.obs["annotation_confidence"]), conf)   # raw untouched
    vals = np.asarray(a.obs["annotation_confidence_calibrated"], dtype=float)
    assert vals.min() >= 0.0 and vals.max() <= 1.0


def test_apply_refuses_to_overwrite_the_raw_column():
    from spatialscribe.analysis import calibration

    a, conf, correct = _adata_with_conf()
    cal = calibration.fit_isotonic(conf, correct)
    with pytest.raises(ValueError):
        calibration.apply(a, cal, src="annotation_confidence", dst="annotation_confidence")


def test_apply_is_a_noop_without_a_calibrator():
    from spatialscribe.analysis import calibration

    a, _, _ = _adata_with_conf()
    assert calibration.apply(a, None) == 0
    assert "annotation_confidence_calibrated" not in a.obs


def test_isotonic_never_reorders_but_may_raise_auc_by_creating_ties():
    """Isotonic is non-STRICTLY increasing: its flat regions merge locally anti-correlated raw scores
    into ties, which moves those pairs from a losing to a neutral AUC contribution. So auc_after can
    exceed auc_before. It still never reorders two cells - a strictly-increasing map would leave AUC
    exactly unchanged. Measured on the real 12k CosMx section: 0.5045 raw -> 0.5674 calibrated.
    """
    from spatialscribe.analysis import calibration

    conf, correct = _miscalibrated()
    cal = calibration.fit_isotonic(conf, correct)
    after = np.clip(cal.predict(conf), 0, 1)

    assert calibration.auc(after, correct) >= calibration.auc(conf, correct)   # ties only help here
    # order is preserved: wherever the calibrated score strictly increases, so did the raw score
    o = np.argsort(conf)
    assert np.all(np.diff(after[o]) >= -1e-9)


def test_abstention_curve_flags_a_gate_that_only_wins_by_discarding_everything():
    from spatialscribe.analysis import calibration

    conf, correct = _miscalibrated()
    curve = calibration.abstention_curve(conf, correct, min_keep=0.10)
    assert curve and {"threshold", "keep_frac", "accuracy", "gain", "useful"} <= set(curve[0])
    # raw confidence anti-ranks correctness, so every threshold that keeps a usable share LOSES accuracy
    assert not any(d["useful"] for d in curve), "an inverted gate must never be reported as useful"


def test_auc_within_groups_unmasks_a_simpson_inversion():
    """Pooled AUC is confounded ACROSS groups. Real case (Atera breast, 50k cells): the pipeline is most
    confident on Endothelial (52.7% accurate) and least confident on Epithelial (98.6% accurate, 58% of
    cells), so pooled AUC inverts to 0.42 while the within-lineage signal is a flat 0.48.

    auc_within rank-normalises confidence inside each group first, removing the across-group offsets, so
    it reports the discrimination that a per-cell gate could actually exploit.
    """
    from spatialscribe.analysis import calibration

    rng = np.random.default_rng(0)
    # group A: HIGH confidence, LOW accuracy. group B: LOW confidence, HIGH accuracy. Within each group
    # confidence is pure noise, so the only signal is the (misleading) across-group offset.
    confA, corrA = rng.uniform(0.6, 0.9, 3000), (rng.random(3000) < 0.30).astype(int)
    confB, corrB = rng.uniform(0.1, 0.4, 3000), (rng.random(3000) < 0.90).astype(int)
    conf = np.concatenate([confA, confB])
    correct = np.concatenate([corrA, corrB])
    groups = np.array(["A"] * 3000 + ["B"] * 3000)

    assert calibration.auc(conf, correct) < 0.25          # pooled: catastrophically inverted
    res = calibration.auc_within(conf, correct, groups)
    assert 0.45 < res["auc_within"] < 0.55                # within: correctly reports "no signal"
    assert res["n_groups"] == 2
    assert set(res["per_group"]) == {"A", "B"}
    assert res["per_group"]["A"]["accuracy"] < res["per_group"]["B"]["accuracy"]
    assert res["per_group"]["A"]["mean_conf"] > res["per_group"]["B"]["mean_conf"]


def test_auc_within_groups_preserves_a_real_within_group_signal():
    from spatialscribe.analysis import calibration

    rng = np.random.default_rng(2)
    conf = rng.uniform(0, 1, 6000)
    correct = (rng.random(6000) < conf).astype(int)       # genuinely informative inside every group
    groups = rng.choice(["A", "B", "C"], 6000)
    assert calibration.auc_within(conf, correct, groups)["auc_within"] > 0.6


def test_auc_within_skips_degenerate_groups():
    from spatialscribe.analysis import calibration

    conf = np.array([0.1, 0.2, 0.3, 0.4, 0.9])
    correct = np.array([0, 1, 0, 1, 1])
    groups = np.array(["A", "A", "A", "A", "tiny"])       # 'tiny' has 1 cell -> dropped, never crashes
    res = calibration.auc_within(conf, correct, groups, min_n=2)
    assert "tiny" not in res["per_group"]
    assert res["n_groups"] == 1


def test_report_compares_against_the_base_rate_null_model():
    """A constant predictor equal to the base rate is already near-perfectly CALIBRATED (ECE ~ 0), so a
    low ece_after proves almost nothing on its own. Brier punishes a constant, so brier_skill (vs that
    null) is what shows the score carries per-cell information. Measured on real sections: skill is only
    0.02 (CosMx), 0.017 (Atera cervical) and 0.001 (Atera breast, i.e. nothing).
    """
    from spatialscribe.analysis import calibration

    conf, correct = _miscalibrated()
    cal = calibration.fit_isotonic(conf, correct)
    rep = calibration.report(conf, correct, cal, baserate=float(correct.mean()))

    # the null model is well-calibrated by construction, and no worse than the calibrated score on ECE
    assert rep["ece_baserate"] < 0.01
    # this fixture has NO discrimination (accuracy is non-monotonic in conf), so isotonic cannot beat it
    assert rep["brier_skill"] < 0.01, "a score with no ranking power must not claim skill over a constant"
    assert rep["brier_baserate"] == pytest.approx(calibration.brier(
        np.full(len(conf), correct.mean()), correct))


def test_brier_skill_is_positive_when_the_score_really_ranks():
    from spatialscribe.analysis import calibration

    # a genuinely informative score: accuracy increases monotonically with confidence
    rng = np.random.default_rng(1)
    conf = rng.uniform(0, 1, 4000)
    correct = (rng.random(4000) < conf).astype(int)
    cal = calibration.fit_isotonic(conf, correct)
    rep = calibration.report(conf, correct, cal, baserate=float(correct.mean()))
    assert rep["auc_after"] > 0.6
    assert rep["brier_skill"] > 0.2          # beats the constant by a wide margin


def test_report_carries_before_and_after():
    from spatialscribe.analysis import calibration

    conf, correct = _miscalibrated()
    cal = calibration.fit_isotonic(conf, correct)
    rep = calibration.report(conf, correct, cal)
    assert rep["ece_after"] < rep["ece_before"]
    assert rep["brier_after"] <= rep["brier_before"] + 1e-9
    assert rep["n"] == len(conf)
    assert "reliability_before" in rep and "reliability_after" in rep
    assert 0.0 <= rep["auc_before"] <= 1.0 and 0.0 <= rep["auc_after"] <= 1.0
    assert rep["abstention_after"]


# --------------------------------------------------------------------------- #
# capability: honest degradation, and calibration ONLY against real labels
# --------------------------------------------------------------------------- #
def test_cap_skips_without_truth_or_reference(processed_adata, ctx):
    from spatialscribe.analysis import capabilities as cap

    res = cap.run(processed_adata, "calibrate_confidence", {}, ctx)
    assert res.ok
    assert res.value["status"].startswith("skipped")
    assert res.value["obs_written"] is None
    assert "annotation_confidence_calibrated" not in processed_adata.obs
    assert "calibration_report" not in processed_adata.uns


def test_cap_requires_annotation_confidence(raw_adata, ctx):
    from spatialscribe.analysis import capabilities as cap

    res = cap.run(raw_adata, "calibrate_confidence", {}, ctx)
    assert not res.ok
    assert res.error["error_type"] == "prerequisite_missing"


def test_cap_with_truth_column_lowers_ece(ctx):
    from spatialscribe.analysis import capabilities as cap

    a, _, _ = _adata_with_conf(n=4000)
    raw = np.asarray(a.obs["annotation_confidence"], dtype=float).copy()
    res = cap.run(a, "calibrate_confidence", {"truth_key": "_gt", "pred_key": "cell_type"}, ctx)
    assert res.ok, res.error
    v = res.value
    assert v["status"] == "ok" and v["source"] == "truth_column"
    assert v["ece_after"] < v["ece_before"]
    assert v["obs_written"] == "annotation_confidence_calibrated"
    assert "annotation_confidence_calibrated" in a.obs
    assert "calibration_report" in a.uns
    assert np.array_equal(np.asarray(a.obs["annotation_confidence"], dtype=float), raw)  # raw untouched
    assert v["note"]                                            # an honest caveat is always attached
    assert "auc_before" in v and "auc_after" in v


def test_cap_gate_threshold_is_chosen_out_of_sample(ctx):
    """The reported gate must be picked on the FIT half and scored on the HELD-OUT half. Picking and
    scoring on the same cells reports a cherry-picked gain that will not survive on new data."""
    from spatialscribe.analysis import capabilities as cap

    a, _, _ = _adata_with_conf(n=4000)
    res = cap.run(a, "calibrate_confidence", {"truth_key": "_gt", "pred_key": "cell_type"}, ctx)
    assert res.ok, res.error
    gate = res.value["gate"]
    if gate is None:                       # a gate that never earns its place is a legitimate outcome
        assert "should be dropped" in res.value["note"]
        return
    assert gate["selected_on"] == "fit half" and gate["measured_on"] == "held-out half"
    assert 0.0 <= gate["keep_frac"] <= 1.0


# --------------------------------------------------------------------------- #
# apply_confidence: the calibrated gate is OPT-IN; default is byte-for-byte unchanged
# --------------------------------------------------------------------------- #
def test_apply_confidence_default_ignores_the_calibrated_column(processed_adata):
    from spatialscribe.analysis import annotate, markers

    ms = markers.for_tissue("melanoma")
    base = annotate.apply_confidence(processed_adata, cluster_key="cell_type", marker_sets=ms)
    base_verdict = np.asarray(processed_adata.obs["annotation_verdict"]).copy()

    processed_adata.obs["annotation_confidence_calibrated"] = 0.99      # would flip every gate
    again = annotate.apply_confidence(processed_adata, cluster_key="cell_type", marker_sets=ms)

    assert np.array_equal(np.asarray(processed_adata.obs["annotation_verdict"]), base_verdict)
    assert again["pct_pass"] == pytest.approx(base["pct_pass"])


def test_apply_confidence_use_calibrated_gates_on_the_calibrated_score(processed_adata):
    from spatialscribe.analysis import annotate, markers

    ms = markers.for_tissue("melanoma")
    base = annotate.apply_confidence(processed_adata, cluster_key="cell_type", marker_sets=ms)
    raw = np.asarray(processed_adata.obs["annotation_confidence"], dtype=float).copy()

    processed_adata.obs["annotation_confidence_calibrated"] = 0.99
    out = annotate.apply_confidence(processed_adata, cluster_key="cell_type", marker_sets=ms,
                                    use_calibrated=True)

    assert out["pct_pass"] > base["pct_pass"]                      # the calibrated score drives the gate
    # ... but the raw heuristic is still reported unchanged, never overwritten by the calibrated one.
    assert np.array_equal(np.asarray(processed_adata.obs["annotation_confidence"], dtype=float), raw)
