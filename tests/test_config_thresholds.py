"""Prove annotate.apply_confidence reads its verdict cutoffs from config, not hardcode.

apply_confidence computes a per-cell confidence, then splits PASS/WARN/FAIL at the
layer5_confidence composite_confidence warn/fail thresholds fetched via config.get. If the
cut were hardcoded, swapping those config values would not move the WARN+FAIL fraction.
Here we monkeypatch config.get: aggressive cutoffs (warn=0.99, fail=0.98) must flag strictly
MORE cells than permissive cutoffs (warn=0.01, fail=0.0) on the very same confidence scores.
"""

from __future__ import annotations

from spatialscribe.analysis import annotate, config, markers


def _patched_get(warn, fail):
    """Wrap config.get: special-case the layer5 warn/fail keys, delegate everything else."""
    real_get = config.get

    def fake_get(section, *keys, default=None, path=None):
        if section == "layer5_confidence" and keys == ("composite_confidence", "warn"):
            return warn
        if section == "layer5_confidence" and keys == ("composite_confidence", "fail"):
            return fail
        return real_get(section, *keys, default=default, path=path)

    return fake_get


def _flagged_fraction(result: dict) -> float:
    """Fraction of cells NOT confidently PASS = WARN (greyed) + FAIL (abstained)."""
    return result["pct_warn"] + result["pct_abstain"]


def test_verdict_cut_is_config_driven(processed_adata, monkeypatch):
    marker_sets = markers.for_tissue("melanoma")

    # Permissive cutoffs: almost nothing should be flagged.
    monkeypatch.setattr(config, "get", _patched_get(warn=0.01, fail=0.0))
    permissive = annotate.apply_confidence(
        processed_adata, cluster_key="cell_type", marker_sets=marker_sets
    )

    # Aggressive cutoffs on the SAME section/confidence scores: nearly everything flagged.
    monkeypatch.setattr(config, "get", _patched_get(warn=0.99, fail=0.98))
    aggressive = annotate.apply_confidence(
        processed_adata, cluster_key="cell_type", marker_sets=marker_sets
    )

    permissive_flagged = _flagged_fraction(permissive)
    aggressive_flagged = _flagged_fraction(aggressive)

    # Config-driven cut: tightening the thresholds must strictly increase the flagged share.
    assert aggressive_flagged > permissive_flagged, (
        f"thresholds not config-driven: permissive={permissive_flagged:.3f} "
        f"aggressive={aggressive_flagged:.3f}"
    )
    # Sanity: the aggressive run should flag a large majority.
    assert aggressive_flagged > 0.5
