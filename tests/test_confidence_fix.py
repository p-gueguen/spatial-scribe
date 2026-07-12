"""Confidence-formula fixes: panel-size-invariant PMP + averaged (not multiplied) soft factors.

Diagnosed on a real Xenium Prime 5K breast section: 68% of cells abstained, dominated by
'low confidence'. Two root causes, both fixed here:
  1. PMP divided marker counts by TOTAL panel counts -> structurally ~0 on 5K/WTA panels,
     a flat ~30% purity tax on ~94% of cells regardless of true purity.
  2. The soft down-weights (purity, coherence, stability, ref posteriors) were MULTIPLIED,
     so a stack of mild penalties compounded (0.7^3 ~= 0.34) and crushed good calls.
"""
from __future__ import annotations

import numpy as np
import pytest

anndata = pytest.importorskip("anndata")


def _adata(genes, counts, labels):
    import pandas as pd

    a = anndata.AnnData(X=np.asarray(counts, dtype="float32"))
    a.var_names = list(genes)
    a.obs_names = [f"c{i}" for i in range(len(labels))]
    a.obs["cell_type"] = pd.Categorical(list(labels))
    return a


def test_pmp_is_panel_size_invariant():
    """PMP must not change when non-marker filler genes are added to the panel (the bug)."""
    from spatialscribe.analysis import purity

    markers = {"T": ["CD3D", "CD3E"], "B": ["MS4A1"]}
    a_small = _adata(["CD3D", "CD3E", "MS4A1"], [[10, 8, 1]], ["T"])
    a_big = _adata(["CD3D", "CD3E", "MS4A1"] + [f"F{i}" for i in range(50)],
                   [[10, 8, 1] + [5] * 50], ["T"])
    purity.pmp(a_small, lineage_markers=markers)
    purity.pmp(a_big, lineage_markers=markers)

    v_small, v_big = float(a_small.obs["pmp"][0]), float(a_big.obs["pmp"][0])
    assert abs(v_small - v_big) < 1e-9                 # panel-size invariant (old metric was NOT)
    # value = assigned-type marker counts (10+8) / all-lineage-marker counts (10+8+1)
    assert abs(v_small - 18 / 19) < 1e-9


def test_pmp_nan_when_no_marker_signal():
    """No lineage-marker transcripts at all -> purity undefined -> NaN (not a false 0)."""
    from spatialscribe.analysis import purity

    a = _adata(["CD3D", "F1", "F2"], [[0, 5, 7]], ["T"])
    purity.pmp(a, lineage_markers={"T": ["CD3D"]})
    assert np.isnan(float(a.obs["pmp"][0]))


def test_combine_soft_averages_not_multiplies():
    """The soft-factor combiner AVERAGES (so mild penalties don't compound), neutral 1.0 if none."""
    from spatialscribe.analysis import annotate

    assert abs(float(annotate._combine_soft([np.array([0.8])] * 3)[0]) - 0.8) < 1e-9   # not 0.512
    assert float(annotate._combine_soft([])) == 1.0
    assert abs(float(annotate._combine_soft([np.array([0.7]), np.array([1.0])])[0]) - 0.85) < 1e-9


def test_apply_confidence_nan_pmp_is_neutral():
    """A NaN PMP (no-marker cell) must be treated as neutral, never corrupting confidence."""
    import pandas as pd

    from spatialscribe.analysis import annotate

    markers = {"T": ["CD3D", "CD3E"], "B": ["MS4A1", "CD79A"]}
    rng = np.random.default_rng(0)
    a = anndata.AnnData(X=rng.poisson(3, size=(30, 4)).astype("float32"))
    a.var_names = ["CD3D", "CD3E", "MS4A1", "CD79A"]
    a.layers["counts"] = a.X.copy()
    a.obs["cell_type"] = pd.Categorical(rng.choice(["T", "B"], 30))
    a.obs["pmp"] = np.where(np.arange(30) < 5, np.nan, 0.9)   # 5 no-marker cells

    annotate.apply_confidence(a, marker_sets=markers)
    conf = np.asarray(a.obs["annotation_confidence"], dtype=float)
    assert not np.isnan(conf).any()


# --------------------------------------------------------------------------- #
# pct_counts_control must measure BACKGROUND, not "everything that is not a gene".
#
# Diagnosed on the real Xenium Prime 5K breast section (2026-07 overnight benchmark): 94.8% of cells
# abstained. Measured per feature type on that section: TRUE negative controls carry 0.000% of a
# cell's counts, while the 3294 'Deprecated Codeword' features - decommissioned GENE probes that
# still decode real transcripts - carry a median 14.96%. The broad control mask counted those, so
# pct_counts_control was ~15%, contam = clip(15/5, 0, 1) = 1 for 86.2% of cells, and
# confidence = posterior * (1-contam) * ceiling went to exactly 0 (frac_conf_zero == frac_contam_1
# == 0.8618) -> abstained as 'ambiguous_mixed'. NOT the count floor (1.1%), NOT panel coverage
# (0 red lineages). Atera WTA / standard Xenium ship no deprecated codewords, hence 0.06% / 0.00%.
# --------------------------------------------------------------------------- #
def _xenium_var(feature_types):
    import pandas as pd

    return pd.DataFrame({"feature_types": list(feature_types)})


def test_deprecated_codewords_are_not_negative_controls():
    """Deprecated/unassigned codewords carry real signal - they are panel non-genes, not background."""
    from spatialscribe.analysis import io

    a = anndata.AnnData(X=np.ones((3, 5), dtype="float32"))
    a.var_names = ["CD3D", "DeprecatedCodeword_1", "UnassignedCodeword_1",
                   "NegControlProbe_1", "NegControlCodeword_1"]
    a.var = _xenium_var(["Gene Expression", "Deprecated Codeword", "Unassigned Codeword",
                         "Negative Control Probe", "Negative Control Codeword"]).set_index(a.var_names)

    broad = io.build_control_mask(a)
    strict = io.build_neg_control_mask(a)
    assert list(broad) == [False, True, True, True, True]      # panel mask: every non-gene
    assert list(strict) == [False, False, False, True, True]   # background: true neg-controls only


def test_neg_control_mask_falls_back_to_names_without_feature_types():
    """CosMx/MERSCOPE exports have no feature_types - the name fallback must exclude deprecated too."""
    from spatialscribe.analysis import io

    a = anndata.AnnData(X=np.ones((2, 4), dtype="float32"))
    a.var_names = ["Snap25", "NegPrb1", "DeprecatedCodeword_7", "BLANK_3"]
    assert list(io.build_neg_control_mask(a)) == [False, True, False, True]


def test_prime5k_like_section_is_not_abstained(tmp_path):
    """End-to-end: a Prime 5K-like section (deprecated codewords hold ~15% of counts) stays usable."""
    from spatialscribe.analysis import annotate, qc

    rng = np.random.default_rng(0)
    n = 60
    genes = ["CD3D", "CD3E", "MS4A1", "CD79A", "F1", "F2", "F3", "F4"]
    names = genes + ["DeprecatedCodeword_1", "NegControlProbe_1"]
    X = np.column_stack([
        rng.poisson(6, size=(n, len(genes))),          # ~48 gene counts / cell
        rng.poisson(9, size=(n, 1)),                   # deprecated codeword: ~15% of the cell
        np.zeros((n, 1)),                              # true negative control: silent
    ]).astype("float32")
    a = anndata.AnnData(X=X)
    a.var_names = names
    a.var = _xenium_var(["Gene Expression"] * len(genes)
                        + ["Deprecated Codeword", "Negative Control Probe"]).set_index(a.var_names)
    a.layers["counts"] = a.X.copy()
    qc.compute_qc(a)                                   # must score background, not deprecated signal

    assert float(np.median(a.obs["pct_counts_control"])) < 1.0

    import pandas as pd
    a.obs["cell_type"] = pd.Categorical(rng.choice(["T", "B"], n))
    res = annotate.apply_confidence(a, marker_sets={"T": ["CD3D", "CD3E"], "B": ["MS4A1", "CD79A"]})
    assert res["pct_abstain"] < 0.5, f"Prime 5K-like section abstained: {res['pct_abstain']:.3f}"
    assert res["usability"] in ("ok", "warn")


def test_a_truly_contaminated_cell_is_still_zeroed():
    """The penalty must survive the fix: >=5% of counts on TRUE negative controls = contaminated."""
    import pandas as pd

    from spatialscribe.analysis import annotate

    rng = np.random.default_rng(0)
    n = 40
    a = anndata.AnnData(X=rng.poisson(6, size=(n, 8)).astype("float32"))
    a.var_names = ["CD3D", "CD3E", "MS4A1", "CD79A", "F1", "F2", "F3", "F4"]
    a.layers["counts"] = a.X.copy()
    a.obs["cell_type"] = pd.Categorical(rng.choice(["T", "B"], n))
    pct = np.full(n, 0.06)                              # healthy background
    pct[0] = 5.5                                        # one genuinely contaminated cell
    a.obs["pct_counts_control"] = pct

    annotate.apply_confidence(a, marker_sets={"T": ["CD3D", "CD3E"], "B": ["MS4A1", "CD79A"]})
    conf = np.asarray(a.obs["annotation_confidence"], dtype=float)
    assert conf[0] == 0.0
    assert conf[1:].mean() > 0.1


def test_usability_gate_flags_an_unusable_annotation():
    """Output-usability gate: 'ran to completion with ok input-QC' must not hide ~95% abstention."""
    from spatialscribe.analysis import annotate

    assert annotate.usability_flag(0.02) == "ok"
    assert annotate.usability_flag(0.30) == "warn"
    assert annotate.usability_flag(0.948) == "fail"     # the real Prime 5K benchmark number
