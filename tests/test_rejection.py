# tests/test_rejection.py
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")


def _seeded(n_extra_cols=True):
    """A small AnnData whose obs seeds one cell per granular condition.

    Layout (10 cells):
      0 PASS  (confident, no reason)
      1 FAIL  low counts        -> too_few_counts
      2 FAIL  few genes         -> too_few_genes
      3 FAIL  crisp_impure      -> mixed_lineages
      4 FAIL  high control %    -> high_background
      5 WARN  low vsi           -> spatial_doublet   (optional col)
      6 WARN  seg fragment      -> poor_segmentation (optional col)
      7 WARN  low coherence     -> spatially_incoherent (optional col)
      8 WARN  nothing special   -> low_margin
      9 FAIL  novel             -> novel_unknown
    """
    n = 10
    a = anndata.AnnData(X=np.ones((n, 3), dtype="float32"))
    a.var_names = ["CD3D", "CD3E", "TRAC"]
    a.obs_names = [f"c{i}" for i in range(n)]

    verdict = np.array(["PASS", "FAIL", "FAIL", "FAIL", "FAIL",
                        "WARN", "WARN", "WARN", "WARN", "FAIL"])
    reason = np.array(["", "low_quality", "low_quality", "ambiguous_mixed", "ambiguous_mixed",
                       "", "", "", "", "novel"])
    a.obs["annotation_verdict"] = pd.Categorical(verdict)
    a.obs["annotation_reason"] = reason
    a.obs["annotation_confidence"] = np.where(verdict == "PASS", 0.9,
                                              np.where(verdict == "WARN", 0.4, 0.1))

    counts = np.full(n, 500.0)
    counts[1] = 4.0            # cell 1: too few counts
    genes = np.full(n, 50.0)
    genes[2] = 3.0             # cell 2: too few genes
    a.obs["total_counts"] = counts
    a.obs["n_genes_by_counts"] = genes

    impure = np.zeros(n, dtype=bool)
    impure[3] = True           # cell 3: mixed lineages
    a.obs["crisp_impure"] = impure

    pct = np.zeros(n)
    pct[4] = 12.0              # cell 4: high control background
    a.obs["pct_counts_control"] = pct

    a.obs["cell_type"] = pd.Categorical(["T cell"] * n)

    if n_extra_cols:
        vsi = np.full(n, 0.95)
        vsi[5] = 0.3           # cell 5: spatial doublet
        a.obs["vsi"] = vsi

        seg = np.array(["ok"] * n, dtype=object)
        seg[6] = "small"       # cell 6: poor segmentation
        a.obs["seg_area_flag"] = seg
        a.obs["nucleus_present"] = np.ones(n, dtype=bool)

        coh = np.full(n, 0.9)
        coh[7] = 0.05          # cell 7: spatially incoherent
        a.obs["spatial_coherence"] = coh
    return a


def test_writes_columns_and_pass_is_blank():
    from spatialscribe.analysis import rejection

    a = _seeded()
    summ = rejection.assign_rejection_reasons(a)
    assert "rejection_reason" in a.obs
    assert "rejection_detail" in a.obs
    # PASS cell (0) gets no reason.
    assert str(a.obs["rejection_reason"].iloc[0]) == ""
    assert a.obs["rejection_detail"].iloc[0] == ""
    # Non-PASS cells (9) are all reasoned.
    assert summ["n_non_pass"] == 9
    assert summ["n_reasoned"] == 9


def test_granular_codes_for_seeded_conditions():
    from spatialscribe.analysis import rejection

    a = _seeded()
    rejection.assign_rejection_reasons(a)
    code = a.obs["rejection_reason"].astype(str).to_numpy()
    assert code[1] == "too_few_counts"
    assert code[2] == "too_few_genes"
    assert code[3] == "mixed_lineages"
    assert code[4] == "high_background"
    assert code[5] == "spatial_doublet"
    assert code[6] == "poor_segmentation"
    assert code[7] == "spatially_incoherent"
    assert code[8] == "low_margin"
    assert code[9] == "novel_unknown"


def test_no_nucleus_is_not_poor_segmentation():
    """A cell missing an in-plane nucleus is EXPECTED in a thin tissue section (the nucleus sits in an
    adjacent cut), not an 'implausible segment'. Only size outliers (fragment/merged) are segmentation
    problems. Regression for over-flagging valid anucleate cells as poor_segmentation."""
    from spatialscribe.analysis import rejection

    a = _seeded()
    nuc = np.ones(a.n_obs, dtype=bool)
    nuc[8] = False                     # cell 8 is a plain WARN with a NORMAL-size segment, no nucleus in-plane
    a.obs["nucleus_present"] = nuc
    rejection.assign_rejection_reasons(a)
    code = a.obs["rejection_reason"].astype(str).to_numpy()
    assert code[8] != "poor_segmentation", "anucleate cell wrongly flagged as an implausible segment"
    assert code[8] == "low_margin"     # it stays a plain low-margin WARN
    assert code[6] == "poor_segmentation"   # a real size outlier is still flagged


def test_details_carry_specific_numbers():
    from spatialscribe.analysis import rejection

    a = _seeded()
    rejection.assign_rejection_reasons(a)
    det = a.obs["rejection_detail"].astype(str).to_numpy()
    assert "4" in det[1] and "transcript" in det[1].lower()
    assert "3" in det[2] and "gene" in det[2].lower()
    assert "no em-dash" not in det[8]
    # No em-dashes anywhere in user-facing strings.
    assert not any("—" in s for s in det)


def test_breakdown_sums_to_non_pass():
    from spatialscribe.analysis import rejection

    a = _seeded()
    rejection.assign_rejection_reasons(a)
    df = rejection.rejection_breakdown(a)
    assert set(["reason", "label", "n_cells", "pct_of_untyped", "description"]) <= set(df.columns)
    assert int(df["n_cells"].sum()) == 9  # equals the number of non-PASS cells
    # Sorted descending.
    assert list(df["n_cells"]) == sorted(df["n_cells"], reverse=True)
    assert abs(df["pct_of_untyped"].sum() - 100.0) < 0.5


def test_robust_when_optional_columns_absent():
    from spatialscribe.analysis import rejection

    a = _seeded(n_extra_cols=False)  # no vsi / seg / coherence columns
    summ = rejection.assign_rejection_reasons(a)
    code = a.obs["rejection_reason"].astype(str).to_numpy()
    # Column-independent reasons still fire.
    assert code[1] == "too_few_counts"
    assert code[2] == "too_few_genes"
    assert code[3] == "mixed_lineages"
    assert code[4] == "high_background"
    assert code[9] == "novel_unknown"
    # Cells that would have needed the missing columns fall through to low_margin (never crash).
    assert code[5] in {"low_margin"}
    assert code[6] in {"low_margin"}
    assert code[7] in {"low_margin"}
    assert summ["n_reasoned"] == 9


def test_works_without_annotation_reason_column():
    """assign_rejection_reasons must still work when the coarse reason is absent (infer)."""
    from spatialscribe.analysis import rejection

    a = _seeded(n_extra_cols=False)
    del a.obs["annotation_reason"]
    rejection.assign_rejection_reasons(a)
    code = a.obs["rejection_reason"].astype(str).to_numpy()
    assert code[0] == ""                        # PASS untouched
    assert code[1] == "too_few_counts"          # inferred from counts
    assert code[3] == "mixed_lineages"          # inferred from crisp_impure


def test_panel_resolvability_warnings():
    from spatialscribe.analysis import rejection

    a = _seeded(n_extra_cols=False)
    a.obs["cell_type"] = pd.Categorical(
        ["NK cell"] * 3 + ["T cell"] * 3 + ["Endothelial"] * 4)
    pc = {
        "coverage": {
            "NK cell": {"status": "red", "n_present": 0, "n_markers": 4},
            "T cell": {"status": "green", "n_present": 5, "n_markers": 7},
            "Endothelial": {"status": "amber", "n_present": 1, "n_markers": 4},
        },
        "confusable_pairs": [{"pair": ["NK cell", "T cell"]}],
    }
    warn = rejection.panel_resolvability_warnings(a, pc)
    kinds = {w.get("cell_type") or tuple(w.get("pair", [])) for w in warn}
    assert "NK cell" in kinds          # red coverage surfaced
    assert "Endothelial" in kinds      # amber coverage surfaced
    assert "T cell" not in kinds       # green not surfaced
    # NK-cell warning reports the affected cell count.
    nk = next(w for w in warn if w.get("cell_type") == "NK cell")
    assert nk["n_cells"] == 3
    assert "NK cell" in nk["message"]


def _typability_fixture():
    """A section where one amber-coverage type IS separable and another is NOT.

    Epithelial/Tumor: 2 on-panel markers (EPCAM, KRT8), strongly expressed -> identifiability
    AUC ~1.0 (the Panel-check tab calls it 'resolved: yes' despite 2/8 coverage).
    Adipocyte: 1 on-panel marker (ADIPOQ) that is pure noise -> AUC ~0.5 (genuinely weak).
    ~80 background genes span the expression range so sc.tl.score_genes finds control bins.
    """
    rng = np.random.default_rng(0)
    per = 100
    types = ["Epithelial/Tumor"] * per + ["Adipocyte"] * per + ["Other"] * per
    n = len(types)
    is_tum = np.array([t == "Epithelial/Tumor" for t in types])
    is_oth = np.array([t == "Other" for t in types])

    markers = ["EPCAM", "KRT8", "ADIPOQ", "CD3D", "CD3E", "CD8A"]
    bg = [f"BG{i}" for i in range(80)]
    var = markers + bg
    X = np.zeros((n, len(var)), dtype="float32")
    # tumour markers: high in tumour, ~0 elsewhere -> strong separation
    for j in (0, 1):
        X[:, j] = rng.poisson(np.where(is_tum, 15.0, 0.2)).astype("float32")
    # adipocyte "marker": uniform noise across ALL cells -> no separation
    X[:, 2] = rng.poisson(3.0, n).astype("float32")
    # Other markers: high in Other -> Other is green + separable
    for j in (3, 4, 5):
        X[:, j] = rng.poisson(np.where(is_oth, 15.0, 0.2)).astype("float32")
    # background genes span the mean-expression range (score_genes control bins)
    for k in range(len(bg)):
        X[:, len(markers) + k] = rng.poisson(rng.uniform(0.5, 20.0), n).astype("float32")

    a = anndata.AnnData(X=X)
    a.var_names = var
    a.obs_names = [f"c{i}" for i in range(n)]
    a.obs["cell_type"] = pd.Categorical(types)
    pc = {
        "coverage": {
            "Epithelial/Tumor": {"status": "amber", "n_present": 2, "n_markers": 8,
                                 "present": ["EPCAM", "KRT8"]},
            "Adipocyte": {"status": "amber", "n_present": 1, "n_markers": 6,
                          "present": ["ADIPOQ"]},
            "Other": {"status": "green", "n_present": 3, "n_markers": 5,
                      "present": ["CD3D", "CD3E", "CD8A"]},
        },
        "confusable_pairs": [],
    }
    return a, pc


def test_panel_warning_defers_to_identifiability_auc():
    """Regression: the Annotate 'weakly resolved' warning must agree with the Panel-check
    'resolved' verdict. A type with few on-panel markers but a HIGH identifiability AUC is
    confidently typable (Panel-check says 'yes'), so the coverage warning must be SUPPRESSED -
    while a genuinely non-separable amber type keeps its warning."""
    from spatialscribe.analysis import rejection, panel_check as _pc

    a, pc = _typability_fixture()

    # The two surfaces must reach the same verdict for each type.
    ms = {t: d["present"] for t, d in pc["coverage"].items()}
    typable = {r["cell_type"] for r in _pc.typability_table(
        a, marker_sets=ms, panel_check_result=pc) if r.get("confidently_typable")}
    assert "Epithelial/Tumor" in typable            # 2/8 markers, but AUC ~1.0
    assert "Adipocyte" not in typable               # noise marker, AUC ~0.5

    warn = rejection.panel_resolvability_warnings(a, pc)
    warned = {w.get("cell_type") for w in warn if w.get("kind") == "coverage"}
    assert "Epithelial/Tumor" not in warned, "separable type wrongly flagged 'weakly resolved'"
    assert "Adipocyte" in warned                    # genuinely weak -> still warned


def test_panel_reasons_assigned_to_cells():
    from spatialscribe.analysis import rejection

    a = _seeded(n_extra_cols=False)
    # Make cell 8 an NK cell with a red panel so it is attributed to the panel, not low_margin.
    ct = a.obs["cell_type"].astype(str).to_numpy()
    ct[8] = "NK cell"
    a.obs["cell_type"] = pd.Categorical(ct)
    pc = {"coverage": {"NK cell": {"status": "red", "n_present": 0, "n_markers": 4},
                       "T cell": {"status": "green", "n_present": 5, "n_markers": 7}},
          "confusable_pairs": []}
    rejection.assign_rejection_reasons(a, panel_check_result=pc)
    assert a.obs["rejection_reason"].astype(str).iloc[8] == "panel_cannot_resolve"
    assert "NK cell" in a.obs["rejection_detail"].iloc[8]
