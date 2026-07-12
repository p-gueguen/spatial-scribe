"""End-to-end pipeline smoke test on a tiny synthetic section (CPU, no API key).

Exercises the analysis chain that does NOT need Claude or the demo download:
io control mask -> QC -> preprocess+cluster -> marker labels -> panel check -> squidpy
neighborhood enrichment. Guards against scanpy/squidpy/decoupler API drift.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")
pytest.importorskip("scanpy")


def _synthetic():
    """Two spatially-separated populations: T cells vs melanocytes, + noise + controls."""
    rng = np.random.default_rng(0)
    t_markers = ["CD3D", "CD3E", "CD2", "CD8A", "TRAC"]
    mel_markers = ["MLANA", "PMEL", "TYR", "DCT", "MITF", "SOX10"]
    noise = [f"NOISE{i}" for i in range(30)]
    controls = ["NegControlProbe_00001", "BLANK_0001"]
    genes = t_markers + mel_markers + noise + controls
    n = 400
    grp = np.array([0] * (n // 2) + [1] * (n - n // 2))  # 0 = T, 1 = melanocyte

    X = rng.poisson(0.3, size=(n, len(genes))).astype("float32")
    gidx = {g: i for i, g in enumerate(genes)}
    for g in t_markers:
        X[grp == 0, gidx[g]] += rng.poisson(8, (grp == 0).sum())
    for g in mel_markers:
        X[grp == 1, gidx[g]] += rng.poisson(8, (grp == 1).sum())

    var = pd.DataFrame(index=genes)
    var["feature_types"] = ["Gene Expression"] * (len(genes) - 2) + ["Negative Control Probe"] * 2
    a = anndata.AnnData(X=X, var=var)
    a.obs_names = [f"cell{i}" for i in range(n)]
    # Two spatial blobs so clustering + neighborhoods have structure.
    coords = np.column_stack([
        rng.normal(np.where(grp == 0, 0, 50), 4),
        rng.normal(0, 4, n),
    ])
    a.obsm["spatial"] = coords
    a.obs["true_group"] = pd.Categorical(np.where(grp == 0, "T", "Mel"))
    return a


def test_end_to_end_cpu(monkeypatch):
    monkeypatch.setenv("SPATIALSCRIBE_FORCE_CPU", "1")
    from spatialscribe.analysis import (
        annotate, cluster, export, io, markers, panel_check, purity, qc, spatial, subcluster,
    )

    a = _synthetic()
    a.var["control"] = io.build_control_mask(a)
    panel = a.var_names[~a.var["control"]].tolist()

    # QC
    qc.compute_qc(a)
    summ = qc.qc_summary(a)
    assert summ["n_cells"] == a.n_obs
    assert "pct_counts_control" in a.obs

    # Panel check with the curated melanoma sets (deterministic, offline)
    pc = panel_check.check_panel(panel, marker_sets=markers.LINEAGE_MARKERS)
    assert pc["n_cell_types"] > 0
    assert pc["coverage"]["Malignant/Melanocyte"]["n_present"] >= 3  # MLANA/DCT/MITF/SOX10...
    assert pc["coverage"]["T cell"]["n_present"] >= 3

    # Cluster (CPU) + marker labels (no LLM)
    cluster.cluster(a, resolution=1.0)
    assert "leiden" in a.obs
    labels = annotate.marker_labels(a, cluster_key="leiden")
    assert set(labels.values()) & {"T cell", "Malignant/Melanocyte"}

    # Consensus without LLM sets cell_type
    df = annotate.consensus_annotate(a, cluster_key="leiden", use_llm=False)
    assert "cell_type" in a.obs
    assert len(df) == a.obs["leiden"].nunique()

    # Per-cell confidence + abstention (Layer 5) + purity metrics (Layer 3)
    pcheck = panel_check.check_panel(panel, marker_sets=markers.LINEAGE_MARKERS)
    summ2 = annotate.apply_confidence(
        a, cluster_key="cell_type", marker_sets=markers.LINEAGE_MARKERS, panel_check_result=pcheck,
    )
    assert {"annotation_confidence", "annotation_verdict", "cell_type_final"} <= set(a.obs.columns)
    assert abs(summ2["pct_pass"] + summ2["pct_warn"] + summ2["pct_abstain"] - 1.0) < 1e-6

    # Full QC funnel headline (Layers 0-6 orchestrated)
    headline = qc.run_funnel(a, cluster_key="cell_type", marker_sets=markers.LINEAGE_MARKERS)
    assert headline["n_cells"] == a.n_obs
    assert "annotatability" in headline
    assert "6_spatial_coherence" in headline["layers_run"]

    assert 0.0 <= purity.crisp_purity(a, markers.LINEAGE_MARKERS) <= 1.0
    assert purity.mecr(a, markers.LINEAGE_MARKERS) >= 0.0

    # Cell states (H6) runs (may be empty on the tiny synthetic panel) + VSI join (H5)
    from spatialscribe.analysis import states

    smat = states.state_by_celltype(a, cluster_key="cell_type")
    assert smat is not None
    import os
    import tempfile

    import pandas as pd

    vp = os.path.join(tempfile.mkdtemp(), "vsi.parquet")
    pd.DataFrame({"cell_id": list(a.obs_names[:50]), "vsi": 0.9,
                  "frac_low_vsi": 0.0, "vsi_n_px": 5}).to_parquet(vp)
    info = qc.apply_ovrlpy_vsi(a, vp)
    assert "vsi" in a.obs and info["threshold"] == 0.7

    # Spatial neighborhood enrichment + co-occurrence + niches
    from spatialscribe.analysis import niches

    a.obs["cell_type"] = a.obs["cell_type"].astype("category")
    res = spatial.nhood_enrichment(a, cluster_key="cell_type")
    assert len(res["categories"]) >= 1
    occ = spatial.co_occurrence(a, cluster_key="cell_type")
    assert "occ" in occ
    ntab = niches.call_niches(a, cluster_key="cell_type", n_niches=4)
    assert "niche" in a.obs and len(ntab) >= 1

    # Subcluster the largest cell type (no LLM)
    top_ct = a.obs["cell_type"].value_counts().idxmax()
    sub, rows = subcluster.subcluster(a, top_ct, use_llm=False)
    assert "subtype" in a.obs and len(rows) >= 1

    # Exports round-trip
    import os
    import tempfile

    d = tempfile.mkdtemp()
    assert export.export_h5ad(a, os.path.join(d, "out.h5ad")).exists()
    scr = export.export_script([{"label": "load", "code": "x = 1"}], os.path.join(d, "run.py"))
    assert scr.exists() and "x = 1" in scr.read_text()


def test_annotation_methods_tool_reference_free_and_with_methods():
    import anndata, numpy as np, pandas as pd
    from spatialscribe.analysis import capabilities as cap

    a = anndata.AnnData(X=np.zeros((6, 3), dtype="float32"))
    a.obs_names = [f"c{i}" for i in range(6)]
    # (a) reference-free: no method columns -> empty methods_run, no crash
    out = cap.run(a, "annotation_methods", {}, cap.RunContext(tissue="melanoma")).value
    assert isinstance(out.get("methods_run"), list) and out["methods_run"] == []

    # (b) with method columns present, the tool reports popV consensus recomputed over the ACTUAL method
    #     labels (vs the final cell_type), not a pre-supplied score. 3 methods -> a trusted ensemble.
    a.obs["cell_type"] = pd.Categorical(["T", "Mel", "T", "B", "B", "T"])
    a.obs["rctd_first_type"] = ["T", "Mel", "T", None, "B", "T"]      # one None -> coverage 5/6
    a.obs["singler_label"] = ["T", "Mel", "Mel", "B", "B", "T"]       # c2 disagrees with cell_type 'T'
    a.obs["scanvi_label"] = ["T", "Mel", "B", "B", "B", "T"]          # c2 disagrees too
    out2 = cap.run(a, "annotation_methods", {}, cap.RunContext(tissue="melanoma")).value
    assert set(out2["methods_run"]) == {"rctd", "singler", "scanvi"}
    assert abs(out2["coverage"]["rctd"] - 5 / 6) < 1e-6       # one None
    assert out2["coverage"]["singler"] == 1.0
    cons = out2["consensus"]
    assert 0.0 <= cons["mean_agreement"] <= 1.0
    assert isinstance(cons["reliability"], dict) and cons["reliability"]     # reliability-bin distribution
    assert cons["pct_trusted_ensemble"] > 0.0                               # >=3 methods voted

    # (c) the consensus reflects REAL cross-method agreement (written to obs by the tool): c0 (all three
    #     agree on 'T') is very_high; c2 (T vs Mel vs B, winner 'T') has only 1/3 agreeing -> low.
    assert str(a.obs["consensus_reliability"].iloc[0]) == "very_high"
    assert int(a.obs["consensus_score"].iloc[2]) == 1
    assert str(a.obs["consensus_reliability"].iloc[2]) == "low"
