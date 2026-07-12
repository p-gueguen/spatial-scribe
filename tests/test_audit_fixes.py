"""Regressions for the app-audit round (bug-classes analogous to the subcluster/typability/count-floor fixes).

- labels.match_cell_type: shared tolerant resolver (case / plural / substring), reused by every cap
  that takes an LLM-supplied cell-type label.
- immune_exclusion: resolves a near-miss label instead of hard-erroring (the copilot passes "T cells").
- segmentation_qc: guards a degenerate cell_area distribution + inclusive boundaries (count-floor twin).
- copilot map_view: obs-producing caps (niches, ...) emit a recolour directive so the chat map updates.
- suggest_count_floor exposes the key "floor" (a `.get("count_floor")` typo silently pinned the floor to 10).
"""
from __future__ import annotations

import anndata as ad
import numpy as np
import pytest

from spatialscribe.analysis import labels, qc


# ---- shared resolver -----------------------------------------------------------------------------
def test_match_cell_type_tolerant():
    labels.demo()   # its own asserts (exact / case / plural / substring / absent / ambiguous)
    cats = ["T cell", "Epithelial/Tumor"]
    assert labels.match_cell_type("T cells", cats) == "T cell"
    assert labels.match_cell_type("tumor", cats) == "Epithelial/Tumor"
    assert labels.match_cell_type("neuron", cats) is None


# ---- segmentation_qc degeneracy ------------------------------------------------------------------
def _adata_with_area(area):
    a = ad.AnnData(X=np.zeros((len(area), 2), dtype="float32"))
    a.obs["cell_area"] = np.asarray(area, dtype=float)
    return a


def test_segmentation_qc_constant_area_flags_nothing():
    out = qc.segmentation_qc(_adata_with_area([10.0] * 200))   # constant -> degenerate span
    assert out["pct_area_outlier"] == 0.0                       # not half the section


def test_segmentation_qc_zero_area_block_flagged_small():
    # a block of cells at 0 area: the 0.5-percentile lands on 0, a strict `< lo` would miss them; the
    # explicit `area <= 0` catches the empty fragments.
    area = np.concatenate([np.zeros(10), np.full(190, 50.0)])
    a = _adata_with_area(area)
    qc.segmentation_qc(a)
    assert (a.obs["seg_area_flag"].astype(str) == "small").sum() >= 10   # the 0-area block is flagged


def test_segmentation_qc_modal_max_not_overflagged():
    # 90% of cells share the top area (a mode at hi): a naive inclusive `area >= hi` would flag them
    # all "large". Strict `> hi` must not.
    area = np.concatenate([np.zeros(20), np.full(180, 50.0)])
    a = _adata_with_area(area)
    qc.segmentation_qc(a)
    assert (a.obs["seg_area_flag"].astype(str) == "large").sum() == 0       # the modal 50s are not "large"


def test_segmentation_qc_warns_on_heavy_outlier_tail():
    area = np.concatenate([np.zeros(20), np.full(180, 50.0)])   # 10% at 0 area -> heavy tail
    out = qc.segmentation_qc(_adata_with_area(area))
    assert out["status"] == "warn" and out["pct_area_outlier"] > 0.05   # dot reflects the tail, not always green


# ---- count-floor key contract (the views.py typo) ------------------------------------------------
def test_suggest_count_floor_uses_floor_key():
    a = ad.AnnData(X=np.zeros((100, 3), dtype="float32"))
    a.obs["total_counts"] = np.arange(100, dtype=float)
    a.uns["n_panel_genes"] = 5101
    sf = qc.suggest_count_floor(a)
    assert "floor" in sf and "count_floor" not in sf     # views.select_cells must read "floor"


# ---- immune_exclusion tolerant labels + map_view emission ----------------------------------------
def test_immune_exclusion_resolves_near_miss_and_errors_helpfully(processed_adata):
    from spatialscribe.analysis import spatial as _sp
    cats = list(processed_adata.obs["cell_type"].astype(str).unique())
    if len(cats) < 2:
        pytest.skip("need >=2 cell types")
    # near-miss wording (lowercase + trailing 's') must resolve, not error
    a, b = cats[0], cats[-1]
    res = _sp.immune_exclusion(processed_adata, a.lower() + "s", b.lower() + "s")
    assert "error" not in res and "zscore" in res
    # a genuinely absent label errors, and the error lists what IS available
    bad = _sp.immune_exclusion(processed_adata, "Klingon", b)
    assert "error" in bad and "not found" in bad["error"]


def test_niches_emits_map_view(processed_adata, ctx):
    from spatialscribe.analysis import capabilities as cap
    res = cap.run(processed_adata, "niches", {}, ctx)
    assert res.ok, res.error
    assert "niche" in processed_adata.obs
    assert any(x.get("kind") == "map_view" and x.get("color_by") == "niche" for x in ctx.artifacts)


# ---- pipeline runs the six-layer QC funnel on EVERY route (not only when annotation produced labels) --
def test_run_funnel_runs_section_layers_without_cell_type():
    # Layers 0-2 (section / segmentation / count-floor) never needed labels; run_funnel must run them
    # even with no cell_type - the enabler for running qc_funnel on the clusters-only route.
    a = ad.AnnData(X=np.ones((60, 4), dtype="float32"))
    a.obs["total_counts"] = np.arange(60, dtype=float) + 5
    a.obs["cell_area"] = np.full(60, 30.0)
    a.uns["n_panel_genes"] = 5101
    out = qc.run_funnel(a)                          # no cell_type in obs
    assert isinstance(out, dict) and out.get("n_cells") == 60   # ran the section layers, no crash


def test_pipeline_runs_qc_funnel_on_every_route():
    # Regression: qc_funnel used to live inside `if cell_type`, so a clusters-only / annotate-failed
    # section got no section QC at all. It must now appear in the stage record regardless.
    import pandas as pd
    from spatialscribe.analysis import capabilities as cap
    from spatialscribe.analysis import io as _io
    from spatialscribe.analysis import pipeline as pl
    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in range(20)] + ["NegControlProbe_00001", "BLANK_0001"]
    X = rng.poisson(1.0, (200, len(genes))).astype("float32")
    var = pd.DataFrame(index=genes)
    var["feature_types"] = ["Gene Expression"] * (len(genes) - 2) + ["Negative Control Probe"] * 2
    a = ad.AnnData(X=X, var=var)
    a.obs_names = [f"c{i}" for i in range(200)]
    a.obsm["spatial"] = rng.normal(0, 10, (200, 2))
    a.var["control"] = _io.build_control_mask(a)
    a.layers["counts"] = a.X.copy()
    rec = pl.run_pipeline(a, cap.RunContext(tissue="melanoma", use_llm=False))
    assert "qc_funnel" in [s["name"] for s in rec["stages"]]   # runs on every route now
