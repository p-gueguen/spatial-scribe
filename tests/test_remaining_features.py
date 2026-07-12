"""Tests for the features added in the 'remaining tasks' pass:

- custom-reference plumbing: RunContext.resolve_reference, reference_transfer capability,
  annotate.celltypist_transfer;
- malignant_concordance (learned callers degrade gracefully off-env);
- self-verify auto-rerun loop (verify.autorerun + self_heal capability);
- the report deliverables: render_analysis_report (self-contained HTML) + 300-DPI PNGs.

All offline / CPU-only (celltypist is in the main env; the subprocess callers degrade to 'skipped').
"""
from __future__ import annotations

import numpy as np


def _labeled_reference(section, n_per: int = 80, k: int = 3, seed: int = 0):
    """A small block-structured labeled reference on the SECTION's genes (so genes overlap 100% and
    a classifier can learn the classes). Returns (AnnData, 'cell_type')."""
    import anndata as ad
    import pandas as pd

    genes = [str(g) for g in section.var_names]
    ng = len(genes)
    rng = np.random.default_rng(seed)
    block = max(1, ng // (k + 1))
    X, lab = [], []
    for ci in range(k):
        base = rng.poisson(0.5, (n_per, ng)).astype("float32")
        base[:, ci * block:(ci + 1) * block] += rng.poisson(9, (n_per, block))
        X.append(base)
        lab += [f"RefType{ci}"] * n_per
    a = ad.AnnData(X=np.vstack(X))
    a.var_names = genes
    a.obs_names = [f"ref{i}" for i in range(a.n_obs)]
    a.layers["counts"] = a.X.copy()
    a.obs["cell_type"] = pd.Categorical(lab)
    return a, "cell_type"


# --------------------------------------------------------------------------- #
# RunContext.resolve_reference
# --------------------------------------------------------------------------- #
def test_resolve_reference_none():
    from spatialscribe.analysis import capabilities as cap
    ref, key, src = cap.RunContext().resolve_reference()
    assert ref is None and key is None and "no reference" in src


def test_resolve_reference_in_memory(processed_adata):
    from spatialscribe.analysis import capabilities as cap
    r, key = _labeled_reference(processed_adata)
    ctx = cap.RunContext(tissue="melanoma", reference=r, ref_label_key=key)
    ref, k, src = ctx.resolve_reference()
    assert ref is r and k == key and "Panel-check step" in src


# --------------------------------------------------------------------------- #
# reference_transfer capability + celltypist_transfer
# --------------------------------------------------------------------------- #
def test_reference_transfer_no_reference(processed_adata, ctx):
    from spatialscribe.analysis import capabilities as cap
    out = cap.run(processed_adata, "reference_transfer", {}, ctx)
    assert out.ok and out.value["status"] == "no_reference"


def test_reference_transfer_with_reference(processed_adata):
    from spatialscribe.analysis import capabilities as cap
    r, key = _labeled_reference(processed_adata)
    ctx = cap.RunContext(tissue="melanoma", reference=r, ref_label_key=key)
    out = cap.run(processed_adata, "reference_transfer", {}, ctx)
    assert out.ok, out.error
    v = out.value
    assert v["status"] == "ok"
    methods = {a["method"]: a.get("status") for a in v["arms"]}
    assert methods.get("celltypist") == "ok"                  # in-env arm ran
    # tacco is now a DECLARED dependency (pyproject extra "transfer", which pixi pulls via
    # extras=["transfer"]), so the documented default transfer engine must actually RUN here.
    # Asserting mere presence let a missing tacco pass silently: capabilities.py catches the
    # ImportError and records status="skipped: ...", the suite stayed green, and the app quietly
    # annotated with the celltypist fallback while the docs claimed TACCO was the default.
    assert methods.get("tacco") == "ok", f"tacco arm did not run: {methods.get('tacco')!r}"
    assert "celltypist_label" in processed_adata.obs


def test_celltypist_transfer_graceful_no_overlap(processed_adata):
    """Disjoint gene names -> the arm skips (never raises)."""
    from spatialscribe.analysis import annotate
    import anndata as ad
    import pandas as pd
    r = ad.AnnData(X=np.random.default_rng(0).poisson(1.0, (60, 10)).astype("float32"))
    r.var_names = [f"ZZZ{i}" for i in range(10)]
    r.obs["cell_type"] = pd.Categorical(["a", "b"] * 30)
    res = annotate.celltypist_transfer(processed_adata, r, "cell_type")
    assert res["status"].startswith("skipped")


def _random_label_reference(section, n=240, k=3, seed=1):
    """A reference on the section's genes but with labels UNRELATED to expression (uniform counts +
    random labels) - a deliberately POOR fit (low depth-matched F1)."""
    import anndata as ad
    import pandas as pd
    genes = [str(g) for g in section.var_names]
    rng = np.random.default_rng(seed)
    a = ad.AnnData(X=rng.poisson(1.0, (n, len(genes))).astype("float32"))
    a.var_names = genes
    a.obs_names = [f"u{i}" for i in range(n)]
    a.layers["counts"] = a.X.copy()
    a.obs["cell_type"] = pd.Categorical(rng.integers(0, k, n).astype(str))
    return a, "cell_type"


# --------------------------------------------------------------------------- #
# reference<->panel single matching-score + recommendation
# --------------------------------------------------------------------------- #
def test_reference_match_score_good(processed_adata):
    from spatialscribe.analysis import reference as _ref
    r, key = _labeled_reference(processed_adata)
    m = _ref.reference_panel_match(r, [str(g) for g in processed_adata.var_names], key)
    assert m["status"] == "ok"
    assert 0.0 <= m["match_score"] <= 1.0
    assert m["match_score"] == m["global"]["match_score"]           # mirrored
    rec = m["recommendation"]
    assert rec["action"] in ("supervised_transfer", "coarsen_then_transfer", "unsupervised_clustering")
    assert rec["verdict"] == m["global"]["verdict"]


def test_reference_match_poor_recommends_clustering(processed_adata):
    from spatialscribe.analysis import reference as _ref
    r, key = _random_label_reference(processed_adata)
    m = _ref.reference_panel_match(r, [str(g) for g in processed_adata.var_names], key)
    # labels unrelated to expression -> poor depth-matched resolvability -> steer away from supervised
    assert m["match_score"] <= 0.6
    assert m["recommendation"]["action"] in ("unsupervised_clustering", "coarsen_then_transfer")


def test_reference_match_cap_surfaces_score(processed_adata):
    from spatialscribe.analysis import capabilities as cap
    r, key = _labeled_reference(processed_adata)
    ctx = cap.RunContext(tissue="melanoma", reference=r, ref_label_key=key)
    out = cap.run(processed_adata, "reference_match", {}, ctx)
    assert out.ok and out.value["match_score"] is not None
    assert "recommendation" in out.value and "alternatives" in out.value


# --------------------------------------------------------------------------- #
# malignant_concordance (graceful degradation off-env)
# --------------------------------------------------------------------------- #
def test_malignant_concordance_non_tumour(processed_adata):
    from spatialscribe.analysis import cnv
    res = cnv.call_malignant_concordance(processed_adata, tissue="kidney")
    assert res["status"].startswith("skipped") and res["tumour_context"] is False


def test_malignant_concordance_tumour_marker_only(processed_adata):
    from spatialscribe.analysis import capabilities as cap
    out = cap.run(processed_adata, "malignant_concordance", {}, cap.RunContext(tissue="melanoma"))
    assert out.ok, out.error
    v = out.value
    assert v["status"] == "ok" and v["tumour_context"] is True
    names = {c["name"]: c["status"] for c in v["callers"]}
    assert names.get("marker score") == "ok"                  # always-on marker path ran
    assert any(s.startswith("skipped") for n, s in names.items() if n != "marker score")  # learned callers off-env
    assert "malignant_score" in processed_adata.obs           # declared produces satisfied


# --------------------------------------------------------------------------- #
# self-verify auto-rerun loop
# --------------------------------------------------------------------------- #
def test_self_heal_runs(processed_adata, ctx):
    from spatialscribe.analysis import capabilities as cap
    out = cap.run(processed_adata, "self_heal", {"max_rounds": 1}, ctx)
    assert out.ok, out.error
    v = out.value
    assert v["status"] == "ok"
    for k in ("initial_failed", "final_failed", "rounds", "advisory_remaining", "n_auto_actions"):
        assert k in v


def test_autorerun_never_raises_without_annotation(raw_adata, ctx):
    """No cell_type -> verify finds nothing to fix -> a clean, empty result (no crash)."""
    from spatialscribe.analysis import verify
    raw_adata.obs["cell_type"] = "A"                          # single trivial label
    res = verify.autorerun(raw_adata, ctx, max_rounds=1)
    assert res["status"] == "ok"


# --------------------------------------------------------------------------- #
# report deliverables
# --------------------------------------------------------------------------- #
def test_render_analysis_report_self_contained(processed_adata, tmp_path):
    from spatialscribe.analysis import export
    p = export.render_analysis_report(processed_adata, action_log=[{"name": "annotate", "params": {}}],
                                      tissue="melanoma", path=tmp_path / "report.html")
    html = p.read_text()
    assert p.exists() and len(html) > 1000
    assert "<img" in html                                     # inline figure(s)
    assert "data:image/png;base64," in html                  # self-contained (no external asset)
    assert "http://" not in html and "https://" not in html  # no external references
    assert "Re-runnable analysis" in html


def test_save_report_pngs_300dpi(processed_adata, tmp_path):
    from spatialscribe.analysis import export
    pngs = export.save_report_pngs(processed_adata, tmp_path, dpi=300)
    assert pngs, "expected at least one figure PNG"
    for p in pngs:
        assert p.exists() and p.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"   # real PNG header


# --------------------------------------------------------------------------- #
# CyteType-style cell-state typing (states.assign_cell_states)
# --------------------------------------------------------------------------- #
def _state_structured_section():
    """A section where the first block of cells is clearly cycling and the second hypoxic."""
    import anndata as ad
    import pandas as pd
    rng = np.random.default_rng(0)
    n = 240
    sg = ["MKI67", "TOP2A", "PCNA", "CCNB1", "CDK1", "BIRC5", "UBE2C",  # cycling
          "VEGFA", "CA9", "HIF1A", "SLC2A1", "NDRG1", "LDHA"]          # hypoxia
    genes = sg + [f"BG{i}" for i in range(40)]
    X = rng.poisson(0.5, (n, len(genes))).astype("float32")
    X[:120, :7] += rng.poisson(10, (120, 7))     # cycling block
    X[120:, 7:13] += rng.poisson(10, (120, 6))   # hypoxia block
    a = ad.AnnData(X=X.copy()); a.var_names = genes
    a.obs_names = [f"c{i}" for i in range(n)]; a.layers["counts"] = X.copy()
    a.obs["cell_type"] = pd.Categorical(["Malignant"] * 120 + ["Endothelial"] * 120)
    return a


def test_assign_cell_states_types_dominant_program():
    from spatialscribe.analysis import states
    a = _state_structured_section()
    res = states.assign_cell_states(a, cluster_key="cell_type")
    assert "cell_state" in a.obs                                      # per-cell label written
    assert {"Cycling", "Hypoxia"} <= set(res["distribution"])         # both programs recovered
    tops = {ct: d["top_state"] for ct, d in res["per_celltype"].items()}
    assert tops["Malignant"] == "Cycling" and tops["Endothelial"] == "Hypoxia"


def test_assign_cell_states_no_signatures_graceful(raw_adata):
    """A panel with no state-program genes -> every cell 'None', obs column still written."""
    from spatialscribe.analysis import states
    raw_adata.obs["cell_type"] = "A"
    res = states.assign_cell_states(raw_adata, cluster_key="cell_type")
    assert "cell_state" in raw_adata.obs
    assert set(res["distribution"]) <= {"None", *res["states"]}


def test_assign_cell_states_capability_produces_key(processed_adata, ctx):
    from spatialscribe.analysis import capabilities as cap
    out = cap.run(processed_adata, "assign_cell_states", {}, ctx)
    assert out.ok, out.error
    assert "cell_state" in processed_adata.obs                        # declared produces satisfied
    assert "distribution" in out.value and "per_celltype" in out.value


# --------------------------------------------------------------------------- #
# Computable per-cell-type typability (panel_check.typability_table) - not gene counts
# --------------------------------------------------------------------------- #
def test_typability_auc_basis_without_reference(processed_adata):
    from spatialscribe.analysis import panel_check as pc, capabilities as cap
    cap.run(processed_adata, "panel_check", {}, cap.RunContext(tissue="melanoma"))
    rows = pc.typability_table(processed_adata, cluster_key="cell_type")
    assert rows and all("confidently_typable" in r for r in rows)
    bases = {r["basis"] for r in rows}
    assert "identifiability_auc" in bases                             # the computed-AUC path is used
    assert "depth_matched_f1" not in bases                            # no reference supplied
    # AUC-based rows are computed evidence ONLY when the type has enough cells: a one-vs-rest AUC on
    # <~50 positives is noise and is flagged weak_evidence (the documented guard, now enforced in code).
    for r in rows:
        if r["basis"] == "identifiability_auc":
            assert "n_cells" in r
            assert r["weak_evidence"] is ((r.get("n_cells") or 0) < 50)


def test_typability_small_n_auc_is_weak_and_not_typable(processed_adata):
    # A handful of cells cannot be "confidently typable" no matter how high the one-vs-rest AUC: the
    # bundled demo's B/Plasma (~30 cells) must be weak_evidence and NOT counted confidently_typable.
    from spatialscribe.analysis import panel_check as pc, capabilities as cap
    cap.run(processed_adata, "panel_check", {}, cap.RunContext(tissue="melanoma"))
    rows = pc.typability_table(processed_adata, cluster_key="cell_type")
    small = [r for r in rows if r["basis"] == "identifiability_auc" and (r.get("n_cells") or 0) < 50]
    assert small, "expected at least one <50-cell AUC row in the demo (B/Plasma)"
    for r in small:
        assert r["weak_evidence"] is True
        assert r["confidently_typable"] is False


def test_typability_zero_cell_coverage_not_typable(processed_adata):
    # A lineage with green marker coverage but 0 assigned cells is a hypothetical, not a confirmed
    # typable type - it must not be counted confidently_typable (was: green coverage => typable).
    from spatialscribe.analysis import panel_check as pc, capabilities as cap
    cap.run(processed_adata, "panel_check", {}, cap.RunContext(tissue="melanoma"))
    rows = pc.typability_table(processed_adata, cluster_key="cell_type")
    for r in rows:
        if (r.get("n_cells") or 0) == 0:
            assert r["confidently_typable"] is False, r


def test_typability_f1_basis_with_reference(processed_adata):
    from spatialscribe.analysis import panel_check as pc, reference as ref, capabilities as cap
    cap.run(processed_adata, "panel_check", {}, cap.RunContext(tissue="melanoma"))
    r, key = _labeled_reference(processed_adata)
    rm = ref.reference_panel_match(r, [str(g) for g in processed_adata.var_names], key)
    rows = pc.typability_table(processed_adata, cluster_key="cell_type", reference_match=rm)
    assert any(x["basis"] == "depth_matched_f1" for x in rows)         # reference -> depth-matched F1 wins


def test_typability_coverage_fallback_is_flagged_weak(processed_adata):
    from spatialscribe.analysis import panel_check as pc, capabilities as cap
    cap.run(processed_adata, "panel_check", {}, cap.RunContext(tissue="melanoma"))
    a = processed_adata.copy()
    del a.obs["cell_type"]                                             # no labels -> AUC impossible
    rows = pc.typability_table(a, cluster_key="cell_type")
    assert rows and all(r["basis"] == "coverage_only" and r["weak_evidence"] for r in rows)


def test_panel_verdict_grounds_on_computed_typability(processed_adata, monkeypatch):
    from spatialscribe.analysis import llm, panel_check as pc, capabilities as cap
    cap.run(processed_adata, "panel_check", {}, cap.RunContext(tissue="melanoma"))
    typ = pc.typability_table(processed_adata, cluster_key="cell_type")
    captured = {}
    monkeypatch.setattr(llm, "complete",
                        lambda system, user, **k: captured.update(system=system, user=user) or "verdict")
    out = llm.panel_verdict(processed_adata.uns["panel_check"], typability=typ)
    assert out == "verdict"
    assert "typability_per_cell_type" in captured["user"]              # the computed table is in the prompt
    assert "confidently_typable" in captured["system"] and "depth_matched_f1" in captured["system"]


# --------------------------------------------------------------------------- #
# SPLIT spillover purification (in-app) + marker dot-plot before/after
# --------------------------------------------------------------------------- #
def test_split_purify_reference_free_fallback(processed_adata, ctx):
    """No SPLIT/rctd envs and no reference -> the reference-free main-env fallback runs (never
    raises), writes the split_corrected layer, and reports it is reference-free."""
    from spatialscribe.analysis import capabilities as cap
    out = cap.run(processed_adata, "split_purify", {}, ctx)
    assert out.ok and out.value["status"] == "ok"
    assert "reference-free" in (out.value.get("note") or "")
    assert "split_corrected" in processed_adata.layers


def test_composition_prior_is_cellspace_and_reflects_the_section():
    """The data-driven TACCO prior is a cell-fraction Series over the reference types, built from the
    section's own CellTypist labels - not the reference's (possibly capped) composition. This is what
    flips the TACCO arm from 0.35 to 0.69 accuracy on a depth-capped reference."""
    import anndata as ad
    import numpy as np
    import pandas as pd

    from spatialscribe.analysis import reference_transfer as rt
    # reference: 3 types present; section is 80% type A (the "tumour-heavy" case the default prior misses)
    ref = ad.AnnData(X=np.ones((6, 4), dtype="float32"))
    ref.obs["cell_type"] = ["A", "A", "B", "B", "C", "C"]
    sec = ad.AnnData(X=np.ones((10, 4), dtype="float32"))
    sec.obs["celltypist_label"] = ["A"] * 7 + ["B", "B", "C"]
    prior = rt.composition_prior(sec, ref, "cell_type")
    assert prior is not None
    assert set(prior.index) == {"A", "B", "C"}
    assert abs(float(prior.sum()) - 1.0) < 1e-9          # a proper distribution
    assert prior["A"] > prior["B"] > prior["C"]          # tracks the section, not the uniform reference
    assert prior["A"] > 0.5                               # 7/10 -> the dominant type dominates the prior
    # no CellTypist column -> None (falls back to TACCO's reference-derived default, honestly)
    assert rt.composition_prior(ad.AnnData(X=np.ones((3, 4), dtype="float32")), ref, "cell_type") is None


def test_split_purify_defaults_to_tacco_not_rctd(processed_adata, monkeypatch):
    """With a reference + R env but NO rctd_python, split_purify must take the TACCO weights path
    (not skip for lack of an rctd-py env). We stub the heavy runner to prove which engine was chosen
    without needing TACCO or R at test time."""
    import anndata as ad

    from spatialscribe.analysis import split
    chosen = {}

    def _fake_purify(adata, reference, ref_label_key, *, engine, rctd_python, rscript, max_cells):
        chosen["engine"] = engine
        chosen["rctd_python"] = rctd_python
        return {"status": "ok", "method": f"{engine}_split", "pct_removed": 0.5}

    monkeypatch.setattr(split, "_purify_split", _fake_purify)
    monkeypatch.setenv("SPATIALSCRIBE_SPLIT_RSCRIPT", __file__)  # any existing path
    monkeypatch.delenv("SPATIALSCRIBE_SPLIT_RCTD_PYTHON", raising=False)
    monkeypatch.delenv("SPATIALSCRIBE_RCTD_PYTHON", raising=False)
    ref = ad.AnnData(X=processed_adata.X[:5].copy())
    ref.obs["cell_type"] = ["A", "B", "A", "B", "A"]
    out = split.split_purify(processed_adata, reference=ref, ref_label_key="cell_type")
    assert out["method"] == "tacco_split"          # default engine, no rctd-py env needed
    assert chosen["engine"] == "tacco" and chosen["rctd_python"] is None


def test_split_purify_rctd_engine_opt_in(processed_adata, monkeypatch, tmp_path):
    """weights_engine='rctd' takes the RCTD path only when an rctd-py env is configured."""
    import anndata as ad

    from spatialscribe.analysis import split
    chosen = {}

    def _fake(adata, reference, ref_label_key, *, engine, rctd_python, rscript, max_cells):
        chosen["engine"] = engine
        return {"status": "ok", "method": f"{engine}_split", "pct_removed": 0.1}

    monkeypatch.setattr(split, "_purify_split", _fake)
    fake_py = tmp_path / "py"; fake_py.write_text("")
    monkeypatch.setenv("SPATIALSCRIBE_SPLIT_RSCRIPT", __file__)
    monkeypatch.setenv("SPATIALSCRIBE_RCTD_PYTHON", str(fake_py))
    ref = ad.AnnData(X=processed_adata.X[:5].copy())
    ref.obs["cell_type"] = ["A", "B", "A", "B", "A"]
    out = split.split_purify(processed_adata, reference=ref, ref_label_key="cell_type",
                             weights_engine="rctd")
    assert out["method"] == "rctd_split" and chosen["engine"] == "rctd"


def test_split_join_aligns_purified_counts(tmp_path):
    """_join_split maps purified_counts.mtx (genes x cells) onto layers['split_corrected'] by name."""
    import anndata as ad
    import pandas as pd
    import scipy.io as sio
    import scipy.sparse as sp
    from spatialscribe.analysis import split
    a = ad.AnnData(X=np.zeros((10, 8), dtype="float32"))
    a.var_names = [f"G{i}" for i in range(8)]; a.obs_names = [f"c{i}" for i in range(10)]
    inp, sd = tmp_path / "inputs", tmp_path / "split"; inp.mkdir(); sd.mkdir()
    cells = [f"c{i}" for i in range(10)]; pgenes = ["G1", "G3", "G5"]
    pd.Series(cells).to_csv(inp / "cells.txt", header=False, index=False)
    pd.Series(pgenes).to_csv(sd / "purified_genes.txt", header=False, index=False)
    sio.mmwrite(str(sd / "purified_counts.mtx"),                     # genes x cells
                sp.csr_matrix(np.arange(30).reshape(3, 10).astype(float)))
    n = split._join_split(a, inp, sd)
    assert n == 10 and "split_corrected" in a.layers
    L = a.layers["split_corrected"].toarray()
    assert L[1, a.var_names.get_loc("G3")] == 11.0                  # gene G3 (row1), cell c1 (col1) = 11
    assert L[0, a.var_names.get_loc("G0")] == 0.0                   # a gene not purified stays 0


def test_marker_dotplot_layer_and_corrected(processed_adata, ctx):
    from spatialscribe.analysis import capabilities as cap, views
    import scipy.sparse as sp
    genes = [str(g) for g in processed_adata.var_names[:3]]
    # a layer-based dot-plot builds
    processed_adata.layers["split_corrected"] = sp.csr_matrix(processed_adata.X)
    assert views.dotplot_view(processed_adata, genes, layer="split_corrected") is not None
    ctx.artifacts.clear()
    out = cap.run(processed_adata, "marker_dotplot", {"genes": genes, "corrected": True}, ctx)
    assert out.ok and out.value.get("compared") is True
    assert len(ctx.artifacts) == 2                                  # raw + purified figures
