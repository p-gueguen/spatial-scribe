"""PLAID-style per-cell scoring of the de-novo NMF gene programs.

PLAID (Pathway-Level Average Intensity Detection) is the fastest single-sample enrichment
primitive: the (standardized) average expression of a gene set per cell. Here each discovered
program's top-loading genes ARE the gene set; the score says how strongly each cell expresses
that program's signature. The correctness claim that matters: cells the NMF *assigned* to a
program must score higher on that program's own signature than the rest of the section, and the
score must separate them well above chance (AUROC). Two well-separated biologies -> two clean,
specific programs.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")
pytest.importorskip("scanpy")


def _two_program_section(n=400, seed=0):
    """T cells vs melanocytes: two disjoint marker sets, each high in half the cells."""
    rng = np.random.default_rng(seed)
    t = ["CD3D", "CD3E", "CD8A", "TRAC"]
    mel = ["MLANA", "PMEL", "TYR", "DCT", "MITF", "SOX10"]
    noise = [f"N{i}" for i in range(30)]
    genes = t + mel + noise
    grp = np.array([0] * (n // 2) + [1] * (n - n // 2))
    X = rng.poisson(0.3, (n, len(genes))).astype("float32")
    gi = {g: i for i, g in enumerate(genes)}
    for g in t:
        X[grp == 0, gi[g]] += rng.poisson(8, (grp == 0).sum())
    for g in mel:
        X[grp == 1, gi[g]] += rng.poisson(8, (grp == 1).sum())
    a = anndata.AnnData(X=X, var=pd.DataFrame(index=genes))
    a.obsm["spatial"] = rng.normal(0, 50, (n, 2))
    return a


def test_score_programs_cells_score_high_on_own_program():
    import scanpy as sc

    from spatialscribe.analysis import programs

    a = _two_program_section()
    sc.pp.normalize_total(a)
    sc.pp.log1p(a)
    programs.discover_programs(a, n_programs=2)

    rows = programs.score_programs(a, top_n=8)

    # writes a cells x programs score matrix
    assert "program_scores" in a.obsm
    assert a.obsm["program_scores"].shape == (a.n_obs, 2)

    # one row per program, each carrying the scoring summary the hover renders
    assert len(rows) == 2
    for r in rows:
        assert {"program", "plaid_in", "plaid_out", "specificity"} <= set(r)
        # a cell assigned to a program over-expresses that program's own signature
        assert r["plaid_in"] > r["plaid_out"]
        # and the score separates the program's cells from the rest well above chance
        assert r["specificity"] > 0.8
        assert 0.0 <= r["specificity"] <= 1.0


def test_score_programs_equals_full_matrix_zscore_mean():
    """Pin the exact numeric contract: each program's per-cell score is the mean of the
    per-gene z-scores over that program's top-loading genes (positive preferred). This is
    the invariant the subset-first optimization must preserve byte-for-byte - z-scoring is
    per-gene independent, so z-scoring only the selected genes is identical to z-scoring the
    whole matrix and then subsetting. Guards the optimization against any silent drift."""
    import scanpy as sc

    from spatialscribe.analysis import programs

    a = _two_program_section()
    sc.pp.normalize_total(a)
    sc.pp.log1p(a)
    programs.discover_programs(a, n_programs=2)

    top_n = 8
    programs.score_programs(a, top_n=top_n)
    got = a.obsm["program_scores"]

    # independent ground truth: z-score EVERY gene across cells, then average each program's
    # top-loading genes - the full-matrix computation the optimization replaces.
    X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
    X = X.astype(float)
    sd = X.std(0)
    Z = (X - X.mean(0)) / np.where(sd == 0, 1.0, sd)
    load = np.asarray(a.varm["program_loadings"])
    expected = np.zeros_like(got)
    for k in range(load.shape[1]):
        order = np.argsort(-load[:, k])
        sig = [int(i) for i in order[:top_n] if load[i, k] > 0] or [int(i) for i in order[:top_n]]
        expected[:, k] = Z[:, sig].mean(1)

    assert got.shape == expected.shape
    assert np.allclose(got, expected, atol=1e-6)


def test_discover_programs_folds_in_scoring():
    """discover_programs returns the PLAID scoring inline so the same table row carries it."""
    import scanpy as sc

    from spatialscribe.analysis import programs

    a = _two_program_section()
    sc.pp.normalize_total(a)
    sc.pp.log1p(a)
    tab = programs.discover_programs(a, n_programs=2)

    assert {"plaid_in", "plaid_out", "specificity"} <= set(tab.columns)
    assert "program_scores" in a.obsm
    # the specificity is a real separation, not a constant
    assert (tab["specificity"] > 0.8).all()


def test_every_default_state_program_clears_the_on_panel_floor_for_a_targeted_panel():
    """Regression: "Antigen presentation" scoped to MHC-II alone left only CIITA on a Xenium Prime 5K
    panel (1 gene), so `states.score_states` silently skipped it on both bundled demos. The processing
    arm (CIITA/NLRC5/TAP1/TAP2/TAPBP/PSMB8/PSMB9/B2M) is what carries the axis on a targeted panel.
    Every default program must land >=2 genes on a panel that has no HLA-D genes."""
    from spatialscribe.analysis import markers as m

    panel = {"CIITA", "NLRC5", "TAP1", "TAP2", "MKI67", "TOP2A", "STAT1", "IRF7", "VEGFA", "HIF1A",
             "HSPA8", "HSPB1", "ZEB1", "SNAI2", "PDCD1", "LAG3", "GZMB", "PRF1", "MMP9", "LOX"}
    for program, genes in m.CELL_STATES.items():
        on_panel = m.present(panel, {program: genes})[program]
        assert len(on_panel) >= 2, f"{program} scores only {on_panel} - below the >=2 floor"


def test_state_programs_resolve_onto_a_mouse_title_case_panel():
    """One human-uppercase set serves both species: `markers._resolve_one` folds case unambiguously.
    The mouse MHC orthologs (H2-Aa, H2-K1) are NOT case-variants of the HLA symbols and are correctly
    dropped, so the antigen program must survive on its conserved machinery genes alone."""
    from spatialscribe.analysis import markers as m

    mouse_panel = {"Ciita", "Cd74", "Nlrc5", "Tap1", "Tap2", "Tapbp", "Psmb8", "Psmb9", "B2m",
                   "H2-Aa", "H2-K1", "Gzmb", "Prf1", "Mmp9", "Lox"}
    ap = m.present(mouse_panel, {"ap": m.CELL_STATES["Antigen presentation"]})["ap"]
    assert len(ap) >= 2 and "H2-Aa" not in ap and "H2-K1" not in ap
    assert m.present(mouse_panel, {"cy": m.CELL_STATES["T-cytotoxicity"]})["cy"] == ["Gzmb", "Prf1"]
