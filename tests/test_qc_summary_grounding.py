"""Regression: qc_summary must never silently zero-default missing QC columns.

The QC "AI verdict" is grounded in qc_summary. When total_counts was absent (verdict runs before the
QC step, or a GPU path didn't persist the columns), qc_summary defaulted to np.zeros -> median genes
0, median counts 0, fraction_empty 100% -> a false "section completely unusable" verdict on a
perfectly good section. It must compute the metrics from the raw counts instead.
"""

from __future__ import annotations


def test_qc_summary_computes_when_columns_absent():
    from spatialscribe.analysis import demo, qc

    a = demo.load_demo().adata
    assert "total_counts" not in a.obs                 # compute_qc has NOT run

    s = qc.qc_summary(a)
    assert s["median_genes_per_cell"] > 0              # computed from raw counts, NOT zero-defaulted
    assert s["median_transcripts_per_cell"] > 0
    assert s["fraction_empty_cells"] < 0.5             # NOT a false 100%-empty
    assert "total_counts" in a.obs                     # now populated for downstream steps


def test_qc_uses_counts_layer_when_X_is_normalized():
    # A pre-normalized section (log1p X) that still carries a raw 'counts' layer - the bundled demo
    # caches, or a section clustered before a QC re-run. QC depth must come from the counts, not from
    # log-space X (which reports median transcripts < median genes - biologically impossible).
    import anndata as ad
    import numpy as np
    import scanpy as sc

    from spatialscribe.analysis import qc

    rng = np.random.default_rng(0)
    counts = rng.poisson(3.0, size=(200, 60)).astype("float32")     # raw counts, deep enough
    a = ad.AnnData(X=counts.copy())
    a.var_names = [f"g{i}" for i in range(60)]
    a.layers["counts"] = counts.copy()
    sc.pp.normalize_total(a)
    sc.pp.log1p(a)                                                  # X is now log1p (uns['log1p'] set)

    s = qc.qc_summary(a)
    raw_tx = float(np.median(counts.sum(1)))
    assert abs(s["median_transcripts_per_cell"] - raw_tx) < 1e-6    # from counts, not log-space X
    assert s["median_transcripts_per_cell"] >= s["median_genes_per_cell"]   # no impossible inversion


def test_merscope_blank_controls_drive_pct_counts_control():
    # MERSCOPE keeps its Blank negative controls in obsm['blank_genes'] (var is panel-only, so the
    # var-based neg_control mask is all-False). compute_qc must derive pct_counts_control from the
    # blanks, not leave it a structural 0 that blinds the contamination QC on MERSCOPE.
    import anndata as ad
    import numpy as np

    from spatialscribe.analysis import qc

    rng = np.random.default_rng(0)
    genes = rng.poisson(5.0, size=(120, 30)).astype("float32")           # ~150 gene counts/cell
    blanks = np.zeros((120, 8), dtype="float32"); blanks[:, 0] = 3.0      # ~3 blank counts/cell
    a = ad.AnnData(X=genes)
    a.var_names = [f"g{i}" for i in range(30)]
    a.var["neg_control"] = False                                          # MERSCOPE: none in var
    a.obsm["blank_genes"] = blanks
    qc.compute_qc(a)

    gene_tot = genes.sum(1)
    expected = 100.0 * 3.0 / (gene_tot + 3.0)                             # % of counts on a control probe
    assert np.allclose(a.obs["pct_counts_control"].to_numpy(), expected, atol=1e-4)
    assert float(np.median(a.obs["pct_counts_control"])) > 0.5           # not a structural 0


def test_qc_summary_matches_after_explicit_compute_qc():
    from spatialscribe.analysis import demo, qc

    a1 = demo.load_demo().adata
    qc.compute_qc(a1)
    ref = qc.qc_summary(a1)                             # the canonical values

    a2 = demo.load_demo().adata                         # no compute_qc -> auto-computed inside qc_summary
    got = qc.qc_summary(a2)
    for k in ("median_genes_per_cell", "median_transcripts_per_cell", "fraction_empty_cells"):
        assert got[k] == ref[k], k                      # identical whether or not compute_qc ran first
