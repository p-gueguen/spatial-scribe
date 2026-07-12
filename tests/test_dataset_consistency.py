"""Cross-dataset / cross-panel consistency workflow.

This suite guards the invariants that keep SpatialScribe's certainty / robustness advice
*self-consistent* across every spatial dataset it supports - Xenium 480-plex and Prime 5K,
CosMx / MERSCOPE, and Atera whole-transcriptome. It is the pytest face of the standalone
``scripts/check_consistency.py`` report (both share ``DATASET_PROFILES`` below).

The failure modes it exists to catch (all found in a July 2026 audit):
  * a panel-size-indexed count floor at Layer 2 but a *fixed* ``counts < 10`` gate at Layer 5
    (the exact non-transferable floor the docs warn against) - inconsistent certainty on 5K/WTA;
  * ``suggest_count_floor``'s ``platform`` arg and the whole ``cross_platform`` YAML profile
    being dead config (never read), so CosMx/MERSCOPE got Xenium thresholds;
  * ``platform`` / panel size not stamped on ``adata.uns``, so downstream QC was platform-blind.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")

SRC = Path(__file__).resolve().parents[1] / "src" / "spatialscribe"

# Mirror of docs/DATASETS.md - every real section the tool has been pointed at, plus the
# expected Layer-2 count-floor MODE. Small targeted panels get a fixed floor; large (>=1000
# gene) panels get a section-relative distributional floor. Add a dataset here -> it is guarded.
# (name, platform, n_panel_genes, expected_floor_mode)
DATASET_PROFILES = [
    ("mouse_brain (rctd-py)", "xenium", 250, "fixed"),
    ("renal 10x demo", "xenium", 405, "fixed"),
    ("breast tumour A", "xenium", 480, "fixed"),
    ("breast tumour merged", "xenium", 541, "fixed"),
    ("skin melanoma demo (Prime 5K)", "xenium", 5101, "distributional"),
    # The atera5k positive-control demo: platform 'atera' on a 5K-SIZED panel (the WTA row below only
    # exercises atera at 18k), so this pins that the floor mode keys off panel size, not platform.
    ("Atera 5K demo (positive control)", "atera", 4935, "distributional"),
    ("Atera breast (WTA)", "atera", 18000, "distributional"),
]


def _synthetic(n=300, n_x_genes=40, n_panel_genes=None, platform="xenium", counts_scale=0.5, seed=0):
    """Tiny AnnData whose *declared* panel size (uns) is decoupled from the # genes in X.

    Real 5K sections in a unit test would be huge; instead we stamp ``uns['n_panel_genes']``
    the way ``io.load`` now does, and assert the pipeline keys its behaviour off that.
    """
    rng = np.random.default_rng(seed)
    X = rng.poisson(counts_scale, size=(n, n_x_genes)).astype("float32")
    a = anndata.AnnData(X=X)
    a.var_names = [f"G{i}" for i in range(n_x_genes)]
    a.obs_names = [f"c{i}" for i in range(n)]
    a.var["control"] = False
    a.uns["platform"] = platform
    if n_panel_genes is not None:
        a.uns["n_panel_genes"] = int(n_panel_genes)
    return a


# --------------------------------------------------------------------------- #
# Invariant 1 - the count-floor MODE tracks the declared panel size, per dataset.
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "name,platform,n_panel,mode", DATASET_PROFILES, ids=[p[0] for p in DATASET_PROFILES]
)
def test_count_floor_mode_matches_dataset_profile(name, platform, n_panel, mode):
    from spatialscribe.analysis import qc

    a = _synthetic(platform=platform, n_panel_genes=n_panel)
    qc.compute_qc(a)
    out = qc.suggest_count_floor(a)  # must read uns['n_panel_genes'], NOT the 40 genes in X
    assert out["mode"] == mode, f"{name}: expected {mode} floor, got {out['mode']}"


# --------------------------------------------------------------------------- #
# Invariant 2 - the Layer-5 low-signal gate is the SAME panel-indexed floor as Layer 2,
# not a divergent hardcoded counts<10. On a large panel it must abstain FEWER cells as
# "low_quality" than a targeted panel would, holding the count distribution fixed.
# --------------------------------------------------------------------------- #
def _confidence_section(n_panel_genes):
    """Every cell: 12 counts over 12 genes. That clears BOTH targeted floors (counts >= 10, genes >= 5)
    yet falls below the RICH-panel gene floor (12 < qc.RICH_PANEL_MIN_GENES = 15). So the SAME cells are
    typable on a 300-gene targeted panel but low_quality on a 5K panel - the gate is panel-indexed."""
    genes = ["CD3D", "CD3E", "TRAC", "MLANA", "SOX10", "DCT"] + [f"N{i}" for i in range(60)]
    gi = {g: i for i, g in enumerate(genes)}
    n = 200
    X = np.zeros((n, len(genes)), dtype="float32")
    for c in range(n):
        block = ("CD3D", "CD3E", "TRAC") if c < n // 2 else ("MLANA", "SOX10", "DCT")
        for g in block:
            X[c, gi[g]] = 1.0          # 3 markers x 1 = 3 counts, 3 genes
        for k in range(9):
            X[c, gi[f"N{k}"]] = 1.0    # +9 noise genes -> total 12 counts, 12 genes
    a = anndata.AnnData(X=X)
    a.var_names = genes
    a.obs_names = [f"c{i}" for i in range(n)]
    a.var["control"] = False
    a.layers["counts"] = X.copy()      # apply_confidence reads the raw gate from here
    a.uns["n_panel_genes"] = int(n_panel_genes)
    a.obs["cell_type"] = pd.Categorical(["T cell"] * (n // 2) + ["Malignant/Melanocyte"] * (n // 2))
    return a


def test_confidence_low_quality_gate_is_panel_indexed():
    from spatialscribe.analysis import annotate

    lineage = {"T cell": ["CD3D", "CD3E", "TRAC"], "Malignant/Melanocyte": ["MLANA", "SOX10", "DCT"]}

    small = _confidence_section(n_panel_genes=300)
    annotate.apply_confidence(small, cluster_key="cell_type", marker_sets=lineage)
    small_lowq = float((small.obs["annotation_reason"] == "low_quality").mean())

    large = _confidence_section(n_panel_genes=5000)
    annotate.apply_confidence(large, cluster_key="cell_type", marker_sets=lineage)
    large_lowq = float((large.obs["annotation_reason"] == "low_quality").mean())

    # 12-gene cells clear the targeted floors (counts >= 10, genes >= 5) -> typable on a 300-gene panel.
    assert small_lowq < 0.1, f"targeted 12-gene cells should not be low_quality, got {small_lowq:.2f}"
    # On a 5K panel the same cells fall below the rich-panel gene floor (12 < RICH_PANEL_MIN_GENES) and
    # abstain: the gate is panel-indexed (a rich panel needs more detected genes to identify a cell).
    # This is the gene arm that kills the breast-5K near-empty artefact clusters; the count-floor arm is
    # covered by test_qc_layers / test_count_floor_shallow.
    assert large_lowq > 0.9, f"5K near-empty (12<15 genes) cells should gate low_quality, got {large_lowq:.2f}"
    assert large_lowq > small_lowq, (
        f"Layer-5 gate is not panel-indexed: 5K abstained {large_lowq:.2f} vs targeted {small_lowq:.2f}"
    )


# --------------------------------------------------------------------------- #
# Invariant 3 - io.load stamps platform + panel size on adata.uns (the root fix that lets
# every downstream, adata-only function be platform/panel aware).
# --------------------------------------------------------------------------- #
def test_io_load_stamps_platform_and_panel_size(tmp_path):
    from spatialscribe.analysis import io

    a = anndata.AnnData(X=np.random.default_rng(0).poisson(1, size=(20, 12)).astype("float32"))
    a.var_names = [f"G{i}" for i in range(12)]
    a.obs_names = [f"c{i}" for i in range(20)]
    p = tmp_path / "sec.h5ad"
    a.write_h5ad(p)

    s = io.load(p)
    assert s.adata.uns.get("platform") == "h5ad"
    assert s.adata.uns.get("n_panel_genes") == len(s.panel_genes)


# --------------------------------------------------------------------------- #
# Invariant 4 - section-QC flagging is platform-aware (wires the cross_platform profile):
# CosMx gets a looser negative-control margin; MERSCOPE FFPE is never hard-failed on
# retention (it legitimately loses most cells).
# --------------------------------------------------------------------------- #
def test_section_qc_flag_is_platform_aware():
    from spatialscribe.analysis import qc

    metrics = {
        "pct_counts_control": 6.0,         # > Xenium fail (5.0)
        "fraction_empty_cells": 0.30,      # > Xenium fail (0.25)
        "median_genes_per_cell": 30,
        "median_transcripts_per_cell": 60,
    }
    xen = qc._flag(metrics, platform="xenium")
    cos = qc._flag(metrics, platform="cosmx")
    mer = qc._flag(metrics, platform="merscope")

    assert xen["pct_counts_control"] == "error"      # Xenium: 6 > 5
    assert cos["pct_counts_control"] != "error"       # CosMx: looser neg-control margin
    assert xen["fraction_empty_cells"] == "error"     # Xenium: 0.30 > 0.25
    assert mer["fraction_empty_cells"] != "error"     # MERSCOPE retention gate: warn at most


# --------------------------------------------------------------------------- #
# Invariant 5 - no dead threshold config: every top-level profile in the YAML is consumed
# by code (read via config.get, or mirrored by a named construct). This is the general
# drift-catcher: add/rename a YAML profile without wiring it -> this fails.
# --------------------------------------------------------------------------- #
# section -> a token whose presence in src/ proves the profile is actually consumed.
_WIRED_TOKENS = {
    "layer0_section": "XENIUM_THRESHOLDS",
    "layer1_segmentation": "def segmentation_qc",
    "layer2_counts": 'config.get("layer2_counts"',
    "layer3_contamination": "def crisp_purity",
    "layer4_panel": "def check_panel",
    "layer5_confidence": 'config.get("layer5_confidence"',
    "layer6_spatial": "def spatial_coherence",
    "consensus_popv": 'config.get("consensus_popv"',
    "cluster_confidence": 'config.get("cluster_confidence"',
    "cross_platform": 'config.get("cross_platform"',
    "abstention_labels": "_ABSTAIN",
    "annotatability_summary": "pct_pass",
    "verify": 'config.get("verify"',
}
_YAML_SECTIONS_EXEMPT = {"meta"}  # pure metadata, not a threshold profile


def test_no_dead_threshold_config():
    import yaml

    from spatialscribe.analysis import config

    with open(config._DEFAULT_YAML) as fh:
        y = yaml.safe_load(fh)
    src_text = "\n".join(p.read_text() for p in SRC.rglob("*.py"))

    for section in y:
        if section in _YAML_SECTIONS_EXEMPT:
            continue
        assert section in _WIRED_TOKENS, (
            f"YAML profile '{section}' has no wiring token - add it to _WIRED_TOKENS and wire "
            f"it into code, or exempt it. Unwired config silently does nothing."
        )
        token = _WIRED_TOKENS[section]
        assert token in src_text, (
            f"YAML profile '{section}' is dead config: its wiring token {token!r} appears "
            f"nowhere in src/. The threshold is documented but never applied."
        )


# --------------------------------------------------------------------------- #
# Invariant 6 - every platform io.load can emit is handled, and the imaging platforms all
# have a cross_platform negative-control scaling entry.
# --------------------------------------------------------------------------- #
def test_supported_platforms_all_covered():
    from spatialscribe.analysis import config

    scaling = config.get("cross_platform", "neg_control_threshold_scaling", default={}) or {}
    assert {"xenium", "cosmx", "merscope"} <= set(scaling), (
        f"cross_platform neg-control scaling is missing an imaging platform: {sorted(scaling)}"
    )
