"""Contamination / purity metrics (QC Layer 3) - thin re-export shim.

The implementations now live in the standalone **spatial-anno-metrics** package
(github.com/p-gueguen/spatial-anno-metrics, ``spatial_anno_metrics.purity``), the single source of
truth for annotation-quality scores. This module keeps the ``spatialscribe.analysis.purity`` import
path + signatures working (``qc.run_funnel``, ``annotate.apply_confidence``, ``cnv``, the app, tests)
and supplies the SpatialScribe default marker sets, so nothing else in the app had to change.

Reference-free (panel markers only): ``crisp_purity``, ``mecr``, ``pmp``.
Reference-based (need a matched reference): ``nmp``, ``ncp``.

Install: editable path dep in ``pixi.toml`` (``spatial-anno-metrics = { path = "../../../spatial-anno-metrics", editable = true }``).
See ``docs/research/cell-annotation-quality-metrics.md`` §3c/§3f for the catalog.
"""
from __future__ import annotations

from spatial_anno_metrics import purity as _p

from . import markers as _m

__all__ = ["crisp_purity", "mecr", "pmp", "nmp", "ncp"]


def crisp_purity(adata, lineage_markers: dict[str, list[str]] | None = None) -> float:
    """CRISP purity (reference-free) - writes ``obs['crisp_impure']``, returns dataset purity."""
    return _p.crisp_purity(adata, lineage_markers or _m.LINEAGE_MARKERS)


def mecr(adata, lineage_markers: dict[str, list[str]] | None = None) -> float:
    """MECR (reference-free) - mean mutually-exclusive co-detection over disjoint lineage pairs."""
    return _p.mecr(adata, lineage_markers or _m.LINEAGE_MARKERS)


def pmp(adata, assigned_label_key: str = "cell_type", lineage_markers: dict | None = None) -> None:
    """Per-cell marker purity (reference-free, panel-size invariant) - writes ``obs['pmp']``."""
    return _p.pmp(adata, lineage_markers or _m.LINEAGE_MARKERS, label_key=assigned_label_key)


def nmp(adata, reference=None, assigned_label_key: str = "cell_type",
        ref_label_key: str = "cell_type") -> dict:
    """Negative-marker proportion (reference-based, guarded) - writes ``obs['nmp']``."""
    return _p.nmp(adata, reference=reference, label_key=assigned_label_key, ref_label_key=ref_label_key)


def ncp(adata, reference=None, ref_label_key: str = "cell_type", coexpr_threshold: float = 0.1,
        max_genes: int = 400, max_pairs: int = 2000) -> dict:
    """Non-coexpression preservation (reference-based, guarded, MECR-style reference pairs)."""
    return _p.ncp(adata, reference=reference, ref_label_key=ref_label_key,
                  coexpr_threshold=coexpr_threshold, max_genes=max_genes, max_pairs=max_pairs)
