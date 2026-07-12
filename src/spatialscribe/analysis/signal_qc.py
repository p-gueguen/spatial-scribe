"""Annotation-independent signal / QC metrics - thin re-export shim.

The implementations now live in the standalone **spatial-anno-metrics** package
(github.com/p-gueguen/spatial-anno-metrics, the single source of truth). This module keeps the
``spatialscribe.analysis.signal_qc`` import path working. See ``pixi.toml`` for the editable dep
and ``docs/research/cell-annotation-quality-metrics.md`` for the catalog.
"""
from spatial_anno_metrics.signal_qc import (  # noqa: F401
    detection_entropy,
    moran_signal,
    run_signal_qc,
    signal_to_noise,
    sparsity,
    tx_per_area,
)

__all__ = [
    "moran_signal", "signal_to_noise", "sparsity", "detection_entropy",
    "tx_per_area", "run_signal_qc",
]
