"""Cell-annotation quality metrics - thin re-export shim.

The metric implementations now live in the standalone **spatial-anno-metrics** package
(github.com/p-gueguen/spatial-anno-metrics), which is the single source of truth. This module
keeps the ``spatialscribe.analysis.eval_metrics`` import path working (used by ``qc.run_funnel``
and the tests) so nothing else in the app had to change.

Install: ``uv pip install -e ../../../spatial-anno-metrics --no-deps`` (declared as an editable
path dep in ``pixi.toml``). See ``docs/research/cell-annotation-quality-metrics.md`` for the catalog.
"""
from spatial_anno_metrics.eval_metrics import (  # noqa: F401
    annotation_quality,
    composition_accuracy,
    conformal_prediction_sets,
    deconvolution_metrics,
    element_centric_similarity,
    external_scores,
    hierarchical_accuracy,
    inter_sample_consistency,
    internal_validity,
    marker_program_fidelity,
    panel_resolvability,
)

__all__ = [
    "internal_validity", "inter_sample_consistency", "marker_program_fidelity",
    "external_scores", "element_centric_similarity", "hierarchical_accuracy",
    "composition_accuracy", "deconvolution_metrics", "panel_resolvability",
    "annotation_quality", "conformal_prediction_sets",
]
