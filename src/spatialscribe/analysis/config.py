"""Tunable QC/annotation thresholds, loaded from annotation_qc_thresholds.yaml.

What it does
------------
Reads the committed draft-threshold YAML so ``qc.py`` / ``annotate.py`` / ``spatial.py``
read pass/warn/fail numbers from one place instead of hardcoding them (the YAML itself
says "ship as tunable config"). Falls back to built-in defaults when the file or PyYAML
is unavailable, so the demo stays offline-deterministic.

How to use it
-------------
>>> from .config import get
>>> get("layer5_confidence", "composite_confidence", "fail")   # -> 0.25

Depends on: PyYAML (optional; falls back if absent).
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any

# docs/research/annotation_qc_thresholds.yaml relative to this file:
#   .../src/spatialscribe/analysis/config.py -> parents[3] == project root.
# NOTE: the docs/ reorg (commit e3db95e) moved this file into docs/research/; the
# path below points at its current location, not the pre-reorg docs/ top level.
_DEFAULT_YAML = (
    Path(__file__).resolve().parents[3] / "docs" / "research" / "annotation_qc_thresholds.yaml"
)

# Built-in fallback (mirrors the doc) used only when the YAML can't be read.
_FALLBACK: dict = {
    "layer1_segmentation": {"cell_area_um2": {"warn_lt": 8, "warn_gt": 400}},
    "layer2_counts": {
        "total_counts_panel_schedule": {"targeted_lt1000_genes": {"warn": 25, "fail": 10}},
        "large_panel_percentile": 0.02,
    },
    "layer5_confidence": {"composite_confidence": {"warn": 0.5, "fail": 0.25}},
    "layer6_spatial": {"spatial_coherence_frac_same_label_neighbors": {"warn": 0.2, "k_neighbors": 15}},
}


@lru_cache(maxsize=8)
def load_thresholds(path: str | None = None) -> dict:
    """Parse the thresholds YAML into a dict; fall back to built-in defaults on any error."""
    p = Path(path) if path else _DEFAULT_YAML
    try:
        import yaml

        with open(p) as fh:
            data = yaml.safe_load(fh)
        return data if isinstance(data, dict) else dict(_FALLBACK)
    except Exception:
        return dict(_FALLBACK)


def get(section: str, *keys: str, default: Any = None, path: str | None = None) -> Any:
    """Nested lookup ``thresholds[section][key1][key2]...``; return ``default`` if any hop misses."""
    node: Any = load_thresholds(path).get(section, {})
    for k in keys:
        if not isinstance(node, dict) or k not in node:
            return default
        node = node[k]
    return node
