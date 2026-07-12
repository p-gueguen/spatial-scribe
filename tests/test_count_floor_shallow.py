"""The distributional count floor must still drop empty cells on a SHALLOW section (regression).

A shallow 5K section can have >= pct% of cells at 0 counts, so the pct-th percentile is 0. A raw floor
of 0 with a strict `counts < floor` removes NOTHING - not even the 0-count segments - so the QC funnel
showed "count floor 0, 0.0% removed". The floor is clamped to >= 1 so the empties are dropped, without
reaching the aggressive fixed <10 floor that over-removes on a 5K panel.
"""
from __future__ import annotations

import numpy as np

from spatialscribe.analysis import qc


def test_shallow_5k_floor_drops_empties_not_zero():
    # Prime-5K-size panel, shallow: 5% zeros + a low-count body -> pct=0.02 raw quantile is 0.
    counts = np.concatenate([np.zeros(500),
                             np.random.default_rng(0).integers(1, 80, 9500)]).astype(float)
    floor, mode = qc.panel_indexed_floor(counts, n_panel_genes=5101)
    assert mode == "distributional"
    assert floor >= 1.0                     # not 0 -> the layer actually removes something
    assert int((counts < floor).sum()) == 500   # exactly the empty (0-count) cells


def test_deep_5k_floor_uses_percentile_unchanged():
    counts = np.random.default_rng(0).integers(50, 2000, 10000).astype(float)
    floor, _ = qc.panel_indexed_floor(counts, n_panel_genes=5101)
    assert floor == float(np.quantile(counts, 0.02))   # deep: the clamp is a no-op, the quantile stands


def test_targeted_panel_uses_fixed_floor():
    counts = np.random.default_rng(0).integers(0, 100, 1000).astype(float)
    _, mode = qc.panel_indexed_floor(counts, n_panel_genes=480)
    assert mode == "fixed"      # < 1000 genes -> a fixed floor, not the distributional percentile
