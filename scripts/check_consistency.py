#!/usr/bin/env python
"""Cross-dataset consistency check - the runnable face of tests/test_dataset_consistency.py.

Verifies that SpatialScribe treats every supported spatial dataset (Xenium 480 / Prime 5K,
CosMx, MERSCOPE, Atera WTA) *consistently* when it advises on cell-typing certainty and
robustness. Catches the class of bug where a panel-size- or platform-dependent threshold is
wired in one layer but hardcoded (or dead config) in another.

Run it:
    pixi run -e main python scripts/check_consistency.py     # or: pixi run check

Exit code 0 = all invariants hold; 1 = at least one inconsistency (prints which). Designed for
CI / pre-commit / a periodic cron - no pytest, no network, no GPU, no API key.

Checks
------
1. count-floor MODE tracks the declared panel size, per dataset profile (targeted -> fixed,
   >=1000 genes -> distributional).
2. the Layer-5 confidence low-signal gate uses the SAME panel-indexed floor as Layer 2
   (a 5K panel abstains fewer cells as "low_quality" than a targeted panel, same counts).
3. section-QC flagging is platform-aware (CosMx looser neg-control; MERSCOPE not hard-failed
   on retention) - i.e. the cross_platform YAML profile is actually consumed.
4. no dead threshold config: every YAML profile is wired into code.
5. every supported platform has a cross_platform negative-control scaling entry.
"""

from __future__ import annotations

import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "src"))
sys.path.insert(0, str(_REPO / "tests"))

# Single source of truth for the dataset roster + wiring tokens (shared with the pytest suite).
from test_dataset_consistency import (  # noqa: E402
    DATASET_PROFILES,
    _WIRED_TOKENS,
    _YAML_SECTIONS_EXEMPT,
    _confidence_section,
    _synthetic,
)

import anndata  # noqa: E402,F401  (ensures the scientific stack is importable before checks run)


class Result:
    def __init__(self, name: str):
        self.name = name
        self.failures: list[str] = []

    def check(self, ok: bool, detail: str) -> None:
        if not ok:
            self.failures.append(detail)

    @property
    def ok(self) -> bool:
        return not self.failures


def check_count_floor_modes() -> Result:
    from spatialscribe.analysis import qc

    r = Result("count-floor mode per dataset profile")
    for name, platform, n_panel, mode in DATASET_PROFILES:
        a = _synthetic(platform=platform, n_panel_genes=n_panel)
        qc.compute_qc(a)
        out = qc.suggest_count_floor(a)
        r.check(out["mode"] == mode,
                f"{name} ({n_panel} genes, {platform}): expected {mode}, got {out['mode']}")
    return r


def check_confidence_gate_panel_indexed() -> Result:
    from spatialscribe.analysis import annotate

    r = Result("Layer-5 confidence gate is panel-indexed")
    lineage = {"T cell": ["CD3D", "CD3E", "TRAC"], "Malignant/Melanocyte": ["MLANA", "SOX10", "DCT"]}

    small = _confidence_section(n_panel_genes=300)
    annotate.apply_confidence(small, cluster_key="cell_type", marker_sets=lineage)
    small_lowq = float((small.obs["annotation_reason"] == "low_quality").mean())

    large = _confidence_section(n_panel_genes=5000)
    annotate.apply_confidence(large, cluster_key="cell_type", marker_sets=lineage)
    large_lowq = float((large.obs["annotation_reason"] == "low_quality").mean())

    r.check(small_lowq > 0.9, f"targeted panel should gate ~all at fixed floor 10, got {small_lowq:.2f}")
    r.check(large_lowq < small_lowq,
            f"5K panel not distinguished: 5K low_quality={large_lowq:.2f} vs targeted={small_lowq:.2f}")
    return r


def check_platform_flagging() -> Result:
    from spatialscribe.analysis import qc

    r = Result("section-QC flagging is platform-aware (cross_platform wired)")
    metrics = {"pct_counts_control": 6.0, "fraction_empty_cells": 0.30,
               "median_genes_per_cell": 30, "median_transcripts_per_cell": 60}
    xen = qc._flag(metrics, platform="xenium")
    cos = qc._flag(metrics, platform="cosmx")
    mer = qc._flag(metrics, platform="merscope")
    r.check(xen["pct_counts_control"] == "error", "Xenium neg-control 6% should be error")
    r.check(cos["pct_counts_control"] != "error", "CosMx neg-control should be looser than Xenium")
    r.check(xen["fraction_empty_cells"] == "error", "Xenium 30% empty should be error")
    r.check(mer["fraction_empty_cells"] != "error", "MERSCOPE retention should not hard-fail")
    return r


def check_no_dead_config() -> Result:
    import yaml

    from spatialscribe.analysis import config

    r = Result("no dead threshold config")
    with open(config._DEFAULT_YAML) as fh:
        y = yaml.safe_load(fh)
    src_text = "\n".join(p.read_text() for p in (_REPO / "src" / "spatialscribe").rglob("*.py"))
    for section in y:
        if section in _YAML_SECTIONS_EXEMPT:
            continue
        if section not in _WIRED_TOKENS:
            r.check(False, f"YAML profile '{section}' has no wiring token (add to _WIRED_TOKENS)")
            continue
        token = _WIRED_TOKENS[section]
        r.check(token in src_text, f"YAML profile '{section}' is dead config (token {token!r} not in src/)")
    return r


def check_platform_coverage() -> Result:
    from spatialscribe.analysis import config

    r = Result("imaging platforms covered by cross_platform scaling")
    scaling = config.get("cross_platform", "neg_control_threshold_scaling", default={}) or {}
    for p in ("xenium", "cosmx", "merscope"):
        r.check(p in scaling, f"platform '{p}' missing from neg_control_threshold_scaling")
    return r


CHECKS = [
    check_count_floor_modes,
    check_confidence_gate_panel_indexed,
    check_platform_flagging,
    check_no_dead_config,
    check_platform_coverage,
]


def main() -> int:
    print("SpatialScribe cross-dataset consistency check")
    print("=" * 60)
    print("\nDataset roster (docs/DATASETS.md):")
    for name, platform, n_panel, mode in DATASET_PROFILES:
        print(f"  - {name:<32} {platform:<8} {n_panel:>6} genes -> {mode} floor")
    print()

    results = [chk() for chk in CHECKS]
    n_fail = 0
    for res in results:
        status = "PASS" if res.ok else "FAIL"
        print(f"[{status}] {res.name}")
        for f in res.failures:
            n_fail += 1
            print(f"         - {f}")
    print("=" * 60)
    if n_fail:
        print(f"{n_fail} inconsistency(ies) found across {len(results)} checks.")
        return 1
    print(f"All {len(results)} consistency checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
