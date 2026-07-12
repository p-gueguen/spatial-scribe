from __future__ import annotations


def test_load_thresholds_reads_yaml():
    from spatialscribe.analysis import config

    th = config.load_thresholds()
    assert isinstance(th, dict)
    # Layer-5 fail cutoff is published in the YAML as 0.25.
    assert config.get("layer5_confidence", "composite_confidence", "fail") == 0.25


def test_get_falls_back_when_missing():
    from spatialscribe.analysis import config

    # Non-existent section returns the caller's default, never raises.
    assert config.get("no_such_section", "x", default=42) == 42
    # A missing file path silently falls back to built-in defaults.
    th = config.load_thresholds("/does/not/exist.yaml")
    assert "layer2_counts" in th
