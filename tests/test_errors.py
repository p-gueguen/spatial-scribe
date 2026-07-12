"""Unit coverage for the structured error contract (errors.py).

Pure: no AnnData, no compute. Asserts the JSON-able dicts the two frontends consume.
"""

from __future__ import annotations

from spatialscribe.analysis.errors import (
    PrerequisiteError,
    CapabilityError,
    to_error_dict,
)
from spatialscribe.analysis.keys import Obs


def test_prerequisite_error_to_dict():
    exc = PrerequisiteError("annotate", [Obs.LEIDEN], "run 'cluster' first")
    d = exc.to_dict()
    assert d["error_type"] == "prerequisite_missing"
    assert d["capability"] == "annotate"
    assert d["missing"] == ["obs:leiden"]
    assert d["hint"] == "run 'cluster' first"
    # the human-readable message names the missing key.
    assert "obs:leiden" in str(exc)


def test_capability_error_to_dict():
    d = CapabilityError("x", "boom").to_dict()
    assert d["error_type"] == "capability_failed"
    assert d["capability"] == "x"
    assert d["message"] == "boom"


def test_to_error_dict_generic_exception():
    d = to_error_dict(ValueError("nope"), "cluster")
    assert d == {
        "error_type": "capability_failed",
        "capability": "cluster",
        "message": "nope",
    }


def test_to_error_dict_preserves_structured_prerequisite():
    exc = PrerequisiteError("a", [Obs.CELL_TYPE], "h")
    d = to_error_dict(exc)
    assert d["error_type"] == "prerequisite_missing"
    assert d["capability"] == "a"
    assert d["missing"] == ["obs:cell_type"]
    assert d["hint"] == "h"
