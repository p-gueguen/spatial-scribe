"""Unit coverage for the data-contract (keys.py) and derived pipeline state (status.py).

Pure-ish tests: they lean on the shared ``raw_adata`` / ``processed_adata`` fixtures from
conftest for real ``.present`` checks, but exercise no heavy compute of their own.
"""

from __future__ import annotations

import pytest

from spatialscribe.analysis import keys
from spatialscribe.analysis import status
from spatialscribe.analysis import capabilities as cap
from spatialscribe.analysis.keys import Key, Obs, Uns


# --------------------------------------------------------------------------- #
# Key.present + __str__ + space validation
# --------------------------------------------------------------------------- #
def test_key_present_absent_on_raw(raw_adata):
    # cell_type / panel_check are absent on the raw section (annotate/panel_check
    # have not run). Requested alone so the processed_adata fixture (which mutates
    # raw_adata in place) cannot contaminate this object.
    assert Obs.CELL_TYPE.present(raw_adata) is False
    assert Uns.PANEL_CHECK.present(raw_adata) is False


def test_key_present_true_on_processed(processed_adata):
    # cell_type set by annotate; panel_check writes its uns entry.
    assert Obs.CELL_TYPE.present(processed_adata) is True
    assert Uns.PANEL_CHECK.present(processed_adata) is True


def test_key_str_is_space_colon_name():
    assert str(Obs.CELL_TYPE) == "obs:cell_type"
    assert str(Obs.LEIDEN) == "obs:leiden"


def test_key_bad_space_raises_value_error():
    with pytest.raises(ValueError):
        Key("not_a_space", "whatever")


# --------------------------------------------------------------------------- #
# all_keys / is_declared
# --------------------------------------------------------------------------- #
def test_all_keys_are_keys_and_unique():
    all_keys = keys.all_keys()
    assert all_keys, "expected at least one declared key"
    assert all(isinstance(k, Key) for k in all_keys)
    pairs = [(k.space, k.name) for k in all_keys]
    assert len(pairs) == len(set(pairs)), "duplicate (space, name) declared in keys.py"


def test_is_declared():
    assert keys.is_declared(Obs.CELL_TYPE) is True
    assert keys.is_declared(Obs.LEIDEN) is True
    assert keys.is_declared(Key("obs", "not_a_real_key")) is False


# --------------------------------------------------------------------------- #
# status.pipeline_status
# --------------------------------------------------------------------------- #
def test_pipeline_status_processed(processed_adata):
    st = status.pipeline_status(processed_adata)
    assert st["cluster"] is True
    assert st["annotate"] is True


def test_pipeline_status_raw(raw_adata):
    st = status.pipeline_status(raw_adata)
    assert st["cluster"] is False
    assert st["annotate"] is False


# --------------------------------------------------------------------------- #
# status.missing_prereqs / can_run / available
# --------------------------------------------------------------------------- #
def test_missing_prereqs_annotate_on_raw(raw_adata):
    missing = status.missing_prereqs(raw_adata, "annotate")
    assert [str(k) for k in missing] == [str(Obs.LEIDEN)]


def test_can_run(raw_adata):
    assert status.can_run(raw_adata, "annotate") is False
    assert status.can_run(raw_adata, "cluster") is True


def test_available(raw_adata):
    avail = status.available(raw_adata)
    for runnable in ("cluster", "compute_qc", "panel_check"):
        assert runnable in avail
    for blocked in ("annotate", "niches"):
        assert blocked not in avail


# --------------------------------------------------------------------------- #
# capabilities.producer_of (the prereq-hint backbone status.py relies on)
# --------------------------------------------------------------------------- #
def test_producer_of():
    assert cap.producer_of(Obs.CELL_TYPE) == "annotate"
    assert cap.producer_of(Obs.LEIDEN) == "cluster"
