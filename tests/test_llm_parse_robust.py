"""Regression: `_parse_json` must degrade to None on a malformed/truncated model reply,
never raise. A raise here propagates out of `annotate_clusters` -> `consensus_annotate` ->
the `annotate` capability and crashes the whole step (the "progress bar stuck at 5/7,
annotate never marked done" bug), instead of falling back to marker-only labels.
"""
from spatialscribe.analysis.llm import _parse_json


def test_parse_json_plain_object():
    assert _parse_json('{"a": 1}') == {"a": 1}


def test_parse_json_extracts_from_prose_and_fences():
    # Model wraps JSON in prose / markdown fences: the object is still recovered.
    assert _parse_json('Sure:\n```json\n{"0": {"label": "T cell"}}\n```') == {"0": {"label": "T cell"}}


def test_parse_json_returns_none_on_truncated_reply():
    # A reply cut off at max_tokens mid-object: the first json.loads fails AND the regex-
    # extracted substring is itself invalid. Must return None, not raise.
    truncated = ('{"0": {"label": "T cell", "confidence": "high", "rationale": "CD3D+ CD2+"}, '
                 '"1": {"label": "B cell", "confid')
    assert _parse_json(truncated) is None


def test_parse_json_returns_none_on_garbage():
    assert _parse_json("no json here at all") is None
