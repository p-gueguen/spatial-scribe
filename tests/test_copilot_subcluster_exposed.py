"""subcluster must be reachable by the copilot (offline guard for the exposure).

Regression: subcluster is a real capability but was NOT copilot_exposed, so a chat request to
"subcluster the T cells" had no valid tool target and both LLMs routed it to self_heal /
cluster_confidence instead. This pins the exposure (the LIVE routing itself is checked by
an internal LLM smoke-check against a real model - a mocked test can't verify routing).
"""
from __future__ import annotations

from spatialscribe.analysis import capabilities as cap


def test_subcluster_is_copilot_exposed():
    tools = {t["name"]: t for t in cap.copilot_tools()}
    assert "subcluster" in tools, "subcluster missing from the copilot toolset"
    schema = tools["subcluster"]["input_schema"]
    # the tool is useless to the copilot without a required, described cell_type argument
    assert "cell_type" in schema["required"]
    assert schema["properties"]["cell_type"].get("description"), "cell_type needs a description to route the arg"


def test_subcluster_description_points_at_self_heal():
    """The disambiguator that keeps whole-annotation 'auto-fix' traffic on self_heal, not subcluster."""
    d = {t["name"]: t for t in cap.copilot_tools()}["subcluster"]["description"].lower()
    assert "self_heal" in d and "subcluster" in d
