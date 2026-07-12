"""Registry conformance - contracts every capability must satisfy, parametrized over the
whole registry so a newly added tool is covered by construction.

The six checks: (1) a complete contract declared with only registered keys; (2) a valid
Anthropic tool schema (copilot-exposed); (3) it runs on the golden fixture and actually
produces every key it declares; (4) it raises a *structured* PrerequisiteError when its
prerequisites are absent; (5) the copilot dispatch and a direct run() are the same path;
(6) its provenance record renders to compilable Python.
"""

from __future__ import annotations

import json

import pytest

from conftest import capability_params
from spatialscribe.analysis import capabilities as cap
from spatialscribe.analysis import export
from spatialscribe.analysis import keys as K

ALL = sorted(cap.REGISTRY)
COPILOT = sorted(cap.copilot_names())
WITH_PREREQS = sorted(n for n, c in cap.REGISTRY.items() if c.requires)


@pytest.mark.parametrize("name", ALL)
def test_contract_declared(name):
    c = cap.REGISTRY[name]
    assert c.name and c.label and c.description and callable(c.fn)
    for k in (*c.requires, *c.produces):
        assert K.is_declared(k), f"{name} references undeclared key {k}"
    assert set(c.required_params) <= set(c.params), f"{name} requires a param it does not declare"


@pytest.mark.parametrize("name", COPILOT)
def test_tool_schema_valid(name):
    t = cap.REGISTRY[name].to_tool_schema()
    assert set(t) == {"name", "description", "input_schema"}
    assert t["name"] == name and t["description"]
    schema = t["input_schema"]
    assert schema["type"] == "object"
    assert set(schema["required"]) <= set(schema["properties"])


@pytest.mark.parametrize("name", ALL)
def test_runs_and_produces(name, processed_adata, ctx):
    params = capability_params(name, processed_adata)
    res = cap.run(processed_adata, name, params, ctx)
    assert res.ok, (name, res.error)
    for k in cap.REGISTRY[name].produces:
        assert k.present(processed_adata), f"{name} declared it produces {k} but did not"
    assert res.record is not None and res.record["name"] == name
    json.dumps(res.value, default=str)  # value must be JSON-able for the copilot


@pytest.mark.parametrize("name", WITH_PREREQS)
def test_prereq_enforced(name, raw_adata, ctx):
    # raw_adata has no leiden/cell_type, so anything with prerequisites must fail structurally.
    res = cap.run(raw_adata, name, capability_params(name, raw_adata), ctx)
    assert not res.ok
    assert res.error["error_type"] == "prerequisite_missing"
    assert res.error["missing"], f"{name} reported no missing keys"
    assert res.error["hint"], f"{name} gave no remediation hint"


def test_frontend_parity(processed_adata, ctx):
    # Both frontends dispatch through cap.run; every copilot-exposed tool resolves to a
    # registry entry with a matching schema name, and running it directly succeeds.
    schema_names = {t["name"] for t in cap.copilot_tools()}
    assert schema_names == cap.copilot_names()
    for name in cap.copilot_names():
        assert name in cap.REGISTRY
        res = cap.run(processed_adata, name, capability_params(name, processed_adata), ctx)
        assert res.ok, (name, res.error)


def test_provenance_roundtrip(processed_adata, ctx, tmp_path):
    records = []
    for nm, p in [("compute_qc", {}), ("cluster", {"resolution": 1.0}), ("niches", {})]:
        r = cap.run(processed_adata, nm, p, ctx)
        assert r.ok, (nm, r.error)
        records.append(r.record)
    out = export.export_script(records, tmp_path / "analysis.py", tissue="melanoma")
    src = out.read_text()
    compile(src, str(out), "exec")               # malformed script -> SyntaxError
    assert "cap.run(adata," in src
    assert 'RunContext(tissue="melanoma")' in src
