# tests/test_load_section.py
"""The copilot can LOAD a section from a path: detect the platform, infer the tissue from the panel,
auto-select the reference, and swap the new section in. Covers the capability contract (what it
stashes on ctx.loaded) and the run_copilot rebind (a load mid-conversation makes the new section the
one the rest of the turn - and later turns - run on). Pure: no LLM, no reference files needed."""
from __future__ import annotations

import numpy as np
import pytest

anndata = pytest.importorskip("anndata")


def _write_section(path, *, organism="Mouse", tissue="Brain", n=60, genes=8):
    """A tiny .h5ad section io.load can open: real genes + obsm['spatial'] + a Xenium panel stamp."""
    rng = np.random.default_rng(0)
    a = anndata.AnnData(X=rng.poisson(0.5, size=(n, genes)).astype("float32"))
    a.var_names = [f"Gene{i}" for i in range(genes)]
    a.obs_names = [f"c{i}" for i in range(n)]
    a.obsm["spatial"] = rng.random((n, 2)) * 100.0
    a.uns["panel"] = {"name": f"Xenium {organism} {tissue}", "organism": organism,
                      "tissue_type": tissue, "n_targets": genes}
    a.uns["panel_name"] = a.uns["panel"]["name"]
    a.write_h5ad(path)
    return str(path)


def test_load_section_capability_infers_tissue_and_stashes_new_section(tmp_path):
    from spatialscribe.analysis import capabilities as cap

    p = _write_section(tmp_path / "mouse_brain.h5ad")
    dummy = anndata.AnnData(X=np.zeros((3, 2), dtype="float32"))   # the "current" section (unrelated)
    ctx = cap.RunContext(tissue="melanoma")

    out = cap.run(dummy, "load_section", {"path": p, "auto_reference": False}, ctx)

    assert out.ok, out.error
    v = out.value
    assert v["status"] == "loaded" and v["n_cells"] == 60
    assert v["tissue"] == "mouse brain"                 # inferred from the panel, NOT the melanoma default
    # The new section is stashed for the app/copilot to swap in - a capability cannot replace adata.
    assert len(ctx.loaded) == 1
    entry = ctx.loaded[-1]
    assert entry["adata"] is not dummy and entry["adata"].n_obs == 60
    assert entry["tissue"] == "mouse brain" and entry["path"] == p


def test_load_section_is_copilot_exposed_and_requires_path():
    from spatialscribe.analysis import capabilities as cap

    assert "load_section" in cap.copilot_names()
    schema = next(t for t in cap.copilot_tools() if t["name"] == "load_section")
    assert "path" in schema["input_schema"]["required"]
    # A missing/empty path is a clean structured error, not a crash.
    err = cap.run(anndata.AnnData(X=np.zeros((2, 2), dtype="float32")), "load_section", {}, cap.RunContext())
    assert not err.ok and "path" in str(err.error).lower()


def test_load_section_auto_reference_degrades_honestly_with_no_registry(tmp_path, monkeypatch):
    # With no SPATIALSCRIBE_REF_* wired, auto-select must NOT invent a reference: it reports
    # no_reference and the strategy routes to marker-based annotation (honest degradation).
    from spatialscribe.analysis import capabilities as cap
    for k in list(__import__("os").environ):
        if k.startswith("SPATIALSCRIBE_REF"):
            monkeypatch.delenv(k, raising=False)

    p = _write_section(tmp_path / "sec.h5ad")
    ctx = cap.RunContext(tissue="melanoma")
    out = cap.run(anndata.AnnData(X=np.zeros((2, 2), dtype="float32")),
                  "load_section", {"path": p, "auto_reference": True}, ctx)
    assert out.ok, out.error
    ref = out.value["reference"]
    assert ref["status"] == "no_reference"                       # nothing fabricated
    assert out.value["recommended_mode"] in ("annotate", "cluster")
    assert ctx.loaded[-1]["reference"] is None                   # no reference registered on the section


def test_load_section_accepts_section_path_alias(tmp_path):
    from spatialscribe.analysis import capabilities as cap

    p = _write_section(tmp_path / "s.h5ad", n=20)
    out = cap.run(anndata.AnnData(X=np.zeros((2, 2), dtype="float32")),
                  "load_section", {"section_path": p, "auto_reference": False}, cap.RunContext())
    assert out.ok, out.error
    assert out.value["n_cells"] == 20


def test_user_reported_markup_flow_loads_via_recovery_and_alias(tmp_path):
    # The exact live failure: the vLLM printed the call as a ```json block that also named the path
    # 'section_path'. Recovery parses the markup; the section_path alias makes load_section load. Both
    # fixes together turn the user's "didn't work" (raw markup echoed) into a real load.
    from spatialscribe.analysis import capabilities as cap
    from spatialscribe.analysis import llm

    p = _write_section(tmp_path / "brain.h5ad", n=33)
    markup = ("I'll load the mouse brain data from the specified path. ```json\n"
              '{ "tool_name": "load_section", "arguments": { "section_path": "%s" } }\n```' % p)
    recovered = llm._recover_tool_calls(markup, {"load_section"})
    assert recovered and recovered[0]["name"] == "load_section"
    args = {**recovered[0]["input"], "auto_reference": False}
    out = cap.run(anndata.AnnData(X=np.zeros((2, 2), dtype="float32")), "load_section", args,
                  cap.RunContext())
    assert out.ok, out.error
    assert out.value["n_cells"] == 33 and out.value["tissue"] == "mouse brain"


def test_run_copilot_rebinds_to_the_loaded_section(tmp_path, monkeypatch):
    # A load mid-conversation swaps the active section: the copilot loop rebinds `adata` so a later
    # tool call in the SAME turn (here describe_sample) runs on the NEWLY loaded section, and the
    # loaded list the caller passed in is populated (that is how the app commits the swap).
    from spatialscribe.agent import tools as agent_tools
    from spatialscribe.analysis import llm

    p = _write_section(tmp_path / "brain.h5ad", n=42)
    seen_cells = {}

    calls = iter([
        # turn 1: ask to load, then describe -> two tool calls in one turn
        {"text": "", "assistant": {"role": "assistant", "content": "ld"},
         "tool_calls": [{"id": "t1", "name": "load_section", "input": {"path": p, "auto_reference": False}},
                        {"id": "t2", "name": "describe_sample", "input": {}}]},
        # turn 2: no tool calls -> final answer
        {"text": "Loaded the mouse brain section (42 cells).", "assistant": {"role": "assistant", "content": "done"},
         "tool_calls": []},
    ])

    def fake_tool_chat(system, messages, tools, max_tokens=1024):
        turn = next(calls)
        # capture what describe_sample saw, to prove the rebind happened
        for m in messages:
            if isinstance(m.get("content"), list):
                for c in m["content"]:
                    if c.get("tool_use_id") == "t2":
                        seen_cells["describe"] = c["content"]
        return turn

    monkeypatch.setattr(llm, "tool_chat", fake_tool_chat)

    old = anndata.AnnData(X=np.zeros((5, 2), dtype="float32"))
    loaded: list = []
    reply = agent_tools.run_copilot(old, "load the brain data and describe it",
                                    use_llm=True, loaded=loaded)
    assert "42" in reply
    assert len(loaded) == 1 and loaded[-1]["adata"].n_obs == 42
    # describe_sample ran AFTER the rebind -> it reported 42 cells (the new section), not 5 (the old).
    assert "42" in seen_cells.get("describe", "")


def test_resolve_server_path_corrects_missing_mount_prefix(tmp_path, monkeypatch):
    """A path that dropped its the cluster mount prefix (a user pasting '/projects/pXXXXX/...' instead of
    '/data/projects/...', or an LLM emitting it) still resolves: _resolve_server_path tries the
    known storage roots and returns the first that exists. Nonexistent -> None (load then errors clearly)."""
    from spatialscribe.analysis import capabilities as cap
    # A fake object storage root: <tmp>/data/projects/internal/section
    real = tmp_path / "srv" / "object storage" / "projects" / "internal" / "section"
    real.mkdir(parents=True)
    monkeypatch.setattr(cap, "_SERVER_ROOTS", (str(tmp_path / "srv" / "object storage"),))

    assert cap._resolve_server_path(str(real)) == str(real)          # already-correct path: returned as-is
    assert cap._resolve_server_path("/projects/internal/section") == str(real)  # prefix dropped -> corrected
    assert cap._resolve_server_path("/projects/internal/nope") is None          # nothing matches -> None
    assert cap._resolve_server_path("") is None
