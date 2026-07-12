"""Endpoint-agnostic LLM backend selection + OpenAI translation (offline, no network)."""
from __future__ import annotations

import pytest

from spatialscribe.analysis import llm


def test_provider_defaults_to_anthropic(monkeypatch):
    monkeypatch.delenv("SPATIALSCRIBE_LLM_BASE_URL", raising=False)
    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-test")
    assert llm.provider() == "anthropic"
    assert llm.available() is True
    assert llm.default_model().startswith("claude")


def test_openai_base_url_takes_precedence(monkeypatch):
    monkeypatch.setenv("SPATIALSCRIBE_LLM_BASE_URL", "http://localhost:8000/v1/")
    monkeypatch.setenv("SPATIALSCRIBE_LLM_MODEL", "Qwen3.6-27B-FP8")
    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-test")  # present, but base URL wins
    assert llm.provider() == "openai"
    assert llm.available() is True
    assert llm.default_model() == "Qwen3.6-27B-FP8"
    assert llm._openai_base() == "http://localhost:8000/v1"  # trailing slash trimmed


def test_available_false_when_nothing_configured(monkeypatch):
    monkeypatch.delenv("SPATIALSCRIBE_LLM_BASE_URL", raising=False)
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    assert llm.available() is False
    assert llm.provider() == "anthropic"


def test_tools_anthropic_to_openai_shape():
    schema = {"type": "object", "properties": {"k": {"type": "integer"}}}
    oa = llm._tools_to_openai([{"name": "niches", "description": "call niches", "input_schema": schema}])
    assert len(oa) == 1
    assert oa[0]["type"] == "function"
    assert oa[0]["function"] == {"name": "niches", "description": "call niches", "parameters": schema}


def _fake_post(content, tool_calls=None):
    """A stand-in for llm._openai_post returning one assistant message (content + structured calls)."""
    def _post(messages, max_tokens, model, json_mode=False, tools=None):
        return {"choices": [{"message": {"role": "assistant", "content": content,
                                         "tool_calls": tool_calls or []}}]}
    return _post


_LS_TOOLS = [{"name": "load_section", "description": "load a section",
              "input_schema": {"type": "object", "properties": {"path": {"type": "string"}}}},
             {"name": "describe_sample", "description": "describe", "input_schema": {"type": "object"}}]


def test_openai_tool_chat_recovers_toolcall_from_fenced_markup(monkeypatch):
    # The the cluster vLLM intermittently emits the call as PROSE (a ```json block) with tool_calls EMPTY.
    # The client must recover it into a structured call so the copilot runs the tool, not echoes markup.
    monkeypatch.setenv("SPATIALSCRIBE_LLM_BASE_URL", "http://x/v1")
    monkeypatch.setenv("SPATIALSCRIBE_LLM_MODEL", "Qwen")
    markup = ('I will load the mouse brain data.\n```json\n'
              '{"tool_name": "load_section", "arguments": {"path": "/data/sec"}}\n```')
    monkeypatch.setattr(llm, "_openai_post", _fake_post(markup, tool_calls=[]))
    turn = llm.tool_chat("SYS", [{"role": "user", "content": "load it"}], _LS_TOOLS)
    assert turn["tool_calls"], "tool call emitted as fenced markup was not recovered"
    assert turn["tool_calls"][0]["name"] == "load_section"
    assert turn["tool_calls"][0]["input"] == {"path": "/data/sec"}
    # replayable as a clean assistant tool_calls turn (so the next round's tool_result ids match)
    asst = turn["assistant"]["_openai"]
    assert asst["tool_calls"][0]["function"]["name"] == "load_section"
    assert asst["tool_calls"][0]["id"] == turn["tool_calls"][0]["id"]


def test_openai_tool_chat_recovers_toolcall_from_tag_markup(monkeypatch):
    monkeypatch.setenv("SPATIALSCRIBE_LLM_BASE_URL", "http://x/v1")
    monkeypatch.setenv("SPATIALSCRIBE_LLM_MODEL", "Qwen")
    markup = 'thinking...\n<tool_call>\n{"name": "describe_sample", "arguments": {}}\n</tool_call>'
    monkeypatch.setattr(llm, "_openai_post", _fake_post(markup, tool_calls=[]))
    turn = llm.tool_chat("SYS", [{"role": "user", "content": "describe"}], _LS_TOOLS)
    assert [c["name"] for c in turn["tool_calls"]] == ["describe_sample"]


def test_openai_tool_chat_does_not_hijack_plain_json_answer(monkeypatch):
    # A normal answer that merely CONTAINS JSON (naming no tool) must stay a text answer, not be
    # mis-parsed into a phantom tool call.
    monkeypatch.setenv("SPATIALSCRIBE_LLM_BASE_URL", "http://x/v1")
    monkeypatch.setenv("SPATIALSCRIBE_LLM_MODEL", "Qwen")
    monkeypatch.setattr(llm, "_openai_post", _fake_post('The section has {"n_cells": 5} cells.', []))
    turn = llm.tool_chat("SYS", [{"role": "user", "content": "how many?"}], _LS_TOOLS)
    assert turn["tool_calls"] == []
    assert "n_cells" in turn["text"]


def test_openai_tool_chat_prefers_structured_calls_when_present(monkeypatch):
    # When the server DOES return structured tool_calls, they win untouched (no markup scan).
    monkeypatch.setenv("SPATIALSCRIBE_LLM_BASE_URL", "http://x/v1")
    monkeypatch.setenv("SPATIALSCRIBE_LLM_MODEL", "Qwen")
    tc = [{"id": "c1", "type": "function",
           "function": {"name": "load_section", "arguments": '{"path": "/p"}'}}]
    monkeypatch.setattr(llm, "_openai_post", _fake_post("", tool_calls=tc))
    turn = llm.tool_chat("SYS", [{"role": "user", "content": "load"}], _LS_TOOLS)
    assert turn["tool_calls"] == [{"id": "c1", "name": "load_section", "input": {"path": "/p"}}]


def test_openai_msgs_translates_tool_results():
    system = "SYS"
    messages = [
        {"role": "user", "content": "where are the T cells?"},
        {"role": "assistant", "_openai": {"role": "assistant", "content": None,
                                          "tool_calls": [{"id": "c1", "type": "function",
                                                          "function": {"name": "niches", "arguments": "{}"}}]}},
        {"role": "user", "content": [{"type": "tool_result", "tool_use_id": "c1",
                                      "content": '{"ok": true}', "is_error": False}]},
    ]
    out = llm._openai_msgs(system, messages)
    assert out[0] == {"role": "system", "content": "SYS"}
    assert out[1]["content"] == "where are the T cells?"
    assert out[2]["role"] == "assistant" and out[2]["tool_calls"][0]["id"] == "c1"
    assert out[3] == {"role": "tool", "tool_call_id": "c1", "content": '{"ok": true}'}


def _fake_openai(content, tool_calls=None):
    """A minimal OpenAI/vLLM chat-completions response with the given assistant content."""
    return {"choices": [{"message": {"content": content, "tool_calls": tool_calls or []}}]}


def test_openai_prose_tool_call_is_recovered_not_leaked(monkeypatch):
    # The the cluster Qwen vLLM sometimes emits its call as PROSE in content (Hermes JSON) instead of a
    # structured tool_calls array. It must be recovered (so the tool actually runs) and the raw
    # <tool_call> markup must never reach the user as an "answer".
    monkeypatch.setattr(llm, "provider", lambda: "openai")
    monkeypatch.setattr(llm, "_openai_post", lambda *a, **k: _fake_openai(
        '<tool_call>{"name": "composition_chart", "arguments": {}}</tool_call>'))
    out = llm.tool_chat("sys", [{"role": "user", "content": "how many T cells?"}],
                        [{"name": "composition_chart", "description": "", "input_schema": {}}])
    assert "<tool_call>" not in out["text"]                                   # never leak raw markup
    assert [c["name"] for c in out["tool_calls"]] == ["composition_chart"]    # recovered -> tool runs
    # the replayed assistant carries a matching structured tool_call so the next round stays valid
    replay = out["assistant"]["_openai"].get("tool_calls") or []
    assert replay and replay[0]["function"]["name"] == "composition_chart"
    assert replay[0]["id"] == out["tool_calls"][0]["id"]


def test_openai_prose_function_tag_tool_call_is_recovered(monkeypatch):
    # The other observed vLLM form: <tool_call><function>NAME</function></tool_call>.
    monkeypatch.setattr(llm, "provider", lambda: "openai")
    monkeypatch.setattr(llm, "_openai_post", lambda *a, **k: _fake_openai(
        "<tool_call><function>niches</function></tool_call>"))
    out = llm.tool_chat("sys", [{"role": "user", "content": "show niches"}],
                        [{"name": "niches", "description": "", "input_schema": {}}])
    assert "<tool_call>" not in out["text"]
    assert [c["name"] for c in out["tool_calls"]] == ["niches"]


def test_provider_override_toggles_between_configured_backends(monkeypatch):
    # The clickable model tag: set_provider flips Anthropic <-> the local vLLM at runtime, but only to
    # a CONFIGURED backend, with $SPATIALSCRIBE_LLM_PROVIDER as the default.
    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-test")
    monkeypatch.setenv("SPATIALSCRIBE_LLM_BASE_URL", "http://vllm:8081/v1")
    monkeypatch.setenv("SPATIALSCRIBE_LLM_MODEL", "Qwen3-X")
    monkeypatch.setenv("SPATIALSCRIBE_LLM_PROVIDER", "anthropic")
    try:
        llm.set_provider(None)                                          # start clean
        assert set(llm.providers_available()) == {"anthropic", "openai"}
        assert llm.provider() == "anthropic" and llm.default_model() == llm._ANTHROPIC_DEFAULT
        assert llm.set_provider("openai") is True
        assert llm.provider() == "openai" and llm.default_model() == "Qwen3-X"
        assert llm.set_provider("bogus") is False and llm.provider() == "openai"   # invalid -> no-op
        monkeypatch.delenv("SPATIALSCRIBE_LLM_BASE_URL")
        llm.set_provider(None)
        assert llm.set_provider("openai") is False                     # not configured -> refuse
    finally:
        llm.set_provider(None)                                          # never leak the override


def test_openai_function_equals_parameter_markup_is_recovered(monkeypatch):
    # The live 2026-07-11 form (the format DRIFTS per vLLM reload): the copilot's `load_section` call
    # arrived as <function=NAME><parameter name="KEY" ...>VALUE</parameter></function> with NO
    # tool_calls field, and the raw <parameter ...> markup was echoed to the user. It must recover the
    # call (with its args) AND strip the markup.
    monkeypatch.setattr(llm, "provider", lambda: "openai")
    monkeypatch.setattr(llm, "_openai_post", lambda *a, **k: _fake_openai(
        'I\'ll load that section.\n<function=load_section>\n'
        '<parameter name="path" schema="string">/data/projects/an internal benchmark/xenium_out/Region_1</parameter>\n'
        '</function>'))
    out = llm.tool_chat("sys", [{"role": "user", "content": "load /data"}],
                        [{"name": "load_section", "description": "", "input_schema": {}}])
    assert [c["name"] for c in out["tool_calls"]] == ["load_section"]
    assert out["tool_calls"][0]["input"]["path"] == "/data/projects/an internal benchmark/xenium_out/Region_1"
    assert "<parameter" not in out["text"] and "<function=" not in out["text"]   # markup never leaks


def test_openai_unparseable_tool_call_markup_is_stripped_not_leaked(monkeypatch):
    # Unrecoverable markup must be stripped (no leak) and leave no tool_calls, so run_copilot falls
    # through to its honest degrade instead of handing the user garbage.
    monkeypatch.setattr(llm, "provider", lambda: "openai")
    monkeypatch.setattr(llm, "_openai_post", lambda *a, **k: _fake_openai(
        "Here you go <tool_call>@@@</tool_call>"))
    out = llm.tool_chat("sys", [{"role": "user", "content": "hi"}],
                        [{"name": "niches", "description": "", "input_schema": {}}])
    assert "<tool_call>" not in out["text"] and "@@@" not in out["text"]
    assert out["tool_calls"] == []


def test_openai_dangling_closing_toolcall_tag_is_not_leaked(monkeypatch):
    # A per-reload Qwen vLLM state emitted a load_section call as LOOSE PROSE with a DANGLING </tool_call>
    # close tag - no opening tag, no JSON, no <function> tags (the exact shape a user reported leaking:
    # "...Region_1__20260306__134108 load_section </tool_call>"). The old strip only removed OPENING
    # <tool_call> tags and complete blocks, so the close tag AND the call-as-prose residue leaked verbatim.
    # When tool-call markup is present but nothing recovers, the turn was a mis-formatted call attempt,
    # not an answer -> degrade to empty text (run_copilot's honest fallback), never leak the residue.
    monkeypatch.setattr(llm, "provider", lambda: "openai")
    monkeypatch.setattr(llm, "_openai_post", lambda *a, **k: _fake_openai(
        "I'll load the spatial section from that path. "
        "/server-side/path/to/Xenium/output/Region_1__20260306__134108 load_section </tool_call>"))
    out = llm.tool_chat("sys", [{"role": "user", "content": "load Region_1"}],
                        [{"name": "load_section", "description": "", "input_schema": {}}])
    assert "</tool_call>" not in out["text"]          # the dangling close tag never reaches the user
    assert "load_section" not in out["text"]           # nor the call-as-prose residue
