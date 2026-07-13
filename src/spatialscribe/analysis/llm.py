"""LLM helpers - grounded annotation, panel/QC verdicts, subtype naming, and the copilot.

Every function passes the *actual computed numbers* (marker genes, coverage tables) into the
prompt and asks the model to reason only over them - never to invent expression. Outputs are
requested as JSON and parsed defensively.

Config - two interchangeable backends, chosen by environment (endpoint-agnostic)
-------------------------------------------------------------------------------
* **Anthropic** (default): ``ANTHROPIC_API_KEY`` [+ ``ANTHROPIC_MODEL``, default Haiku].
* **OpenAI-compatible** (any ``/v1`` server - a self-hosted vLLM/TGI/Ollama endpoint, or
  OpenAI itself): set ``SPATIALSCRIBE_LLM_BASE_URL`` (e.g. ``http://<your-vllm-host>:8000/v1``),
  ``SPATIALSCRIBE_LLM_MODEL`` (e.g. ``local-llm``), and optionally
  ``SPATIALSCRIBE_LLM_API_KEY`` (default ``dummy`` - vLLM ignores it). When the base URL is
  set it takes precedence over Anthropic. No extra SDK: we POST ``/chat/completions`` with httpx.
"""

from __future__ import annotations

import json
import os
import re

# Anthropic default: Haiku (cheap + fast for the interactive copilot / testing). Override with
# ANTHROPIC_MODEL=claude-sonnet-5 for the higher-quality model when it matters.
_ANTHROPIC_DEFAULT = os.environ.get("ANTHROPIC_MODEL", "claude-haiku-4-5")


def _openai_base() -> str | None:
    """OpenAI-compatible base URL (``.../v1``) if one is configured, else None (-> Anthropic)."""
    b = os.environ.get("SPATIALSCRIBE_LLM_BASE_URL")
    return b.rstrip("/") if b else None


# Runtime provider override (set via set_provider / the UI's clickable model tag), so the operator can
# flip Anthropic <-> the local vLLM WITHOUT a redeploy - e.g. when the vLLM's tool-calling drifts.
# None = fall back to $SPATIALSCRIBE_LLM_PROVIDER, then the base-URL heuristic.
_PROVIDER_OVERRIDE: str | None = None


def _configured(p: str) -> bool:
    """Is provider ``p`` usable in this process (its credentials/endpoint are set)?"""
    return bool(os.environ.get("ANTHROPIC_API_KEY")) if p == "anthropic" else bool(_openai_base())


def providers_available() -> list[str]:
    """The providers configured in this process, so the UI can offer a toggle only when >1 exists."""
    return [p for p in ("anthropic", "openai") if _configured(p)]


def provider() -> str:
    """Active backend. Precedence: runtime override -> ``$SPATIALSCRIBE_LLM_PROVIDER`` -> the base-URL
    heuristic (``'openai'`` when an OpenAI-compatible base URL is set, else ``'anthropic'``)."""
    if _PROVIDER_OVERRIDE in ("anthropic", "openai"):
        return _PROVIDER_OVERRIDE
    env = (os.environ.get("SPATIALSCRIBE_LLM_PROVIDER") or "").strip().lower()
    if env in ("anthropic", "openai"):
        return env
    return "openai" if _openai_base() else "anthropic"


def set_provider(p: str | None) -> bool:
    """Set the runtime provider override. Only switches to a CONFIGURED provider (else no-op -> False);
    ``None`` clears the override (back to env/heuristic). Process-wide, resets on restart."""
    global _PROVIDER_OVERRIDE
    if p is None:
        _PROVIDER_OVERRIDE = None
        return True
    p = str(p).strip().lower()
    if p in ("anthropic", "openai") and _configured(p):
        _PROVIDER_OVERRIDE = p
        return True
    return False


def available() -> bool:
    """True when the ACTIVE provider is configured - the single gate the app uses for LLM features."""
    return _configured(provider())


def default_model() -> str:
    """The model id for the active backend (Anthropic default, or the configured OpenAI model)."""
    if provider() == "openai":
        return os.environ.get("SPATIALSCRIBE_LLM_MODEL", "")
    return _ANTHROPIC_DEFAULT


# Back-compat module attribute (some callers read llm.DEFAULT_MODEL). Anthropic default; the
# OpenAI path resolves its model from SPATIALSCRIBE_LLM_MODEL at call time.
DEFAULT_MODEL = _ANTHROPIC_DEFAULT


def get_client():
    """Return an Anthropic client (raises a clear error if the key is missing)."""
    import anthropic

    if not os.environ.get("ANTHROPIC_API_KEY"):
        raise RuntimeError("ANTHROPIC_API_KEY is not set - export it before using the copilot.")
    return anthropic.Anthropic()


def _openai_post(messages: list, max_tokens: int, model: str | None,
                 json_mode: bool = False, tools: list | None = None) -> dict:
    """POST an OpenAI-compatible ``/chat/completions`` and return the raw JSON dict.

    Sends vLLM's ``chat_template_kwargs={"enable_thinking": false}`` - reasoning models otherwise
    spend the token budget on a hidden reasoning trace and truncate the real answer. If the server
    rejects that extension (a strict OpenAI endpoint -> HTTP 400) we retry without it, so the layer
    stays endpoint-agnostic. ``json_mode`` requests a strict JSON object; ``tools`` (OpenAI schema)
    enables function-calling.
    """
    import httpx

    base = _openai_base()
    if not base:
        raise RuntimeError("SPATIALSCRIBE_LLM_BASE_URL is not set.")
    model = model or os.environ.get("SPATIALSCRIBE_LLM_MODEL")
    if not model:
        raise RuntimeError("SPATIALSCRIBE_LLM_MODEL is not set (required for the OpenAI backend).")
    key = os.environ.get("SPATIALSCRIBE_LLM_API_KEY", "dummy")
    body: dict = {"model": model, "messages": messages, "max_tokens": max_tokens, "temperature": 0}
    if json_mode:
        body["response_format"] = {"type": "json_object"}
    if tools:
        body["tools"] = tools
        body["tool_choice"] = "auto"
    url = f"{base}/chat/completions"
    headers = {"Authorization": f"Bearer {key}"}
    with httpx.Client(timeout=180) as c:
        r = c.post(url, json={**body, "chat_template_kwargs": {"enable_thinking": False}},
                   headers=headers)
        if r.status_code == 400:  # strict server rejected the vLLM reasoning-off extension
            r = c.post(url, json=body, headers=headers)
        r.raise_for_status()
        return r.json()


def _openai_text(data: dict) -> str:
    """Assistant text from an OpenAI chat response, stripping any stray ``<think>`` reasoning trace."""
    try:
        txt = data["choices"][0]["message"].get("content") or ""
    except (KeyError, IndexError, TypeError):
        return ""
    return re.sub(r"<think>.*?</think>", "", txt, flags=re.DOTALL).strip()


def complete(system: str, user: str, max_tokens: int = 1024, model: str | None = None,
             json_mode: bool = False) -> str:
    """Single-turn completion; returns the assistant text. Dispatches to the active backend.

    ``json_mode`` asks OpenAI-compatible servers for a strict JSON object (``response_format``) -
    important for reasoning models like local-llm that otherwise wrap or truncate JSON; it is a no-op on
    Anthropic (the prompts already request JSON and ``_parse_json`` is defensive either way).

    Cost (Anthropic): the system prompt is sent as a cache-control 'ephemeral' block so prompt
    caching charges the large, reused system prompt at the cheap cache-read rate on repeat calls.
    """
    if provider() == "openai":
        data = _openai_post(
            [{"role": "system", "content": system}, {"role": "user", "content": user}],
            max_tokens=max_tokens, model=model, json_mode=json_mode)
        return _openai_text(data)

    client = get_client()
    model = model or _ANTHROPIC_DEFAULT
    try:
        msg = client.messages.create(
            model=model, max_tokens=max_tokens,
            system=[{"type": "text", "text": system, "cache_control": {"type": "ephemeral"}}],
            messages=[{"role": "user", "content": user}],
        )
    except Exception:
        msg = client.messages.create(
            model=model, max_tokens=max_tokens, system=system,
            messages=[{"role": "user", "content": user}],
        )
    return "".join(block.text for block in msg.content if getattr(block, "type", None) == "text")


def tool_chat(system: str, messages: list, tools: list, max_tokens: int = 1024,
              model: str | None = None) -> dict:
    """One turn of a tool-use chat, normalized across backends, for the copilot.

    ``messages`` and ``tools`` are in **Anthropic** shape (the registry's native format). For the
    OpenAI backend they are translated on the way out and the reply is normalized back. Returns
    ``{"text": str, "tool_calls": [{"id", "name", "input"}], "assistant": <native asst message>}``
    so :func:`spatialscribe.agent.tools.run_copilot` stays provider-agnostic.
    """
    if provider() == "openai":
        return _openai_tool_chat(system, messages, tools, max_tokens, model)

    client = get_client()
    msg = client.messages.create(model=model or _ANTHROPIC_DEFAULT, max_tokens=max_tokens,
                                 system=system, tools=tools, messages=messages)
    text = "".join(b.text for b in msg.content if getattr(b, "type", None) == "text")
    calls = [{"id": b.id, "name": b.name, "input": b.input}
             for b in msg.content if getattr(b, "type", None) == "tool_use"]
    return {"text": text, "tool_calls": calls, "assistant": {"role": "assistant", "content": msg.content}}


def _tools_to_openai(tools: list) -> list:
    """Anthropic tool schema ``[{name, description, input_schema}]`` -> OpenAI ``functions``."""
    return [{"type": "function",
             "function": {"name": t["name"], "description": t.get("description", ""),
                          "parameters": t.get("input_schema", {"type": "object", "properties": {}})}}
            for t in tools]


def _openai_msgs(system: str, messages: list) -> list:
    """Anthropic-shaped copilot messages -> OpenAI messages (system + user/assistant/tool turns).

    The copilot threads two message shapes: an assistant turn we stored as ``{"role":"assistant",
    "_openai": <native>}`` and tool results as ``{"role":"user","content":[{type:tool_result,...}]}``.
    Plain string-content user turns pass through unchanged.
    """
    out: list = [{"role": "system", "content": system}]
    for m in messages:
        if m.get("_openai") is not None:            # a normalized assistant turn from a prior round
            out.append(m["_openai"])
            continue
        content = m.get("content")
        if isinstance(content, list):               # tool_result blocks -> OpenAI 'tool' messages
            for blk in content:
                if isinstance(blk, dict) and blk.get("type") == "tool_result":
                    out.append({"role": "tool", "tool_call_id": blk["tool_use_id"],
                                "content": blk["content"] if isinstance(blk["content"], str)
                                else json.dumps(blk["content"], default=str)})
        else:
            out.append({"role": m.get("role", "user"), "content": content})
    return out


_TOOLCALL_TAG = re.compile(r"<tool_call>\s*(.+?)\s*</tool_call>", re.DOTALL)
_TOOLCALL_FENCE = re.compile(r"```(?:json|tool_call|python|tool_code)?\s*(.+?)\s*```", re.DOTALL)
_FUNCTION_TAG = re.compile(r"<function>\s*([\w.-]+)\s*</function>", re.DOTALL)
# vLLM/local-llm "pythonic" tool markup: <function=NAME> <parameter name="KEY" ...>VALUE</parameter> ... </function>
# (also the <parameter=KEY>VALUE</parameter> attribute-less variant). vLLM leaves this in the content
# when its tool-call parser is off/mismatched, with NO structured tool_calls field.
_FUNCTION_EQ = re.compile(r"<function=([\w.-]+)\s*>(.*?)(?=</function>|<function=|$)", re.DOTALL)
_PARAM_TAG = re.compile(r"<parameter(?:\s+name=[\"']?([\w.-]+)[\"']?[^>]*|=([\w.-]+))\s*>(.*?)</parameter>", re.DOTALL)


def _recover_tool_calls(text: str, valid_names: set) -> list:
    """Recover tool calls a model emitted as PROSE MARKUP instead of the OpenAI ``tool_calls`` field.

    The vLLM's function-calling silently stops emitting structured ``tool_calls`` at some
    (serving-config-dependent, per-reload) schema size and instead prints the call as text -
    ``<tool_call>{...}</tool_call>``, a ```json fenced ``{"tool_name"/"name": ..., "arguments": ...}``
    block, or ``<function>NAME</function>`` - which the client would otherwise surface to the user as
    raw markup (see ``the vLLM markup-as-prose fallback``). Scans those shapes and returns
    ``[{"id","name","input"}]`` for any candidate that NAMES A REAL TOOL (``valid_names`` guards
    against mis-reading a legitimate JSON answer as a call). Never raises.
    """
    if not text or ("{" not in text and "<function>" not in text and "<function=" not in text):
        return []
    # Candidate JSON payloads, best-signal first: explicit <tool_call> tags, then fenced code blocks,
    # then (only if neither matched) the outermost bare {...} object.
    cands = _TOOLCALL_TAG.findall(text) + _TOOLCALL_FENCE.findall(text)
    if not cands and "{" in text:
        m = re.search(r"\{.*\}", text, re.DOTALL)
        if m:
            cands.append(m.group(0))
    out: list = []
    for raw in cands:
        try:
            obj = json.loads(raw)
        except Exception:
            continue
        for d in (obj if isinstance(obj, list) else [obj]):
            if not isinstance(d, dict):
                continue
            fn = d.get("function")
            if isinstance(fn, dict):            # {"function": {"name":..., "arguments":...}}
                d = fn
            name = d.get("name") or d.get("tool_name") or (fn if isinstance(fn, str) else None)
            args = d.get("arguments")
            if args is None:
                args = d.get("parameters", d.get("input", {}))
            if isinstance(args, str):
                try:
                    args = json.loads(args)
                except Exception:
                    args = {}
            if name in valid_names and isinstance(args, dict):
                out.append({"id": f"recovered-{len(out)}", "name": name, "input": args})
    # <function>NAME</function> form (local-llm tool-code): the args are a trailing JSON object, if any.
    if not out:
        for m in _FUNCTION_TAG.finditer(text):
            name = m.group(1)
            if name not in valid_names:
                continue
            am = re.search(r"\{.*\}", text[m.end():], re.DOTALL)
            args = {}
            if am:
                try:
                    args = json.loads(am.group(0))
                except Exception:
                    args = {}
            out.append({"id": f"recovered-{len(out)}", "name": name,
                        "input": args if isinstance(args, dict) else {}})
    # <function=NAME><parameter name="KEY">VALUE</parameter>...</function> form (the local-llm vLLM's
    # actual leak: seen live emitting a load_section call as this markup with no tool_calls field).
    if not out:
        for fm in _FUNCTION_EQ.finditer(text):
            name = fm.group(1)
            if name not in valid_names:
                continue
            args = {}
            for pm in _PARAM_TAG.finditer(fm.group(2)):
                key = pm.group(1) or pm.group(2)
                if key:
                    args[key] = pm.group(3).strip()
            out.append({"id": f"recovered-{len(out)}", "name": name, "input": args})
    return out


def _openai_tool_chat(system: str, messages: list, tools: list, max_tokens: int,
                      model: str | None) -> dict:
    data = _openai_post(_openai_msgs(system, messages), max_tokens=max_tokens, model=model,
                        tools=_tools_to_openai(tools))
    choice = (data.get("choices") or [{}])[0]
    asst = choice.get("message", {}) or {}
    text = re.sub(r"<think>.*?</think>", "", asst.get("content") or "", flags=re.DOTALL).strip()
    calls = []
    for tc in asst.get("tool_calls") or []:
        fn = tc.get("function", {})
        try:
            args = json.loads(fn.get("arguments") or "{}")
        except Exception:
            args = {}
        calls.append({"id": tc.get("id"), "name": fn.get("name"), "input": args})
    # Recovery: the server returned NO structured call but printed one as markup in the content (a
    # per-reload vLLM failure). Parse it back into a structured call so the tool actually RUNS.
    # A per-reload vLLM state can emit the call as loose prose with a DANGLING </tool_call> close tag
    # (no opening tag, no JSON, no <function> tags) - "... load_section </tool_call>". Note it BEFORE the
    # strips: any tool_call marker means the turn was a tool-call ATTEMPT, so if nothing recovers below we
    # must not leak the residue as if it were an answer.
    had_toolcall_markup = bool(re.search(r"</?tool_call>", text))
    recovered_from_prose = False
    if not calls:
        recovered = _recover_tool_calls(text, {t["name"] for t in tools})
        if recovered:
            calls, recovered_from_prose = recovered, True
    # Never surface raw <tool_call> markup to the user: strip whole blocks AND any truncated
    # (max_tokens) open tag. Unrecoverable markup -> empty text -> run_copilot's honest degrade.
    text = re.sub(r"<tool_call>.*$", "", _TOOLCALL_TAG.sub("", text), flags=re.DOTALL).strip()
    # Same for the <function=...>/<parameter ...> markup (recovered or not): strip a whole function
    # block plus any stray function/parameter tags so partial/unrecoverable markup never reaches the user.
    text = re.sub(r"<function=.*$", "", text, flags=re.DOTALL)
    text = re.sub(r"</?(?:function|parameter)[^>]*>", "", text).strip()
    # The block/open-tag strips above miss a DANGLING </tool_call> (close tag, no opening) and leave the
    # call-as-prose residue ("... load_section"). When the turn carried tool-call markup but recovered NO
    # call, it was a mis-formatted call attempt, not an answer: blank the residue so it degrades honestly.
    if had_toolcall_markup and not calls:
        text = ""
    if recovered_from_prose:
        # Rebuild a CLEAN assistant turn (structured tool_calls, markup gone) so the next round's
        # history replay and tool_result id-matching work as if the server had emitted it structurally.
        text = ""
        asst = {"role": "assistant", "content": None, "tool_calls": [
            {"id": c["id"], "type": "function",
             "function": {"name": c["name"], "arguments": json.dumps(c["input"])}}
            for c in calls]}
    # keep the (possibly rewritten) assistant message so the next round can replay it as history
    return {"text": text, "tool_calls": calls, "assistant": {"role": "assistant", "_openai": asst}}


def _parse_json(text: str):
    """Best-effort JSON extraction from a model reply; returns None if unrecoverable.

    Must NEVER raise: a truncated/malformed reply (e.g. the model hit ``max_tokens``
    mid-object) has to degrade to None so annotation falls back to marker-only labels,
    instead of crashing the whole ``annotate`` capability.
    """
    try:
        return json.loads(text)
    except Exception:
        pass
    m = re.search(r"\{.*\}|\[.*\]", text, re.DOTALL)
    if not m:
        return None
    try:
        return json.loads(m.group(0))
    except Exception:
        return None


# --------------------------------------------------------------------------- #
# Annotation
# --------------------------------------------------------------------------- #
_ANNOTATE_SYS = (
    "You are a spatial-transcriptomics annotation expert. You are given the top marker "
    "genes for each cluster from a targeted imaging panel. Assign the most likely cell "
    "type to each cluster. Reason ONLY from the provided markers - do not assume genes "
    "that are not listed are absent or present. Tissue context: {context}. "
    'Reply as JSON: {{"<cluster>": {{"label": str, "confidence": "high|medium|low", '
    '"rationale": str}}}}.'
)


# Closed-vocabulary clause. The label is constrained AT GENERATION TIME, never snapped afterwards:
# string-nearest is not biology-nearest ("pDC" is near nothing in a 10-lineage set and would snap to a
# wrong lineage). The escape hatch is mandatory - a section whose panel genuinely misses a lineage must
# be able to say so instead of being forced into a confident wrong label.
_ANNOTATE_VOCAB = (
    ' The "label" MUST be copied EXACTLY from this list: {allowed}.'
    ' If none of them fits the markers, use "{novel}". Never invent a label outside the list.'
)
NOVEL_LABEL = "Novel / unknown"


def annotate_clusters(markers_per_cluster: dict[str, list[str]],
                      context: str = "human skin / melanoma",
                      allowed_labels: list[str] | None = None) -> dict:
    """Claude-as-annotator: cluster -> {label, confidence, rationale} from markers.

    ``allowed_labels`` closes the label vocabulary in the PROMPT (plus the ``NOVEL_LABEL`` escape
    hatch). Callers must still validate the reply - a model can ignore an instruction - which
    :func:`annotate.consensus_annotate` does by falling back to the marker label.
    """
    sys = _ANNOTATE_SYS.format(context=context)
    if allowed_labels:
        allowed = list(dict.fromkeys([*allowed_labels, NOVEL_LABEL]))
        sys += _ANNOTATE_VOCAB.format(allowed=", ".join(f'"{a}"' for a in allowed), novel=NOVEL_LABEL)
    user = "Clusters and their top markers:\n" + "\n".join(
        f"- Cluster {k}: {', '.join(v)}" for k, v in markers_per_cluster.items()
    )
    # Budget scales with cluster count (label+confidence+rationale is ~120 tokens/cluster);
    # too small a cap truncates the JSON mid-object and the whole reply is lost. Floor 1500,
    # ~160/cluster, cap 8000 (Haiku's output ceiling).
    max_tokens = int(min(8000, max(1500, 200 * len(markers_per_cluster) + 400)))
    out = _parse_json(complete(sys, user, max_tokens=max_tokens, json_mode=True))
    return out or {}


def name_subtypes(parent_type: str, markers_per_subcluster: dict[str, list[str]],
                  context: str = "human skin / melanoma") -> dict:
    """Name subclusters of one parent cell type (H2), e.g. T cell -> CD4/CD8/Treg."""
    sys = (
        f"You are subtyping cells already annotated as '{parent_type}' in {context}. "
        "Given the top markers per subcluster, name the most likely subtype. Reason only "
        'from the markers. Reply JSON: {"<subcluster>": {"label": str, '
        '"confidence": "high|medium|low", "rationale": str}}.'
    )
    user = "\n".join(f"- Subcluster {k}: {', '.join(v)}" for k, v in markers_per_subcluster.items())
    # Scale the token budget with the number of subclusters (as annotate_clusters does). A fixed small
    # cap truncates the JSON mid-object for many subclusters, so _parse_json returns None and every
    # subtype falls back to "subcluster N" with confidence "n/a" - the observed all-n/a bug.
    max_tokens = int(min(8000, max(1200, 200 * len(markers_per_subcluster) + 400)))
    return _parse_json(complete(sys, user, max_tokens=max_tokens, json_mode=True)) or {}


# --------------------------------------------------------------------------- #
# De-novo program naming - a plain-language label per NMF program
# --------------------------------------------------------------------------- #
_PROGRAM_SYS = (
    "You are naming de-novo gene expression programs discovered by NMF in {context}. For each "
    "program you are given its top loading genes (highest-weighted first). Give a short plain-language "
    "name (under ~4 words) for the biological program those genes suggest (e.g. 'interferon response', "
    "'cell cycle', 'extracellular matrix', 'hypoxia'), or 'mixed / unclear' when there is no coherent "
    "theme. Reason ONLY from the listed genes - do not invent genes. Write in plain prose using "
    'hyphens, never em-dashes. Reply ONLY as JSON: {{"<program>": {{"label": str, "confidence": '
    '"high|medium|low", "rationale": str}}}}.'
)


def name_programs(top_genes_per_program: dict[str, list[str]], context: str = "human tumour") -> dict:
    """Name each de-novo NMF program from its top loading genes.

    ``top_genes_per_program`` = ``{"Program 0": [genes...], ...}``. Returns
    ``{"Program k": {label, confidence, rationale}}`` (empty ``{}`` on an unrecoverable reply / no key)
    so relabeling degrades cleanly to "Program k"."""
    user = "\n".join(f"- {k}: {', '.join(v)}" for k, v in top_genes_per_program.items())
    max_tokens = int(min(6000, max(1200, 160 * len(top_genes_per_program) + 400)))
    out = _parse_json(complete(_PROGRAM_SYS.format(context=context), user,
                               max_tokens=max_tokens, json_mode=True))
    return out if isinstance(out, dict) else {}


# --------------------------------------------------------------------------- #
# Panel-adequacy verdict (H3)
# --------------------------------------------------------------------------- #
_PANEL_SYS = (
    "You are advising a wet-lab scientist on the limits of their targeted spatial panel. "
    "For each cell type you are given a COMPUTED discriminability score and a `confidently_typable` "
    "flag - THIS is what your judgment must rest on, NOT how many markers are on the panel. The "
    "score's `basis` is one of: `depth_matched_f1` (a single-cell reference was provided - the "
    "fraction of that type's cells correctly recovered when the reference is thinned to THIS "
    "section's real sequencing depth; the gold standard, panel-size aware), `identifiability_auc` "
    "(one-vs-rest ROC-AUC of the type's marker score on this section; ~1 = separable, ~0.5 = not), "
    "or `coverage_only` (marker PRESENCE only - WEAK evidence, use it only when no computed score "
    "exists). Rules: a type with its markers present but a low F1/AUC is NOT confidently typable - "
    "say so; call out every type flagged `confidently_typable=false`; prefer the `depth_matched_f1` "
    "basis when present. Also note pairs the panel cannot separate (consider merging). Ground every "
    "statement in the provided numbers; a marker being on-panel is necessary but not sufficient. Do "
    "not invent genes or cell types. Write in plain prose using hyphens, never em-dashes."
)


def panel_verdict(panel_check_result: dict, typability: list | None = None) -> str:
    """Turn a ``panel_check.check_panel`` result (+ the computed per-type typability table) into a
    grounded plain-language verdict. When ``typability`` is given, the verdict is grounded in the
    COMPUTED discriminability (depth-matched F1 / identifiability AUC), not marker presence."""
    cov = {ct: {"n_present": d["n_present"], "n_markers": d["n_markers"],
                "status": d["status"], "missing": d["missing"][:8]}
           for ct, d in panel_check_result.get("coverage", {}).items()}
    payload = {
        "typability_per_cell_type": typability or [],
        "coverage": cov,
        "confusable_pairs": panel_check_result.get("confusable_pairs", []),
        "merge_groups": panel_check_result.get("merge_groups", []),
    }
    return complete(_PANEL_SYS, json.dumps(payload, indent=2, default=str), max_tokens=900)


# --------------------------------------------------------------------------- #
# Cell-state naming (CyteType-style) - the dominant functional program per group
# --------------------------------------------------------------------------- #
_STATE_SYS = (
    "You are annotating the functional STATE of cells (orthogonal to their lineage identity) in "
    "{context}. For each group you are given its top marker genes and its mean z-scores for canonical "
    "state programs (cycling, interferon/ISG, hypoxia, stress/heat-shock, EMT, T-exhaustion, "
    "T-cytotoxicity, ECM remodeling, and antigen presentation). Name the single dominant state in plain "
    "language - e.g. 'cycling', 'interferon-activated', 'hypoxic', 'stressed', 'EMT-like', 'exhausted', "
    "'cytotoxic', 'ECM-remodeling', 'antigen-presenting', or 'no dominant program (resting)'. Reason "
    "ONLY from the provided scores and markers - a high positive z-score means the "
    "program is elevated in that group; near-zero means it is not. Do not invent genes. Write in plain "
    'prose using hyphens, never em-dashes. Reply ONLY as JSON: {{"<group>": {{"state": str, '
    '"confidence": "high|medium|low", "rationale": str}}}}.'
)


def name_cell_states(per_group: dict, context: str = "human tumour") -> dict:
    """Name each group's dominant functional state from its marker genes + program z-scores.

    ``per_group`` = ``{group: {"markers": [...], "state_scores": {program: z}}}``. Returns
    ``{group: {state, confidence, rationale}}`` (empty on an unrecoverable reply). Prompt-cached."""
    max_tokens = int(min(6000, max(1200, 160 * len(per_group) + 400)))
    out = _parse_json(complete(_STATE_SYS.format(context=context),
                               json.dumps(per_group, indent=2, default=str),
                               max_tokens=max_tokens, json_mode=True))
    return out if isinstance(out, dict) else {}


# --------------------------------------------------------------------------- #
# QC verdict - plain-language "can I trust this section's signal?"
# --------------------------------------------------------------------------- #
_QC_SYS = (
    "You are advising a wet-lab scientist on whether their spatial section's signal is good enough "
    "to trust the cell-type annotation. You are given the six-layer QC funnel headline (per-layer "
    "flags, the panel-indexed count floor and how much it removes, segmentation, spatial coherence, "
    "and the pass/tentative/abstained annotatability split) and, if present, the top reasons cells "
    "could not be confidently typed. Write a short plain-language verdict (3-5 sentences): is the "
    "section trustworthy overall, which layer is the weakest link, and one concrete next step (e.g. "
    "raise the count floor, exclude a low-quality region, or coarsen a confusable label). Ground "
    "every statement in the numbers provided - do not invent metrics. Write in plain prose using "
    "hyphens, never em-dashes."
)


def qc_verdict(funnel_headline: dict, rejection_breakdown: list | None = None,
               section_metrics: dict | None = None) -> str:
    """Turn the QC funnel headline (+ optional rejection breakdown) into a grounded plain-language
    verdict on whether the section's signal can be trusted for annotation."""
    payload = {
        "funnel": funnel_headline or {},
        "top_rejection_reasons": (rejection_breakdown or [])[:6],
        "section_metrics": section_metrics or {},
    }
    return complete(_QC_SYS, json.dumps(payload, indent=2, default=str), max_tokens=700)


# --------------------------------------------------------------------------- #
# Free-text tissue -> cell-type marker panel (any tissue / organism)
# --------------------------------------------------------------------------- #
_MARKER_PANEL_SYS = (
    "You are a spatial-transcriptomics expert assembling a cell-type marker panel for one tissue. "
    "Given the tissue/organism context and the genes actually on the targeted panel, return the "
    "major cell types expected in that tissue and, for each, its canonical marker genes. PREFER "
    "genes from the provided panel list; a few canonical markers not in the list are acceptable. "
    "Use official gene symbols exactly as a reference database would write them. Cover the 6-15 "
    "major cell types, not rare states. Do not invent gene symbols. "
    "If the user provides an explicit list of cell types, return markers for EXACTLY those types "
    "and use their names verbatim as the JSON keys (do not add, drop, merge, or rename any). "
    'Reply ONLY as JSON: {"<cell type>": ["GENE1", "GENE2", ...]}.'
)


def generate_marker_panel(tissue: str, panel_genes: list[str], max_types: int = 15,
                          cell_types: list[str] | None = None) -> dict:
    """Claude proposes ``{cell_type: [markers]}`` grounded to the panel.

    Without ``cell_types`` it proposes the major types for a free-text ``tissue``. With
    ``cell_types`` (the section's CURRENT annotated categories) it returns markers for exactly
    those types, keyed verbatim - so panel-check reflects the labels actually on the section.

    The panel list is stride-sampled to at most 1200 genes for the prompt (a 5K/WTA panel would
    otherwise blow the context, and an alphabetical head would bias the sample); the caller filters
    the reply to on-panel genes regardless. Never raises for a bad reply - returns ``{}`` so the
    resolver falls back to a curated set.
    """
    genes = [str(g) for g in panel_genes]
    step = max(1, len(genes) // 1200)
    shown = genes[::step][:1200]
    cts = [str(c) for c in (cell_types or []) if str(c).strip()]
    limit = len(cts) or max_types
    user = (f"Tissue / organism context: {tissue}\n\n"
            + (f"Return markers for EXACTLY these cell types (verbatim keys): "
               f"{', '.join(cts)}\n\n" if cts else "")
            + f"Panel genes ({len(genes)} total"
            + f"{f'; {len(shown)} sampled across the panel' if len(genes) > len(shown) else ''}):\n"
            + ", ".join(shown))
    out = _parse_json(complete(_MARKER_PANEL_SYS, user, max_tokens=2000, json_mode=True))
    if not isinstance(out, dict):
        return {}
    clean: dict[str, list[str]] = {}
    for ct, gl in list(out.items())[:limit]:
        if isinstance(gl, list):
            clean[str(ct)] = [str(g) for g in gl if isinstance(g, str)]
    return clean


# --------------------------------------------------------------------------- #
# Grouping cell-type NAMES into biologically related categories (marker dot-plot rows)
# --------------------------------------------------------------------------- #
_GROUP_SYS = (
    "You are organizing the cell types found in a {context} spatial section into a small number of "
    "biologically related categories, so related types sit next to each other in a plot (e.g. CD4 T "
    "and CD8 T under one lymphoid category). You are given the EXACT list of cell-type names present. "
    "Group EVERY given name into exactly one category; use the names VERBATIM (do not rename, split, "
    "merge, add, or invent a type). Choose 4-8 short category names that fit the tissue. This is a "
    "NAME-grouping task only - do not comment on genes or expression. "
    'Reply ONLY as JSON: {{"<category>": ["<cell type>", ...], ...}}.'
)

# Grouping only depends on the NAME list + tissue, so cache per (tissue, sorted names). Small dict.
_GROUP_CACHE: dict[tuple, dict[str, list[str]]] = {}


def group_cell_types(cell_types, tissue: str = "") -> dict[str, list[str]]:
    """Group cell-type NAMES into biologically related categories via the LLM (marker dot-plot rows).

    Returns ``{category: [cell_type, ...]}``. GROUNDED: every returned cell type is one of the inputs
    (strays the model invents are dropped) and every input lands in exactly ONE category (anything the
    model left unassigned goes to ``"Other"``), so the result is total over the inputs. The LLM only
    groups NAMES - it never picks genes or numbers. Returns ``{}`` when no LLM backend is configured,
    so callers fall back to the deterministic keyword grouping. Cached per ``(tissue, sorted types)``.
    """
    cts = [str(c) for c in (cell_types or []) if str(c).strip()]
    if not cts or not available():
        return {}
    key = ((tissue or "").strip().lower(), tuple(sorted(cts)))
    if key in _GROUP_CACHE:
        return _GROUP_CACHE[key]
    want = {c.lower(): c for c in cts}                       # re-map to the requested names verbatim
    user = "Cell types present (group all of these, verbatim):\n" + "\n".join(f"- {c}" for c in cts)
    out = _parse_json(complete(_GROUP_SYS.format(context=tissue or "this"), user,
                               max_tokens=int(min(4000, 60 * len(cts) + 400)), json_mode=True))
    groups: dict[str, list[str]] = {}
    assigned: set[str] = set()
    if isinstance(out, dict):
        for cat, members in out.items():
            if not isinstance(members, list):
                continue
            keep = []
            for m in members:
                orig = want.get(str(m).strip().lower())      # only real inputs, once each
                if orig is not None and orig not in assigned:
                    keep.append(orig)
                    assigned.add(orig)
            if keep:
                groups[str(cat)] = keep
    leftover = [c for c in cts if c not in assigned]          # totality: every input lands somewhere
    if leftover:
        groups.setdefault("Other", []).extend(leftover)
    _GROUP_CACHE[key] = groups
    return groups
