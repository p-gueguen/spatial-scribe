"""The copilot (H4): a Claude tool-use loop over the capability registry.

The tool schema, the dispatch, and the structured errors all come from
:mod:`spatialscribe.analysis.capabilities` - the single source of truth the wizard also
uses. Claude answers a plain-language question by calling whitelisted capabilities that run
the *real* analysis and return computed numbers (grounded, never invented). Each successful
call is appended to the shared action log (same provenance record the wizard emits), so
copilot-driven analyses land in the re-runnable export too.
"""

from __future__ import annotations

import json


# Friendly progress labels for the tools that take a visible moment, so the copilot stream shows what
# it is DOING instead of a static "thinking …" (a load + reference select can be ~10s). Others fall
# back to the tool name.
_TOOL_STATUS = {
    "load_section": "loading the section …", "auto_select_reference": "selecting a reference …",
    "reference_transfer": "transferring labels …", "annotation_strategy": "planning annotation …",
    "annotate": "annotating cells …", "cluster": "clustering …", "niches": "finding niches …",
    "panel_check": "checking the panel …", "qc_funnel": "running QC …", "subcluster": "subclustering …",
    "malignant_concordance": "calling malignant cells …", "split_purify": "purifying spillover …",
    "neighborhood_enrichment": "scoring neighborhoods …", "assign_cell_states": "typing cell states …",
}


def run_copilot(adata, prompt: str, tissue: str = "melanoma", max_turns: int = 6,
                action_log: list | None = None, use_llm: bool = False,
                artifacts: list | None = None, reference=None, ref_label_key: str | None = None,
                loaded: list | None = None, on_status=None) -> str:
    """Answer a plain-language question by letting Claude call capability tools.

    Parameters
    ----------
    action_log: optional session action log; successful tool calls append their provenance
        record ``{name, params}`` here so the export captures copilot-driven steps too.
    artifacts: optional list; view-producing capabilities (``highlight_cells``) append figure
        artifacts ``{kind, engine, title, fig}`` here so the app can render the copilot's views.
    loaded: optional list; the ``load_section`` capability appends the section it loaded here (the
        newest is at the end). The copilot rebinds to it for the rest of the turn (so a same-turn
        "load X and annotate it" annotates the NEW section), and the caller swaps it onto the session.
    """
    from ..analysis import capabilities as cap
    from ..analysis import llm

    ctx = cap.RunContext(tissue=tissue, use_llm=use_llm,
                         artifacts=artifacts if artifacts is not None else [],
                         reference=reference, ref_label_key=ref_label_key,
                         loaded=loaded if loaded is not None else [])
    tools = cap.copilot_tools()

    cats = list(adata.obs["cell_type"].astype("category").cat.categories) if "cell_type" in adata.obs else []
    system = (
        "You are SpatialScribe's copilot for a wet-lab scientist analyzing an imaging-based "
        f"spatial ({tissue}) section. Behave strictly:\n"
        "1. ALWAYS answer by calling a tool - the tools run the real analysis and return computed "
        "numbers. Never state a number, fraction, cell type, or gene fact that did not come from a "
        "tool result; if no tool can answer, say so plainly instead of speculating.\n"
        "2. When the user asks to SEE / show / highlight / color / plot / 'where are' something, call "
        "the matching VIEW tool (highlight_cells, show_spatial, marker_dotplot, expression_violin, "
        "composition_chart) so a plot is drawn - do not merely describe it in words.\n"
        "3. Report a tool's own honest note verbatim when it has one (e.g. 'no cells were flagged as "
        "low-quality; showing the lowest-count decile') - do not smooth it over.\n"
        "4. On an 'error_type':'prerequisite_missing' error, tell the user which step to run (its "
        "'hint'), do not guess.\n"
        "Be terse - a sentence or two plus the plot. "
        f"Available cell types: {cats}."
    )
    # Provider-agnostic tool-use loop: llm.tool_chat normalizes Anthropic tool-use and
    # OpenAI-compatible function-calling to the same {text, tool_calls, assistant} shape, so the same
    # loop drives Claude or an internal vLLM endpoint. (A server that returns no tool_calls just
    # yields a text answer, which we return - a graceful degrade rather than a crash.)
    messages: list[dict] = [{"role": "user", "content": prompt}]
    for _ in range(max_turns):
        turn = llm.tool_chat(system, messages, tools, max_tokens=1024)
        if not turn["tool_calls"]:
            return turn["text"] or (
                "I could not turn that into a concrete analysis - try a narrower question.")
        messages.append(turn["assistant"])
        results = []
        for call in turn["tool_calls"]:
            if on_status is not None:            # stream what we're about to run (see copilot_stream)
                try:
                    on_status(_TOOL_STATUS.get(call["name"], f"running {str(call['name']).replace('_', ' ')} …"))
                except Exception:
                    pass
            res = cap.run(adata, call["name"], call["input"] or {}, ctx)
            # A load_section call replaces the active section: rebind so the rest of THIS turn's tool
            # calls (and later turns) run against the newly loaded adata, not the stale one.
            if ctx.loaded and ctx.loaded[-1]["adata"] is not adata:
                adata = ctx.loaded[-1]["adata"]
            payload = res.error if res.error is not None else res.value
            if res.ok and res.record is not None and action_log is not None:
                action_log.append(res.record)
            results.append({
                "type": "tool_result", "tool_use_id": call["id"],
                "content": json.dumps(payload, default=str),
                "is_error": res.error is not None,
            })
        messages.append({"role": "user", "content": results})
    return "I ran several analyses but could not converge on a concise answer - try a narrower question."
