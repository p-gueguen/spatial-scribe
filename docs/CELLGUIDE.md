# CellGuide markers — what it is, where it's used, and how it could skill up an LLM

`analysis/cellguide.py` grounds cell-type labels in the **Cell Ontology** and serves markers from the
CZI **CellGuide** corpus. It complements the curated `markers.LINEAGE_MARKERS` (which stay the
benchmarked default) and replaces "ask the LLM to invent marker genes" with a deterministic,
corpus-backed lookup. Companion to [ANNOTATION_SCHEME.md](ANNOTATION_SCHEME.md).

## The data source

Undocumented static JSON, versioned by snapshot. Treat it as **best-effort, never load-bearing** (the
provider degrades to curated markers on any failure).

| path | content |
|---|---|
| `/latest_snapshot_identifier` | the current snapshot id (e.g. `1764612212`); everything else lives under `/<snapshot>/` |
| `/<snap>/celltype_metadata.json` | 1186 CL types → `{name, id, clDescription, synonyms}` — the label→CL resolver |
| `/<snap>/canonical_marker_genes/CL_xxxx.json` | literature-curated `{tissue, symbol, ...}`. **Human symbols, no organism field.** |
| `/<snap>/computational_marker_genes/CL_xxxx.json` | census-derived `{marker_score, specificity, symbol, groupby_dims{organism, tissue}}`. `marker_score` = the UI **"Effect Size"**. |

Two traps baked into the code (both cost a review cycle): rank the computational markers on the
**tissue-agnostic aggregate only** (mixing tissue-specific records surfaces census non-markers like
SKAP1 over CD3D), and re-rank the **canonical** set (which has no intrinsic order) by that same effect
size so flagship markers lead (fibroblast → DCN/COL1A2/PDGFRA, not an alphabetical head). The bare
host returns `200` with an empty body for wrong paths, and responses are gzip that `urllib` will not
auto-decompress — so a 200 does not mean the path is valid.

## The API

```python
from spatialscribe.analysis import cellguide as cg
cg.resolve_cl("Stromal/CAF")                 # -> "CL:0000057" (fibroblast); None if no ontology match
cg.markers_for_label("T cell", organism="Homo sapiens")
#   -> {"cl_id": "CL:0000084", "source": "canonical", "markers": ["CD3D","IL7R","CD3G",...]}
cg.markers_for_label("T cell", organism="mouse")     # -> source "computational", ["Cd3g","Cd3d",...]
```

- `resolve_cl` matches exact (case/punctuation-insensitive) names + synonyms, composite aliases
  (`Stromal/CAF`→fibroblast), and a `"<label> cell"` suffix so bare lineages ("Mast", "Myeloid")
  resolve. Never fuzzy-snaps — a wrong ontology anchor is worse than none.
- `markers_for_label`: **canonical** (human, ≥3 symbols, effect-ordered) → else **computational**
  (effect-ordered; mouse uses this with mouse-case symbols) → `source:"none"`.
- Offline-safe: snapshot-keyed disk cache under `/data/.../cache/cellguide`; `<cache>/.latest` lets
  an outage reuse the last snapshot; the ~1.3 MB computational file is distilled to the top-20
  `[symbol, score]` per organism before caching.

## Where it is used today (and where it is not)

CellGuide is wired into `markers.markers_for_types`, and reached from **two Panel-check endpoints**:
`GET /api/{sid}/panel_report` and `GET /api/{sid}/panel_verdict` (via
`_apply_llm_markers_for_categories`). That helper grounds **only the NOVEL labels** — those not in the
curated tissue dictionary — and **merges** them onto the curated baseline; it never replaces curated
markers for known lineages (which are benchmarked at 8/10 vs 6/10 typable).

**So on the stock demo it does not fire** — every demo label is already a curated breast lineage, so
`novel` is empty and the panel endpoints return before touching CellGuide. It activates once the user
introduces a label the curated dicts never anticipated: **click-to-subcluster** (e.g. "CD8 T"),
**rename a cell type**, or **adopt reference-transferred labels** (CellTypist / RCTD). Then those
labels get ontology-grounded markers instead of LLM-invented ones.

The unique thing CellGuide adds — a real `CL:xxxx` id + description per label — is currently consumed
only internally (for marker selection); it is **not yet surfaced in the UI**.

## Extension: CellGuide as an LLM tool ("skill up" a weaker model)

The most interesting use is turning the lookup into a **tool the copilot LLM calls directly**, rather
than a function only the panel-check path invokes. The copilot (`agent.run_copilot`) is a tool-use
loop over `capabilities.copilot_tools()` — every `Capability` with `copilot_exposed=True` becomes an
Anthropic/OpenAI tool schema the model can call. Adding one CellGuide capability would let the model
*look up* grounded markers and the CL id instead of relying on its parametric knowledge:

```python
# sketch — an additive capability, no regression surface
def _cap_lookup_markers(adata, ctx, cell_type, organism="human", **_):
    from . import cellguide
    return cellguide.markers_for_label(cell_type, organism=organism)  # {cl_id, source, markers}

Capability(
    name="lookup_markers", copilot_exposed=True,
    label="Look up canonical markers for a cell type",
    description="Cell-Ontology-grounded marker genes (CellGuide): canonical > computational, "
                "human or mouse. Use this instead of recalling markers from memory.",
    params={"type": "object", "properties": {
        "cell_type": {"type": "string"}, "organism": {"enum": ["human", "mouse"]}},
        "required": ["cell_type"]},
    fn=_cap_lookup_markers,
)
```

**Why this helps, concretely.** The copilot defaults to Claude, but it is endpoint-agnostic and can
run against a smaller local model (any OpenAI-compatible `/v1` server via `SPATIALSCRIBE_LLM_BASE_URL`).
A smaller model has weaker biological recall, so a deterministic marker/ontology lookup is
RAG-for-annotation: instead of guessing "what are T-cell markers?" from memory it **calls the list**, so
its answers are grounded in the CZI corpus and carry a CL id. That is the "skill up" — the tool
supplies the knowledge the small model lacks, and the grounding is auditable (canonical vs
computational, with a snapshot id) rather than a hallucination.

**Honest caveats before building it:**
- Some local-llm vLLM builds intermittently emit **no `tool_calls`** and degrade to a text answer (seen on
  the live deploy; see the grounding-eval note in progress). A lookup tool only helps a backend that
  reliably tool-calls; verify the deployed local-llm build does before relying on it.
- Panel-check now **augments** curated with CellGuide (union, curated first), it does not replace it.
  `_apply_llm_markers_for_categories` unions CellGuide's canonical markers onto the curated baseline
  for EVERY annotated category, so under-covered curated lineages (Mast 4, NK 4, Endothelial 4 curated
  markers) gain the literature markers the curated dicts never named and more land on the panel ->
  higher identifiability / more typable types. This is a deliberate reversal of the earlier
  "lookup-only" stance (requested after users saw good panels read "not typable"): safe because
  CellGuide is a **static, versioned** DB (deterministic across runs, unlike an LLM) and curated
  markers still LEAD the list (their benchmarked flagship genes rank first, never discarded); only a
  label CellGuide cannot resolve falls through to the LLM. It does shift the panel-adequacy numbers vs
  the curated-only baseline - that is the intent (fill the missing markers), not a regression.
- It is additive and offline-degrading by construction (the provider returns `source:"none"` when the
  endpoint is unreachable), so it never makes the copilot worse than today.

Not built yet — this section documents the design so it can be picked up deliberately.
