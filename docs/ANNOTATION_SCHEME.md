# The cell-type annotation scheme

What the app assigns, who decides it, and what the evidence says it *should* be.
Companion to [annotation-method-selection.md](research/annotation-method-selection.md) (which tool is
best) and [confidence-calibration.md](research/confidence-calibration.md) (what the confidence score
is worth).

## What a `cell_type` is

The primary label lives in `obs['cell_type']` and is assigned **per leiden cluster, not per cell**
(`annotate.consensus_annotate`). Two signals are computed:

| signal | how | label space |
|---|---|---|
| `marker_label` | argmax of per-cluster mean marker scores | a key of the curated per-tissue dict |
| `claude_label` | an LLM naming the cluster from its top markers | the **same** closed vocabulary |

The curated dict is chosen by tissue: `LINEAGE_MARKERS` (skin/melanoma, 10 lineages),
`EPITHELIAL_LINEAGES` (breast / solid tumour, 10), `BRAIN_MARKERS` (9). Granularity is deliberately
**coarse and hierarchical**: `cell_type` stays at the lineage level, while `subtype` (de-novo
subclustering), `cell_state` (orthogonal programs) and `program` (NMF) live in separate, opt-in
columns. There is **no Cell Ontology / CL vocabulary**, on purpose - `consensus.py` uses exact-label
matching, which a ~10-label space does not need an ontology graph to do.

## Three invariants (each fixed a real bug; each pinned by a test)

### 1. The LLM names; it does not invent, and it does not overrule

The allowed lineages plus a `Novel / unknown` escape hatch are constrained **in the prompt**
(`llm.annotate_clusters(allowed_labels=...)`), and the reply is validated on the way back by exact,
case/whitespace-tolerant matching. Outcomes:

| the model says | result |
|---|---|
| the marker argmax | used |
| `Novel / unknown` | used - an honest "no lineage on this panel fits" |
| something off-vocabulary (`pDC`) | discarded, marker argmax used, `off_vocabulary=True` |
| a **different valid lineage** | `conflict=True`: marker argmax kept, `obs['label_conflict']` set, cells capped at WARN |

There is **no fuzzy/nearest-string snapping**, deliberately. `pDC` and `RBC` are near nothing in a
ten-lineage space and would snap to a *wrong* lineage, converting an honest "I don't know" into a
confident error.

Why the marker argmax wins a genuine conflict: on the 2026-07 five-section benchmark (two different
LLMs naming the *same* leiden clusters) **44 of 59** cluster-label disagreements were naming-only -
which the closed vocabulary now deletes outright - and the remaining **15 were real biology that
nothing adjudicates**. Agreement with the CellTypist proxy (0.457 / 0.510) sits at or below its own
**~0.48 majority-class baseline**, so that proxy cannot say who is right. Absent evidence that the
LLM beats the argmax, a conflict is *flagged*, not silently resolved. Measured live on Prime 5K:
Qwen named 2 of 9 clusters `Myeloid`/`Endothelial` where the markers said `T cell`.

Closing the vocabulary makes a **missing lineage a real cost**: audit the tissue dict first. Adding a
lineage is free on panels that cannot see it - `markers.present()` drops any lineage whose markers are
all off-panel.

### 2. An abstention is not a cell type

`cell_type_final` overlays six pseudo-labels (`Unassigned: low quality`, `Ambiguous: mixed`,
`Unresolvable: panel`, `Uncertain: low confidence`, `Ambiguous: label conflict`, `Novel / unknown`)
into the *same* `pd.Categorical` as real lineages. `annotate` owns the one vocabulary
(`ABSTENTION_LABELS`, matched **exactly, never by prefix** - `Novel epithelial subtype X` is a real
label) and the one grouping column (`annotation_key`).

| consumer | treatment |
|---|---|
| composition (report **and** interactive chart), state heatmap, niche **names** | abstained cells **excluded**, excluded share stated |
| spatial map, neighbour graph, niche **features** | abstained cells **kept**, collapsed into one grey `Not assigned` |

Abstained cells stay in the neighbour graph because they occupy real tissue: deleting them would
silently rewire every neighbourhood. They therefore remain in the permutation null, which
`nhood_enrichment` discloses as `pct_abstained_in_graph`.

**Two abstained numbers, both correct, different questions:**

- `composition_table`'s fraction = every cell whose final label is not a lineage (includes the honest
  `Novel / unknown` and the leiden-NaN `Unassigned`). **15.5%** on a real Prime 5K section.
- `apply_confidence`'s `pct_abstain` = the confidence-verdict FAIL rate. **10.8%** on the same section.

### 3. Scope

The closed vocabulary constrains the **marker + LLM path only**. `consensus_annotate(method_label_cols=…)`
(the RCTD / SingleR / scANVI / panhumanpy consensus vote) overwrites `cell_type` with the winning
*reference* label, which carries the reference's own controlled vocabulary - usually finer than these
~10 lineages, and deliberately not coerced into them.

## What *should* be used

**Coarse lineage is the floor, not the ceiling, and depth chooses the ceiling.** Panel size and
per-cell depth gate label granularity: a K=29 atlas is 93% resolvable on a 5K panel while K=50 is only
54%; Atera WTA at ~1131 counts/cell resolves 10/10 lineages and 98% of a K=50 scheme. Measure with
**depth-matched per-class F1** (`eval_metrics.panel_resolvability`), never one-vs-rest AUC, which
inflates at full depth and deflates to ~0.5 on shallow sections.

**When a tissue-matched reference exists, route to reference transfer.** Reference choice dominates
method choice (see [annotation-method-selection.md](research/annotation-method-selection.md)): RCTD
doublet-mode with a matched, ideally single-nucleus reference, SingleR as a fast cross-check. Do not
ship a zero-shot foundation model. The routing ladder already exists -
`reference.plan_annotation_strategy` scores the reference-panel match, **coarsens** confusable types,
tries a better reference, then recommends supervised transfer vs de-novo clustering. The default path
only lands on marker + LLM because no reference has been loaded.

**Do not** adopt a Cell Ontology for a ten-label space; **do not** wire RCTD/SingleR/scANVI "to improve
accuracy" measured against the CellTypist proxy, which scores at or below its own majority-class
baseline and therefore cannot adjudicate anything; and **do not** cap the primary column at coarse
lineage on deep panels.

## Where the code lives

| concern | symbol |
|---|---|
| abstention vocabulary + grouping column | `annotate.ABSTENTION_LABELS`, `is_abstention`, `annotation_key`, `collapse_abstention`, `typed_mask` |
| closed vocabulary | `llm.annotate_clusters(allowed_labels=…)`, `annotate._canonicalise` |
| conflict flag | `obs['label_conflict']`, `apply_confidence` -> `pct_label_conflict` |
| composition denominator | `export.composition_table` (used by the PNG *and* `views.composition_view`) |
| tests | `tests/test_annotation_vocabulary.py`, `tests/test_abstention_grouping.py` |
