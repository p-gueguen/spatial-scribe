---
name: spatialscribe
description: Drive SpatialScribe - a spatial-transcriptomics cell-type annotation copilot - end-to-end on your own imaging-based spatial data, compute, and Anthropic key. Use when annotating a Xenium, CosMx, MERSCOPE, Atera, or .h5ad section; running the panel-adequacy check ("which cell types can my panel resolve?"), the QC funnel, or a reference match; deciding supervised annotation vs unsupervised clustering from the reference-to-panel match; clustering and consensus cell-type annotation with per-cell confidence and abstention; malignant/tumour calling, spatial niches, neighborhood enrichment, or gene programs; self-verifying an annotation; or exporting an annotated .h5ad + a self-contained HTML report + a re-runnable script. Also use to launch the SpatialScribe web app or to drive its grounded Claude copilot in plain language. Triggers on "SpatialScribe", "annotate my Xenium/CosMx section", "panel check", "spatial annotation QC", "reference match", "should I cluster or annotate", "run_copilot".
---

# SpatialScribe

Take a raw imaging-based spatial-transcriptomics section (Xenium / CosMx / MERSCOPE / Atera / `.h5ad`) to QC'd, clustered, annotated cell types with honest per-cell confidence, spatial insight, and a shareable report - in the user's own environment, on their own compute, with their own Anthropic key. **Data never leaves their machine.**

This skill teaches you to **drive** SpatialScribe's engine; it does not reimplement it. All logic lives in the installed `spatialscribe` package (`src/spatialscribe/analysis/`) behind ONE capability registry that both the app and the copilot use. Your job is to orchestrate that registry, pick the run surface, and stay grounded in its computed numbers.

## When to use this skill

- "Annotate this Xenium/CosMx run" / "type the cells" -> the headless flow (`scripts/run.py`).
- "Which cell types can my panel resolve?" / "is my panel adequate?" -> `panel_check` (+ `reference_match`).
- "Should I annotate from this reference or just cluster?" -> `annotation_strategy` (the gate).
- "Where are the T cells?" / "color the tissue by CD8A" / one-off exploration -> the copilot (`scripts/ask.py`).
- "Let me explore this visually" -> the web app (`scripts/launch_app.sh`).

## Install (user's own environment)

```bash
bash scripts/install.sh            # clone + pixi install -e main (or pip -e .) + verify
export ANTHROPIC_API_KEY=sk-...    # never commit it; ANTHROPIC_MODEL=claude-sonnet-5 optional (default Haiku)
```
Runs CPU-only out of the box (`export SPATIALSCRIBE_FORCE_CPU=1` to force it). GPU is optional (install `rapids-singlecell` on a CUDA node; `backend.py` auto-detects, else CPU).

## The one rule: drive the capability registry

`src/spatialscribe/analysis/capabilities.py` is the single source of truth. Never hand-roll scanpy / QC / annotation code here - orchestrate the registry:

```python
from spatialscribe.analysis import capabilities as cap
ctx = cap.RunContext(tissue="human breast", use_llm=True, reference=ref, ref_label_key=key)
res = cap.run(adata, "panel_check", {}, ctx)                 # -> RunResult(ok, value, record, error)
res = cap.ensure(adata, "cluster", ctx, params={"resolution": 0.5})   # compute-if-absent
```

- Every capability declares `requires` / `produces` keys; `run()` checks prerequisites and returns a structured `prerequisite_missing` error whose `hint` names the step to run first. Follow the hint; never guess.
- Collect each `res.record` into an `action_log` list - it becomes the re-runnable export, and copilot-driven steps append to it too.
- Full catalog (names, requires/produces, params, copilot-exposed) is in [references/capabilities.md](references/capabilities.md).

## The flow and its gate

```
load -> compute_qc -> panel_check -> reference_match -> annotation_strategy (GATE)
     -> cluster -> { reference_transfer | annotate | clusters-only }
     -> niches -> malignant (tumour-gated) -> split (opt-in) -> self_heal -> qc_funnel -> export
```

The whole spine runs as ONE degrade-graceful function, `analysis.pipeline.run_pipeline(adata, ctx,
opts)` - each stage records an honest `ok | skipped | failed` status, the run never aborts on an
optional arm, and BOTH `scripts/run.py` and the app's "Run full analysis" background job drive it.

**`annotation_strategy`** is the decision the user asked about most: it self-verifies the reference<->panel match and auto-reruns a remediation ladder (coarsen the labels the panel cannot separate -> reselect a better-overlap tissue-matched reference), then recommends `reference_transfer` (supervised), `annotate` (marker-based, the no-reference default), or `cluster` (de-novo, when the reference genuinely cannot resolve the panel). It **never** routes to clustering on a single low metric. Full ladder + rationale in [references/flow.md](references/flow.md).

Run the FULL pipeline headless (obeys the gate; runs niches + malignant + self-heal too; exports an
annotated `.h5ad` + HTML report + re-runnable script). `--tumour/--no-tumour` sets the malignant gate:
```bash
python scripts/run.py --path <run_dir_or.h5ad> --tissue "human breast" [--reference ref.h5ad] \
                      [--tumour|--no-tumour] [--split] [--rctd] --out results/
python scripts/run.py --demo --out results/          # synthetic melanoma, no data needed
```

## The copilot: run_copilot

For plain-language, exploratory, or "show me" requests, use the grounded copilot instead of the full flow:

```bash
python scripts/ask.py --demo "which cell types can this panel resolve?"
python scripts/ask.py --path <dir> --tissue "human breast" "are the T cells excluded from the tumour?"
```
It is a whitelisted tool-use loop over the capability registry (NOT arbitrary exec): it runs real capabilities and answers only from their computed results; view tools stash plotly figures for optional PNG export. Needs `ANTHROPIC_API_KEY`.

## The interactive app

```bash
bash scripts/launch_app.sh        # React SPA + FastAPI (single-origin) -> http://localhost:8000
```
Right surface when the user wants to lasso regions, click-to-subcluster, or browse - not for scripted batch runs.

## Grounding and honesty (do not violate)

- **Never invent numbers or labels.** Report only figures a capability returned. Gene PRESENCE on a panel is not DETECTABILITY.
- **Per-cell confidence is a heuristic** driving abstention; it is NOT calibrated. Abstention (`Unassigned`/`Ambiguous`/`Unresolvable`/`Uncertain`/`Novel`) is honest, never a confident wrong label.
- **Reference choice matters more than method choice**, and a wrong-tissue reference can pass overlap yet annotate nonsense - match the tissue first. **Panel size gates label granularity** - coarsen for small panels.

## Gotchas (each has cost hours)

- The count floor is **panel-indexed** (`qc.panel_indexed_floor`), never a fixed `<10` (that drops ~60% of a Xenium Prime 5K section).
- Use **`maxRank=150`** for Xenium signature scoring (the 1500 default dilutes to noise).
- Platform + panel size live on **`adata.uns`** (`platform`, `n_panel_genes`) - read them, never re-detect from `X`.
- **`ovrlpy` runs only as a subprocess** in its own env (clashing polars/umap-learn) - never `import ovrlpy` in the main interpreter.
- The learned malignant callers (infercnvpy CNV, Cancer-Finder) are isolated subprocesses configured by env vars (`SPATIALSCRIBE_CNV_*`, `CANCERFINDER_*`) and **degrade to `skipped`** when absent - the skill works fully without them.
- Data + references are injected via env vars with **empty public defaults**: `io.load()` takes any platform dir or `.h5ad`; tissue references resolve from `SPATIALSCRIBE_REF_*` / `SPATIALSCRIBE_REFERENCE_REGISTRY` or an explicit `--reference`. No storage root is hardcoded; there is no database - state is one section + an `action_log`.

## References (read on demand)

- [references/capabilities.md](references/capabilities.md) - the full capability catalog (auto-generated from the registry).
- [references/flow.md](references/flow.md) - the flow spine, the annotation_strategy gate ladder, reference-choice-matters.
