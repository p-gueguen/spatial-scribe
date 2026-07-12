# The SpatialScribe flow and the supervised-vs-clustering gate

```
load -> compute_qc -> panel_check -> reference_match -> annotation_strategy (GATE)
     -> cluster -> { reference_transfer | annotate | clusters-only }
     -> niches -> malignant (tumour-gated) -> split (opt-in) -> self_heal -> qc_funnel -> export
```
(`qc_funnel` - the six-layer section QC + per-cell confidence + section AQI - runs on EVERY route,
not only when annotation produced labels.)

`analysis.pipeline.run_pipeline(adata, ctx, opts)` runs this whole spine as ONE degrade-graceful
function (each stage records an honest `ok | skipped:<why> | failed:<why>` status and the run never
aborts on an optional arm). It is the single engine that BOTH the headless CLI (`scripts/run.py`) and
the app's "Run full analysis" background job drive, so the autonomous run and the interactive app can
never diverge. `self_heal` is the enhanced re-run loop: a failed marker-check subclusters a
heterogeneous type OR abstains a mislabelled one, so a type can leave the failing set (never a silent
confident wrong label).

## Steps

1. **load** - `io.load(path)` auto-detects the platform (Xenium / CosMx / MERSCOPE / Atera) or reads a `.h5ad`, returning a `SpatialSample` whose `.adata` has `obsm['spatial']`, `var['control']`, and `uns['platform']`/`uns['n_panel_genes']` stamped. `demo.load_demo().adata` is the synthetic melanoma section.
2. **compute_qc** - per-cell counts / genes / %-control and a section summary. The count floor is **panel-indexed** (`qc.panel_indexed_floor`), never a fixed `<10`.
3. **panel_check** - which cell types the panel CAN and CANNOT resolve, confusable pairs, and merge groups. Gene presence is necessary, not sufficient.
4. **reference_match** - depth-matched per-type resolvability of a single-cell reference against the panel (`eval_metrics.panel_resolvability`, which supersedes one-vs-rest AUC). Returns overlap, per-type F1, `frac_resolvable`, and a good/fair/poor verdict.
5. **annotation_strategy (the gate)** - self-verifies the reference<->panel match and auto-reruns a remediation ladder, then recommends the annotation mode. Writes `uns['annotation_route']` with `recommended_mode` + the full ladder.

## The gate ladder (why it is a ladder, not a single threshold)

`frac_resolvable` is **granularity-dependent** - a reference with many fine subtypes scores low even when it is an excellent *coarse-lineage* reference. So the gate never routes to clustering on a single low metric. It climbs an ordered ladder, re-scoring at each rung, and stops as soon as the verdict reaches good/fair:

1. **score** the reference as given.
2. **coarsen** - if fair/poor with confusable types, merge the labels the panel cannot separate (from `confused_with`/`confused_frac`) and re-score. This is the highest-yield move: panel size gates label granularity (a K=50 atlas may be 54% resolvable on a 5K panel but 93% at K=29 coarse).
3. **reselect** - if still poor, try the next best-overlap, tissue-matched, *available* reference in the registry (`choose_reference` + `load_reference`) and re-score.
4. **route**:
   - reference resolves the panel (good/fair, possibly coarsened) -> **`reference_transfer`** (supervised label transfer; keep per-cell abstention).
   - no reference supplied -> **`annotate`** (marker-based supervised; `panel_check` governs which types resolve).
   - reference genuinely uninformative - still poor after coarsening + reselecting, or overlap below the classifier floor (<~0.2 / insufficient_overlap) -> **`cluster`** (de-novo Leiden + marker-based naming). Forcing supervised labels here would invent cell types the data cannot support.

## Reference choice matters more than method choice

A wrong-tissue reference can post high overlap + F1 (the resolvability check is reference-INTERNAL cross-validation) yet annotate the section as nonsense. Match the tissue first; prefer single-nucleus references for imaging panels. Pair the internal metric with marker sanity on the actual section.

## Honesty

- Per-cell confidence is a **heuristic** that drives abstention; it is **not calibrated**. Abstention (`Unassigned`/`Ambiguous`/`Unresolvable`/`Uncertain`/`Novel`) is honest - never a confident wrong label.
- Report only figures a capability computed. Never invent numbers or labels.
