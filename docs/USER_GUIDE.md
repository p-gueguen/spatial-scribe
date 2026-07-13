# User Guide

SpatialScribe walks you from a raw spatial section to annotated cell types and spatial
insight, in plain language. The app has three panes: the **wizard** (left), the
**interactive canvas** (center), and the **Claude copilot** (right).

## The workflow (left rail)

1. **Load** - click *Xenium 5K (shallow)* for the bundled public 10x Prime 5K breast
   Xenium section (a 100k-cell processed subset, so it opens instantly), or load your own
   Xenium folder / `.h5ad`. *Load synthetic melanoma* gives an offline toy section.

2. **Panel check** - before any annotation, SpatialScribe tells you which cell types your
   panel can and cannot resolve (🟢 resolvable / 🟠 weak / 🔴 cannot resolve), and which
   pairs are indistinguishable. With a Claude key it writes a plain-language verdict. This
   is the honest guardrail: it never offers a label the panel can't support.

3. **QC** - section quality metrics vs the/10x thresholds, plus **region QC**:
   box-select-select any area on the map to get QC just for that region, and exclude bad regions
   (folds, necrosis) from the analysis.

4. **Cluster** - Leiden clustering (GPU-accelerated when available); pick a resolution.

5. **Annotate** - multi-method annotation (marker scoring + optional Claude), reconciled to
   a consensus. Every cell gets a **confidence** and, when the data can't support a call,
   an honest **abstention** (`Unassigned` / `Ambiguous` / `Unresolvable` / `Uncertain` /
   `Novel`) instead of a made-up label. The **annotatability headline** shows what fraction
   is confidently typed vs abstained. Expand **cell-type × cell-state** to see which
   lineages are cycling / interferon-activated / hypoxic. **Click a cell type on the map,
   box-select it, and subcluster** it into subtypes (Claude names them with a key set).

6. **Spatial + niches** - ask *are the T cells excluded from the tumor?* (neighborhood
   enrichment z-score + verdict), see the full neighborhood-enrichment heatmap, and call
   **TME niches** (tumor core / immune-excluded margin / …).

7. **Report** - export the annotated `.h5ad`, and a **re-runnable Python script** built from
   every action you took, so the whole analysis reproduces without the app.

## The copilot (right pane)

With `ANTHROPIC_API_KEY` set, ask questions in plain English. Claude runs the *real*
analysis (neighborhood enrichment, niches, panel summary, composition) and answers grounded
in the computed numbers - it never invents values. Try:

- "Are the T cells excluded from the tumor?"
- "Which cell types can't this panel resolve?"
- "What TME niches are in this section?"
- "Summarize the composition and QC."

## Trust

- Every step shows the code it ran; the report reproduces it.
- The copilot only calls whitelisted analysis tools on your loaded data.
- Your tissue data stays local; only small computed summaries go to the Claude API.
- Marker *presence* on a panel is necessary, not sufficient (probes drop out) - coverage
  numbers are always shown next to any verdict.
