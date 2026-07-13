# QC Metrics for Cell-Type Annotation in Imaging-Based Spatial Transcriptomics

**Status:** design reference / evidence baseline · **Scope:** Xenium-first (MERFISH/Vizgen and CosMx/NanoString as secondary comparison) · **Audience:** SpatialScribe developers building the annotation pass/warn/fail criteria (`analysis/qc.py`, `analysis/panel_check.py`, and the future `analysis/annotate.py`).

> **Why this document exists.** SpatialScribe annotates cells for wet-lab scientists who cannot sanity-check the call themselves. That raises the bar: every label the app shows must carry a **confidence** and, where the data cannot support a call, an **honest warning** rather than a confident guess. Imaging-based spatial transcriptomics (iST) makes this unusually hard because the measurement is **sparse** (a targeted panel, not the whole transcriptome; 100–5,000 genes), **noisy** (segmentation errors, spatial doublets, transcript diffusion / spillover), and **has no ground truth**. This document is the evidence base from which we set the concrete pass/warn/fail thresholds. It catalogs the metric families, states what each one actually measures, gives Xenium-first evidence and typical value ranges, and proposes **provisional** cutoffs clearly flagged as starting points to tune, not settled truth.

**Companion artifacts:** [`annotation_qc_thresholds.yaml`](annotation_qc_thresholds.yaml) — the machine-readable version of every draft threshold, ready to load into the QC code · [`research-evidence-appendix.md`](research-evidence-appendix.md) — the verified findings + full source list from the adversarial research pass that grounds this doc.

---

## Table of contents

1. [The problem: why annotation QC is hard in iST](#1-the-problem)
2. [The layered QC model (six gates)](#2-the-layered-qc-model)
3. [Layer 0 — Run / section acquisition QC](#layer-0)
4. [Layer 1 — Segmentation QC](#layer-1)
5. [Layer 2 — Per-cell count & expression QC](#layer-2)
6. [Layer 3 — Contamination / spillover / signal-integrity QC](#layer-3)
7. [Layer 4 — Panel adequacy (resolvability)](#layer-4)
8. [Layer 5 — Annotation-confidence QC](#layer-5)
9. [Layer 6 — Spatial-context QC of labels](#layer-6)
10. [Confidence scores, abstention & the "cannot type this cell" warning](#7-confidence-scores)
11. [Solutions for sparsity, spatial doublets & transcript diffusion](#8-solutions)
12. [Consolidated draft pass / warn / fail criteria](#9-consolidated-criteria)
13. [How this maps to SpatialScribe code](#10-code-mapping)
14. [Cross-platform QC comparison (secondary: MERSCOPE, CosMx)](#14-cross-platform)
15. [Caveats, open problems & anti-patterns](#11-caveats)
16. [Glossary](#12-glossary)
17. [References](#13-references)

---

<a name="1-the-problem"></a>
## 1. The problem: why annotation QC is hard in iST

Sequencing-based scRNA-seq annotates whole-transcriptome profiles of physically dissociated cells. Imaging-based spatial transcriptomics (Xenium, MERFISH/MERSCOPE, CosMx) is different on three axes that each degrade annotation and each needs its own QC:

| Axis | scRNA-seq | Imaging spatial (iST) | Consequence for annotation |
|---|---|---|---|
| **Depth** | whole transcriptome, 1,000s–10,000s genes/cell | **targeted panel** (100–5,000 genes); typically **10–100s** genes and **50–500 transcripts** per cell | Many cell types have **few or no discriminating markers on-panel**; scores computed over a handful of counts are high-variance. |
| **Cell boundaries** | physical dissociation (1 cell = 1 barcode, plus ambient RNA / doublet artifacts) | **in-situ segmentation** of an image; nucleus ± expansion or membrane stain | Segmentation errors create **empty cells**, **merged doublets**, and **truncated** cells; the "cell" is a hypothesis, not a fact. |
| **Signal purity** | ambient RNA in droplet | **transcript diffusion / optical spillover** between adjacent cells; **vertical (z) overlap** in a 2D projection | A cell's counts are a **mixture** of its own transcripts and its neighbors' → false co-expression → mis-annotation toward the locally dominant type. |

**The core failure mode this document guards against:** a sparse, contaminated count vector is confidently assigned to *some* cell type because every classifier returns *a* label. Without QC, "T cell" and "artifact of tumor spillover onto a stromal fragment" look identical in the output. The job of annotation QC is to (a) attach a **calibrated confidence** to each call and (b) **abstain** — return `Unassigned` / `Ambiguous` / `Unresolvable` — when the evidence cannot support any call, with a reason the scientist can act on.

**No ground truth.** There is no gold-standard label for a Xenium cell. Every metric here is therefore either (i) an **internal-consistency** check (does this label cohere with markers, neighbors, and the reference?) or (ii) a **contamination / artifact** check. Treat all thresholds as **distributional heuristics**, per 10x's own guidance that "there are no universal thresholds … no single metric determines success or failure." QC here is a **triage funnel**, not a pass/fail oracle.

---

<a name="2-the-layered-qc-model"></a>
## 2. The layered QC model (six gates)

Annotation quality is the *last* thing you can measure but the *last* thing to fail — most bad annotations are caused by upstream artifacts. We therefore model QC as an ordered funnel where each layer both (a) filters/flags cells and (b) feeds a signal into the final per-cell confidence. A cell that survives all gates with strong signals gets a confident label; a cell that trips a gate gets flagged, down-weighted, or abstained.

```
 ┌─────────────────────────────────────────────────────────────────────┐
 │ LAYER 0  Run / section acquisition QC        (is the SECTION usable?)│  ── decode Q, assigned frac, neg-control rate
 ├─────────────────────────────────────────────────────────────────────┤
 │ LAYER 1  Segmentation QC                      (is this a real CELL?) │  ── cell area, nucleus present, empty-cell frac
 ├─────────────────────────────────────────────────────────────────────┤
 │ LAYER 2  Per-cell count & expression QC       (enough SIGNAL?)       │  ── counts/cell, genes/cell, effective depth
 ├─────────────────────────────────────────────────────────────────────┤
 │ LAYER 3  Contamination / spillover / integrity(is the signal PURE?)  │  ── neg-control %, VSI, doublet score, MECR/CRISP
 ├─────────────────────────────────────────────────────────────────────┤
 │ LAYER 4  Panel adequacy (resolvability)       (CAN this type exist?) │  ── marker coverage, confusable pairs (per type)
 ├─────────────────────────────────────────────────────────────────────┤
 │ LAYER 5  Annotation-confidence QC             (is the LABEL trusted?)│  ── posterior prob, margin, entropy, consensus
 ├─────────────────────────────────────────────────────────────────────┤
 │ LAYER 6  Spatial-context QC of the label      (does it FIT context?) │  ── spatial coherence / PAS, neighborhood enrichment
 └─────────────────────────────────────────────────────────────────────┘
                                │
                                ▼
                 Per-cell verdict:  PASS  ·  WARN  ·  FAIL(abstain)
                 + confidence [0,1] + machine-readable reason code
```

Layers 0–3 are **cell-quality** gates (largely label-independent). Layer 4 is a **panel-capability** gate (dataset-level, per cell type). Layers 5–6 are **label-quality** gates (per annotation arm). The final verdict combines them (see [§7](#7-confidence-scores) and [§9](#9-consolidated-criteria)).

**Cost tiers** (what data each layer needs — matters for the app, which starts from a cell×gene AnnData + coords):

- **Cheap (cell×gene + labels + coords only):** empty-cell fraction, counts/genes per cell, neg-control %, MECR, CRISP purity, marker-score margin, spatial coherence / PAS, consensus agreement.
- **Needs the reference kept:** RCTD/deconvolution weights & spot_class, label-transfer posterior & margin, pseudobulk-reference correlation.
- **Needs molecule-level data (`transcripts.parquet`) or the z-stack:** ovrlpy VSI, NPMI molecular coherence, Baysor/proseg re-segmentation, per-probe Q-score audits.

---

<a name="layer-0"></a>
## 3. Layer 0 — Run / section acquisition QC

This is the section-level gate that decides whether annotation should be attempted **at all**. It is already implemented for the run/section level in `XENIUM_QC_THRESHOLDS.md` and mirrored in [`qc.py`](../../src/spatialscribe/analysis/qc.py) (`XENIUM_THRESHOLDS`). Reproduced here for completeness and to anchor the funnel; these are 10x-published cutoffs, not our invention.

| Metric | Direction | Warn | Fail/Error | What it means for annotation |
|---|---|---|---|---|
| `fraction_transcripts_decoded_q20` | lower worse | < 0.70 | < 0.50 | Poor decoding → the counts themselves are unreliable; annotate with caution flag. **Custom-panel-only variant** (`fraction_custom_transcripts_decoded_q20`, WARN ≤0.60/ERR ≤0.50) is computed over *custom add-on probes only* — a low value flags specific failing probes, not the whole sample (see [caveat](#11-caveats)). |
| `fraction_transcripts_assigned` | lower worse | < 0.50 | < 0.30 | Few transcripts land in cells → segmentation or diffusion problem; expect noisy annotation. |
| `fraction_empty_cells` | higher worse | > 0.15 | > 0.25 | Over-segmentation; empty "cells" will be forced into a label unless gated at Layer 1. |
| `median_genes_per_cell` | lower worse | < 20 | < 10 | Sensitivity floor; below this, per-cell annotation is unreliable and should lean on neighborhood/deconvolution. |
| `median_transcripts_per_cell` | lower worse | < 50 | < 25 | As above. |
| `negative_control_probe_counts_per_control_per_cell` | higher worse | > 0.02 | > 0.05 | Non-specific binding / background floor; high values inflate false co-expression. |

**Draft rule.** If a section fails ≥1 Layer-0 metric at the **error** level, annotate but stamp the whole section `low-confidence` and surface the failing metric prominently. Two or more error-level fails → recommend the scientist not trust cell-level calls (aggregate/region-level only).

**Cross-platform note.** MERSCOPE and CosMx do not expose identical section metrics. Analogues: MERSCOPE reports per-FOV transcript density and a **blank/negative-barcode misidentification rate**; CosMx (AtoMx) reports **SystemControl / NegPrb** counts and per-FOV QC. Normalize all three to the same three signals — decode quality, assignment fraction, negative-control rate — before applying the funnel.

---

<a name="layer-1"></a>
## 4. Layer 1 — Segmentation QC

**Why it comes before everything label-related.** The single largest driver of mis-annotation in iST is bad cell boundaries: a merged doublet reads as a co-expressing "hybrid," a truncated cell loses its markers, and an empty fragment gets a spurious label. Segmentation QC decides *whether the object being annotated is a plausible single cell.*

### 4.1 What Xenium (and others) actually do

- **Xenium** — early Xenium Onboard Analysis (XOA <2.0) used **DAPI nucleus segmentation + a fixed ~15 µm expansion**. XOA ≥2.0 uses **multi-modal segmentation** driven by interior/boundary stains (nucleus + cell-boundary + interior RNA morphology), which reduces the expansion artifact. Know which one produced your data — expansion-based objects have systematically more spillover.
- **MERSCOPE** — Cellpose-based on DAPI + poly-T / cell-boundary stains, per-z-plane then 3D-stitched.
- **CosMx** — nuclear (Morphology markers) + membrane-marker-guided segmentation.
- **Third-party re-segmentation** (often materially better): **Baysor** (molecule-density-aware), **proseg** (probabilistic, fast), **BIDCell** (self-supervised deep), **Cellpose/Mesmer**, **ComSeg**, **segger**, and orchestration via **sopa**. If the app ever re-segments, its QC metrics change but the *families* below stay the same.

### 4.2 Segmentation QC metrics

| Metric | What it measures | Typical Xenium range | Draft warn / fail (per cell unless noted) |
|---|---|---|---|
| **Cell area** (µm²) | Physical plausibility of the segment | median ≈ 30–200 µm²; long right tail | **WARN** < 8 or > 400 µm²; **FAIL** < 4 (fragment) or implausibly large merged blob (dataset-relative: > 99.5th pct). |
| **Nucleus:cell area ratio** | Over-expansion (ratio→0) vs no cytoplasm (ratio→1) | ~0.2–0.6 | Flag extremes; ratio ≈ 1 with high counts can be a merged nucleus doublet. |
| **Nucleus present** (bool) | Cell has a detected nucleus vs cytoplasm-only fragment | most cells True | Cytoplasm-only + low counts → strong `Unassigned` candidate. |
| **`fraction_empty_cells`** (section) | Over-segmentation rate | < 0.10 good | as Layer 0 (WARN 0.15 / FAIL 0.25). |
| **Transcript-in-cell fraction** (section) | Under-segmentation / diffusion | > 0.50 good | as Layer 0 `fraction_transcripts_assigned`. |
| **Cell density / cells-per-area** | Regional over/under-segmentation | tissue-dependent | Flag FOVs/regions in the extreme tails vs the section. |
| **Segmentation-method disagreement** (optional) | Re-segment with a 2nd method; fraction of cells whose label/counts change materially | — | High disagreement region → down-weight annotation confidence there. |

**Draft rule.** Cells with area below the fragment floor **or** empty (counts < 1) are **FAIL → not annotated** (label `Unassigned:low-quality-segment`). Cells in the implausibly-large tail with multi-lineage marker co-expression are candidate **doublets** → Layer 3.

### 4.3 Peer-reviewed best practice & the expansion problem (Salas et al., *Nat Methods* 2025)

- **Recommended segmentation:** Cellpose (nuclei model) **+ Baysor read-based refinement**; fall back to segmentation-free methods (SSAM, Points2Regions) when standard segmentation fails.
- **Validate segmentation with `NMP` (negative marker purity) + reads-assigned**, and compare competing segmentations with **ARI** and median / 5th-percentile reads-per-cell (see [§6.3](#63-dataset-level-co-expression-contamination-metrics-cheap-no-reference-needed)).
- **The expansion problem, quantified.** In their reference dataset, nuclei have a mean radius of **~5.06 µm**, the *ideal* cell expansion is only **~5.64 µm** (cell-type dependent), and transcripts **more than ~10.71 µm from the cell centroid** correlate with **domain-specific background** rather than the cell's own signal. This is direct, cited evidence that **fixed large expansions pull in neighbor/background transcripts** — the mechanistic root of spillover-driven mis-annotation. **10x itself acted on this: XOA reduced the default nucleus expansion from ~15 µm (v1.x) to ~5 µm (v2.0+)** — so "default expansion" is *version-dependent*; check which XOA version produced your data. If your data used a large fixed expansion, treat peripheral transcripts with suspicion and lean harder on the contamination gates (Layer 3).

### 4.4 Segmentation error is a first-order axis — but classification survives it better than state analysis (Mitchel et al., *Nat Genet* 2025)

A reassuring-but-important nuance from the Baysor authors' segmentation-error study ("Impact and correction of segmentation errors in spatial transcriptomics"):

- Segmentation-driven transcript **misassignment can dominate downstream results** and **confounds most single-cell analyses of cellular *state*** — differential expression, neighbor-influence, and ligand–receptor inference are **highly sensitive** to it.
- **Cell-type *classification* itself is comparatively robust** to segmentation error. **Practical implication for SpatialScribe:** the *lineage-level annotation* the app shows is the most defensible output; **cell-*state* / subtype / DE / spatial-interaction claims built on the same segments deserve stronger warnings.** Bake this asymmetry into confidence: coarse labels can be `PASS` where a fine state call on the same cell is only `WARN`.
- **cellAdmix** (matrix factorization of local molecular neighborhoods) identifies and isolates molecular admixtures — a contamination correction analogous to scRNA-seq doublet filtering — and is a candidate cleanup step upstream of annotation.

---

<a name="layer-2"></a>
## 5. Layer 2 — Per-cell count & expression QC

The sensitivity gate: does this cell carry **enough signal** to support any annotation? This is where iST sparsity bites hardest.

| Metric | What it measures | Xenium-first evidence / range | Draft warn / fail |
|---|---|---|---|
| **`total_counts`** (transcripts/cell) | Effective library size | median 50–500 depending on panel/tissue; T cells & neutrophils sit low | **FAIL** < 10 counts (repo `minCounts=10`, keeps low-RNA immune cells that spillover hits hardest); **WARN** < 25. |
| **`n_genes_by_counts`** (genes/cell) | Feature breadth | median 20–150 | **FAIL** < 5 genes; **WARN** < 10. |
| **Effective depth vs panel size** | Counts relative to what the panel *could* yield | ratio to section median | Cell in bottom decile of its section → down-weight. |
| **`pct_counts_control`** (per cell) | Fraction of a cell's counts that are control probes | should be ≈ 0 | **WARN** > 2%; **FAIL** > 5% (repo `max_pct_control=5`). |
| **Max single-gene fraction** | One gene dominating the cell (e.g. a 16S/rRNA or a spillover MALAT1) | — | Flag cells where 1 gene > 50% of counts (often artifact or extreme specialization). |

**Peer-reviewed anchor (Salas et al., *Nat Methods* 2025).** Across a large Xenium panel benchmark they report a mean of **186.6 reads/cell**, **76.8% of reads assigned to cells**, and — critically — they **exclude cells with < 10 assigned reads**, which removed only **0.21%** of cells. This independently validates the repo's `minCounts = 10` floor as a low, biology-preserving cutoff rather than an aggressive filter.

> **⚠ The count floor is NOT transferable across panels — this bites SpatialScribe's own demo dataset.** The very same `< 10` floor that removes **0.21%** of cells on a *targeted* Xenium panel removes **~60% of cells on the Xenium Prime 5K panel** (Marconato/SPLIT team, *Nat Methods* 2026), because a 5,000-plex panel spreads the same capture budget across ~10× more genes → far fewer counts per gene per cell. SpatialScribe's public demo **is** the Prime 5K skin melanoma dataset, so a naive `minCounts=10` would either nuke most cells or (if lowered) admit near-empty ones. **Rule: index the count floor to panel size** — a **panel-size-indexed threshold schedule** (small targeted → Prime 5K → whole-transcriptome), or a **distributional floor** (e.g. section-relative percentile), not one fixed number. This is the clearest example of why every threshold here ships as tunable config.

**Why not just raise `minCounts`?** Because real immune cells (T, NK, neutrophils, pDC) are genuinely low-RNA and are exactly the cells spillover corrupts. A high count floor silently deletes the biology you care about. Keep the floor low (10) but **route low-count cells through stronger confidence gating** (Layer 5) and **hierarchical/coarse-only annotation** (see [§8](#8-solutions)), rather than excluding them.

**Sparsity-aware scoring gotcha (Xenium-specific).** Signature/marker scoring (`UCell`, `AddModuleScore`, `score_genes`) must use a **small `maxRank` (≈150)** on Xenium — the scRNA-seq default (1500) ranks past the ~<500 detected genes into zero-count noise and dilutes every score toward zero. This is a *QC-relevant* setting: a mis-set `maxRank` makes every marker score look weak and every cell look "un-typeable" for the wrong reason.

---

<a name="layer-3"></a>
## 6. Layer 3 — Contamination / spillover / signal-integrity QC

The purity gate. A cell can pass Layers 0–2 (real segment, plenty of counts) and still be **un-annotatable** because its counts are a *mixture* of itself and its neighbors. This layer is the heart of iST-specific QC and the direct cause of the "false co-expression → mis-annotation" failure mode.

### 6.1 Background / negative-control contamination

- **Negative-control probe rate** (per cell and per section): 10x's `negative_control_probe_counts_per_control_per_cell` (WARN 0.02 / FAIL 0.05). Exact 10x definition: **`(Q20+ negative-control-probe counts) / (# negative control probes) / (# cells)`** — a standardized per-control-per-cell rate; a **high** value flags nonspecific/off-target probe binding, poor assay conditions, or low transcript content. 10x also derives a per-cell **"estimated false positives" = (NCP counts per control per cell) × (# target genes)** — a directly interpretable "how many of this cell's counts are probably spurious" number. MERFISH uses **blank barcodes** → a misidentification / false-discovery rate; CosMx uses **NegPrb / SystemControl** counts. All estimate the **non-specific background floor**. A cell whose real signal is not clearly above this floor cannot be confidently typed.
- **Per-cell control burden** feeds directly into confidence: high control % → lower confidence, candidate `Unassigned:noisy`.

> **Terminology trap (do not conflate).** "Negative **control probe** rate" (above) is 10x's *background-binding* metric. It is a completely different thing from **NCP = negative *co-expression* purity** ([§6.3](#63-dataset-level-co-expression-contamination-metrics-cheap-no-reference-needed)), a *reference-anchored specificity* metric from Salas et al. 2025. Both abbreviate awkwardly; keep them separate in code and prose. (Our own research pass flagged three sources that conflated them.)

### 6.2 Spatial doublets & vertical overlap (segmentation-adjacent contamination)

- **ovrlpy Vertical Signal Integrity (VSI)** — for imaging platforms where transcripts have a z-coordinate but cells are segmented in a 2D projection, ovrlpy models the vertical coherence of the transcript signal and returns a per-pixel / per-cell **integrity score**; low integrity flags **vertically overlapping cells, tissue folds, and doublets** that 2D segmentation cannot separate. It needs molecule-level data with z (`transcripts.parquet`) and runs as a subprocess (see the repo's `subprocesses/ovrlpy/`). This is SpatialScribe's designated **confident-cell** filter (H5).
- **3D coherence score** (Salas et al., *Nat Methods* 2025) — the peer-reviewed formalization of the same idea: the **cosine similarity between the latent-space vector fields of the top vs bottom half** of the section along z. Regions with **coherence < 0.2** indicate potentially overlapping cells (mixed-source signal); in their reference Xenium dataset **~1.8% of cells** were flagged this way. Use it (or VSI) as the vertical-doublet gate; the fraction flagged is itself a section-level QC number.
- **RCTD doublet mode** (`spot_class` ∈ {`singlet`, `doublet_certain`, `doublet_uncertain`, `reject`}) — when a reference is available, RCTD explicitly models each cell as up to two cell types and reports whether a **single** type explains it. `doublet_*`/`reject` cells are, by construction, low-confidence annotations. This is the single most useful reference-based purity signal for iST.
  - **Reframing for imaging data (Marconato/SPLIT team, *Nat Methods* 2026):** because imaging ST has **no droplet co-encapsulation**, an RCTD "doublet" is almost never two real cells in one barcode — it is **spillover contamination** from neighbors. Contamination can be strong enough that **RCTD swaps the primary and the contaminating cell type**, i.e. mis-labels the cell as its neighbor. So treat **RCTD secondary weights and the singlet/doublet call as *contamination / annotation-confidence* flags**, not as biological doublet biology. **SPLIT** exploits exactly these weights (segmentation-agnostic, post-hoc) to purify the cell and recover its true identity.

### 6.3 Dataset-level co-expression contamination metrics (cheap, no reference needed)

These quantify how often *mutually-exclusive* lineage markers appear in the same cell — a direct readout of spillover:

- **MECR — Mutually Exclusive Co-expression Rate** (Hartman & Satija, 2024): for a pair of markers from **disjoint lineages**, MECR = `#(both > 0) / #(either > 0)`. Averaged over lineage pairs it is a **dataset-level, panel-arm-independent** contamination score. Lower is cleaner. Used as the platform-comparison and purification-benchmark yardstick.
- **CRISP purity** (Center for Spatial OMICs): the **per-cell** extension of MECR — a marker-positive cell is *impure* if it detects markers from ≥2 disjoint lineages; `purity = 1 − N_impure / N_marker+`. Cell-resolved, so it can gate individual cells. (In our internal CosMx benchmark, CRISP purity rose **0.27 → 0.99 after SPLIT**, independently confirming a ~270× MECR improvement — evidence that these metrics track real contamination and that purification works.)
- **PMP — Positive Marker Purity**: per-cell fraction of a cell's transcripts that come from its *assigned* type's marker set. Needs per-type marker sets; complements CRISP.
- **NCP — Negative Co-expression Purity** (Salas et al., *Nat Methods* 2025): the **percentage of gene pairs that are *not* co-expressed in a matched scRNA-seq reference and remain non-co-expressed in situ**. Scale 0–1, **higher = better specificity** (near 0 = many reference-exclusive pairs now co-occur → contamination). All commercial platforms in that study showed **mean NCP > 0.8**, making **≈0.8 a reasonable dataset-level "good specificity" line**. NCP is *reference-anchored* (needs a matched scRNA-seq dataset) — complementary to MECR, which needs only disjoint-lineage marker pairs.
- **NMP — Negative Marker Purity** (Salas et al., *Nat Methods* 2025): the **percentage of reads from "negative markers" (genes a cell type should NOT express) found in cells assigned to that type**. It is the most direct readout of **segmentation-driven mis-assignment / spillover into annotations** and is the metric they recommend for *comparing segmentation strategies*. Lower = cleaner. Per-cell-type and dataset-level.

### 6.4 Molecular-neighborhood coherence (needs molecule data)

- **NPMI molecular coherence** (TRACER): normalized pointwise mutual information over molecule-level transcript neighborhoods — measures whether co-located transcripts are the ones that *should* co-occur in a real cell. More sensitive than cell×gene co-expression but requires `transcripts.parquet`. Defer unless molecule data is loaded.

### 6.5 Contamination-mitigation tools (not metrics, but they set the achievable purity)

- **SPLIT** (Spatial purification of Layered Intracellular Transcripts) — post-RCTD removal of spatial contamination; the app's canonical purifier. Improves CRISP/MECR dramatically (above) but is a **tradeoff**: spatial coherence can rise or fall by region (a genuine protein-gating tradeoff, not free lunch).
- Molecule-aware re-segmentation (Baysor/proseg) reduces spillover at the source.

**Draft rule.** Per-cell purity signals combine into a **contamination flag**: (`pct_counts_control` high) OR (RCTD `spot_class` ∈ {doublet_uncertain, reject}) OR (CRISP-impure) OR (ovrlpy VSI below section threshold) → **WARN**, and the annotation confidence is capped (see [§7](#7-confidence-scores)). Multiple purity fails → **FAIL → `Ambiguous:mixed`**.

---

<a name="layer-4"></a>
## 7. Layer 4 — Panel adequacy (resolvability)

**This is the layer most other pipelines skip, and it is the honest answer to "why can't you type this cell?"** A targeted panel cannot resolve cell types whose discriminating markers are not on it — *no amount of counts fixes a missing gene.* This gate is **dataset-level and per-cell-type**, and it must run **before** annotation so the app never offers a label the panel cannot support. It is implemented in [`panel_check.py`](../../src/spatialscribe/analysis/panel_check.py) (feature H3).

| Metric | What it measures | Rule |
|---|---|---|
| **Marker coverage** per cell type | Fraction of a type's canonical markers present on the panel | traffic light: **≥3 present = green**, **1–2 = amber**, **0 = red**. Red types **cannot** be called; amber types get a low confidence ceiling. |
| **Pairwise discriminability** | For each type pair, is there ≥1 *private* on-panel marker either way? | Pairs with **no private on-panel marker** are flagged `cannot separate` → they must be **merged** into one reported label. |
| **Confusability merge-groups** | Connected components of the "cannot separate" graph | Report the merged group (e.g. "CD4 T / CD8 T — not separable on this panel") instead of an arbitrary split. |

**Proof case (do not hide it).** On the Xenium 5K panel the melanocyte markers `PMEL/TYR/TYRP1` are **absent** — only `MLANA/MITF/DCT/SOX10` are present. A naive melanocyte score is therefore under-powered, and the honest output is a **panel-adequacy warning**, not a confident melanocyte call built on 4 genes. The app must *surface* this gap, never patch the marker list to conceal it.

**Gene presence ≠ detectability.** A gene on the panel can still drop out (low expression, failing probe). Panel adequacy is **necessary, not sufficient** — always show the coverage numbers alongside the verdict and let empirical detection (Layer 2) further down-weight.

**Draft rule → warning taxonomy.** Layer 4 produces the `Unresolvable-by-panel` abstention reason: if the best-supported label belongs to a red-coverage type or a confusable merge-group, the app reports the **coarsest resolvable label** plus an explicit "panel cannot resolve <subtype> here" warning.

---

<a name="layer-5"></a>
## 8. Layer 5 — Annotation-confidence QC

The label-quality gate: given that the cell is real, has signal, is reasonably pure, and its type is resolvable — **how much do we trust the specific label?** This is where the per-cell confidence score is primarily built. Three independent evidence streams, then a consensus.

### 8.1 Reference-based confidence (needs a reference)

| Method | Native confidence signal | Notes |
|---|---|---|
| **RCTD / spacexr** | per-type **weights** + `spot_class`; **first-vs-second weight ratio** | Doublet-mode `spot_class` is both a purity (Layer 3) and confidence signal. The **weight margin** (top1 − top2) is a strong continuous confidence. |
| **CellTypist** | posterior **probability** per class; `majority_voting` over-clustering | Probability + top1–top2 **margin**; low prob → abstain. |
| **SingleR** | **delta** (score gap to next label) + `pruned.labels` (NA = low-confidence) | `pruned.labels == NA` is a built-in abstention. |
| **Cell2location / Tangram / Stereoscope** | posterior cell-type proportions | Spot/deconvolution-oriented; use the dominant-fraction and its margin. |
| **Pseudobulk-reference correlation** | Pearson/Spearman of the cell (or cluster) profile vs each reference type's pseudobulk | Simple, robust; low max-correlation → out-of-distribution / novel. |

### 8.2 Marker-based confidence (no reference needed)

- **Marker-set score** (`UCell` `maxRank≈150`, `AddModuleScore`, `score_genes`) per candidate type.
- **Score margin** = top1 − top2 marker score. **The margin, not the absolute score, is the confidence.** A high top-1 score with an equally high top-2 is *ambiguous*, not confident.
- **Positive-fraction** = fraction of the type's markers actually detected in the cell (ties to PMP).

### 8.3 Consensus across methods

- **Agreement rate**: fraction of annotation methods (marker-score / CellTypist / reference-transfer / Claude-LLM verdict) that agree on the label. SpatialScribe's design already reconciles marker scoring + CellTypist + Claude into a **consensus with flagged disagreements** — the disagreement flag *is* a confidence signal.
- **Cluster-vs-cell coherence**: does the cell's individual call match its Leiden cluster's majority call? Persistent minorities inside a cluster are either rare subtypes or errors → flag.
- **Annotation stability under transcript subsampling** (evidence-backed, cheap): randomly drop a fraction of a cell's transcripts, re-annotate, and measure the **label-flip rate**. A cell whose call flips under mild subsampling is inherently low-confidence. Motivation: in MERFISH, dropping **~20% of transcripts changed the cell-type label of 10–15% of cells**, and dropping 40 gene species flipped ~10% (eLife 2025) — detection completeness *directly* degrades annotation stability, so subsampling robustness is a principled per-cell confidence proxy that needs no reference.

### 8.4 Turning signals into a calibrated confidence

Raw scores are **not** probabilities and are **not comparable across cell types**. Two principled options (see [§8/solutions](#8-solutions) for the full recommendation):

- **Posterior summaries**: max class probability, **margin** (top1−top2), and **entropy** of the class distribution. Entropy near its max (≈ log K) → the classifier is guessing.
- **Calibration**: temperature/Platt scaling of the scores against a held-out or pseudobulk-labeled set so that "0.9 confidence" means "~90% correct."
- **Conformal prediction** (recommended): produces a **set of plausible labels** with a guaranteed coverage level; abstain when the set has >1 label (ambiguous) or 0 labels (novel). This is the cleanest formalization of "warn when we cannot type the cell."

**Draft rule.** Per-cell annotation confidence ∈ [0,1] is the **minimum** of (calibrated posterior confidence) and (1 − contamination penalty) and (panel-coverage ceiling). Below **0.5 → WARN** (show label greyed / "tentative"); below **0.25 → FAIL → abstain** with a reason code.

---

<a name="layer-6"></a>
## 9. Layer 6 — Spatial-context QC of the label

Uniquely available in spatial data and blind to non-spatial metrics (ARI, silhouette): does the label make **spatial sense**? A lone "hepatocyte" surrounded by T cells is probably a contamination/segmentation error.

| Metric | Definition | Draft use |
|---|---|---|
| **Spatial coherence / neighbor-label agreement** (STEAM) | Mean fraction of a cell's *k* (≈15) spatial nearest neighbors sharing its label | Per-cell continuous confidence input. |
| **PAS — Proportion of Abnormal Spots** | % of cells with **< 20%** same-label neighbors | Dataset/arm-level annotation-noise summary; high PAS = spatially incoherent labels. |
| **Neighborhood enrichment / co-occurrence** (squidpy `nhood_enrichment`, `co_occurrence`) | Whether label pairs are spatially enriched/depleted vs a permutation null | Sanity check on the *annotation as a whole*; biologically impossible adjacencies flag systematic errors. |
| **spARI** (spatially-aware ARI) | ARI that accounts for spatial arrangement of labels | Needs a ground-truth spatial-domain partition; benchmark-only. |

**Caveat — spatial coherence is a double-edged sword.** Rare-but-real infiltrating cells (a single CD8 T in a tumor nest) *legitimately* have low same-label neighbor fractions. Use spatial coherence to **down-weight**, never to hard-delete, and always combine with marker evidence. It penalizes exactly the rare biology that is often most interesting. (In our benchmark, spatial-coherence gains were a genuine tradeoff — helping some samples, hurting others.)

---

<a name="7-confidence-scores"></a>
## 10. Confidence scores, abstention & the "cannot type this cell" warning

*(This section directly addresses the requirement that the app provide confidence scores and warn when cell types cannot be typed, given sparse and noisy spatial data.)*

### 10.1 Design principle: every cell gets a label **and** a reason

The output contract for each cell is:

```
{ label, confidence ∈ [0,1], verdict ∈ {PASS, WARN, FAIL}, reason_code, evidence:{...} }
```

`FAIL` cells are **not** silently dropped and **not** given a made-up type — they get an **abstention label** with a machine- and human-readable reason. The four abstention classes map cleanly onto the layer that produced them:

| Abstention label | Trigger (layer) | Plain-language message to the scientist |
|---|---|---|
| `Unassigned: low quality` | Layer 1–2 (segmentation fail / too few counts) | "This object is too small / has too few transcripts to identify." |
| `Ambiguous: mixed` | Layer 3 (doublet / contamination / low VSI) | "This cell's signal is a mix of neighboring cell types (spatial doublet or spillover)." |
| `Unresolvable: panel` | Layer 4 (red coverage / confusable pair) | "Your panel does not contain the genes needed to tell <A> from <B> here — reported as <coarse group>." |
| `Uncertain: low confidence` | Layer 5–6 (no decisive marker/reference signal, or spatially incoherent) | "No cell type scored clearly above the others for this cell." |
| `Novel / unknown` | Layer 5 (low max reference similarity, conformal empty set) | "This cell does not match any reference type — possibly a state/type not in the reference." |

This is the single most important product decision: **a warning is more useful than a confident wrong label** to a scientist who cannot check it.

### 10.2 How the confidence number is built (recommended)

A per-cell confidence that is honest under sparsity and noise is a **penalized posterior**:

```
confidence = clip01(
    calibrated_posterior_conf            # Layer 5: temperature-scaled max prob (or conformal 1−ambiguity)
  * (1 − contamination_penalty)          # Layer 3: f(pct_control, VSI, RCTD spot_class, CRISP purity)
  * panel_coverage_ceiling               # Layer 4: green=1.0, amber≈0.6, red→abstain
  * spatial_coherence_weight             # Layer 6: mild, ∈[0.8,1.0]; never zero out rare cells
)
low_signal_gate:  if counts < panel-indexed floor (qc.panel_indexed_floor, NOT a fixed 10) or genes < 15 (RICH_PANEL_MIN_GENES on >=1000-gene panels; 5 on targeted) → verdict=FAIL, label=Unassigned:low quality
```

Key properties:
- **Margin-based, not score-based** — top1 vs top2 separation drives the posterior; a high score that is not *separated* from the runner-up yields low confidence.
- **Multiplicative penalties** — any one bad layer can pull confidence down; you cannot buy back a contaminated cell with a strong marker score.
- **Calibrated** — the number means something (see 10.3). Uncalibrated softmax/UCell scores must **not** be shown as "confidence."

### 10.3 Making the confidence trustworthy (calibration)

- **Temperature / Platt scaling** on a held-out labeled subset (or pseudobulk-anchored labels) so displayed confidence ≈ empirical accuracy. Report a **reliability diagram** in the QC section of the HTML export.
- **Conformal prediction** for coverage guarantees: pick a target error α; output the **label set** whose calibrated scores exceed the conformal threshold. Set size 1 → confident; >1 → `Ambiguous` with the candidate list; 0 → `Novel`. This is the most defensible framing for a self-serve tool because the guarantee holds without assuming the classifier is well-specified.
- **Do not invent numbers.** Per the project rule, confidence and verdicts must be grounded strictly in the computed tables (gene presence ≠ detectability). The LLM copilot narrates the numbers; it never manufactures a confidence.

### 10.4 What to *show* the scientist

- Color cells by **verdict** (pass/warn/fail) as a first-class spatial layer, and by **confidence** as a continuous layer — so a scientist sees at a glance *where* the annotation is trustworthy.
- For every reported cell type, show the **panel-adequacy badge** (green/amber/red) next to it.
- A **section-level "annotatability" summary**: % cells confidently typed, % warned, % abstained, and the top reasons — this is the honest headline number.

---

<a name="8-solutions"></a>
## 11. Solutions for sparsity, spatial doublets & transcript diffusion

*(Concrete, implementable mitigations for the three iST pathologies, mapped to where they plug into SpatialScribe.)*

### 11.1 Sparsity (targeted panel → few genes/cell, and missing markers)

| Problem | Solution | Where it lives |
|---|---|---|
| Few detected genes/cell → per-cell scores high-variance | **Hierarchical annotation**: call the **coarse lineage** first (well-powered), attempt **subtypes only where panel coverage + counts allow**; otherwise report the coarse label + "subtype not resolvable." | `annotate.py` (new), gated by `panel_check.py` |
| Single-gene scores noisy | **Module/signature scores** (UCell `maxRank=150`) over marker *sets*, not single genes; borrow strength across the whole set. | `markers.py`, `annotate.py` |
| Missing discriminating markers | **Panel-adequacy gate (Layer 4)** — never offer a label the panel can't support; emit `Unresolvable:panel`. | `panel_check.py` (H3) |
| Per-cell too sparse to call at all | **Reference deconvolution (RCTD)** pools panel-wide signal and yields calibrated weights + doublet class; or **neighborhood-aware smoothing** (borrow from spatial neighbors) with a coherence guard so rare cells aren't blurred away. | RCTD/rctd-py; `spatial.py` |
| Low-count immune cells wrongly deleted | Keep `minCounts=10` (not 20) and route low-count cells through **stronger confidence gating**, not exclusion. | `qc.py` `DEFAULT_FILTER` |

### 11.2 Spatial doublets (merged / vertically-overlapping cells)

| Problem | Solution | Where it lives |
|---|---|---|
| 2D segmentation merges stacked cells | **ovrlpy VSI** to flag low-integrity (overlapping) cells → `Ambiguous:mixed`. | `subprocesses/ovrlpy/` (H5) |
| A cell explained by two types | **RCTD doublet mode** `spot_class`; report `doublet_*`/`reject` as ambiguous, optionally show the **two** most likely types. | rctd-py |
| Systematic co-expression of exclusive lineages | **MECR / CRISP purity** to quantify and, after **SPLIT** purification, verify the improvement. | metrics module; `split-purification` |
| Under-segmentation at the source | Optional **re-segmentation** (Baysor/proseg via sopa) then re-QC. | (stretch) |

### 11.3 Transcript diffusion / optical spillover

| Problem | Solution | Where it lives |
|---|---|---|
| Neighbor transcripts leak into a cell | **Negative-control monitoring** (per-cell + per-section) as the background-floor reference; require real signal clearly above floor. | `qc.py` |
| Spillover inflates false markers | **SPLIT** post-RCTD purification; **nucleus-restricted** counting where boundaries are unreliable; molecule-aware assignment (Baysor). | `split-purification`, io |
| Dense regions worst-hit | Compute contamination metrics **per region** (dense tumor nests vs sparse stroma) and down-weight annotation confidence in high-diffusion regions; expose via H1 region-QC. | `qc.py` `region_qc` (H1) |
| Can't tell spillover from real co-expression | **Panel-arm-independent MECR** as the arbiter; and prefer **marker margin** over absolute score (spillover raises absolutes but rarely the margin of the *correct* type). | metrics; `annotate.py` |
| Over-large expansion pulls in neighbor transcripts | **Tune the expansion distance to biology** (Salas et al.: ideal ≈5.64 µm vs nuclear radius ≈5.06 µm; transcripts >10.71 µm from centroid are mostly background) and **validate with NMP** — re-segment or nucleus-restrict rather than trusting a fixed 15 µm halo. | io / (re-seg, stretch) |

> **Quantified evidence (Salas et al., *Nat Methods* 2025):** diffusion is real and measurable — peripheral transcripts (>~10.71 µm from centroid) carry background rather than cell-own signal, and **~1.8% of cells** additionally showed z-axis mixed-source signal (3D-coherence < 0.2). These are the numbers behind the "your cell's counts are a mixture" warning.

### 11.4 The general recipe

> **Abstain gracefully, aggregate when you can't resolve individually, and always show the reason.** Where a single cell is too sparse/noisy to type, fall back to (a) a coarser label, (b) a neighborhood/region-level composition estimate, or (c) an explicit abstention — never a confident guess. Deconvolution and hierarchical calling recover most of the signal that per-cell hard-calling throws away.

---

<a name="9-consolidated-criteria"></a>
## 12. Consolidated draft pass / warn / fail criteria

**These are starting points to tune on real data, not settled cutoffs.** Direction: "low" = smaller is worse; "high" = larger is worse. Section-level metrics gate the whole run; per-cell metrics gate individual annotations. The machine-readable copy is [`annotation_qc_thresholds.yaml`](annotation_qc_thresholds.yaml).

### 12.1 Section-level (Layer 0) — 10x-published, do not invent

| Metric | Dir | Warn | Fail |
|---|---|---|---|
| `fraction_transcripts_decoded_q20` | low | 0.70 | 0.50 |
| `fraction_transcripts_assigned` | low | 0.50 | 0.30 |
| `fraction_empty_cells` | high | 0.15 | 0.25 |
| `median_genes_per_cell` | low | 20 | 10 |
| `median_transcripts_per_cell` | low | 50 | 25 |
| `negative_control_probe_counts_per_control_per_cell` | high | 0.02 | 0.05 |

### 12.2 Per-cell quality (Layers 1–3) — provisional

| Metric | Dir | Warn | Fail | Basis |
|---|---|---|---|---|
| `total_counts` | low | 25 | 10 | repo `minCounts=10` |
| `n_genes_by_counts` | low | 10 | 5 | repo `DEFAULT_FILTER` |
| `pct_counts_control` | high | 2% | 5% | repo `max_pct_control=5` / 10x neg-control scale |
| `cell_area_um2` | — | <8 or >400 | <4 or >99.5th pct | plausibility, dataset-relative |
| `ovrlpy_VSI` | low | section 25th pct | section 10th pct | relative; tune per section |
| `crisp_purity` (per cell) | low | impure (≥2 lineages) | — | binary CRISP definition |
| `rctd_spot_class` | categorical | `doublet_uncertain` | `reject` | RCTD doublet mode |

### 12.3 Per-cell annotation confidence (Layers 4–6) — provisional

| Signal | Dir | Warn | Fail/abstain | Basis |
|---|---|---|---|---|
| calibrated `confidence` (composite, §10.2) | low | 0.50 | 0.25 | design default; calibrate to accuracy |
| marker score **margin** (top1−top2) | low | small (arm-relative) | ≈0 | ambiguity → abstain |
| reference posterior max prob | low | 0.50 | 0.30 | CellTypist-style |
| conformal label-set size | high | 2 | 0 (novel) or >3 | coverage-guaranteed abstention |
| panel marker coverage (per type) | low | amber (1–2) | red (0) | H3 traffic light |
| spatial coherence (frac same-label neighbors) | low | <0.20 (PAS) | — (down-weight only, never delete) | STEAM/PAS |
| consensus agreement across methods | low | split | strong disagreement | multi-method reconcile |

### 12.4 Section-level annotatability summary (report headline)

- **% cells confidently typed** (verdict PASS) — target > 60–70% on a good FFPE Xenium section (tune).
- **% warned**, **% abstained**, and **top-3 abstention reasons**.
- **Panel-adequacy roster**: which requested cell types are green / amber / red on this panel.

---

<a name="10-code-mapping"></a>
## 13. How this maps to SpatialScribe code

| Layer / concept | Module | Status |
|---|---|---|
| Layer 0 section thresholds | `analysis/qc.py` (`XENIUM_THRESHOLDS`, `qc_summary`, `_flag`) | ✅ implemented |
| Layer 1–2 per-cell filters | `analysis/qc.py` (`compute_qc`, `DEFAULT_FILTER`, `apply_filter`) | ✅ implemented (extend with cell-area, single-gene-dominance) |
| Region QC (H1) | `analysis/qc.py` (`region_qc`) | ✅ implemented |
| Layer 3 contamination metrics (MECR/CRISP/PMP) | **new** `analysis/purity.py` (proposed) | ⬜ to build — port from `an internal benchmark script` |
| Layer 3 VSI | `subprocesses/ovrlpy/run_ovrlpy.py` + join | ◻ scaffolded (H5) |
| Layer 4 panel adequacy | `analysis/panel_check.py` (`check_panel`) | ✅ implemented |
| Layer 5 annotation + confidence | **new** `analysis/annotate.py` (multi-method + consensus + calibrated confidence + abstention) | ⬜ to build — **this document is its QC spec; [`annotation-method-selection.md`](annotation-method-selection.md) is its method-choice spec** (SingleR/RCTD defaults, reference > method, no zero-shot FMs) |
| Layer 6 spatial coherence / PAS | `analysis/spatial.py` (add `spatial_coherence`, `pas`) | ◻ extend |
| Confidence → UI (verdict/reason colors, badges, headline) | `app/` + `llm.py` verdict narration | ◻ extend |
| Machine-readable thresholds | [`annotation_qc_thresholds.yaml`](annotation_qc_thresholds.yaml) | ✅ this deliverable |

**Build order suggestion:** (1) extend `qc.py` with cell-area + single-gene-dominance; (2) add `purity.py` (MECR/CRISP/PMP — cheap, cell×gene only); (3) build `annotate.py` around the penalized-posterior confidence + the five abstention classes; (4) wire the verdict/confidence/reason into the spatial canvas and the HTML report; (5) ride-along ovrlpy VSI and RCTD `spot_class` when molecule data / a reference is present.

---

<a name="14-cross-platform"></a>
## 14. Cross-platform QC comparison (secondary: MERSCOPE, CosMx)

Xenium is the primary target, but SpatialScribe also ingests CosMx and Atera data, and the QC gates must be **platform-aware** — the same numeric cutoff means different things per platform. Verified cross-platform findings (2025 benchmarks):

| QC axis | Finding | Implication for gates |
|---|---|---|
| **Specificity / FDR** | On-target fraction & false-discovery ranking: **Xenium > MERSCOPE > CosMx**. Xenium consistently lowest FDR / highest on-target; CosMx highest FDR. After total-count normalization, Xenium 5K showed a lower negative-control proportion than CosMx 6K (which also shows stronger *spatial aggregation* of negative-control signal). | The negative-control / NCP thresholds should be **looser (more forgiving) for CosMx**, tighter for Xenium; a CosMx cell needs a bigger margin over background to be trusted. |
| **Segmentation precision** | On densely packed cells, **CosMx & Xenium ≈ 0.90 vs MERSCOPE ≈ 0.83** (p < 0.01), against >31,400 manually annotated cells (dense/sparse/elongated). | Segmentation-derived confidence penalties should be **larger for MERSCOPE in dense tissue**. |
| **Per-cell QC retention** | Standard filtering keeps **~92–96% of Xenium & CosMx cells** but only **~28% and ~3% of MERSCOPE cells** in two FFPE TMAs (i.e. 70–97% loss). | The "acceptable % cells retained" gate must be **platform- and assay-aware** — a 60% retention that's alarming on Xenium can be *normal* on MERSCOPE FFPE. Don't hard-fail a section on retention alone. |
| **Detection completeness → label stability** | In MERFISH/MERSCOPE, losing ~20% of transcripts flips **10–15% of cell labels**; 40 dropped gene species flip ~10%. | Lower-completeness platforms need **stronger annotation-stability gating** (see the subsampling metric, [§8.3](#layer-5)) and more conservative abstention. |
| **Sensitivity** | Xenium/CosMx/MERSCOPE differ in per-cell counts, but a specific "Xenium 2× CosMx / 20× MERSCOPE sensitivity ceiling" claim **was refuted in verification** and is *not* asserted here. | Do not hard-code cross-platform sensitivity ratios; calibrate the count floor per dataset. |

**Bottom line:** ship a small **per-platform threshold profile** (Xenium / MERSCOPE / CosMx) layered on top of the panel-size indexing. The *metric families* are identical across platforms; only the cutoffs move.

---

<a name="11-caveats"></a>
## 15. Caveats, open problems & anti-patterns

- **No universal thresholds (10x's own words, and the field's).** Everything scales with tissue (fresh vs FFPE), panel size, and goal (discovery vs validation). Ship thresholds as **tunable config**, expose them, and default to the distributional/relative variants where possible. The **Spatial Touchstone / SpatialQM** effort (Plummer et al., *Nat Biotechnol* 2025) — the field's reference multi-platform standardization project (six tissues, multiple sites, centralized sectioning) — deliberately provides **tissue-specific baselines and an adaptable metrics toolkit rather than universal pass/fail cutoffs**. Treat our numbers as **local defaults calibrated against a tissue-matched baseline**, not constants.
- **Custom-panel Q20 trap.** `fraction_custom_transcripts_decoded_q20` is computed over **custom add-on probes only**. A low value flags *specific failing custom probes* (verify per-probe from `transcripts.parquet`, `qv≥20`), **not** a sample-quality problem — the predesigned panel and overall data are usually fine. Do not tell a client their sample failed based on this metric alone.
- **Gene presence ≠ detectability.** Panel adequacy is necessary, not sufficient; a listed gene can drop out. Always pair Layer-4 coverage with Layer-2 empirical detection.
- **Spatial coherence penalizes rare biology.** Infiltrating/rare cells legitimately have incoherent neighborhoods. Use it to down-weight, never to delete.
- **Protein-anchored / reference-anchored metrics are not independent.** If you validate an RNA annotation with a protein gate or a reference-derived silver label, that metric favors the reference/protein arm by construction. Keep at least one **method-independent** check (MECR/CRISP, or a least-circular pseudospot ground truth) in the battery. Use the **resting** lineage cluster as a pair proxy; reserve **cycling** clusters for the proliferation axis (a cycling-CD8 cluster used as an "NK" proxy silently flips results).
- **Min-max scaling inflates tiny gaps.** When comparing a few methods, min-max normalization turns a negligible raw difference into a big visual bar — annotate raw values.
- **Calibration drift.** A confidence calibrated on one tissue/panel does not transfer; re-calibrate per dataset, or use conformal prediction (distribution-free).
- **Reference/ground-truth is imperfect.** Matched snRNA-seq (the usual contamination ground truth) captures **nuclear, not cytoplasmic** transcripts, so "precise contamination quantification" is aspirational; NCP/NMP/RCTD are all reference-conditioned and inherit its blind spots.
- **Numbers are chemistry-version- and era-specific.** The cross-platform rankings (Xenium > MERSCOPE > CosMx) and most magnitudes come from ~2023-era panels/chemistries; MERFISH 2.0, updated Xenium Prime 5K, and CosMx WTx/6K can shift them. The 10x default expansion already moved 15 µm → 5 µm. Re-baseline against a **tissue- and chemistry-matched** reference (e.g. a Spatial Touchstone profile), not against these constants.
- **Segmentation precision is one metric on one scenario;** don't over-generalize the 0.90-vs-0.83 numbers. Retention percentages also varied internally across TMAs.
- **Open problems (also flagged by the research pass):** (i) **no validated numeric cutoffs for reference-based deconvolution confidence** (RCTD singlet/doublet boundaries, secondary-weight fractions) as a per-cell annotation gate — an explicit gap to close with our own data; (ii) **spatial-coherence and marker-specificity thresholds remain largely qualitative** in the literature; (iii) no consensus **panel-size-indexed threshold schedule** exists — we have to build it; (iv) doublet vs rare-real-hybrid is genuinely ambiguous without orthogonal (protein) evidence; (v) novel-type / out-of-distribution detection in targeted panels is under-studied.

---

<a name="12-glossary"></a>
## 16. Glossary

- **iST** — imaging-based spatial transcriptomics (Xenium, MERFISH/MERSCOPE, CosMx).
- **VSI** — Vertical Signal Integrity (ovrlpy); per-cell overlap/doublet score from z-coherence.
- **MECR** — Mutually Exclusive Co-expression Rate (dataset-level contamination).
- **CRISP purity** — per-cell version of MECR; a cell is impure if it detects ≥2 disjoint lineages.
- **PMP** — Positive Marker Purity (per-cell fraction of transcripts from the assigned type's markers).
- **PAS** — Proportion of Abnormal Spots (% cells with <20% same-label spatial neighbors).
- **spot_class** — RCTD doublet-mode category (`singlet`/`doublet_certain`/`doublet_uncertain`/`reject`).
- **Margin** — top1 − top2 score/probability; the true confidence signal, more than the absolute.
- **Conformal prediction** — distribution-free method yielding label *sets* with coverage guarantees; empty set → novel, >1 → ambiguous.
- **Abstention** — returning `Unassigned`/`Ambiguous`/`Unresolvable`/`Uncertain`/`Novel` instead of a forced label.

---

<a name="13-references"></a>
## 17. References

### 17.1 Load-bearing primary sources (verified in the research pass)

1. **Salas S.M., Kuemmerle L.B., Mattsson-Langseth C., et al.** "Optimizing Xenium In Situ data utility by quality assessment and best-practice analysis workflows." *Nature Methods* 22(4):813–823 (2025). DOI [10.1038/s41592-025-02617-2](https://doi.org/10.1038/s41592-025-02617-2) · open access [PMC11978515](https://pmc.ncbi.nlm.nih.gov/articles/PMC11978515/). — **The definitive Xenium QC best-practices paper.** Source for **NCP** (negative co-expression purity, mean >0.8), **NMP** (negative marker purity, for segmentation quality), the **3D-coherence** overlap score (<0.2 → ~1.8% cells), the **transcript-diffusion distances** (nuclear radius ~5.06 µm, ideal expansion ~5.64 µm, background beyond ~10.71 µm), per-cell read stats (186.6 mean, 76.8% assigned, exclude <10 reads = 0.21%), gene detection efficiency, and the Cellpose+Baysor best-practice pipeline.
2. **Plummer J., Cook D.P., Martelotto L.G., et al.** "Standardized metrics for assessment and reproducibility of imaging-based spatial transcriptomics datasets." *Nature Biotechnology* (2025). DOI [10.1038/s41587-025-02811-9](https://doi.org/10.1038/s41587-025-02811-9). — **Spatial Touchstone / SpatialQM**: a **7-axis technical QC taxonomy** (reproducibility, sensitivity, dynamic range, signal-to-noise, false-discovery rate, **cell-type annotation**, and congruence with single-cell profiling), an open-source tool (**SpatialQM**, also does reference-based annotation transfer/imputation), standardized SOPs, and a 254-profile public repository (GEO **GSE277080**). Provides **tissue-specific baselines rather than universal thresholds** (the basis for "calibrate locally"). This is a concrete tool-based annotation-QC option alongside squidpy / sopa / RCTD / SPLIT.
4. **Mitchel J., et al.** "Impact and correction of segmentation errors in spatial transcriptomics." *Nature Genetics* (2025). DOI [10.1038/s41588-025-02497-4](https://doi.org/10.1038/s41588-025-02497-4). — Segmentation error is first-order; cell-type *classification* is comparatively robust while *state* analyses (DE, neighbor-influence, ligand–receptor) are highly sensitive; introduces **cellAdmix** admixture correction.
5. **Marconato L., et al. (SPLIT).** "SPLIT: Spatial purification of layered intracellular transcripts." *Nature Methods* (2026). DOI [10.1038/s41592-026-03089-8](https://doi.org/10.1038/s41592-026-03089-8) · bioRxiv [2025.04.23.649965](https://www.biorxiv.org/content/10.1101/2025.04.23.649965v1.full) · code [github.com/bdsc-tds/SPLIT](https://github.com/bdsc-tds/SPLIT). — Source for the **~60% Prime-5K QC-fail** finding and the RCTD-weight-based post-hoc purification (segmentation-agnostic).
6. **Cross-platform benchmarks (secondary platforms):** segmentation precision, per-cell retention & FDR ranking — *Nat Commun* [s41467-025-64990-y](https://www.nature.com/articles/s41467-025-64990-y); negative-control cross-platform — *Nat Commun* [s41467-025-64292-3](https://www.nature.com/articles/s41467-025-64292-3); MERFISH detection-completeness → label-stability — *eLife* [reviewed-preprint 105149](https://elifesciences.org/reviewed-preprints/105149).
3. **10x Genomics** — Xenium Onboard Analysis QC metrics, definitions & troubleshooting (decode Q20, transcripts assigned, empty cells, genes/transcripts per cell, negative-control probe rate). [Metric calculations](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/advanced/metric-calculations) · [Analysis-summary troubleshooting](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/analysis/analysis-summary-troubleshooting) · [Negative-control-probe metric KB](https://kb.10xgenomics.com/hc/en-us/articles/18385764969613).

### 17.2 Methods, metrics & tools

- **RCTD** (Robust Cell Type Decomposition), doublet mode & `spot_class` — Cable D.M. et al., *Nat Biotechnol* 2022.
- **MECR** — Mutually Exclusive Co-expression Rate — Hartman & Satija, 2024.
- **CRISP** purity — Center for Spatial OMICs. · **PMP** — positive marker purity.
- **STEAM** — spatial coherence / PAS for annotation evaluation. · **TRACER** — NPMI molecular-neighborhood coherence. · **spARI** — spatially-aware ARI.
- **ovrlpy** — Vertical Signal Integrity for imaging spatial doublets/overlap.
- **SPLIT** — Spatial purification of layered intracellular transcripts (post-RCTD).
- Segmentation: **Baysor** (Petukhov et al., *Nat Biotechnol* 2022), **proseg**, **BIDCell**, **sopa**, **segger**, **ComSeg**, **Cellpose/Mesmer**, **SSAM**, **Points2Regions**.
- **squidpy** (neighborhood enrichment, co-occurrence) — Palla et al., *Nat Methods* 2022.
- Reference-based annotation with native confidence signals: **CellTypist**, **SingleR**, **Cell2location**, **Tangram**, **SpaGE** (imputation).
- Platform refs: **Xenium** (Janesick et al.), **MERFISH/MERSCOPE** (Chen/Moffitt et al.), **CosMx** (He et al.).

### 17.3 internal sources

- `the 10x section-level QC thresholds` — the section-level 10x thresholds mirrored in `qc.py`.
- internal CosMx SPLIT benchmark: `an internal benchmark script` (CRISP 0.27→0.99 post-SPLIT; MECR ~270×).
- `an internal QC-metric note` memory — the cheap-vs-needs-extra-data metric-feasibility taxonomy.

> **Verification status.** The numeric claims above were checked by an adversarial-verification research pass (6 search angles, 24 sources fetched, 113 claims extracted → **25 verified: 23 confirmed, 2 refuted-and-excluded**, 12 synthesized). The **2 refuted** claims are deliberately **not** asserted here: (a) an "NCP ≈ 0.8 warn/fail threshold" *misattributed to the SPLIT paper* — NCP and its >0.8 line come only from Salas et al. (PMC11978515), which this doc cites correctly; and (b) a "Xenium ~2× CosMx / ~20× MERSCOPE sensitivity ceiling" claim, which failed verification — so no cross-platform sensitivity *ratio* is hard-coded. Separately, the NCP-terminology trap (negative *control probe* rate ≠ negative *co-expression* purity) is kept explicit throughout. Only one headline number — the **~60% Prime-5K QC-fail** — rested on a split (2-1) vote; treat it as strong-but-single-source. Citation venues/years for §17.2 items are from working knowledge + the metric taxonomy and should be link-verified before external publication.
</invoke>
