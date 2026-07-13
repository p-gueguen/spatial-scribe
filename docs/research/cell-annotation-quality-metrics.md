# Cell Annotation Quality Metrics — a consolidated reference

**Status:** authoritative design reference (the metrics catalog) · **Scope:** cell-type annotation quality for scRNA-seq **and** imaging spatial transcriptomics (Xenium / MERSCOPE / CosMx) · **Audience:** anyone at a genomics facility evaluating "is this cell-type annotation any good?" — in SpatialScribe, in a custom report, or in a benchmark.

> **Why this document exists, and how it relates to the others.** We already have [`cell-annotation-qc.md`](cell-annotation-qc.md), which is excellent — but it is a **spatial pipeline QC funnel** (six ordered gates: segmentation → counts → contamination → panel adequacy → confidence → spatial coherence) scoped Xenium-first and oriented toward *per-cell gating and abstention inside the SpatialScribe pipeline*. It does **not** consolidate two whole families of metric: (a) **reference-based / external-validation** metrics used when a ground truth or trusted reference exists (accuracy, F1, ARI, NMI, kappa, biological fidelity), and (b) **internal cluster-validity + inter-sample-consistency** metrics used when there is no ground truth (silhouette variants, neighborhood purity, medoid/Ward concordance, inter-sample consistency — the family the Carmona lab's [**scTypeEval**](https://github.com/carmonalab/scTypeEval) packages). This document is the **catalog**: every metric worth using to *evaluate* an annotation, organized by what evidence you have, what each measures, whether it needs a reference, and where it lives (or should live) in code. Use this to *choose a battery*; use `cell-annotation-qc.md` for the spatial per-cell funnel and [`annotation_qc_thresholds.yaml`](annotation_qc_thresholds.yaml) for the tunable numbers; use [`annotation-method-selection.md`](annotation-method-selection.md) to choose *which annotator to run*.

**The one thing to internalize:** in real the cluster projects there is **almost never a ground truth**. So the **internal / reference-free** metrics (§3) are the default regime, and the **external** metrics (§2) apply only to benchmarks with labels. Do not reach for accuracy/F1 on a delivery — you have no truth to compute them against.

---

## Table of contents

1. [Decision guide — which metrics to compute](#1-decision-guide)
2. [External / reference-based metrics (ground truth or trusted reference)](#2-external)
3. [Internal / reference-free metrics (the default regime)](#3-internal)
   - 3a. [Cluster-validity / separation (scTypeEval + classics)](#3a-cluster-validity)
   - 3b. [Inter-sample consistency — ISC (scTypeEval)](#3b-isc)
   - 3c. [Marker / contamination purity (reference-free)](#3c-purity)
   - 3d. [Per-cell confidence / uncertainty](#3d-confidence)
   - 3e. [Spatial-context metrics (spatial only)](#3e-spatial)
   - 3f. [Reference-anchored specificity (needs a matched profile, not labels)](#3f-anchored)
   - 3g. [Marker-program fidelity (needs marker sets, not truth)](#3g-marker-fidelity)
   - 3h. [Signal-quality QC (annotation-independent; SpatialQM)](#3h-signal-qc)
4. [What we standardize on — default metric batteries](#4-batteries)
5. [Master metric table (quick reference)](#5-master-table)
6. [Implementation status in SpatialScribe (`eval_metrics.py`)](#6-implementation)
7. [scTypeEval — what it adds and how to fold it in](#7-sctypeeval)
8. [Cross-cutting rules & anti-patterns](#8-rules)
9. [References](#9-references)

---

<a name="1-decision-guide"></a>
## 1. Decision guide — which metrics to compute

Pick the column that matches what you actually have. This is the whole document in one table; the rest is detail.

| You have… | Regime | Compute (in priority order) |
|---|---|---|
| **Ground-truth labels** (manual gold standard, or a benchmark dataset) | External (§2) | **macro-F1** + **balanced accuracy** + **ARI** + per-class F1; add **AvgBIO / program concordance** because accuracy ≠ biological fidelity. Never rank methods on plain accuracy alone. |
| **A trusted reference atlas** but no per-cell truth | External-ish (§2.4) + Internal | Label-transfer native confidence (SingleR `delta`/`pruned`, CellTypist prob, RCTD `spot_class`/weight margin) as a quality proxy; **NCP/NMP** (§3f) if the reference is matched scRNA/snRNA. |
| **Biological replicates** (≥2 samples per condition) but no truth | Internal (§3b) | **Inter-sample consistency (ISC)** — the scTypeEval idea: a real cell type looks consistent across samples. Pair with cluster-validity (§3a). |
| **Just labels + an embedding / expression matrix** (the common case) | Internal (§3a, §3d) | **Silhouette + 2label-silhouette + neighborhood purity** (scTypeEval); **marker-score margin + entropy**; **subsampling label-stability**. |
| **Marker gene sets** (curated per lineage) | Internal (§3c) | **MECR + CRISP purity + PMP + on-target marker specificity** — reference-free contamination readouts. |
| **Spatial coordinates** | Internal (§3e) | **Spatial coherence / PAS**, neighborhood enrichment sanity. Down-weight, never delete (rare infiltrating cells are legitimately incoherent). |
| **A targeted panel** (Xenium/CosMx) | Internal (§3c, §3f) + panel gate | Everything above **plus panel adequacy** (can the panel even resolve the type — [`cell-annotation-qc.md` §7](cell-annotation-qc.md)). |

> Most SpatialScribe deliveries land in rows 4–7 simultaneously (labels + embedding + markers + coords, no truth). The **default spatial battery** in [§4](#4-batteries) is exactly that combination.

---

<a name="2-external"></a>
## 2. External / reference-based metrics (ground truth or trusted reference)

Use these **only** when you have labels you trust as truth — a manual gold standard, a held-out benchmark, or (with caveats) a reference atlas. They answer "how close is the predicted labeling to the reference labeling?"

### 2.1 Per-cell agreement

| Metric | Definition | Range / dir | When to use / gotcha |
|---|---|---|---|
| **Accuracy** | fraction of cells whose predicted label == truth | 0–1, higher better | Dominated by the majority class; **misleading under class imbalance** (a "call everything tumor" classifier scores high). Almost never report alone. |
| **Balanced accuracy** | mean of per-class recall | 0–1, higher better | Corrects majority-class dominance; the right "overall accuracy" number when types are imbalanced (they always are). |
| **Cohen's κ** | agreement corrected for chance | −1…1, higher better | Good single-number agreement that accounts for label prevalence; κ>0.8 strong, 0.6–0.8 moderate. |

### 2.2 Per-class quality (precision / recall / F1)

- **Precision**(t) = TP / (TP+FP) — of cells called t, how many really are t.
- **Recall**(t) = TP / (TP+FN) — of true-t cells, how many were found.
- **F1**(t) = harmonic mean of the two.
- **macro-F1** = unweighted mean of per-class F1 → **weights rare types equally** (the metric that exposes "great on tumor, useless on the 200 T cells"). **weighted-F1** weights by class size (closer to accuracy).

**Rule:** for cell typing, **macro-F1 + per-class F1 table** is the headline, because the failure that matters (dropping a rare immune type) is invisible in accuracy and weighted-F1. This is exactly the RCTD-drops-T-cells failure the SPLIT comparison surfaced.

### 2.3 Partition agreement (label sets as clusterings)

| Metric | Measures | Range | Notes |
|---|---|---|---|
| **ARI** (Adjusted Rand Index) | pairwise co-assignment agreement, chance-corrected | −1…1 (0 = random, 1 = identical) | The standard partition-agreement metric; robust, label-name-independent (doesn't need matched label vocab). |
| **NMI / AMI** (Normalized / Adjusted Mutual Information) | shared information between partitions | 0–1 | NMI is intuitive; **AMI** is the chance-corrected version — prefer AMI when cluster counts differ. |
| **ECS** (Element-Centric Similarity) | cell-level label consistency vs a reference labeling, via an affinity-matrix formulation (Gates & Ahn 2019; α=0.9) | 0–1 | **Resolution-agnostic** — designed to compare labelings with *different numbers of clusters*, where ARI/NMI are unfair. Zhu 2026's headline structure-metric; per-cell, so it also localizes disagreement. Add it whenever prediction/truth granularity differ. |
| **V-measure** (homogeneity + completeness) | homogeneity = each cluster is one type; completeness = each type in one cluster | 0–1 | Decomposes *why* a partition disagrees — over- vs under-splitting. |
| **spARI** (spatially-aware ARI) | ARI that credits spatial arrangement | — | Spatial benchmarks only; needs a ground-truth spatial-domain partition. |

Granularity trap: ARI/NMI between a fine prediction and a coarse truth (or vice-versa) is low **by construction**, not because the annotation is wrong. Match granularity, or evaluate hierarchically.

### 2.4 Reference native-confidence as a quality proxy (no hard truth, trusted reference)

When you ran a reference-based annotator, its own posterior is a *soft* quality signal (not ground truth, but informative):

- **SingleR** `delta` (score gap to next label) and `pruned.labels` (NA = abstained).
- **CellTypist** posterior probability + top1−top2 **margin**; `majority_voting`.
- **RCTD** `spot_class` (`singlet`/`doublet_*`/`reject`) + first-vs-second **weight margin**.
- **Pseudobulk-reference correlation**: Pearson/Spearman of the cell (or cluster) profile vs each reference type; low max-correlation → out-of-distribution / novel.

These also appear in the spatial funnel ([`cell-annotation-qc.md` §8.1](cell-annotation-qc.md)); listed here so the external-regime picture is complete.

### 2.5 Biological fidelity (accuracy ≠ fidelity)

The load-bearing finding of the Zhu et al. 2026 spatial benchmark (this paper — `sources/2026.06.16.732716v1.full.pdf`; also in [`annotation-method-selection.md`](annotation-method-selection.md)): **high classification accuracy does not guarantee preserved biology** — rankings by fidelity frequently *differ* from accuracy rankings, so accuracy alone "may substantially overestimate" real quality. Complement accuracy with:

- **AvgBIO** (Average Biological Conservation, from scIB) — a summary of structure-preservation metrics, each oriented larger-is-better and scaled to [0,1], then averaged. **Exact form used in Zhu 2026:** `AvgBIO = mean(ARI, NMI, ASW_scaled)` with `ASW_scaled = (ASW+1)/2` (ASW = Average Silhouette Width on the labeled embedding). Reduces over-reliance on any single metric.
- **Marker-gene overlap** — overlap between each method's data-derived per-type DE genes (top 10 / 20 / 30) and the reference cell-type marker genes: do the labels reproduce the *known* markers?
- **ssGSEA program concordance + AUC-ROC on enrichment** — score each cell against ground-truth marker/pathway sets with ssGSEA, then per type compute the **AUC-ROC** of that enrichment score ranking predicted-positive above predicted-negative cells. Directly asks: does a method's label assignment track the cells' actual marker-program enrichment? (Reference-free given marker sets — see [§3g](#3g-marker-fidelity).)
- **Cohen's d (score separability)** — standardized effect size of how cleanly a method's binary type call separates target vs non-target cells on the marker-enrichment score.

**Rule:** in any benchmark, report at least one fidelity metric next to accuracy/F1.

### 2.6 Hierarchy-, composition- & robustness-aware external metrics

The Zhu 2026 framework adds three axes plain accuracy misses:

- **Hierarchical annotation accuracy** (the "accuracy tree") — accuracy scored *along the cell-type lineage hierarchy*, so a call that gets the **major class right but the subtype wrong** earns partial credit instead of a flat miss. This is the principled fix for the granularity trap (§8): a coarse-vs-fine mismatch should not be penalized like a lineage error. Pair with **subtype-level accuracy** (accuracy restricted to fine subtypes) to see *where* resolution breaks down.
- **Cell-type composition / proportion accuracy** — compare predicted vs true **cell-type proportions**. Two settings: (i) **hard labels → global proportions** per section/condition (`composition_accuracy`: L1 / Pearson / JSD on the length-K proportion vectors) — catches methods that get per-cell labels roughly right but distort the population mix (the paper tracks this across a kidney injury–repair time course); (ii) **soft per-spot proportion matrices** (spots × cell-types), i.e. **deconvolution** output — scored with the **OpenProblems `task_spatial_decomposition`** metrics: **R²** (`r2_score(true, pred, multioutput='uniform_average')`, per cell type averaged; ↑) and **JSD** (mean over cell types of `scipy.jensenshannon(true[:,k], pred[:,k])`, axis=0; ↓), plus **RMSE** and per-spot JSD. Use (ii) to score a deconvolver's proportion output (RCTD/TACCO/cell2location weights) against known mixtures / pseudospots (`eval_metrics.deconvolution_metrics`).
- **Robustness to reference volume** — re-run with the reference **down-sampled** (e.g. 25/50/75%) and track accuracy / macro-F1 degradation. A method whose accuracy collapses as the reference shrinks is fragile in the common low-reference regime. Distinct from **query** subsampling ([§3d](#3d-confidence)), which needs no reference.

*(Practical fourth axis in the paper — **scalability**: wall-clock time + peak memory. Not a quality metric, but the reason GraphST OOM'd on 1.37 M cells; track it when method choice must scale.)*

---

<a name="3-internal"></a>
## 3. Internal / reference-free metrics (the default regime)

No ground truth. These test whether a labeling is **self-consistent** — with the embedding, the markers, the neighbors, and across samples. This is what you compute on a real delivery.

<a name="3a-cluster-validity"></a>
### 3a. Cluster-validity / separation (scTypeEval + classics)

Treat each annotated cell type as a cluster and ask: is it **cohesive** (cells of a type are similar) and **separated** (distinct from other types)? Computed in an embedding (PCA/latent) or in a pseudobulk/expression space.

| Metric | Definition | Range / dir | Source |
|---|---|---|---|
| **Silhouette** | per cell, (b−a)/max(a,b) where a = mean within-type distance, b = mean nearest-other-type distance; averaged per type | −1…1, higher better | classic (Rousseeuw 1987); **scTypeEval** |
| **2label-silhouette** | silhouette variant comparing a cell's own type vs **all others pooled** (one-vs-rest) — more stable with many types | −1…1, higher | **scTypeEval** |
| **Neighborhood purity** | fraction of a cell's *k* nearest neighbors (in the embedding) sharing its label | 0–1, higher | **scTypeEval** (`NeighborhoodPurity`); classic kNN purity |
| **Orbital-medoid** | fraction of a type's cells closer to their **own type medoid** than to any other type's medoid | 0–1, higher | **scTypeEval** |
| **Ward-PropMatch** | proportion of a type that falls in its **dominant Ward cluster** (agreement between the labels and an unsupervised Ward clustering) | 0–1, higher | **scTypeEval** |
| **Average-similarity** | within-type similarity relative to between-type similarity | higher | **scTypeEval** |
| **Davies–Bouldin** | mean over types of (within-scatter / between-separation) to the worst-confused neighbor | ≥0, **lower** better | classic; candidate, not in scTypeEval |
| **Calinski–Harabasz** | between/within variance ratio (F-like) | higher | classic; candidate |
| **Graph modularity / connectivity** | modularity of the label partition on the kNN graph; connectivity penalty for split types | higher modularity better | classic; candidate |

**scTypeEval's key framing** (its benchmark over 31 datasets): use a **local** consistency metric to catch **over-partitioning** (silhouette in a reciprocal-classification match space) and a **global** consistency metric to catch **signal degradation** (2label-silhouette in pseudobulk-cosine space), then combine them with a **geometric-mean integrated score**. The important idea for us: *neighborhood purity / silhouette computed on the labels are a direct, reference-free "is this type real?" signal* — under-used in the cluster spatial reports, and cheap.

Interpreting on sparse panels: silhouette on a raw sparse Xenium matrix is noisy — compute it on the **PCA/latent embedding** (or on pseudobulk profiles for a per-type read), never on raw counts, and remember the `maxRank=150` sparsity caveat when scores feed in.

<a name="3b-isc"></a>
### 3b. Inter-sample consistency — ISC (scTypeEval)

**The distinctive scTypeEval contribution.** If you have biological replicates, a *correctly* annotated cell type should have a **reproducible profile across samples**; an over-split or mislabeled type will not. ISC turns replicates into a reference-free quality signal:

- **Pseudobulk distance** across samples per type — aggregate each (sample, type) to a pseudobulk profile, measure within-type across-sample distance (euclidean / **cosine** / pearson); tight = consistent.
- **Single-cell distribution distance** — Wasserstein distance between per-type single-cell distributions across samples (no aggregation).
- **Reciprocal classifiers** — train on sample A, predict sample B and vice-versa; consistent types are recovered reciprocally.

This is genuinely additive to `cell-annotation-qc.md` (which has no multi-sample notion) and to the cluster practice: when a report has replicates, ISC answers "is this cluster a real cell type or a batch/over-clustering artifact?" without any atlas.

<a name="3c-purity"></a>
### 3c. Marker / contamination purity (reference-free)

Do a cell's transcripts belong to its assigned type, or are they a spillover mixture? These need only **marker gene sets** (no reference atlas). Detailed in [`cell-annotation-qc.md` §6](cell-annotation-qc.md); summarized here because they are core annotation-quality metrics, and several are already implemented in [`compare_rctd_tacco.py`](../../scripts/split/compare_rctd_tacco.py).

| Metric | Definition | Dir | Implemented |
|---|---|---|---|
| **MECR** (Mutually Exclusive Co-expression Rate) | over disjoint-lineage marker pairs, `#(both>0)/#(either>0)`, averaged | lower cleaner | `purity.py` (mecr); SPLIT compare |
| **CRISP purity** | per-cell: impure if it detects markers of ≥2 disjoint lineages | higher cleaner | `purity.py` (`crisp_purity`) |
| **PMP** (Positive Marker Purity) | per-cell fraction of counts from the **assigned** type's markers | higher cleaner | `purity.py` (`pmp`) |
| **On-target marker specificity** | on-target fraction of a lineage's marker expression: on/(on+off) mean, bounded [0,1] | higher cleaner | `compare_rctd_tacco.py` (`_specificity`) |
| **On-target signal retention** | fraction of a lineage's own-marker counts kept after a purification step | higher = less over-purified | `compare_rctd_tacco.py` |

These are the metrics the RCTD-vs-TACCO-into-SPLIT comparison uses, and they are the reference-free way to *quantify contamination* — the dominant iST annotation-quality problem.

<a name="3d-confidence"></a>
### 3d. Per-cell confidence / uncertainty

Per-cell "how sure are we about *this* label?" Detailed in [`cell-annotation-qc.md` §8 + §10](cell-annotation-qc.md); the metrics:

- **Marker-score margin** = top1 − top2 marker score. **The margin, not the absolute score, is the confidence.**
- **Entropy** of the class-score distribution (near log K → the classifier is guessing).
- **Reference posterior max-prob** (CellTypist-style) — where a reference exists.
- **Conformal label-set size** — distribution-free: set size 1 = confident, >1 = ambiguous, 0 = novel. The cleanest formalization of "abstain." **Implemented:** `eval_metrics.conformal_prediction_sets` (split-conformal LAC on any annotator's probability matrix; marginal + **class-conditional/Mondrian** — the latter is required or rare types are silently under-covered — + optional lineage-hierarchy collapse à la [scConform](https://github.com/ccb-hms/scConform)). Coverage holds under calibration↔query exchangeability; on native spatial calibrated against a dissociated reference, treat nominal coverage as approximate (recalibrate on a small platform-matched labeled subset, e.g. protein-confirmed cells) — the set SIZE stays a valid relative uncertainty.
- **Annotation stability under subsampling** — drop ~20% of a cell's transcripts, re-annotate, measure label-flip rate; reference-free, evidence-backed (MERFISH: ~20% dropout flips 10–15% of labels). Implemented in `annotate.py` (`annotation_stability`).

<a name="3e-spatial"></a>
### 3e. Spatial-context metrics (spatial only)

Does the label make **spatial sense**? Blind to non-spatial metrics. Detailed in [`cell-annotation-qc.md` §9](cell-annotation-qc.md):

- **Spatial coherence** = mean fraction of a cell's *k* (~15) spatial neighbors sharing its label (`spatial.py` `spatial_coherence`).
- **PAS** (Proportion of Abnormal Spots) = % cells with <20% same-label neighbors — a dataset-level annotation-noise summary.
- **Neighborhood enrichment / co-occurrence** (squidpy) vs a permutation null — biologically impossible adjacencies flag systematic errors.

**Caveat (load-bearing):** rare infiltrating cells are legitimately spatially incoherent. Use these to **down-weight**, never to hard-delete.

<a name="3f-anchored"></a>
### 3f. Reference-anchored specificity (needs a matched profile, not per-cell labels)

Between purely internal and fully external: these need a **matched scRNA/snRNA reference profile** but not per-cell ground truth. From Salas et al. 2025:

- **NCP** (Negative Co-expression Purity) — % of gene pairs **not** co-expressed in the reference that stay non-co-expressed in situ. 0–1, higher = better specificity; commercial platforms mean **>0.8**.
- **NMP** (Negative Marker Purity) — % of reads from "negative markers" (genes a type should NOT express) found in cells assigned to that type. Lower cleaner; the most direct readout of segmentation-driven mis-assignment, and the recommended metric for **comparing segmentations**.

> Terminology trap (kept from `cell-annotation-qc.md`): **NCP** (negative *co-expression* purity, reference-anchored) ≠ the 10x **negative *control-probe* rate** (background binding). Different things; keep separate.

<a name="3g-marker-fidelity"></a>
### 3g. Marker-program fidelity (needs marker sets, not ground truth)

A powerful family the Zhu 2026 benchmark formalizes: given curated marker/program sets, you can measure whether the **labels track the biology without per-cell ground truth** — so these run on real deliveries, not just benchmarks. They are the §2.5 fidelity metrics in their reference-free form:

- **ssGSEA enrichment** (GSVA) — score each cell against each type's marker/pathway set.
- **AUC-ROC on enrichment** — per type, how well the method's predicted-positive cells are ranked above the rest by that type's enrichment score (1 = perfect ranking, 0.5 = chance). A clean, threshold-free "do the cells labeled *T cell* actually carry the highest T-cell-program score?" This **generalizes the on-target specificity** in [`compare_rctd_tacco.py`](../../scripts/split/compare_rctd_tacco.py) into a proper ranking metric.
- **Cohen's d** — standardized target-vs-non-target separation on the enrichment score (effect size, not just AUC).
- **Marker-gene overlap** — overlap of a label's data-derived top-N DE genes with its curated markers.

Use these when you have marker sets but no truth (the usual the cluster case), as the **fidelity** complement to the **purity** metrics in [§3c](#3c-purity): purity asks "is the cell contaminated?", fidelity asks "does the label match the cell's program?".

<a name="3h-signal-qc"></a>
### 3h. Signal-quality QC (annotation-independent; SpatialQM)

Upstream of annotation quality is *signal* quality: does the section carry real, structured biology above background? These are **label-independent** - they score the data, not the labels - and are ported from the **SpatialQM** R package (Center for Spatial OMICs; NOT the identically-named Plummer/Spatial-Touchstone tool - see §9). Implemented in `analysis/signal_qc.py`, **squidpy-compatible** and fast:

- **Moran's I gene-vs-control** (`getMorans`) - the flagship. Real genes are **spatially autocorrelated** (Moran's I well above 0); negative-control probes are spatially random (I ≈ 0). The **separation** between the gene and control I distributions (and the fraction of genes above the controls' 95th percentile) is a reference-free spatial signal-quality score. Computed via squidpy's `sq.gr.spatial_autocorr` on the `sq.gr.spatial_neighbors` graph - parallel, scales to 1e6 cells.
- **Signal-to-noise** (`getMeanSignalRatio` / `getMaxRatio`) - per-cell log2(mean gene expr / mean negative-control expr) and log2(max gene / control mean); a cell whose signal doesn't clear the control floor is untypeable.
- **Sparsity** (`getSparsity`) - zero fraction of the count matrix (dataset + per cell).
- **Detection entropy** (`getEntropy`) - per-cell normalized Shannon entropy over detected genes; low → one gene dominates (rRNA / probe artifact / extreme specialization) - the continuous version of the "max single-gene fraction" flag in [`cell-annotation-qc.md` §5](cell-annotation-qc.md).
- **Transcript density per area** (`getTxPerArea`) - counts / cell area (needs `cell_area`).

These gate whether annotation should be trusted *at all* on a section (they sit alongside Layer 0-3 of the spatial funnel); run them before the label-quality metrics above. `signal_qc.run_signal_qc(adata)` returns the headline.

---

<a name="4-batteries"></a>
## 4. What we standardize on — default metric batteries

The authoritative recommendation: don't compute everything; compute the right battery for the scenario. Each battery is ordered (compute the first few always; add the rest if cheap/available).

**A. Spatial delivery, no reference (the typical Xenium/CosMx report):**
1. **Panel adequacy** (can the panel resolve the type at all — [`panel_check.py`](../../src/spatialscribe/analysis/panel_check.py))
2. **MECR + CRISP purity + PMP** (reference-free contamination)
3. **Marker-score margin + entropy** (per-cell confidence)
4. **Spatial coherence / PAS** (down-weight only)
5. **Subsampling label-stability**
6. *(if ≥2 replicates)* **Inter-sample consistency (ISC)**
7. *(reporting)* section **annotatability headline**: % PASS / WARN / abstain + top reasons

**B. Spatial with a matched snRNA reference:** battery A **plus** RCTD `spot_class` + weight-margin, **NCP/NMP**, label-transfer posterior margin. (This is the highest-quality regime; reference choice matters more than method — see [`annotation-method-selection.md`](annotation-method-selection.md).)

**C. scRNA-seq, no ground truth:** **silhouette + 2label-silhouette + neighborhood purity** (scTypeEval internal) + **ISC across samples** + **marker-score margin** + **AUC-ROC on marker enrichment** (§3g fidelity). Add Davies–Bouldin / Ward-PropMatch to catch over-clustering.

**D. Benchmark with ground-truth labels (the Zhu 2026 template):** **macro-F1 + per-class F1 + balanced accuracy** (classification) + **ARI/AMI + ECS** (partition, ECS when granularity differs) + **AvgBIO** (structure conservation) + **≥1 fidelity metric** (AUC-ROC on marker enrichment / marker-gene overlap) + **hierarchical accuracy** (so subtype-vs-lineage errors are scored differently) + **composition/proportion accuracy** if the question is compositional/dynamic; report **scalability** (time/mem) when method choice must scale. Never rank on accuracy alone; never on a single metric.

**Three cross-cutting rules (apply to every battery):**
1. **Margin > absolute score** — a high score not separated from the runner-up is ambiguous, not confident.
2. **Keep ≥1 method-independent metric** — if you validate an RNA annotation with a reference-derived silver label, that metric favors the reference arm by construction; MECR/CRISP (marker-only) is the independent check.
3. **Calibrate/normalize locally** — thresholds don't transfer across tissue/panel/platform (10x's own guidance); prefer distributional/relative variants; min-max scaling across a few methods inflates tiny gaps (annotate raw values).

---

<a name="5-master-table"></a>
## 5. Master metric table (quick reference)

Regime: **E** = external (needs truth) · **I** = internal (reference-free) · **A** = reference-anchored (needs a matched profile, not labels). Needs: L=labels, R=reference atlas, M=marker sets, S=biological replicates (samples), XY=spatial coords, Z=molecule/z data.

| Metric | Regime | Needs | Measures | Dir | In code |
|---|---|---|---|---|---|
| Accuracy / balanced accuracy | E | L+truth | per-cell correctness | high | `eval_metrics.py` |
| Precision / recall / **macro-F1** (+ micro=acc, weighted) | E | L+truth | per-class quality; rare-type sensitivity | high | `eval_metrics.py` |
| **Hierarchical accuracy** · subtype-level accuracy | E | L+truth+hierarchy | partial credit along the lineage tree | high | `eval_metrics.py` |
| **ARI** / AMI / NMI / V-measure | E | L+truth | partition agreement | high | `eval_metrics.py` |
| **ECS** (Element-Centric Similarity) | E | L+truth | resolution-agnostic, cell-level partition agreement | high | `eval_metrics.py` |
| **AvgBIO** = mean(ARI,NMI,ASW_scaled) · **ASW** | E | L+truth+embed | biological-structure conservation | high | — |
| **Composition / proportion accuracy** | E | L+truth | predicted vs true cell-type proportions (hard-label global) | high | `eval_metrics.py` |
| **Deconvolution props** (R²/JSD/RMSE, OpenProblems) | E | proportion matrices | soft per-spot proportions vs known mixtures | R²↑/JSD↓ | `eval_metrics.py` |
| **Robustness to reference volume** | E | L+truth+R | accuracy decay as reference shrinks | stable | — |
| Scalability (runtime · peak memory) | — | — | practical cost, not quality | lower | bench scripts |
| SingleR delta · CellTypist prob · RCTD spot_class | E/I | R | native confidence proxy | high | `methods.py` join |
| **Silhouette / 2label-silhouette** | I | L+embedding | cohesion vs separation | high | `eval_metrics.py` |
| **Neighborhood purity** | I | L+embedding | kNN label agreement | high | `eval_metrics.py` |
| Orbital-medoid · Ward-PropMatch · Avg-similarity | I | L+embedding | cluster–label concordance | high | `eval_metrics.py` |
| Davies–Bouldin · Calinski–Harabasz | I | L+embedding | separation | DB low / CH high | — |
| **Inter-sample consistency (ISC)** | I | L+S | cross-replicate reproducibility | high | `eval_metrics.py` |
| **MECR** | I | M | cross-lineage co-expression | low | `purity.py`, SPLIT compare |
| **CRISP purity** | I | M | per-cell lineage purity | high | `purity.py` |
| **PMP** | I | M | own-marker transcript fraction | high | `purity.py` |
| On-target specificity / retention | I | M | marker specificity after purification | high | `compare_rctd_tacco.py` |
| **AUC-ROC on marker enrichment** · Cohen's d | I | L+M | do labels track the marker-program enrichment | high | `eval_metrics.py` |
| **Marker-gene overlap** (top-N DE vs markers) | I | L+M | labels reproduce known markers | high | — |
| Marker-score **margin** · entropy | I | L+M | per-cell ambiguity | margin high / entropy low | `annotate.py` |
| Conformal set size | I | L+scores | abstention (ambiguous/novel) | =1 | `eval_metrics.py` |
| Subsampling label-stability | I | counts | detection-robustness | low flip | `annotate.py` |
| Spatial coherence / **PAS** | I | L+XY | spatial label agreement | high / PAS low | `spatial.py` |
| Neighborhood enrichment sanity | I | L+XY | impossible adjacencies | — | `spatial.py` |
| **Moran's I gene-vs-control** | I | XY+control | spatial signal vs background noise | high separation | `signal_qc.py` |
| Signal-to-noise · sparsity · detection entropy · tx/area | I | counts (±control/area) | label-independent signal quality | mixed | `signal_qc.py` |
| **NCP** / **NMP** | A | R+M | reference-anchored specificity | NCP high / NMP low | — |
| Panel adequacy (coverage/confusable) | I | M+panel | can the panel resolve the type | green | `panel_check.py` |
| ovrlpy VSI / 3D-coherence | I | Z | vertical doublets/overlap | high | `qc.py` join (H5) |

---

<a name="6-implementation"></a>
## 6. Implementation status in SpatialScribe

> **Single source of truth.** The metric battery (`eval_metrics.py` + `signal_qc.py`) now lives in the standalone package **[spatial-anno-metrics](https://github.com/p-gueguen/spatial-anno-metrics)** (`pip install spatial-anno-metrics`, MIT). In SpatialScribe, `analysis/eval_metrics.py` and `analysis/signal_qc.py` are **thin re-export shims** over that package (declared as an editable path dep in `pixi.toml`), so the `spatialscribe.analysis.*` import paths and the funnel wiring are unchanged. The function inventory below is the package's; edit metrics in the package, not the shims.

**Contamination + confidence + spatial** (reference-free): `purity.py` (MECR, CRISP, PMP, NMP/NCP guards), `annotate.py` (marker margin, entropy path, subsampling stability, penalized-posterior confidence + abstention), `panel_check.py` (panel adequacy + identifiability AUC), `spatial.py` (spatial coherence, PAS, neighborhood sanity), and `scripts/split/compare_rctd_tacco.py` (on-target specificity, marker purity, on-target retention, MECR panel).

**Metric battery** — `spatial_anno_metrics.eval_metrics` (**TDD**; re-exported by `analysis/eval_metrics.py`). Covers exactly the families this doc exposed:

```python
# reference-free internal validity (scTypeEval family) on an embedding + labels
internal_validity(adata, label_key="cell_type", embedding="X_pca", k=30) -> dict
    # {silhouette, silhouette_2label, neighborhood_purity, orbital_medoid,
    #  ward_propmatch, avg_similarity, integrated (geometric mean)}
inter_sample_consistency(adata, label_key, sample_key, metric="cosine") -> dict   # ISC
marker_program_fidelity(adata, label_key, marker_sets) -> dict     # §3g: per-type AUC-ROC + Cohen's d
external_scores(pred, truth) -> dict     # balanced_acc, macro/weighted/per-class F1, ARI/AMI/NMI, kappa, ECS
element_centric_similarity(a, b, alpha=0.9) -> float               # exact vectorized ECS (Gates & Ahn)
hierarchical_accuracy(pred, truth, hierarchy, partial=0.5) -> dict # partial credit along the lineage tree
composition_accuracy(pred, truth) -> dict                         # hard-label global proportions: L1 / Pearson / JSD
deconvolution_metrics(true_prop, pred_prop) -> dict               # soft per-spot proportion matrices: R2 / JSD / RMSE (OpenProblems)
annotation_quality(adata, label_key, marker_sets=None, sample_key=None) -> dict   # reference-free battery
```

`internal_validity` runs the O(N²) silhouettes on a `subsample` for scalability and falls back to a quick PCA when no embedding is present; `external_scores` is a thin `sklearn.metrics` wrapper plus the exact ECS closed form. The reference-free battery (`annotation_quality`) is **wired into `qc.run_funnel`** (guarded, non-breaking) so every labeled section's funnel headline carries `annotation_quality.internal_validity` + `marker_fidelity`; `external_scores` / `hierarchical_accuracy` / `composition_accuracy` / `deconvolution_metrics` are for benchmark scripts with ground truth. Tested against well-separated vs shuffled fixtures (clean labeling → purity≈1, high integrated; shuffled → near-chance).

**Signal-quality QC** — `spatial_anno_metrics.signal_qc` (**TDD**; re-exported by `analysis/signal_qc.py`; §3h) ports the label-independent SpatialQM metrics, squidpy-compatible: `moran_signal` (Moran's I gene-vs-control via `sq.gr.spatial_autocorr`), `signal_to_noise`, `sparsity`, `detection_entropy`, `tx_per_area`, and `run_signal_qc` (headline). Standalone (run before annotation to gate whether a section is trustworthy at all); not auto-wired into the funnel because Moran's I over the full panel is heavier than the per-cell metrics.

> **Two scTypeEval implementations — pick the right one.** `eval_metrics.internal_validity` is an **embedding-space approximation** built for a fast, always-on `qc.run_funnel` headline: it computes the family on the PCA embedding (subsampled) and does **not** implement scTypeEval's dissimilarity spaces or its `recip_classif:Match` local metric, and its `orbital_medoid` / `avg_similarity` / `ward_propmatch` use classic embedding-space definitions rather than the exact `IntClusVal.R` ones. A **byte-exact port** (dissimilarity spaces + `recip_classif:Match` + exact IntClusVal formulas + per-sample pseudobulk, TDD-verified) lives in an internal byte-exact port (not shipped) — use it when you need scTypeEval's recommended local(recip-match)/global(pseudobulk-cosine) combo or paper parity; use `eval_metrics.py` for the production funnel. Empirically the two **rank annotations/normalizations identically but differ in absolute scale**, so the approximation is fine for a QC headline. See an internal normalization benchmark (not shipped) — the **normalization × annotation-quality benchmark** on the 100k breast Xenium 5K section, which anchors to this catalog and reports both implementations plus the implemented contamination/confidence/spatial families across 9 normalizations. Headline finding: normalization is **~invariant** for the contamination metrics (§3c, they read raw counts) but has a **strong** effect on internal separability (§3a/b), the abstention rate (§3d), and above all the clustering the annotation is built on.

---

<a name="7-sctypeeval"></a>
## 7. scTypeEval — what it adds and how to fold it in

[**scTypeEval**](https://github.com/carmonalab/scTypeEval) (Cancer Systems Immunology Lab / Carmona) is an R package for **reference-free** assessment of cell-type annotations, benchmarked internally across **31 datasets**. Inputs: count matrix + metadata / Seurat / SCE; feature selection by HVGs (scran), marker genes, or a custom list; strategies incl. pseudobulk + distance (euclidean/cosine/pearson), single-cell Wasserstein, reciprocal classifiers, optional PCA.

**Metrics it defines** (all reference-free, higher = more consistent): `Silhouette`, `2label_silhouette`, `NeighborhoodPurity`, `ward_PropMatch`, `Orbital_medoid`, `Average_similarity`. **Recommended combo from its benchmark:** *local consistency* (silhouette in reciprocal-classification match space → catches **over-partitioning**) + *global consistency* (2label-silhouette in pseudobulk-cosine space → catches **signal degradation**) + a geometric-mean **integrated score**.

**Why fold it in:** it fills the two families `cell-annotation-qc.md` lacks — general internal cluster-validity applied to *annotation* quality, and the **inter-sample-consistency** idea (use replicates as a truth-free quality signal). Both are directly useful for the cluster reports that have labels + replicates but no atlas.

**How to add (two options):**
1. **Port the metric definitions** into `analysis/eval_metrics.py` (§6) — they are standard distance/clustering computations (silhouette, kNN purity, medoid/Ward concordance); no new heavy dependency, stays in the Python stack. Recommended.
2. **Call the R package** as a subprocess (same pattern as the annotation subprocesses) if we want its exact reciprocal-classifier / Wasserstein internals. Heavier; only if we need byte-parity with the paper.

**Caveats / attribution:** scTypeEval is a **GitHub tool with an internal 31-dataset self-benchmark**; I did not find a peer-reviewed preprint (as of Jul 2026) — cite the repo, flag as `[tool / self-benchmark]`, and treat the specific recommended-combo claims as developer-benchmarked, not independently validated. It is **scRNA-oriented**; on targeted spatial panels, compute its metrics on the PCA/latent embedding (not raw sparse counts) and expect the sparsity caveats (§3a) to apply.

---

<a name="8-rules"></a>
## 8. Cross-cutting rules & anti-patterns

- **No ground truth is the default.** Reach for §3, not §2, on real deliveries. Report external metrics only on benchmarks.
- **No single metric is an oracle.** Every real evaluation is a *battery* (§4). Accuracy alone hides rare-type failure; silhouette alone can be gamed by over-merging; MECR alone can be gamed by deleting transcripts (pair it with on-target retention — the SPLIT comparison does exactly this).
- **Accuracy ≠ biological fidelity** (Zhu 2026): always pair accuracy/F1 with a fidelity metric in benchmarks.
- **Margin > absolute score** for confidence.
- **Keep one method-independent metric** in every battery (marker-only MECR/CRISP) to avoid circularity with reference-anchored ones.
- **Spatial coherence & similar metrics penalize rare biology** — down-weight, never delete.
- **Thresholds are local** — tune per tissue/panel/platform; prefer distributional variants; annotate raw values (don't min-max a 3-method bar chart).
- **Granularity must match** for partition metrics (ARI/NMI) — fine-vs-coarse is low by construction, not by error. Fixes: use **ECS** (resolution-agnostic) for the partition score and **hierarchical accuracy** (§2.6) so a right-lineage/wrong-subtype call earns partial credit instead of a flat miss.
- **Internal metrics are consistency, not correctness** — a confidently-wrong-but-self-consistent annotation (e.g. a whole cluster mislabeled) scores well on silhouette/ISC. Internal metrics catch *incoherence*, not systematic bias; only a reference or orthogonal modality (protein) catches the latter.

---

<a name="9-references"></a>
## 9. References

**This catalog builds on the research docs** — see them for the spatial-funnel detail and verified sources:
- [`cell-annotation-qc.md`](cell-annotation-qc.md) — the six-layer spatial QC funnel; source for Layers 0–6, contamination metrics (MECR/CRISP/PMP/NCP/NMP), confidence/abstention, spatial coherence, and all verified primary sources (Salas *Nat Methods* 2025; Plummer/Spatial-Touchstone *Nat Biotechnol* 2025; Mitchel *Nat Genet* 2025; SPLIT *Nat Methods* 2026).
- [`annotation-method-selection.md`](annotation-method-selection.md) — which annotator to run; source for the accuracy≠fidelity finding and the 20-method spatial benchmark (Zhu et al., bioRxiv 2026; AvgBIO / program concordance).
- [`annotation_qc_thresholds.yaml`](annotation_qc_thresholds.yaml) — the tunable pass/warn/fail numbers consumed by the code.

**scTypeEval (the metrics folded in here):**
- Carmona lab, `scTypeEval` — reference-free annotation-quality assessment (silhouette / 2label-silhouette / neighborhood purity / Orbital-medoid / Ward-PropMatch / average-similarity; inter-sample consistency; local/global/integrated framing). GitHub: [github.com/carmonalab/scTypeEval](https://github.com/carmonalab/scTypeEval) `[tool / 31-dataset self-benchmark]`. Lab: [github.com/carmonalab](https://github.com/carmonalab).

**Classic metric definitions (external + internal validity):**
- **ARI** — Hubert & Arabie, *J. Classification* 1985. **AMI/NMI** — Vinh, Epps & Bailey, *JMLR* 2010.
- **Silhouette** — Rousseeuw, *J. Comput. Appl. Math.* 1987. **Davies–Bouldin** 1979; **Calinski–Harabasz** 1974.
- **F1 / precision / recall / balanced accuracy / Cohen's κ** — standard classification metrics (`sklearn.metrics`).
- **spARI** — spatially-aware ARI (spatial-domain benchmarks).
- **Deconvolution proportion metrics** — **OpenProblems `task_spatial_decomposition`** ([github.com/openproblems-bio/task_spatial_decomposition](https://github.com/openproblems-bio/task_spatial_decomposition)): `r2` = `r2_score(uniform_average)` over cell types; `jsd` = mean per-cell-type `scipy.spatial.distance.jensenshannon` (axis=0). Ported into `eval_metrics.deconvolution_metrics` (+ RMSE, per-spot JSD).

**Benchmark-derived metrics — Zhu et al. 2026 (the paper read for this update):**
- Zhu, Hu, … Meltzer & Zhou, "Benchmarking cell type annotation in spatial transcriptomics: resolving cellular hierarchies, biological fidelity, and dynamic cell states," bioRxiv 2026-06-16, DOI [10.64898/2026.06.16.732716](https://www.biorxiv.org/content/10.64898/2026.06.16.732716v1.full) · local `sources/2026.06.16.732716v1.full.pdf`. Source for the **ECS**, **AvgBIO = mean(ARI,NMI,ASW_scaled)**, **hierarchical/subtype accuracy**, **composition/proportion accuracy**, **reference-volume robustness**, **AUC-ROC-on-enrichment / Cohen's d / marker-gene-overlap** fidelity metrics folded in here, and the accuracy≠fidelity finding.
- **ECS** — Gates & Ahn, "Element-centric clustering comparison," *Sci. Rep.* 2019 (resolution-agnostic, affinity-matrix formulation).
- **AvgBIO / ASW-based conservation** — Luecken et al., "Benchmarking atlas-level data integration" (scIB), *Nat. Methods* 2022.
- **ssGSEA** — Barbie et al., *Nature* 2009; **GSVA** — Hänzelmann et al., *BMC Bioinformatics* 2013. **Cohen's d** — Cohen 1988.

**Signal-quality QC (SpatialQM, §3h — the metrics ported into `signal_qc.py`):**
- **SpatialQM** — Center for Spatial OMICs, [github.com/Center-for-Spatial-OMICs/SpatialQM](https://github.com/Center-for-Spatial-OMICs/SpatialQM). R package; source for Moran's I gene-vs-control (`getMorans`), signal-to-noise (`getMeanSignalRatio` / `getMaxRatio`), sparsity, detection entropy, tx-per-area, and its own `getMECR`. `[tool]`
- ⚠ **Name collision:** this is **NOT** the *other* "SpatialQM" (Plummer et al., Spatial Touchstone, *Nat Biotechnol* 2025, GEO GSE277080) cited in `cell-annotation-qc.md` §17 for tissue-specific baselines. Two different tools, same acronym — keep them distinct.
- Moran's I — Moran 1950; computed here via squidpy `sq.gr.spatial_autocorr` (Palla et al., *Nat Methods* 2022).

**Spatial / contamination metric sources** (full detail + DOIs in [`cell-annotation-qc.md` §17](cell-annotation-qc.md)): Salas et al. *Nat Methods* 2025 (NCP, NMP, 3D-coherence); Hartman & Satija 2024 (MECR); Center for Spatial OMICs (CRISP); Cable et al. *Nat Biotechnol* 2022 (RCTD); Marconato et al. *Nat Methods* 2026 (SPLIT); STEAM (spatial coherence / PAS).

> **Provenance note.** The spatial/contamination/confidence metrics and their evidence are carried over (with cross-links, not duplicated) from `cell-annotation-qc.md`, which was built from an adversarial-verification research pass. The **new** material here — the external-validation family (§2), the internal cluster-validity + ISC family (§3a–3b), the marker-program fidelity family (§3g), the decision guide (§1), the batteries (§4), the scTypeEval integration (§7), and the ECS / AvgBIO-form / hierarchical / composition / reference-volume-robustness metrics (§2.5–2.6) — is synthesized from standard metric definitions, the scTypeEval repository, and a **direct reading of the Zhu et al. 2026 benchmark Methods** (`sources/2026.06.16.732716v1.full.pdf`: exact ECS/ASW/AvgBIO/AUC-ROC/Cohen's-d formulas transcribed from its Methods section). Classic-metric citation years are working knowledge and should be link-verified before any external publication.
