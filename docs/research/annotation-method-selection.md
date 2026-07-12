# Which tool best annotates Xenium? — a method-selection evidence base

**Status:** design reference / evidence baseline · **Scope:** cell-type annotation *methods* for 10x Xenium (imaging spatial), preprint-weighted, last ~12 months · **Companion to** [`cell-annotation-qc.md`](cell-annotation-qc.md) (which covers *QC of* annotations; this covers *which annotator to run*). · **Feeds:** the future `analysis/annotate.py` method stack.

> **Provenance.** Built from a deep-research adversarial-verification pass (6 angles, 20 sources, 94 claims → **23 confirmed / 2 refuted → 9 synthesized findings**) plus a structured OpenAlex primary-literature sweep. Full evidence + sources: [`annotation-method-research-appendix.md`](annotation-method-research-appendix.md). Preprints are flagged **[preprint]**; developer self-benchmarks are flagged **[self-bench]** (structurally favor the introduced method — weight accordingly).

---

## TL;DR — the four things that actually matter

1. **There is no universal "best" annotator for Xenium.** The single most robust, cross-cutting result is that **reference choice dominates method choice** — a good matched (ideally single-nucleus) reference moves accuracy more than swapping annotators does. *Invest in the reference before agonizing over the method.*
2. **Accuracy ≠ biological fidelity.** The most comprehensive recent benchmark shows that **high classification accuracy does not guarantee** preserved cellular relationships or coherent downstream pathway/enrichment analysis. Ranking methods on accuracy alone is misleading for a tool whose output feeds biology.
3. **Among battle-tested tools, SingleR and RCTD are the safe defaults.** SingleR wins the only dedicated peer-reviewed Xenium benchmark; RCTD (in doublet mode, + snRNA-seq reference, + SPLIT purification) is the leading recipe for Xenium's contamination/spillover problems.
4. **Foundation models are not (yet) the answer for Xenium.** No foundation model has been shown to win head-to-head on native Xenium cell annotation, and in **zero-shot** (how a self-serve app would use them, without fine-tuning) they can be *outperformed by simpler methods*. Don't bet the app on them.

---

## The evidence, by paradigm

### A. Classical reference label-transfer — the only dedicated, peer-reviewed Xenium benchmark

**Cheng, Jin, Smyth & Chen, *BMC Bioinformatics* 2025** (10.1186/s12859-025-06044-0; WEHI/Smyth lab). 5 reference tools on two HER2+ breast-cancer Xenium samples, paired 10x-Flex snRNA-seq reference, manual marker annotation as ground truth:

| Rank | Method | Consistency (2 samples) | Speed |
|---|---|---|---|
| **1** | **SingleR** | **73.5% / 65%** | **fastest (~25 min)** |
| 2 | **RCTD** | 69.9% / 55% | ~2 h |
| 3–5 | Azimuth · scPred · scmapCell | all **< 60%** | scPred ~3 h |

> **Verdict:** "SingleR was the best performing reference-based tool for the Xenium platform — fast, accurate, easy to use." **Scope caveat:** ONE tissue, 2 samples, only 5 R-based tools — it did **not** test cell2location, Tangram, CellTypist, Seurat TransferData, Symphony, scANVI, or any FM/spatial-context method. "SingleR is best" holds only within that tested set.

### B. Broadest recent multi-method benchmark — different winners (full text now read)

**Zhu, Hu, … Meltzer & Zhou (Vanderbilt), "Benchmarking cell type annotation in spatial transcriptomics: resolving cellular hierarchies, biological fidelity, and dynamic cell states," bioRxiv 2026-06-16 [preprint]** (DOI 10.64898/2026.06.16.732716). **20 methods** (17 reference-based + reference-free/spatial; incl. SingleR, Seurat, RCTD, scANVI, Tangram, cell2location, DestVI, CARD, TACCO, Spatial-ID, GraphST, SpaGCN, BANKSY, and FMs scGPT/scCello/Nicheformer/SToFM), across **4 datasets incl. Xenium** (MERFISH mouse spinal cord · Open-ST human metastatic lymph node · Stereo-seq · **Xenium mouse kidney injury–repair**). **Independent** (methods are not the authors' own). 6-dimension eval: classification (acc/macro-F1/precision/recall), subtype-level, hierarchical, marker-gene concordance, scalability (time+peak mem), downstream biological interpretability (ssGSEA/ARI/AMI/AvgBIO).

- **Aggregate top-3:** **scANVI, Seurat, TACCO** — but context-dependent. Per-scenario: **scANVI** best for fine-grained **subtypes**; **Spatial-ID** highest overall/dominant-type accuracy + ARI (but weak macro-F1 → misses rare types); **RCTD** highest **recall**; **TACCO/Tangram** highest **precision**; **Seurat** highest **AvgBIO** (0.451, biological-structure preservation).
- **Xenium-specific (now confirmed — I read the full PDF):** the Xenium dataset is **mouse kidney ischemia-reperfusion, 12 sections / 6 stages, 300-gene panel, ~1.37 M cells, matched snRNA-seq reference**. **TACCO consistently ranked top** here (Macro-F1/accuracy/precision/recall/ARI across the injury–repair course); DestVI and Seurat also recovered the injury biology; **CARD failed on rare injury states**; **GraphST ran out of memory** on 1.37 M cells; accuracy dropped most during **acute injury**. Note this is a deliberately *hard, dynamic-state* Xenium case, not a clean tumor atlas — generalize cautiously.
- **FMs are more nuanced than "they lose":** zero-shot / minimal adaptation underperforms, but **partial fine-tuning gives marked gains** — **scGPT (partial fine-tuning) had the *highest* accuracy (0.672)** on the lymph-node TME (scCello 0.632); compartment-specific winners (SToFM→tumor, Nicheformer→CAFs, SingleR→fibroblasts/plasma). Head-tuned FMs recover broad organization but with diffuse niche boundaries.

> **Two load-bearing conceptual findings:** (1) **high accuracy ≠ biological fidelity** — rankings by *program concordance* (do labels retain expected transcriptional programs?) frequently **differ** from accuracy rankings, so reported accuracy "may substantially overestimate" real quality; (2) **robustness to domain shift may matter more than the underlying statistical framework** — and closed-label-space methods force novel/stage-specific states into existing categories (→ use abstention, see the QC doc's `Novel/unknown`).

### C. Contamination / sparsity / spillover — the leading recipe

**Bilous et al. (SPLIT), *Nature Methods* 2026** (10.1038/s41592-026-03089-8; bioRxiv 2025.04.23.649965). 40+ breast/lung tumor sections:

- **RCTD is the annotation engine**, "highly effective," its **doublet mode** explicitly modeling "mixed signals arising from segmentation errors, overlapping cells or transcript spillover."
- **snRNA-seq references markedly improve annotation** and sharpen diffusion quantification.
- **SPLIT** (post-RCTD, spatially-aware purification) gives "the best cell-type separation and highest cosine similarity with Chromium snRNA-seq," recovering signal (e.g. T-cell exhaustion near tumor) otherwise buried in contamination; **ProSeg segmentation + SPLIT** is the best overall combination.

**Caveat:** the RCTD endorsement is a single-group preference ("we and others have found"), and SPLIT's "best" is a **[self-bench]** — the independent BMC benchmark ranks SingleR *above* RCTD. So: use RCTD **when you need doublet-mode contamination handling + SPLIT**; use SingleR when you want the fastest accurate single-label call.

### D. Deconvolution (spot-based context) — reference stability is the story

**Spotless benchmark, *eLife* 2024** (Saeys lab; 11 methods). Spot-based (Visium/simulated), *not* native single-cell Xenium — corroborating context only:

- **RCTD and cell2location are the top two** across all metrics (RMSE/AUPR/JSD); then SpatialDWLS, stereoscope, MuSiC.
- **"The choice of reference dataset has the largest impact"**; methods modeling technical variability (**cell2location, RCTD**) are **most stable across references** (cell2location #1 in cross-reference stability). ← this is the origin of TL;DR #1.

### E. Foundation models — not consistently better, and untested on native Xenium

- **scEval, *Advanced Science* 2026** (Liu et al.; 10 FMs × 8 tasks): "single-cell FMs may **not** consistently excel [over] task-specific methods." For annotation the best FMs are **CellPLM, scGPT, Geneformer** (comparable); top overall by performance+accessibility: scGPT, Geneformer, CellFM. **Major caveat:** scEval is dissociated scRNA-seq — **no FM shown to win on native Xenium**.
- **Kedzierska et al., *Genome Biology* 2025** (zero-shot FM evaluation): scGPT & Geneformer in **zero-shot can be outperformed by simpler methods** and face reliability issues. A self-serve app uses FMs zero-shot → this is the operative regime.
- **Fine-tuning is the deciding factor (Zhu et al. 2026, §B):** in that independent spatial benchmark, **zero-shot / head-tuned FMs underperformed**, but **partial fine-tuning let scGPT reach the top accuracy (0.672) on one dataset** (Open-ST lymph node — *not* the Xenium dataset, where TACCO won). So the honest statement is: **no FM has been shown to win on native Xenium, and any FM win requires per-dataset fine-tuning** — infeasible for a zero-shot self-serve app.
- **Refuted & excluded:** a claim that fine-tuned **Nicheformer** systematically tops other FMs on spatial tasks failed verification (vote 1-2) — not asserted.

### F. Newer specialist methods (promising, but self-benchmarked)

| Method | Claim | Criterion it targets | Flag |
|---|---|---|---|
| **STAMapper** (heterogeneous GNN, *Genome Biology* Oct 2025) | Best accuracy on **75/81 spatial datasets**, beats scANVI/RCTD/Tangram; robust on **<200-gene panels** (51.6% vs scANVI 34.4% at 0.2 downsampling) | Sparsity robustness | **[self-bench]**; Xenium **not** in the main 81-dataset comparison (MERFISH/seqFISH/osmFISH/STARmap/CosMx) |
| **RankMap** (rank-based mapping, bioRxiv Mar 2026) | Competitive-or-superior accuracy vs SingleR/Azimuth/RCTD with **consistently lower runtime** on Xenium/MERFISH/Stereo-seq | Scalability to 1e5–1e6 | **[self-bench]** (same first author as the SingleR/BMC benchmark — effectively "SingleR-style rank mapping, faster") |
| **Marker-gene scoring** (Princeton, bioRxiv Jan 2026) | Competitive/better than RCTD/cell2location **for RARE cell types** in simulated Xenium pseudo-spots (0.72 vs RCTD 0.63) | Rare-type sensitivity | The broader "marker scoring beats deconvolution universally" claim was **REFUTED** — advantage is rare-type-specific only |

---

## Winners per criterion (the direct answer)

| Criterion | Best-supported choice | Evidence & confidence |
|---|---|---|
| **(a) Accuracy vs ground truth** | **SingleR** (peer-reviewed, Xenium breast) · **TACCO** (Xenium mouse-kidney, confirmed) · **scANVI / Seurat / TACCO** (aggregate across platforms); scANVI best for subtypes, Spatial-ID best for dominant-type accuracy | BMC 2025 (high) · Zhu et al. bioRxiv 2026 (high, full text read) |
| **(b) Robustness to sparsity / segmentation contamination / spillover** | **RCTD doublet-mode + snRNA-seq reference + SPLIT** (contamination) · **STAMapper** (sparse <200-gene panels) | Nat Methods 2026 (high) · Genome Biol 2025 [self-bench] (medium) |
| **(c) Confidence / uncertainty calibration** | **OPEN — no method demonstrated best-calibrated on Xenium.** Field only warns "accuracy ≠ downstream fidelity." Use QC-layer confidence (see [`cell-annotation-qc.md`](cell-annotation-qc.md) §10) + abstention. | open question (see appendix) |
| **(d) Scalability to 1e5–1e6 cells** | **SingleR** (fast, proven) · **RankMap** (competitive accuracy, lowest runtime) | BMC 2025 (high) · bioRxiv 2026 [self-bench] (medium) |

---

## Recommendation for SpatialScribe (`annotate.py`)

1. **Default engine: RCTD (doublet mode) with a good matched reference**, because it *simultaneously* gives the annotation **and** the contamination/confidence signals the QC layer needs (`spot_class`, secondary weights → Layer 3/5 of the QC doc), and it pairs natively with **SPLIT** (already in the app's toolkit). Offer **SingleR** as the fast alternative / cross-check (it's the peer-reviewed Xenium winner and the natural consensus partner).
2. **Reference quality is the highest-leverage knob.** Bake in a reference-selection step (prefer **single-nucleus**, tissue-matched); surface reference provenance to the user. This beats any method swap.
3. **Multi-method consensus** (marker-score + RCTD/SingleR + CellTypist + Claude) is the right design — the *disagreement* is a confidence signal, and different methods win in different regimes (marker scoring for rare types; RCTD for contamination; SingleR for speed).
4. **Do not ship a foundation-model annotator as the default.** Zero-shot FMs underperform simpler methods and are unvalidated on Xenium. Revisit only if a fine-tuned, Xenium-validated FM benchmark appears.
5. **Report confidence, but treat calibration as unsolved.** No annotator is demonstrably well-calibrated on Xenium; lean on the QC doc's abstention framework rather than trusting any single method's native probability.

---

## Caveats & open questions (don't over-read this)

- **Almost every Xenium-specific "winner" is either single-tissue (SingleR/RCTD = one HER2+ breast study) or a developer self-benchmark (STAMapper, RankMap, SPLIT).** The big multi-method rankings either exclude Xenium (STAMapper) or pool it with other platforms (the 20-method preprint). Treat all rankings as directional.
- **No dedicated calibration/uncertainty benchmark exists for Xenium annotation.** (Open.)
- **FMs have never been run head-to-head against SingleR/RCTD/scANVI on native single-cell Xenium.** (Open.)
- **Spatial-context annotators (TACCO, STELLAR, CellCharter, Banksy-assisted)** are under-tested for contamination-robustness and 1e6-cell scalability specifically. (Open.)
- **Refuted claims excluded:** (i) marker scoring universally beats deconvolution (holds only for rare types); (ii) Nicheformer tops other FMs on spatial tasks.

---

## References (preprint-weighted; [P]=preprint, [SB]=self-benchmark)

- **BMC Bioinformatics 2025** — Cheng et al., dedicated Xenium annotation benchmark (SingleR #1). [10.1186/s12859-025-06044-0](https://doi.org/10.1186/s12859-025-06044-0)
- **bioRxiv 2026-06-16** [P] — 20-method spatial annotation benchmark (scANVI/Seurat/TACCO; accuracy≠fidelity). [biorxiv 10.64898/2026.06.16.732716](https://www.biorxiv.org/content/10.64898/2026.06.16.732716v1.full)
- **Nature Methods 2026** — Bilous et al., SPLIT / RCTD-doublet + snRNA recipe. [10.1038/s41592-026-03089-8](https://doi.org/10.1038/s41592-026-03089-8) · bioRxiv [2025.04.23.649965](https://www.biorxiv.org/content/10.1101/2025.04.23.649965v1)
- **eLife 2024** — Spotless deconvolution benchmark (RCTD & cell2location top; reference dominates). [elifesciences 88431](https://elifesciences.org/reviewed-preprints/88431v1)
- **Advanced Science 2026** — scEval FM benchmark (FMs not consistently better). [10.1002/advs.202514490](https://advanced.onlinelibrary.wiley.com/doi/10.1002/advs.202514490)
- **Genome Biology 2025** — Kedzierska et al., zero-shot FM limitations. [10.1186/s13059-025-03574-x](https://doi.org/10.1186/s13059-025-03574-x)
- **Genome Biology 2025** [SB] — STAMapper GNN (sparse-panel robustness). [10.1186/s13059-025-03773-6](https://doi.org/10.1186/s13059-025-03773-6)
- **bioRxiv 2026-03** [P][SB] — RankMap (scalability). [biorxiv 10.64898/2026.03.01.708931](https://www.biorxiv.org/content/10.64898/2026.03.01.708931v1)
- **bioRxiv 2026-01** [P] — Princeton marker-scoring for rare types. [biorxiv 10.64898/2026.01.13.699379](https://www.biorxiv.org/content/10.64898/2026.01.13.699379v1.full.pdf)

Full source list + verbatim evidence: [`annotation-method-research-appendix.md`](annotation-method-research-appendix.md).
