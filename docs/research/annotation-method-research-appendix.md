# Research evidence appendix — best tool to annotate Xenium

> Auto-generated from the deep-research adversarial-verification pass that grounds [`annotation-method-selection.md`](annotation-method-selection.md). Do not hand-edit; re-run the research to refresh.

**Pass stats:** 6 angles · 20 sources fetched · 94 claims extracted → **25 verified (23 confirmed, 2 refuted)** → 9 findings.

## Summary

There is no single universal "best" cell-type annotation method for 10x Xenium; the strongest cross-cutting result is that reference choice dominates method choice, and that accuracy alone does not guarantee biological fidelity. Among classic reference-based tools, the only dedicated peer-reviewed Xenium benchmark (Cheng/Smyth et al., BMC Bioinformatics 2025, HER2+ breast cancer) ranks SingleR #1 and RCTD #2, with Azimuth/scPred/scmapCell trailing; SingleR is also the fastest. The most comprehensive and recent spatial benchmark (20 methods across four platforms including Xenium, bioRxiv June 2026) instead puts scANVI, Seurat, and TACCO at the top overall, but stresses the ranking is context-dependent and that high classification accuracy does not preserve global cellular relationships or downstream pathway coherence. For handling Xenium's core weaknesses (sparsity, segmentation contamination, transcript spillover), the leading recipe is RCTD (doublet mode) paired with a single-nucleus RNA-seq reference and post-hoc purification via SPLIT (Nature Methods 2026), which markedly improves specificity and recovers signals otherwise obscured by contamination. Foundation models (scGPT, Geneformer, CellPLM, CellFM, Nicheformer) do not consistently beat task-specific methods and have not been shown to win head-to-head on native Xenium annotation; developer self-benchmarks for newer graph/rank methods (STAMapper, RankMap) report strong accuracy and speed but lack independent validation.

## Verified findings (ranked)

### 1. [HIGH · 5 claims, all CONFIRMED (three 3-0, two 2-1)]

On the only dedicated peer-reviewed head-to-head Xenium annotation benchmark, SingleR is the best reference-based method (fast, accurate, closely matching manual annotation) and RCTD is second; Azimuth, scPred, and scmapCell all fall below 60% consistency.

> Cheng/Jin/Smyth et al. (BMC Bioinformatics 26:22, Jan 2025) benchmarked 5 reference-based tools on two 10x Xenium HER2+ breast-cancer samples with a paired 10x Flex snRNA-seq reference, using manual marker-based annotation as ground truth. SingleR ranked 1st (73.49% / 65% consistency across samples), RCTD 2nd (69.92% / 55%); Azimuth, scPred, scmapCell all <60%. SingleR was also fastest (~25 min vs RCTD ~2 h, scPred ~3 h). Verbatim conclusion: 'SingleR was the best performing reference-based cell type annotation tool for the Xenium platform, being fast, accurate and easy to use.' Key scope caveat: single tissue, 2 samples, only 5 R-based tools (NO cell2location, Tangram, CellTypist, Seurat TransferData, scANVI, or spatial-context methods), and manual annotation as imperfect ground truth.

Sources: [1](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-025-06044-0) · [2](https://link.springer.com/article/10.1186/s12859-025-06044-0) · [3](https://pmc.ncbi.nlm.nih.gov/articles/PMC11744978/)

### 2. [HIGH · 3 claims, all 3-0 CONFIRMED]

In the most comprehensive and recent spatial annotation benchmark (20 methods, reference-based + reference-free, across four platforms including Xenium), scANVI, Seurat, and TACCO consistently ranked among the top methods overall, but relative advantages were context-dependent and high accuracy did NOT guarantee biological fidelity.

> bioRxiv preprint (2026-06-16, DOI prefix 10.64898 = legitimate new bioRxiv/openRxiv prefix), 'Benchmarking cell type annotation in spatial transcriptomics: resolving cellular hierarchies, biological fidelity, and dynamic cell states.' Benchmarks 20 state-of-the-art methods across four ST datasets (MERFISH, Xenium, Open-ST, Stereo-seq) with expert-curated labels as ground truth, using accuracy/F1 plus structure-aware metrics. Verbatim: 'scANVI, Seurat, and TACCO consistently ranked among the top-performing methods, although their relative advantages were context dependent.' Critically: 'high classification accuracy did not necessarily correspond to preservation of global cellular relationships or biologically coherent downstream pathway and gene-set enrichment analyses' — i.e., accuracy-only rankings mislead for downstream/calibration use. Caveat: preprint (not peer-reviewed); full text was 403-blocked so the exact method count and whether the scANVI/Seurat/TACCO ranking holds Xenium-specifically (vs pooled across platforms) rest on search snippets; one snippet showed Spatial-ID slightly edging scANVI on raw accuracy in one config.

Sources: [1](https://www.biorxiv.org/content/10.64898/2026.06.16.732716v1.full)

### 3. [HIGH · 5 claims (Nature Methods + bioRxiv versions), CONFIRMED (four 3-0, one 2-1)]

For handling Xenium's core weaknesses (sparsity, segmentation contamination, transcript spillover), the leading approach is RCTD in doublet mode plus a single-nucleus RNA-seq reference, with post-hoc purification via SPLIT giving the best cell-type separation and highest agreement with matched snRNA-seq.

> Bilous et al., 'Resolving sensitivity, specificity and signal contamination in Xenium spatial transcriptomics' (Nature Methods, 2026-04-30; bioRxiv 2025.04.23.649965), on 40+ breast/lung tumor sections. They use RCTD as their annotation engine, describing it as 'highly effective for annotating Xenium data' with a doublet mode that 'helps account for mixed signals arising from segmentation errors, overlapping cells or transcript spillover.' They show snRNA-seq references 'markedly improve cell type annotation and enable more precise quantification of diffusion.' They introduce SPLIT (Spatial Purification of Layered Intracellular Transcripts), a post-RCTD, spatially-aware transcript-purification step: 'SPLIT achieves the best cell-type separation and the highest cosine similarity with Chromium snRNA-seq profiles,' recovering signals (e.g., T-cell exhaustion near malignant cells) otherwise obscured by contamination; combining the ProSeg segmenter with SPLIT yields the greatest overall improvement. CAVEATS: the RCTD endorsement is a single-group preference ('we and others have found'), NOT an independent ranking — the independent BMC benchmark ranks SingleR above RCTD; and SPLIT's 'best' claims are from its own method-introduction paper (self-benchmark), not third-party replication.

Sources: [1](https://www.nature.com/articles/s41592-026-03089-8) · [2](https://www.biorxiv.org/content/10.1101/2025.04.23.649965v1)

### 4. [HIGH · multiple sources (eLife 2-1 + Nature Methods 3-0 corroboration)]

Reference dataset choice, not method choice, is the dominant determinant of spatial annotation/deconvolution quality; methods that model technical variability (cell2location, RCTD) are most stable across references/platforms.

> The Spotless benchmark (eLife 88431 / PMC11126312, Saeys lab; 11 deconvolution methods) states 'the choice of the reference dataset has...the largest impact' on predictions and that 'methods that accounted for technical variability in their models, such as cell2location and RCTD, were more stable to changes in the reference dataset.' cell2location ranked #1 in cross-reference stability. Independently reinforced by the Xenium-specific SPLIT paper showing snRNA-seq reference markedly improves annotation. This is the strongest cross-cutting theme: where reference choice dominates, method choice is secondary. Caveat: Spotless is spot-based/simulated Visium + seqFISH+/STARMap, not native single-cell Xenium — the reference-dominance principle transfers but the specific method ranking is spot-deconvolution, not Xenium single-cell annotation.

Sources: [1](https://elifesciences.org/reviewed-preprints/88431v1) · [2](https://www.biorxiv.org/content/10.1101/2025.04.23.649965v1)

### 5. [HIGH · 1 claim, 3-0 CONFIRMED]

In spot-based spatial deconvolution benchmarks, RCTD and cell2location are the top two performers across all metrics (RMSE, AUPR, JSD), followed by SpatialDWLS, stereoscope, and MuSiC.

> Spotless benchmark of 11 methods (cell2location, DestVI, DSTG, RCTD, SpatialDWLS, SPOTlight, stereoscope, STRIDE, MuSiC, Seurat, Tangram): 'RCTD and cell2location were the top two performers across all metrics, followed by SpatialDWLS, stereoscope, and MuSiC.' Peer-reviewed (published eLife 2024). RELEVANCE CAVEAT: this ranking is scoped to the silver (simulated Visium spot) standards and is SPOT-based deconvolution, NOT single-cell Xenium annotation — so it is corroborating context for RCTD/cell2location strength, not a native Xenium ranking.

Sources: [1](https://elifesciences.org/reviewed-preprints/88431v1)

### 6. [HIGH · 2 claims, both 3-0 CONFIRMED]

Single-cell foundation models do NOT consistently outperform task-specific methods; for the cell-type annotation task the best FMs are CellPLM, scGPT, and Geneformer (comparable to each other), with scGPT/Geneformer/CellFM the top overall by combined performance and accessibility.

> Liu et al., 'Evaluating the Utilities of Foundation Models in Single-Cell Data Analysis' (Advanced Science 2026, the published scEval benchmark of 10 FMs over 8 tasks): 'single-cell FMs may not consistently excel than task-specific methods in all tasks, which challenges the necessity of developing foundation models.' For annotation: 'CellPLM, scGPT, and Geneformer had comparable performance...better than the other models.' Overall top FMs (performance + accessibility): scGPT, Geneformer, CellFM. MAJOR RELEVANCE CAVEAT: scEval is dissociated scRNA-seq, NOT Xenium/imaging spatial — no FM was shown to win head-to-head on native Xenium annotation. The one spatial-FM head-to-head claim (Nicheformer topping Geneformer/scGPT/UCE on spatial tasks) was REFUTED in verification. Nicheformer (Nature Methods 2025) is a legitimate transformer FM whose downstream tasks include spatial cell-type/niche/region label prediction on imaging spatial data, positioning it as a spatial-FM annotator, but its superiority over other FMs is unproven.

Sources: [1](https://advanced.onlinelibrary.wiley.com/doi/10.1002/advs.202514490)

### 7. [MEDIUM · 2 claims, both 3-0 CONFIRMED (but developer self-benchmark)]

STAMapper, a heterogeneous graph neural network transferring scRNA-seq labels to single-cell spatial data, achieves the best accuracy on 75/81 spatial datasets (significantly beating scANVI, RCTD, Tangram) and is notably robust to sparse panels (<200 genes) and downsampling.

> Chen et al., Genome Biology (Oct 2025): STAMapper 'achieves the best performance on 75 out of 81 datasets...in accuracy,' significantly more accurate than scANVI (p=2.2e-14), RCTD (p=1.3e-27), Tangram (p=1.3e-36). Robustness: on <200-gene panels at 0.2 downsampling, STAMapper median accuracy 51.6% vs scANVI 34.4% (2nd best) — directly relevant to Xenium sparsity. CAVEATS: (1) this is a developer SELF-benchmark (zhanglabtools), inherently favorable, no independent validation; (2) the 81 datasets are MERFISH/seqFISH/osmFISH/STARmap/CosMx/Slide-tags — 10x Xenium is NOT in the main benchmark, appearing only as one Discussion validation (93.91% accuracy). Treat as promising but not independently confirmed for Xenium.

Sources: [1](https://link.springer.com/article/10.1186/s13059-025-03773-6)

### 8. [MEDIUM · 1 claim CONFIRMED 3-0; the broader universal claim was REFUTED]

Marker-gene signature scoring can be competitive with or better than complex deconvolution methods specifically for RARE cell types in simulated Xenium-derived pseudo-spots, with RCTD second and cell2location/UCDBase third-fourth.

> Sun/Pritykin et al. (Princeton, bioRxiv Jan 2026): on Visium-like spots simulated by aggregating Xenium breast-cancer single cells, all methods did well for abundant types (Pearson 0.74-0.93); for rare types marker-gene scoring had the highest median correlation (0.72), RCTD 2nd (0.63), cell2location/UCDBase 3rd-4th. IMPORTANT: the broader framing that simple marker scoring universally beats RCTD/cell2location/CARD/DestVI/SPOTlight was REFUTED (vote 1-2) — the advantage is specific to rare cell types in simulated pseudo-spots, and abundant-type performance is comparable across methods. Also this is simulated spot data, not native single-cell Xenium annotation.

Sources: [1](https://www.biorxiv.org/content/10.64898/2026.01.13.699379v1.full.pdf)

### 9. [MEDIUM · 2 claims, both 3-0 CONFIRMED (developer self-benchmark)]

RankMap, a new rank-based reference-mapping annotator, achieves competitive-or-superior accuracy versus SingleR, Azimuth, and RCTD while consistently reducing runtime, benchmarked on Xenium/MERFISH/Stereo-seq — addressing the scalability criterion.

> Cheng et al. (Duke-NUS, bioRxiv March 2026; GitHub jinming-cheng/RankMap + Bioconductor package): benchmarked on five ST datasets (Xenium, MERFISH, Stereo-seq) + two single-cell datasets vs SingleR/Azimuth/RCTD. 'RankMap achieved competitive or superior annotation accuracy while consistently reducing runtime...particularly for large spatial datasets.' CAVEAT: developer self-benchmark, not independent; 'competitive' signals accuracy parity with runtime/scalability as the real advantage. Relevant primarily to the scalability-to-large-datasets criterion, not as an independently validated accuracy winner.

Sources: [1](https://www.biorxiv.org/content/10.64898/2026.03.01.708931v1)

## Caveats (verbatim)

SCOPE OF SOURCES: Almost every Xenium-specific "winner" comes from either (a) a single-tissue benchmark — the SingleR/RCTD result is from ONE tissue (HER2+ breast cancer), 2 samples, only 5 R-based tools, with manual annotation (not orthogonal ground truth) as reference — or (b) developer self-benchmarks (STAMapper, RankMap, SPLIT), which structurally favor the introduced method and lack independent replication. The peer-reviewed dedicated Xenium annotation benchmark (BMC/Smyth 2025) did NOT test cell2location, Tangram, CellTypist, Seurat TransferData, Symphony, scANVI, or any foundation-model / spatial-context method, so "SingleR is best" holds only within that tested panel. PLATFORM MISMATCH: The two big head-to-head rankings that DO include many methods are either not Xenium-specific (STAMapper's 81 datasets exclude Xenium from the main comparison) or pool Xenium with MERFISH/Open-ST/Stereo-seq (the 20-method June-2026 preprint). Deconvolution benchmarks (Spotless, Princeton) are Visium spot-based or simulated pseudo-spots, not native single-cell Xenium annotation. FOUNDATION MODELS: the scEval FM benchmark is dissociated scRNA-seq, not spatial; no FM has been shown to win head-to-head on native Xenium cell-level annotation. PREPRINT/RECENCY: several key sources are unreviewed 2026 preprints (the 20-method benchmark, Princeton marker-scoring, RankMap); the 20-method preprint's full text was 403-blocked, so exact method count and Xenium-specific ranking rest on search snippets. The DOI prefix 10.64898 is the legitimate new bioRxiv/openRxiv 2026 prefix (verified), not fabricated. TWO REFUTED CLAIMS to keep in mind: (1) that simple marker-gene scoring universally beats complex deconvolution (only holds for rare types in one simulation), and (2) that Nicheformer tops all foundation models on spatial tasks. CALIBRATION GAP: no source cleanly identifies a best-calibrated-uncertainty method for Xenium; the strongest calibration-relevant result is a caution that accuracy does not track biological fidelity. PANEL/TISSUE + REFERENCE DEPENDENCE: results are repeatedly panel/tissue-dependent, and reference choice (especially matched snRNA-seq) recurrently dominates method choice.

## Open questions

- How do foundation-model annotators (scGPT, Geneformer, CellPLM, CellFM, Nicheformer) perform head-to-head against SingleR/RCTD/scANVI on NATIVE single-cell 10x Xenium annotation? No source provides this direct comparison — FM benchmarks are scRNA-seq, and the dedicated Xenium benchmark omitted FMs entirely.
- Which method is genuinely best-CALIBRATED in its confidence/uncertainty on Xenium data? No dedicated calibration/uncertainty benchmark surfaced; the field currently only warns that accuracy does not equal downstream biological fidelity.
- Does the June-2026 20-method benchmark's scANVI/Seurat/TACCO top ranking hold XENIUM-SPECIFICALLY, or is it an artifact of pooling across MERFISH/Open-ST/Stereo-seq? Full text was inaccessible, so the per-platform breakdown for Xenium is unconfirmed.
- How do spatial-context-aware annotators (TACCO, STELLAR, CellCharter, Banksy-assisted) compare to reference label-transfer specifically on robustness to segmentation contamination/spillover AND on scalability to 1e5-1e6 cells? Contamination handling is currently best-documented only for the RCTD-doublet + snRNA + SPLIT recipe, not for graph/spatial methods, and no independent large-scale (1e6-cell) runtime/accuracy comparison exists.

## Refuted & excluded (NOT asserted)

- A simple marker gene signature enrichment approach (Scanpy's default gene-set score_genes function, equivalent to Seurat AddModuleScore) performs competitively with or better than complex reference-based deconvolution methods (RCTD, cell2location, CARD, DestVI, SPOTlight, UniCell Deconvolve), especially for rare cell types, and does not require a matched scRNA-seq reference.
- For spatial cell-type/niche annotation, both the fine-tuned Nicheformer and a linear-probing model on its embedding systematically outperform other foundation-model annotators (Geneformer, scGPT, UCE) and embedding baselines (scVI, PCA) - a concrete head-to-head ranking placing Nicheformer on top among foundation models on spatial tasks.

## All sources fetched

- https://www.nature.com/articles/s41467-025-64990-y
- https://advanced.onlinelibrary.wiley.com/doi/10.1002/advs.202514490
- https://elifesciences.org/reviewed-preprints/88431v1
- https://advanced.onlinelibrary.wiley.com/doi/10.1002/advs.202518949
- https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-025-06044-0
- https://www.biorxiv.org/content/10.64898/2026.06.16.732716v1.full
- https://www.biorxiv.org/content/10.1101/2025.08.12.669903.full.pdf
- https://www.biorxiv.org/content/10.64898/2026.01.13.699379v1.full.pdf
- https://www.nature.com/articles/s41592-026-03089-8
- https://academic.oup.com/bib/article/25/4/bbae250/7682297
- https://link.springer.com/article/10.1186/s12859-025-06044-0
- https://www.nature.com/articles/s41592-025-02814-z
- https://www.biorxiv.org/content/10.1101/2025.02.05.636714v1
- https://link.springer.com/article/10.1186/s13059-025-03773-6
- https://www.nature.com/articles/s41592-022-01651-8
- https://www.biorxiv.org/content/10.1101/2025.04.23.649965v1
- https://www.biorxiv.org/content/10.64898/2026.03.01.708931v1
- https://pmc.ncbi.nlm.nih.gov/articles/PMC11744978/
- https://www.biorxiv.org/content/10.1101/2025.02.05.636714.full.pdf
- https://www.nature.com/articles/s41587-025-02811-9
