# Research evidence appendix — SRT cell-annotation QC

> Auto-generated from the deep-research adversarial-verification pass that grounds [`cell-annotation-qc.md`](cell-annotation-qc.md). Do not hand-edit; re-run the research to refresh.

**Pass stats:** 6 search angles · 24 sources fetched · 113 claims extracted → **25 verified (23 confirmed, 2 refuted)** → 12 synthesized findings.

## Summary

Across peer-reviewed 2025-2026 primary literature and vendor documentation, a defensible evidence baseline exists for automated cell-type-annotation QC in imaging-based spatial transcriptomics, with 10x Xenium the best-characterized platform. The QC chain is now well supported at every link: (1) segmentation quality - the default rigid ~15 um nuclear-mask expansion is suboptimal and transcripts beyond ~10.71 um from the centroid correlate more with background than with the cell's own type, so cell-type-aware or restrained expansion is a segmentation-QC criterion, and 10x itself has since dropped the default expansion to 5 um; (2) transcript/count QC - a low-count fail floor of <10 transcripts/cell excludes only ~0.21% of cells in targeted Xenium panels but ~60% of cells in the large Xenium Prime 5K panel, so the threshold must be calibrated per panel; (3) negative-control/FDR - vendor-defined per-cell metrics (Q20 decoding cutoff, negative-control-probe counts per control per cell, estimated false positives) plus the reference-based Negative Co-expression Purity (NCP >0.8 = high specificity) give interpretable contamination/false-discovery indicators, with cross-platform specificity ranking Xenium > MERSCOPE > CosMx; (4) spillover/contamination - segmentation-driven transcript misassignment is a first-order effect that can dominate downstream results, is quantifiable against matched snRNA-seq, and is severe enough that RCTD can swap primary and contaminating cell types, making RCTD secondary weights and singlet/doublet calls annotation-confidence QC signals; and (5) correction tools exist (matrix factorization / cellAdmix analogous to doublet filtering, SPLIT post-hoc purification, and the standardized SpatialQM framework with a 7-axis QC taxonomy). The strongest gap is that no single source provides a validated pass/warn/fail lookup table for annotation confidence; concrete thresholds must be assembled across studies, and several numbers are panel-, tissue-, or chemistry-version-specific rather than universal.

## Verified findings (ranked)

### 1. [HIGH · vote 3-0 (merged from two unanimous claims)]

Transcript/count QC: exclude Xenium cells with <10 assigned high-quality transcripts as a draft FAIL floor. This permissive cutoff removes only ~0.21% of cells across 25 Xenium datasets (mean ~186.6 reads/cell), and is the default used in the squidpy Xenium tutorial (min_counts=10).

Sources: [1](https://www.nature.com/articles/s41592-025-02617-2) · [2](https://pmc.ncbi.nlm.nih.gov/articles/PMC11978515/)

### 2. [MEDIUM · vote 2-1 (split)]

Count-based QC thresholds must be calibrated per panel: the large Xenium Prime 5K panel fails ~60% of cells on low transcript counts, versus ~0.21% for a targeted panel. Large/whole-transcriptome-scale panels have markedly lower per-cell sensitivity, so a single fixed count floor is not transferable across panel sizes.

Sources: [1](https://www.nature.com/articles/s41592-026-03089-8)

### 3. [HIGH · vote 3-0 (merged from two unanimous claims)]

Segmentation QC: the default rigid ~15 um nuclear-mask radius expansion is suboptimal because different cell types have different optimal expansion distances, and transcripts located on average >~10.71 um from the cell centroid correlate more with domain-specific background (contamination) than with nuclear cell-type-specific signatures. Draft criterion: prefer cell-type-aware or restrained expansion (or probabilistic assignment); flag transcripts beyond ~10.71 um as likely spillover. 10x itself has since reduced the XOA default expansion from 15 um (v1.x) to 5 um (v2.0+).

Sources: [1](https://www.nature.com/articles/s41592-025-02617-2) · [2](https://pmc.ncbi.nlm.nih.gov/articles/PMC11978515/)

### 4. [HIGH · vote 3-0 (merged from three unanimous claims)]

Segmentation errors are a first-order (not negligible) QC axis: transcript misassignment can be striking enough to dominate downstream results and confounds most single-cell spatial analyses of cellular state (differential expression, neighbor-influence, ligand-receptor). Cell-type classification itself is relatively robust, but context-dependent analyses are highly sensitive, so accurate segmentation is an upstream prerequisite for trustworthy annotation QC. Matrix factorization of local molecular neighborhoods (e.g. cellAdmix) can identify and isolate molecular admixtures, reducing their impact analogously to doublet filtering in scRNA-seq.

Sources: [1](https://www.nature.com/articles/s41588-025-02497-4)

### 5. [HIGH · vote 3-0 (merged from four unanimous claims)]

Vendor-defined Xenium negative-control / false-discovery QC metrics (interpretable per-cell contamination indicators): (a) high-quality decoded transcripts use a fixed on-instrument decoding cutoff of Q-Score >=20 (Q20) for both genes and negative controls, and the cell-feature matrix is filtered to Q20; (b) 'Negative control probe counts per control per cell' = (Q20+ NCP counts) / (# negative control probes) / (# cells), a standardized rate where a HIGH value flags nonspecific/off-target probe binding, poor assay conditions, or low transcript content; (c) per-cell 'estimated false positives' = (NCP counts per control per cell) x (# target genes).

Sources: [1](https://kb.10xgenomics.com/hc/en-us/articles/18385764969613-What-is-the-Xenium-Negative-control-probe-counts-per-control-per-cell-metric) · [2](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/advanced/metric-calculations)

### 6. [HIGH · vote 3-0 (merged from three unanimous claims)]

Reference-based contamination/annotation-confidence QC: Xenium contains substantial transcript spillover between neighboring cells, matched single-nucleus RNA-seq enables precise quantification of that contamination, and contamination can be strong enough that RCTD deconvolution swaps the primary and secondary (contaminating) cell types. In imaging ST (no droplet co-encapsulation) RCTD secondary weights and singlet/doublet calls are contamination signals rather than biological doublets, so they should be treated as annotation-confidence QC flags. SPLIT (Spatial Purification of Layered Intracellular Transcripts) is a segmentation-agnostic post-hoc purification method that uses these deconvolution weights to resolve mixed signals and improve per-cell purity and cell-type resolution.

Sources: [1](https://www.nature.com/articles/s41592-026-03089-8)

### 7. [HIGH · vote 3-0]

Negative Co-expression Purity (NCP) is a reference-based spillover/specificity metric: the fraction of gene pairs that are non-co-expressed in a scRNA-seq reference and stay non-co-expressed in situ. Draft threshold: NCP > 0.8 (close to 1) indicates high specificity. All tested SRT platforms met mean NCP > 0.8, and Xenium was consistently higher than CosMx (though not the top performer overall).

Sources: [1](https://pmc.ncbi.nlm.nih.gov/articles/PMC11978515/)

### 8. [HIGH · vote 3-0 (merged from two unanimous claims)]

Cross-platform specificity/FDR ranking (for platform-aware negative-control thresholds): Xenium > MERSCOPE > CosMx. Xenium consistently showed the lowest false-discovery rate and highest on-target fraction, CosMx the highest FDR and lowest on-target fraction across most cancer/tissue types. In a separate subcellular-platform benchmark, after normalization by total counts Xenium 5K showed a lower proportion of negative-control signals than CosMx 6K (with stronger spatial aggregation of negative-control signal in CosMx). Negative-control rate is thus a discriminating QC metric between platforms.

Sources: [1](https://www.nature.com/articles/s41467-025-64990-y) · [2](https://www.nature.com/articles/s41467-025-64292-3)

### 9. [HIGH · vote 3-0]

Cross-platform segmentation-precision benchmark (concrete segmentation-quality QC anchor): on densely packed cells, CosMx and Xenium achieved segmentation precision 0.90 versus MERSCOPE 0.83 (p < 0.01), evaluated against >31,400 manually annotated cells across dense/sparse/elongated scenarios.

Sources: [1](https://www.nature.com/articles/s41467-025-64990-y)

### 10. [HIGH · vote 3-0]

Per-cell QC pass rates are strongly platform- and sample-quality-dependent: standard filtering retains ~92-96% of Xenium and CosMx cells but only ~28% and ~3% of MERSCOPE cells in two FFPE tissue microarrays (i.e. >70-97% loss). MERSCOPE/MERFISH FFPE data can lose the majority of cells to filtering, so an acceptable-retention QC gate must be platform- and assay-aware.

Sources: [1](https://www.nature.com/articles/s41467-025-64990-y)

### 11. [HIGH · vote 3-0]

Transcript-detection completeness is a QC metric for annotation reliability: in MERFISH/MERSCOPE, losing ~20% of transcripts changes the cluster (cell-type) labels of 10-15% of cells, and dropping 40 transcript species (dropped-image loss) changed labels of ~10% of cells. Detection completeness / dropout therefore directly degrades annotation stability.

Sources: [1](https://elifesciences.org/reviewed-preprints/105149)

### 12. [HIGH · vote 3-0 (merged from two unanimous claims)]

A standardized, operationalizable QC framework exists for annotation-QC baselining: the Spatial Touchstone project defines a 7-axis technical QC taxonomy (reproducibility, sensitivity, dynamic range, signal-to-noise ratio, false-discovery rate, cell-type annotation, and congruence with single-cell profiling) and ships an open-source tool, SpatialQM, that computes these metrics and performs reference-based imputation/transfer of cell annotations, alongside standardized SOPs (STSOPs) and a public 254-profile repository (GEO GSE277080). This positions SpatialQM as a concrete tool-based annotation-QC option alongside squidpy/sopa/RCTD/SPLIT.

Sources: [1](https://www.nature.com/articles/s41587-025-02811-9)

## Caveats (verbatim from the pass)

No single source provides a validated pass/warn/fail lookup table for cell-type-annotation confidence; the concrete thresholds here are assembled across studies and several are panel-, tissue-, or chemistry-version-specific rather than universal. (1) Numeric values are context-bound: the ~10.71 um radial cutoff, the ~60% Xenium 5K QC-fail rate, the <10 transcripts/cell floor, and the 0.21% below-floor figure all come from specific studies/datasets and should be treated as draft anchors, not constants; 10x has already moved the default segmentation expansion from 15 um to 5 um, so 'default expansion' guidance is version-dependent. (2) Panel/platform dependence is strong: the same count floor that removes 0.21% of cells on a targeted Xenium panel removes ~60% on Xenium Prime 5K, and MERSCOPE FFPE can lose 70-97% of cells; thresholds must be calibrated per panel and per assay. (3) Cross-platform rankings (Xenium > MERSCOPE > CosMx for FDR/on-target/specificity) are relative rankings on 2023-era panel/chemistry versions, not numeric cutoffs, and newer chemistries (MERFISH 2.0 in 2026, updated Xenium/CosMx panels) could shift magnitudes. (4) Reference-based metrics have limits: matched snRNA-seq captures nuclear but not cytoplasmic transcripts, making it an imperfect contamination ground truth ('precise' is aspirational). (5) Vote confidence: one finding (panel-size sensitivity / ~60% 5K fail) rested on a 2-1 split; all others were 3-0. (6) Two related claims were refuted in verification and are excluded: an NCP~0.8 warn/fail threshold misattributed to the SPLIT paper (valid only when sourced to the PMC11978515 best-practices paper), and a Xenium per-cell sensitivity-ceiling claim (1-2). (7) Segmentation benchmark precision is one of three metrics on one morphological scenario; retention percentages for CosMx showed some internal inconsistency (~83% aggregate vs ~95-96% per-TMA). (8) Vendor metric definitions (10x) are authoritative for Xenium but do not carry vendor-published numeric pass/warn/fail cutoffs; interpretation ('high NCP = red flag') is qualitative.

## Open questions

- What are validated numeric cutoffs for reference-based deconvolution confidence as a per-cell annotation-QC gate - specifically RCTD singlet-vs-doublet class boundaries, secondary-weight fractions, and spot-class thresholds - that reliably separate true annotation from contamination in Xenium/MERSCOPE/CosMx?
- Are there established quantitative pass/warn/fail thresholds for spatial coherence of labels (e.g. neighborhood label-consistency, Moran's I on annotations) and for marker-gene specificity as a per-cell annotation-confidence metric, or do these remain qualitative in the current literature?
- How do the numeric thresholds and platform rankings transfer to the newest chemistries and panels (Vizgen MERFISH 2.0, updated Xenium Prime 5K, CosMx WTx/6K), given that most benchmarks used 2023-era versions?
- What are the recommended per-panel low-count and negative-control cutoffs across the sensitivity spectrum (small targeted -> Xenium Prime 5K -> whole-transcriptome), i.e. a panel-size-indexed threshold schedule rather than a single fixed floor?

## Refuted & excluded (NOT asserted in the doc)

- Specificity of imaging-based spatial transcriptomics platforms can be quantified with a negative-control-probe-based metric (NCP, range 0-1); across all tested SRT technologies mean specificity was high (NCP > 0.8), with Xenium slightly lower than other commercial platforms but consistently higher than CosMx, giving a candidate warn/fail threshold of NCP ~0.8 for negative-control/specificity QC.
- Across FFPE tissues, 10x Xenium delivers the highest per-cell gene diversity and detection sensitivity among the three imaging platforms, with roughly 2-fold higher counts than CosMx and up to ~20-fold higher median expression than MERSCOPE on shared genes -- setting a sensitivity ceiling that per-cell transcript/gene QC thresholds should be calibrated against per platform.

## All sources fetched

- https://www.nature.com/articles/s41592-025-02617-2
- https://www.nature.com/articles/s41587-025-02811-9
- https://kb.10xgenomics.com/hc/en-us/articles/18385764969613-What-is-the-Xenium-Negative-control-probe-counts-per-control-per-cell-metric
- https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/advanced/metric-calculations
- https://www.nature.com/articles/s41467-025-64990-y
- https://elifesciences.org/reviewed-preprints/105149
- https://www.nature.com/articles/s41588-025-02497-4
- https://www.nature.com/articles/s41592-026-03089-8
- https://pmc.ncbi.nlm.nih.gov/articles/PMC11978515/
- https://www.nature.com/articles/s41467-025-64292-3
- https://www.10xgenomics.com/support/software/xenium-onboard-analysis/1.9/analysis/xoa-output-analysis-summary
- https://www.nature.com/articles/s41587-021-01044-w
- https://www.biorxiv.org/content/10.1101/2025.04.23.649965v1.full
- https://github.com/bdsc-tds/SPLIT
- https://bioconductor.org/packages//release/bioc/vignettes/spacexr/inst/doc/rctd-tutorial.html
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11167053/
- https://www.nature.com/articles/s41592-021-01358-2
- https://scanpy.readthedocs.io/en/1.10.x/tutorials/spatial/basic-analysis.html
- https://pmc.ncbi.nlm.nih.gov/articles/PMC11744978/
- https://academic.oup.com/bioinformatics/article/40/8/btae458/7720780
- https://evo-byte.com/imaging-based-spatial-transcriptomics-preprocessing-and-qc/
- https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/algorithms-overview/xoa-algorithms
- https://elifesciences.org/reviewed-preprints/107070
- https://bioconductor.org/books/release/OSTA/pages/img-quality-control.html
