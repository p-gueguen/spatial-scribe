# Deep research reports

The evidence base behind SpatialScribe's QC and annotation logic - literature reviews and
adversarially-verified research passes, plus the design docs derived from them. These
ground the code's thresholds and method choices, so nothing in the pipeline is an ad-hoc
default. User-facing guides live one level up in [`../`](../).

## Design references (grounded in the research)

- **[cell-annotation-quality-metrics.md](cell-annotation-quality-metrics.md)** - the
  **authoritative catalog** of metrics for evaluating cell-type annotation quality, across
  scRNA-seq **and** spatial. Organized by what evidence you have: **external / reference-based**
  (accuracy, macro-F1, ARI/NMI, biological fidelity), **internal / reference-free** (silhouette
  variants, neighborhood purity, inter-sample consistency, MECR/CRISP/PMP, confidence margins,
  spatial coherence), and **reference-anchored** (NCP/NMP). Includes a decision guide, default
  metric **batteries** per scenario, a master table, and the **scTypeEval** integration. Start
  here to *choose which metrics to compute*; the doc below is the spatial per-cell funnel.
- **[cell-annotation-qc.md](cell-annotation-qc.md)** - QC metrics for cell-type annotation
  in imaging spatial transcriptomics. The six-layer model (segmentation → counts →
  contamination → panel adequacy → annotation confidence → spatial coherence), the per-cell
  **confidence + abstention** framework, and mitigations for sparsity / spatial doublets /
  diffusion. This is the spec `analysis/annotate.py` and `qc.py` implement.
- **[confidence-calibration.md](confidence-calibration.md)** - what `annotation_confidence` is
  actually worth, measured on three sections with independent expert ground truth. Isotonic
  calibration makes the number honest but mostly just learns the section's base rate (Brier skill
  +0.001 to +0.020); pooled AUC is Simpson-confounded and inverted on a real section; within-lineage
  AUC is 0.48-0.57 everywhere; a calibrator does not transfer between sections. Read this before
  claiming anything about confidence, abstention, or ECE.
- **[annotation-method-selection.md](annotation-method-selection.md)** - which annotator to
  actually run. Reference choice > method choice; SingleR / RCTD-doublet + snRNA reference +
  SPLIT as the Xenium defaults; foundation models do not win zero-shot on native Xenium.
- **[annotation_qc_thresholds.yaml](annotation_qc_thresholds.yaml)** - machine-readable
  pass/warn/fail thresholds (provisional, tunable) consumed by the QC code.

## Verified-findings appendices (sources + claims)

- **[research-evidence-appendix.md](research-evidence-appendix.md)** - verified findings +
  full source list behind `cell-annotation-qc.md`.
- **[annotation-method-research-appendix.md](annotation-method-research-appendix.md)** -
  verified findings + sources behind `annotation-method-selection.md`.

## [sources/](sources/)

Source PDFs cited by the reports. (These are large binaries; consider referencing the
biorxiv/journal URLs instead of committing them if the repo needs to stay lean.)
