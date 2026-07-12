# Extensions & roadmap

Where SpatialScribe stands against the "genuinely missing rungs" for a melanoma/Xenium lab,
and the tools that fill them.

## The six gaps

| # | Gap | Status |
|---|-----|--------|
| 1 | "Which cells are the tumor?" (CNV-based malignant calling) | ✅ **built** - `cnv.py`: dependency-free `malignant_score` (marker proxy, always on) + `call_malignant_cnv` (optional `infercnvpy`, chr-arm loss e.g. chr3 in UM, normal-reference-aware). Wired into the copilot (`malignant_score` / `malignant_concordance`). |
| 2 | "What are the TME niches?" (domain caller) | ✅ **built** - `niches.py`: neighbor-composition niches (no deps) + optional `novae` zero-shot backend. |
| 3 | Reference-anchored annotation (own the melanoma-lab atlas) | ✅ **built** - `reference_transfer.py` (CellTypist in-env primary; TACCO OT optional), CPU; runs when a reference `.h5ad` is supplied. |
| 4 | Spillover correction (not just flagging) | ✅ **built** - ovrlpy flags (H5) + `split_purify` (SPLIT spillover purification, reference-free fallback), copilot-exposed; the marker dot-plot shows raw vs SPLIT-purified. |
| 5 | De-novo spatial programs (beyond fixed marker/state lists) | ✅ **built** - `programs.py`: NMF gene programs (cNMF-lite, no deps), optional spatial smoothing. Wired into the copilot (`discover_programs` / `name_programs`). |
| 6 | Calibrated annotation uncertainty ("how much should I trust this label?") | ⬜ roadmap - deepen the (built) Layer-5 abstention from heuristic thresholds to a per-cell **composition posterior + Fisher/Cramer-Rao separability bound** off the differentiable `rctd-py`, plus a **reference-fit** ("wrong atlas for this tissue") signal. See Tier 2 below + [spatial-as-code-dialogue.md](spatial-as-code-dialogue.md). |

## Tier 1 - added now (CPU-friendly, one-click, high value)

- **squidpy neighborhood enrichment + co-occurrence** (`spatial.py`) - the literal answer to
  "are T cells excluded from the tumor?"; wired into the copilot (`immune_exclusion`,
  `neighborhood_enrichment`). Zero install. **Validated on real data: z = -140.8.**
- **novae zero-shot niches** (`niches.py`, `method="novae"`) - scverse graph foundation model
  pretrained on ~30M iST cells; optional (heavy torch dep). The default composition caller
  needs no install and ships in the demo. Cache the HF weights for the reproducibility gate.
- **TACCO reference transfer** (`reference_transfer.py`) - CPU-native OT annotation transfer
  from a reference atlas; the compositional output doubles as an ambiguity/confidence flag
  that feeds the same abstention machinery as the marker path.

## Tier 2 - roadmap: calibrated abstention from differentiable RCTD

The single most transferable idea from [spatial-as-code-dialogue.md](spatial-as-code-dialogue.md),
and the one that lands squarely on SpatialScribe's differentiator (the Layer-5 confidence/abstention
layer). Because `rctd-py` ports RCTD to PyTorch, the deconvolution is now **differentiable and
GPU-native** - which puts principled per-cell uncertainty within reach:

- **Per-cell composition posterior (backs `annotate.py`).** The per-cell solver is a constrained
  IRWLS (a weighted MLE). The Hessian of the negative log-likelihood at the MLE is the observed
  **Fisher information**; its inverse is the asymptotic covariance of the composition estimate. That
  turns "type A 0.6 / type B 0.4" into a calibrated posterior plus a **Cramer-Rao bound** on how well
  *any* method could separate those types given this cell's UMI depth and the reference's
  separability. This is the principled backbone for the abstention classes that today run on
  thresholds ([annotation_qc_thresholds.yaml](research/annotation_qc_thresholds.yaml)): "these two
  types are information-non-separable here" maps almost 1:1 onto `Ambiguous` / `Unresolvable` /
  `Uncertain`.
- **Reference-fit / misspecification signal (the high-value one).** A restatement of SpatialScribe's
  own "reference choice > method choice" thesis and the 5K-panel melanocyte-gap gotcha. A
  misspecified reference (e.g. a tumor slide deconvolved against a normal-only atlas) leaves a
  fingerprint in the deconvolution residuals; a differentiable, probabilistic RCTD could surface a
  calibrated "this atlas does not fit this tissue" flag instead of returning a confident wrong label.
  Exactly the abstention a wet-lab user needs when they supply a mismatched reference.
- **Panel-adequacy as design, not just diagnosis (extends H3).** With the same per-cell Fisher
  information you can score a candidate gene panel by expected reduction in composition uncertainty
  for the cell types of interest - flipping H3 from *diagnosing* an inadequate panel (the melanocyte
  gap) to *recommending* an information-optimal one.

**Two load-bearing caveats (do not ship without them):**
- **Simplex-boundary fragility.** Weights live on a simplex; the Gaussian/Laplace approximation
  breaks as a weight -> 0, which is *exactly* the singlet/doublet frontier. Use a constraint-aware
  posterior (logistic-normal or Dirichlet on the simplex, or a projected covariance), not a naive
  Laplace.
- **Calibration != correction.** A well-calibrated posterior under a *misspecified* reference is
  worse than an honest point estimate, because it looks trustworthy. The reference-fit signal above
  is the guard: never surface confidence without also checking that the reference fits.

Speculative - a hypothesis from the dialogue, not a benchmarked result. Roadmap, not built.

## Install notes

- `novae` and `tacco` are **optional** extras (not in the core env) - they pull heavier deps
  (novae needs PyTorch). Add them in a dedicated env when needed; the modules import-guard and
  raise a clear message if absent.
- CNV (`infercnvpy`) and SPLIT live behind existing the cluster skills (`insitucnv-analysis`,
  `split-purification`) and would be wired as optional steps.

## Annotation-QC layers (see `research/cell-annotation-qc.md`)

Layer 0-2 (section + per-cell QC) ✅, Layer 3 purity (CRISP/MECR) ✅ + ovrlpy VSI ◻ (subprocess),
Layer 4 panel adequacy ✅, Layer 5 confidence/abstention ✅, Layer 6 spatial coherence ◻ (nhood
present; per-cell coherence/PAS is a small extension).
