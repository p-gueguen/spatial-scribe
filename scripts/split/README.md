# SPLIT reference-path purification (residual contamination)

The RCTD -> SPLIT recipe from the annotation-method-selection doc: deconvolve with rctd-py, then
purify transcript spillover with `SPLIT::rctd_free_purify` (annotation-agnostic, no R RCTD object
needed). Demonstrated on the 10k demo section: **median library size 76 -> 35 (54.5% spillover
removed)**, T-cell x Epithelial co-expression (MECR) 0.7% -> 0%.

1. `rctd_full.py <section.h5ad> <reference.h5ad> <ref_key> <outdir>` (rctd-py env) - runs RCTD
   doublet mode and writes the plain SPLIT inputs (counts.mtx genes x cells, weights.csv cells x
   types full-mode row-normalized, primary.csv, reference.csv types x genes, genes/cells txt).
2. `Rscript run_split.R <indir> <outdir>` (R 4.5 with SPLIT >=0.3.0) - `rctd_free_purify` with
   `DO_remove_residual_contamination=TRUE`, `belonging_threshold=0.5`. Writes purified_counts.mtx +
   cell_meta.csv. GOTCHA: `primary_cell_type` MUST keep cell-id names (SPLIT indexes it by name;
   an unnamed vector -> NA -> `reference[NA,]` subscript-out-of-bounds).
3. `make_split_dotplots.py <section.h5ad> <indir> <split_outdir> <report.html>` (main env) -
   before/after marker dotplots + library-size + MECR spillover metrics.

SPLIT is R-only (bdsc-tds/SPLIT); the rest is Python. Bridge is plain files (mtx/csv).

## Annotation-agnostic: TACCO weights -> the same SPLIT (v0.2.0+)

SPLIT >=0.2.0 takes deconvolution weights from *any* annotator (its own
`convert_rctd_result_to_purify_input` just fills the `deconvolution_weights` / `reference` /
`primary_cell_type` slots). `tacco_weights.py` is the TACCO analogue of `rctd_full.py`: it runs
`tc.tl.annotate` (compositional OT weights = the deconvolution matrix), computes the reference
profile with the **identical** formula, and writes the same 6 plain files, so the SAME
`run_split.R` purifies unchanged.

- `tacco_weights.py <section.h5ad> <reference.h5ad> <ref_key> <outdir> [--match-cells cells.txt]`
  (tacco-env) - `--match-cells split_track/inputs/cells.txt` restricts output to the
  RCTD-surviving cells so the head-to-head runs on an identical cell population.
- `Rscript run_split.R <inputs_tacco> <split_out_tacco>` - unchanged (annotation-agnostic).
- `compare_rctd_tacco.py <rctd_in> <rctd_out> <tacco_in> <tacco_out> <report.html>` (main env) -
  side-by-side report on the shared cells: KPI tiles, type coverage, library + genes/cell
  reduction, cross-lineage MECR (per-pair + mean), marker specificity (on-target fraction),
  marker purity, on-target retention, per-method **raw-vs-purified marker dotplots**, status.

### Result (10k demo, same 7,231 RCTD-surviving cells, same reference + SPLIT params)

Only the weight source differs (RCTD vs TACCO):

| metric | RCTD -> SPLIT | TACCO -> SPLIT |
|---|---|---|
| reference types kept | **5** (drops T/NK/Mast/Adipocyte/Basal via UMI+gene-list filter) | **10** (all) |
| median library removed | 54.5% | 61.6% |
| MECR T x Stromal (purified) | 1.18% | **0.07%** |
| MECR B/Plasma x Epithelial (purified) | 0.06% | 0.00% |
| on-target signal kept (Epi/Mye/B) | 100/97/97% | 100/99/96% |
| on-target signal kept (T cell) | - (type dropped) | 85% (n=778) |

**Takeaway:** TACCO weights purify at least as well as RCTD on the shared lineages (matching
on-target retention) and strictly better on the sparse immune compartment RCTD discards - RCTD
has no T-cell profile, so T-cell transcripts stay misassigned (T x Stromal 1.18%), whereas TACCO
keeps T cells and separates them cleanly (0.07%) without eroding real signal. The extra ~7%
library removal is genuine cross-lineage spillover reassigned into the rare types RCTD ignored,
not deleted signal. Consistent with the method-selection evidence (RCTD filters low-count immune
cells; retain-rare-types annotators do better there). Report:
`/data/spatial-scribe/split_track/rctd_vs_tacco_split.html`.

### Replication on the real full section (100k Prime 5K breast, 72,185 RCTD-assigned cells)

`rctd_result_to_split.py` extracts SPLIT inputs from an already-computed rctd-py result h5ad
(no GPU re-run); TACCO ran on the same 72,185 cells. The finding replicates at ~10x the cells:

| metric | RCTD -> SPLIT | TACCO -> SPLIT |
|---|---|---|
| reference types kept | 5 | 10 |
| median library removed | 54.3% | 61.8% |
| MECR T x Stromal (purified) | 1.14% | **0.05%** |
| on-target signal kept (Epi/Mye/B) | 100/98/94% | 100/99/96% |
| on-target signal kept (T cell) | - (type dropped) | 86% (n=8,137) |

Same conclusion, now on a robust 8k-T-cell population. Reports on object storage:
`https://your-server.example.com/projects/an internal benchmark/Analyses/spatial-scribe/rctd_vs_tacco_split.html`
(10k demo) and `.../rctd_vs_tacco_split_100k.html` (real 100k). Extra recipe step:

- `rctd_result_to_split.py <rctd_result.h5ad> <reference.h5ad> <ref_key> <section.h5ad> <outdir>`
  (main env) - reuse a precomputed rctd-py result instead of re-running RCTD, then `run_split.R`.
