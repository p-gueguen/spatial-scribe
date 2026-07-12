#!/usr/bin/env Rscript
# SPLIT residual-contamination purification driven from rctd-py slot exports (rctd_free_purify,
# the annotation-agnostic path). Reads the plain inputs written by rctd_full.py and writes the
# purified counts + cell_meta. SPLIT 0.3.0, R 4.5.
suppressMessages({library(Matrix); library(SPLIT)})
args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]; outdir <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

counts <- as(Matrix::readMM(file.path(indir, "counts.mtx")), "CsparseMatrix")  # genes x cells
genes  <- readLines(file.path(indir, "genes.txt"))
cells  <- readLines(file.path(indir, "cells.txt"))
rownames(counts) <- genes; colnames(counts) <- cells

W <- as.matrix(read.csv(file.path(indir, "weights.csv"), row.names = 1, check.names = FALSE))  # cells x types
prim <- read.csv(file.path(indir, "primary.csv"))
# MUST keep cell-id names: rctd_free_purify does primary_cell_type[shared_cells] (indexes BY NAME);
# an unnamed vector returns NA for every cell -> reference[NA,] subscript-out-of-bounds.
primary <- setNames(as.character(prim$primary[match(cells, prim$cell_id)]), cells)
ref <- as.matrix(read.csv(file.path(indir, "reference.csv"), row.names = 1, check.names = FALSE))  # types x genes

cat("counts", nrow(counts), "x", ncol(counts), "| weights", nrow(W), "x", ncol(W),
    "| reference", nrow(ref), "x", ncol(ref), "\n")

res <- SPLIT::rctd_free_purify(
  counts = counts,
  deconvolution_weights = W,
  reference = ref,
  primary_cell_type = primary,
  DO_run_in_chunks = TRUE, chunk_size = 5000,
  DO_remove_residual_contamination = TRUE,   # the residual-contamination step (Notion request)
  belonging_threshold = 0.5
)

pc <- as(res$purified_counts, "CsparseMatrix")
Matrix::writeMM(pc, file.path(outdir, "purified_counts.mtx"))
write.csv(res$cell_meta, file.path(outdir, "cell_meta.csv"), row.names = FALSE)
writeLines(rownames(pc), file.path(outdir, "purified_genes.txt"))
writeLines(colnames(pc), file.path(outdir, "purified_cells.txt"))

libs_raw <- Matrix::colSums(counts)
libs_pur <- Matrix::colSums(pc)
cat("SPLIT done:", nrow(pc), "genes x", ncol(pc), "cells\n")
cat("median library size raw ->", median(libs_raw), " purified ->", median(libs_pur),
    sprintf(" (%.1f%% removed)\n", 100 * (1 - median(libs_pur) / median(libs_raw))))
if ("purification_status" %in% names(res$cell_meta))
  print(table(res$cell_meta$purification_status))
