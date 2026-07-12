"""Write the SPLIT::rctd_free_purify inputs from a TACCO optimal-transport transfer, no RCTD.

SPLIT's ``rctd_free_purify`` is annotation-agnostic: it needs a cells x types weight matrix, a
primary type per cell, and a types x genes reference profile - it does NOT care that RCTD produced
them. TACCO's ``tc.tl.annotate`` returns exactly a cells x types compositional matrix, so we can feed
SPLIT from TACCO and drop the second (RCTD) deconvolution entirely. The default annotation engine
already ran TACCO at the Annotate step, so this reuses that work rather than paying for RCTD again.

Usage: python tacco_split_inputs.py <section.h5ad> <reference.h5ad> <ref_label_key> <outdir>
Writes into <outdir> the SAME files rctd_full.py writes (so run_split.R is unchanged):
  counts.mtx (genes x cells), weights.csv (cells x types, rows sum to 1), primary.csv (cell_id,
  primary), reference.csv (types x genes, rows sum to 1), genes.txt, cells.txt.

The section is expected to carry RAW COUNTS in ``layers['counts']`` (SPLIT purifies counts, and TACCO
itself needs integer counts); the loader restores that view. Runs in the main env (tacco is the
declared ``transfer`` extra), so unlike rctd_full.py it needs no separate rctd-py env.
"""
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp


def _counts_view(a):
    """The raw integer-count matrix TACCO and SPLIT both need (never the log1p'd X)."""
    X = a.layers["counts"] if "counts" in a.layers else a.X
    return sp.csr_matrix(X)


def main(sec_h5, ref_h5, ref_key, outdir):
    os.makedirs(outdir, exist_ok=True)
    import tacco as tc

    sec = ad.read_h5ad(sec_h5)
    ref = ad.read_h5ad(ref_h5)

    # TACCO annotate on raw counts (a writable float32 copy - scanpy >= 1.12 divides in place, and a
    # DataFrame-backed / read-only X raises "output array is read-only"; see reference_transfer._as_counts).
    sec_counts = _counts_view(sec)
    a = ad.AnnData(np.asarray(sec_counts.todense(), dtype=np.float32) if sec_counts.shape[0] < 60000
                   else sec_counts.astype(np.float32).copy())
    a.obs_names = sec.obs_names.astype(str)
    a.var_names = sec.var_names.astype(str)
    rc = _counts_view(ref)
    r = ad.AnnData(rc.astype(np.float32).copy())
    r.obs_names = ref.obs_names.astype(str)
    r.var_names = ref.var_names.astype(str)
    r.obs[ref_key] = ref.obs[ref_key].astype(str).to_numpy()

    tc.tl.annotate(a, r, annotation_key=ref_key, result_key="ref_transfer")
    comp = a.obsm["ref_transfer"]
    names = list(comp.columns) if hasattr(comp, "columns") else [str(t) for t in sorted(r.obs[ref_key].unique())]
    W = np.asarray(comp, dtype=float)
    W = W / np.clip(W.sum(1, keepdims=True), 1e-12, None)   # rows sum to 1, same as rctd_full
    ids = a.obs_names.astype(str).to_numpy()
    primary = [names[int(i)] for i in W.argmax(1)]

    # counts genes x cells (all cells - TACCO drops none, unlike RCTD's pixel_mask)
    counts = _counts_view(sec).T.tocsr()

    # reference profile types x genes: per ref cell library-normalized, mean per type, row-renormalized.
    # Built over the REFERENCE's own genes, so its columns are ref.var_names (NOT the section's) - SPLIT
    # intersects reference.csv genes with counts.mtx genes by name, so a mismatched label set is fine but
    # a mismatched LENGTH is a shape error.
    refX = _counts_view(ref).astype(float)
    lib = np.asarray(refX.sum(1)).ravel(); lib[lib == 0] = 1.0
    refn = refX.multiply((1.0 / lib)[:, None]).tocsr()
    labels = ref.obs[ref_key].astype(str).to_numpy()
    P = np.vstack([np.asarray(refn[labels == t].mean(0)).ravel() for t in names])
    P = P / np.clip(P.sum(1, keepdims=True), 1e-12, None)
    ref_genes = list(map(str, ref.var_names))

    genes = list(map(str, sec.var_names))
    sio.mmwrite(os.path.join(outdir, "counts.mtx"), counts)
    pd.DataFrame(W, index=ids, columns=names).to_csv(os.path.join(outdir, "weights.csv"))
    pd.DataFrame({"cell_id": ids, "primary": primary}).to_csv(os.path.join(outdir, "primary.csv"), index=False)
    pd.DataFrame(P, index=names, columns=ref_genes).to_csv(os.path.join(outdir, "reference.csv"))
    pd.Series(genes).to_csv(os.path.join(outdir, "genes.txt"), index=False, header=False)
    pd.Series(ids).to_csv(os.path.join(outdir, "cells.txt"), index=False, header=False)
    print(f"saved TACCO SPLIT inputs to {outdir} | cells: {len(ids)} genes: {len(genes)} types: {len(names)}")
    print("primary dist:", pd.Series(primary).value_counts().to_dict())


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
