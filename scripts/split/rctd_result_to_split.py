"""Extract SPLIT inputs from an ALREADY-COMPUTED rctd-py result h5ad (no RCTD re-run).

Same output contract as rctd_full.py, but reads a stored rctd-py result
(``obsm['rctd_weights']`` full-mode cells x types, ``uns['rctd_cell_type_names']``,
``obs['rctd_first_type']``/``rctd_spot_class``) instead of running RCTD. Use it to drive SPLIT
from a big section whose RCTD was run once on GPU. Keeps only the RCTD-assigned cells (valid,
non-NaN weights + a real primary type; drops 'filtered'/'reject').

Usage: python rctd_result_to_split.py <rctd_result.h5ad> <reference.h5ad> <ref_key> <section.h5ad> <outdir>
Writes: counts.mtx, weights.csv, primary.csv, reference.csv, genes.txt, cells.txt (as rctd_full.py).
"""
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp


def _reference_profile(ref, ref_key, type_order):
    """types x genes gene-frequency profile (IDENTICAL formula to rctd_full.py / tacco_weights.py)."""
    refX = sp.csr_matrix(ref.X).astype(float)
    lib = np.asarray(refX.sum(1)).ravel(); lib[lib == 0] = 1.0
    refn = refX.multiply((1.0 / lib)[:, None]).tocsr()
    labels = ref.obs[ref_key].astype(str).to_numpy()
    P = np.vstack([np.asarray(refn[labels == t].mean(0)).ravel() for t in type_order])
    return P / np.clip(P.sum(1, keepdims=True), 1e-12, None)


def main(rctd_h5, ref_h5, ref_key, sec_h5, outdir):
    os.makedirs(outdir, exist_ok=True)
    r = ad.read_h5ad(rctd_h5)
    ref = ad.read_h5ad(ref_h5)
    sec = ad.read_h5ad(sec_h5)
    assert (r.obs_names == sec.obs_names).all(), "rctd result and section obs_names differ"

    names = list(map(str, r.uns["rctd_cell_type_names"]))
    W_all = np.asarray(r.obsm["rctd_weights"], dtype=float)          # cells x types (full mode)
    first = r.obs["rctd_first_type"].astype(str).to_numpy()
    valid = (~np.isnan(W_all).any(1)) & np.isin(first, names)         # drop filtered/reject
    ids = np.asarray(sec.obs_names)[valid].astype(str)

    W = W_all[valid]
    W = W / np.clip(W.sum(1, keepdims=True), 1e-12, None)
    primary = first[valid]
    counts = sp.csr_matrix(sec[valid].X).T.tocsr()                   # genes x cells
    P = _reference_profile(ref, ref_key, names)

    genes = list(map(str, sec.var_names))
    sio.mmwrite(os.path.join(outdir, "counts.mtx"), counts)
    pd.DataFrame(W, index=ids, columns=names).to_csv(os.path.join(outdir, "weights.csv"))
    pd.DataFrame({"cell_id": ids, "primary": primary}).to_csv(os.path.join(outdir, "primary.csv"), index=False)
    pd.DataFrame(P, index=names, columns=genes).to_csv(os.path.join(outdir, "reference.csv"))
    pd.Series(genes).to_csv(os.path.join(outdir, "genes.txt"), index=False, header=False)
    pd.Series(ids).to_csv(os.path.join(outdir, "cells.txt"), index=False, header=False)
    print(f"saved RCTD SPLIT inputs to {outdir} | assigned cells: {len(ids)} (of {sec.n_obs}) "
          f"genes: {len(genes)} types: {len(names)}")
    print("primary dist:", pd.Series(primary).value_counts().to_dict())


if __name__ == "__main__":
    main(*sys.argv[1:6])
