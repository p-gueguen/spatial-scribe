"""Re-run rctd-py on a section keeping the FULL weights matrix + per-type reference profile,
and write the plain inputs SPLIT::rctd_free_purify needs (counts, weights, primary, reference).

Usage: python rctd_full.py <section.h5ad> <reference.h5ad> <ref_label_key> <outdir>
Writes into <outdir>: counts.mtx (genes x cells, RCTD-surviving), weights.csv (cells x types,
row-normalized full-mode), primary.csv (cell_id, primary_type), reference.csv (types x genes,
each row sums to 1), genes.txt, cells.txt.
"""
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp


def main(sec_h5, ref_h5, ref_key, outdir):
    os.makedirs(outdir, exist_ok=True)
    sec = ad.read_h5ad(sec_h5)
    ref = ad.read_h5ad(ref_h5)
    from rctd import Reference, RCTDConfig, run_rctd

    R = Reference(ref, cell_type_col=ref_key)
    res = run_rctd(sec, R, mode="doublet", config=RCTDConfig(device="auto", UMI_min=20), batch_size=5000)
    names = list(res.cell_type_names)
    mask = np.asarray(res.pixel_mask).astype(bool)
    ids = np.asarray(sec.obs_names)[mask].astype(str)

    # Full-mode weights (cells x types), row-normalized to sum 1
    W = np.asarray(res.weights, dtype=float)
    W = W / np.clip(W.sum(1, keepdims=True), 1e-12, None)
    primary = [names[int(i)] for i in np.asarray(res.first_type)]

    # counts genes x cells (RCTD-surviving cells only)
    counts = sp.csr_matrix(sec[mask].X).T.tocsr()   # genes x cells

    # reference profile: types x genes, each ref cell normalized to sum 1 then mean per type,
    # then each type row renormalized to sum 1 (matches rctd cell_type_info gene-frequency profile)
    refX = sp.csr_matrix(ref.X).astype(float)
    lib = np.asarray(refX.sum(1)).ravel(); lib[lib == 0] = 1.0
    refn = refX.multiply((1.0 / lib)[:, None]).tocsr()
    labels = ref.obs[ref_key].astype(str).to_numpy()
    P = np.vstack([np.asarray(refn[labels == t].mean(0)).ravel() for t in names])  # types x genes
    P = P / np.clip(P.sum(1, keepdims=True), 1e-12, None)

    genes = list(map(str, sec.var_names))
    sio.mmwrite(os.path.join(outdir, "counts.mtx"), counts)
    pd.DataFrame(W, index=ids, columns=names).to_csv(os.path.join(outdir, "weights.csv"))
    pd.DataFrame({"cell_id": ids, "primary": primary}).to_csv(os.path.join(outdir, "primary.csv"), index=False)
    pd.DataFrame(P, index=names, columns=genes).to_csv(os.path.join(outdir, "reference.csv"))
    pd.Series(genes).to_csv(os.path.join(outdir, "genes.txt"), index=False, header=False)
    pd.Series(ids).to_csv(os.path.join(outdir, "cells.txt"), index=False, header=False)
    print(f"saved SPLIT inputs to {outdir} | surviving cells: {len(ids)} genes: {len(genes)} types: {len(names)}")
    print("spot_class dist:", pd.Series(primary).value_counts().to_dict())


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
