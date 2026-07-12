"""TACCO annotation -> the SAME plain SPLIT inputs that rctd_full.py writes, so the exact same
`run_split.R` (annotation-agnostic `rctd_free_purify` / SPLIT >=0.2.0 `purify`) can purify from
TACCO deconvolution weights instead of RCTD's. This is the TACCO analogue of SPLIT's own
`convert_rctd_result_to_purify_input`: it fills the three annotation-agnostic slots
(deconvolution_weights, reference, primary_cell_type) from `tc.tl.annotate`'s compositional output.

TACCO (`tc.tl.annotate`) returns a per-cell composition over the reference cell types
(rows sum to 1) - that IS the deconvolution-weights matrix SPLIT wants. The reference profile
(types x genes) is computed IDENTICALLY to rctd_full.py so that, cell-for-cell, the ONLY thing
that differs between the two SPLIT runs is the weight source (RCTD vs TACCO).

Usage:
  python tacco_weights.py <section.h5ad> <reference.h5ad> <ref_label_key> <outdir> [--match-cells cells.txt]

`--match-cells` restricts the output to a given cell-id list (e.g. the RCTD-surviving
split_track/inputs/cells.txt) so the head-to-head SPLIT comparison runs on an identical cell
population. TACCO still annotates the FULL section first (its natural global-OT mode); only the
written outputs are subset.

Writes into <outdir> (same contract as rctd_full.py): counts.mtx (genes x cells), weights.csv
(cells x types, rows sum to 1), primary.csv (cell_id, primary), reference.csv (types x genes,
each row sums to 1), genes.txt, cells.txt.
"""
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp


def _reference_profile(ref, ref_key, type_order):
    """types x genes gene-frequency profile: each ref cell library-normalised, mean per type,
    each type row renormalised to sum 1. IDENTICAL formula to rctd_full.py (so the SPLIT
    reference is the same object regardless of which annotator produced the weights)."""
    refX = sp.csr_matrix(ref.X).astype(float)
    lib = np.asarray(refX.sum(1)).ravel(); lib[lib == 0] = 1.0
    refn = refX.multiply((1.0 / lib)[:, None]).tocsr()
    labels = ref.obs[ref_key].astype(str).to_numpy()
    P = np.vstack([np.asarray(refn[labels == t].mean(0)).ravel() for t in type_order])
    return P / np.clip(P.sum(1, keepdims=True), 1e-12, None)


def main(sec_h5, ref_h5, ref_key, outdir, match_cells=None):
    os.makedirs(outdir, exist_ok=True)
    sec = ad.read_h5ad(sec_h5)
    ref = ad.read_h5ad(ref_h5)

    import tacco as tc  # heavy import kept after IO

    # Compositional annotation over the reference cell types -> cells x types (rows sum to 1).
    tc.tl.annotate(sec, ref, annotation_key=ref_key, result_key="tacco")
    comp = sec.obsm["tacco"]
    names = list(comp.columns) if hasattr(comp, "columns") else \
        sorted(ref.obs[ref_key].astype(str).unique())
    W_full = np.asarray(comp, dtype=float)
    W_full = W_full / np.clip(W_full.sum(1, keepdims=True), 1e-12, None)
    primary_full = [names[i] for i in W_full.argmax(1)]

    all_ids = np.asarray(sec.obs_names).astype(str)
    if match_cells:
        keep = pd.read_csv(match_cells, header=None)[0].astype(str).tolist()
        pos = pd.Index(all_ids).get_indexer(keep)
        if (pos < 0).any():
            missing = int((pos < 0).sum())
            raise SystemExit(f"{missing} match-cells ids not found in the section")
        sel = pos
        ids = np.asarray(keep)
    else:
        sel = np.arange(sec.n_obs)
        ids = all_ids

    W = W_full[sel]
    primary = [primary_full[i] for i in sel]
    counts = sp.csr_matrix(sec[sel].X).T.tocsr()          # genes x cells
    P = _reference_profile(ref, ref_key, names)           # types x genes (same formula as RCTD)

    genes = list(map(str, sec.var_names))
    sio.mmwrite(os.path.join(outdir, "counts.mtx"), counts)
    pd.DataFrame(W, index=ids, columns=names).to_csv(os.path.join(outdir, "weights.csv"))
    pd.DataFrame({"cell_id": ids, "primary": primary}).to_csv(os.path.join(outdir, "primary.csv"), index=False)
    pd.DataFrame(P, index=names, columns=genes).to_csv(os.path.join(outdir, "reference.csv"))
    pd.Series(genes).to_csv(os.path.join(outdir, "genes.txt"), index=False, header=False)
    pd.Series(ids).to_csv(os.path.join(outdir, "cells.txt"), index=False, header=False)
    print(f"saved TACCO SPLIT inputs to {outdir} | cells: {len(ids)} (of {sec.n_obs}) "
          f"genes: {len(genes)} types: {len(names)}")
    print("primary dist:", pd.Series(primary).value_counts().to_dict())


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    mc = None
    if "--match-cells" in sys.argv:
        mc = sys.argv[sys.argv.index("--match-cells") + 1]
    main(args[0], args[1], args[2], args[3], match_cells=mc)
