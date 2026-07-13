"""Cancer-Finder malignant-cell probability - runs in an ISOLATED torch env.

Cancer-Finder (Patchouli-M/SequencingCancerFinder) is a domain-adaptation (VREx) classifier that
labels each cell malignant/normal from its transcriptome. It needs torch + the SequencingCancerFinder
repo importable + a pretrained checkpoint, none of which belong in the main SpatialScribe env, so it
runs as a SEPARATE process and writes a per-cell ``{cell_id, cancerfinder_prob}`` parquet that the
main app joins back (see ``cancerfinder.call_cancerfinder``). Same isolation pattern as ovrlpy / the
annotation methods / the infercnvpy CNV subprocess.

Env: a python with torch + scanpy + anndata + the CF repo on the path (``models/``, ``utils/``).
Run:
    <cf_python> subprocesses/cancerfinder/run_cancerfinder.py --h5ad section.h5ad --out cf.parquet \
        --repo <SequencingCancerFinder> --ckpt <sc_pretrain_article.pkl> [--threshold 0.5]

Note (from the internal Atera benchmark): Cancer-Finder is over-sensitive on single-cell Xenium -
prefer the ranking (probability / AUROC) over a hard 0.5 threshold, or raise the threshold.
"""
from __future__ import annotations

import argparse
import sys


def main() -> int:
    ap = argparse.ArgumentParser(description="Cancer-Finder per-cell malignant probability (isolated env).")
    ap.add_argument("--h5ad", required=True, help="section h5ad; X (or layer 'counts') = counts")
    ap.add_argument("--out", required=True, help="output {cell_id, cancerfinder_prob} parquet")
    ap.add_argument("--repo", required=True, help="SequencingCancerFinder repo (models/ + utils/)")
    ap.add_argument("--ckpt", required=True, help="pretrained checkpoint (e.g. sc_pretrain_article.pkl)")
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--max-cells", type=int, default=0, help="uniform subsample cap (0 = all cells)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    sys.path.insert(0, args.repo)
    import anndata as ad
    import numpy as np
    import pandas as pd
    import torch
    from models import model
    from utils import args_utils, opt_utils

    # The vendored opt_utils.normalize_matrix_counts calls `anndata.AnnData(raw_df, ...)` on a
    # DataFrame. Under pandas copy-on-write that AnnData's X is a READ-ONLY view, and scanpy >= 1.12's
    # normalize_total divides IN PLACE, so the runner dies with "ValueError: output array is read-only"
    # before a single cell is scored. We cannot shim `opt_utils.anndata.AnnData`: the module also does
    # `isinstance(obj, anndata.AnnData)`, which needs the real type. So replace the one function, keeping
    # its maths line for line and changing only that X is a writable float32 copy.
    # Verified: probabilities match the reference run to max|delta| = 1.8e-07 (corr 1.000).
    import scanpy as _sc

    def _normalize_matrix_counts_writable(raw_df, HVG_list, target_sum=10000):
        gene_list = pd.DataFrame(HVG_list).set_index(0)
        raw_df = raw_df.copy()
        raw_df["sum"] = raw_df.sum(axis=1)
        raw_df = raw_df.sort_values(by="sum", ascending=False).drop(columns="sum")
        raw_df = raw_df[~raw_df.index.duplicated()].T
        X = np.array(raw_df.values, dtype=np.float32, order="C")   # the only change: writable
        adata = ad.AnnData(X, raw_df.index.to_frame(), raw_df.columns.to_frame())
        _sc.pp.normalize_total(adata, target_sum=target_sum)
        _sc.pp.log1p(adata)
        out = pd.DataFrame(adata.X, index=adata.obs.index, columns=adata.var.index)
        sel = pd.merge(gene_list, out.T, how="left", left_index=True, right_index=True)
        return sel.fillna(0.0).T[gene_list.index]

    opt_utils.normalize_matrix_counts = _normalize_matrix_counts_writable

    a = ad.read_h5ad(args.h5ad)
    if "counts" in a.layers:
        a.X = a.layers["counts"].copy()
    if args.max_cells and a.n_obs > args.max_cells:
        rng = np.random.default_rng(args.seed)
        a = a[np.sort(rng.choice(a.n_obs, args.max_cells, replace=False))].copy()

    # Cancer-Finder consumes a genes x cells matrix (the transpose of the usual cells x genes).
    adt = a.T
    adt.obs_names = a.var_names
    adt.var_names = a.obs_names

    cf_args = args_utils.create_args_from_infering(matrix=adt, ckp=args.ckpt, threshold=args.threshold)
    algo = model.VREx(cf_args)
    # weights_only=True: the checkpoint is a state dict (+ an HVG list), no arbitrary pickled objects.
    algo.load_state_dict(torch.load(args.ckpt, map_location="cpu", weights_only=True)["model_dict"])
    algo.eval()

    cells: list = []
    probs: list = []
    with torch.no_grad():
        for input_data, input_loader in opt_utils.InferLoaders(cf_args):
            for data in input_loader:
                out = torch.softmax(algo.predict(data[0].float()), axis=1)[:, 1]
            cells += list(input_data.index)
            probs += [p.item() for p in out]

    out_df = pd.DataFrame({"cell_id": [str(c) for c in cells], "cancerfinder_prob": probs})
    import pyarrow as pa
    import pyarrow.parquet as pq
    pq.write_table(pa.Table.from_pandas(out_df, preserve_index=False), args.out)

    p = np.asarray(probs, dtype=float)
    print(f"cancerfinder: {len(cells)} cells, mean prob {p.mean():.3f}, "
          f"frac>{args.threshold} {(p > args.threshold).mean():.3f}", file=sys.stderr)
    print(args.out)   # convention: parquet path is the last stdout line
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
