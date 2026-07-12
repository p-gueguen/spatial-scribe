#!/usr/bin/env python3
"""Build the bundled public demo (data/demo_public.h5ad) from a larger processed section.

Stratified-subsamples cells (all cell types kept, >=40 per type) and preserves the annotation,
spatial coordinates, neighbour graph, and per-cluster markers, so the app opens fully-analyzed.

Usage:  python scripts/make_public_demo.py <source_processed.h5ad> <out.h5ad> [target_cells]
"""
import sys
import numpy as np
import anndata as ad


def main() -> int:
    src, dst = sys.argv[1], sys.argv[2]
    target = int(sys.argv[3]) if len(sys.argv) > 3 else 20000
    a = ad.read_h5ad(src)
    rng = np.random.default_rng(0)  # fixed seed -> reproducible
    ct = a.obs["cell_type"].astype(str).values
    keep = []
    for t in np.unique(ct):
        ids = np.where(ct == t)[0]
        n = min(len(ids), max(40, round(target * len(ids) / len(ct))))
        keep.append(rng.choice(ids, size=n, replace=False))
    idx = np.sort(np.concatenate(keep))
    sub = a[idx].copy()  # anndata slices X, layers, obsm, obsp consistently; uns carried as-is
    sub.write_h5ad(dst, compression="gzip")
    print(f"{a.shape} -> {sub.shape} written to {dst}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
