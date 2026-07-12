"""CNV-based malignant calling via infercnvpy - runs in an ISOLATED env (cnv_env).

infercnvpy 0.6.0 (and its numba/llvmlite stack) lives only in the ``cnv_env`` conda
env and its shared objects will not load in the main SpatialScribe interpreter, so this runs as
a SEPARATE process and writes a per-cell ``{cell_id, cnv_score, is_malignant}`` parquet that the
main app joins back (see ``cnv.call_malignant_cnv``). Never import infercnvpy in the main env.

Env (required for the import to resolve):
    PYTHONNOUSERSITE=1
    LD_LIBRARY_PATH=/data
    (and avoid SLURM node a compute node - CXXABI mismatch)
Run:
    .../cnv_env/bin/python subprocesses/cnv/run_infercnv.py \
        --h5ad section.h5ad --out cnv.parquet --gtf genes.gtf \
        --cell-type-key cell_type --reference "Stromal/CAF,Endothelial,T cell,Myeloid,B/Plasma"

Method (the cluster insitucnv-analysis skill, Jensen et al. InSituCNV):
  1. Genomic coordinates from a GTF (manual parse - infercnvpy's genomic_position_from_gtf
     crashes on a Categorical chromosome column).
  2. **Neighbor-based smoothing (k=100)** - essential on sparse spatial data; raw per-cell counts
     are too sparse for direct CNV inference.
  3. ``cnv.tl.infercnv`` with a non-tumor reference (window 100 / step 10).
  4. Per-cell CNV burden = RMS of the X_cnv row (the published-figure convention), NOT the
     per-cluster ``cnv.tl.cnv_score``. Malignant = burden above the reference's 95th percentile.

CNV quality scales with panel size: >=5K Prime is strong (AUROC 0.96-0.98); a 480-gene panel is
unreliable. This is intended for Prime-5K / WTA sections.
"""
from __future__ import annotations

import argparse
import sys


def gene_positions_from_gtf(gtf: str, symbols) -> dict:
    """Map gene symbol -> (``chrN``, start, end) by parsing GTF gene rows directly.

    Avoids ``infercnvpy.io.genomic_position_from_gtf`` (it does ``"chr" + var["chromosome"]``,
    which crashes when the chromosome column parses as a pandas Categorical). ~2s for a WTA panel.
    """
    want = set(symbols)
    keep = {str(c) for c in list(range(1, 23)) + ["X"]}
    pos: dict[str, tuple] = {}
    with open(gtf) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            p = line.split("\t")
            if len(p) < 9 or p[2] != "gene":
                continue
            i = p[8].find('gene_name "')
            if i < 0:
                continue
            nm = p[8][i + 11:p[8].find('"', i + 11)]
            if nm in want and nm not in pos and p[0] in keep:
                pos[nm] = ("chr" + p[0], int(p[3]), int(p[4]))
    return pos


def _smooth_for_cnv(adata, n_neighbors: int) -> None:
    """Jensen InSituCNV neighbor-based smoothing: weighted mean of normalized counts over k-NN."""
    import numpy as np
    import scanpy as sc
    from scipy.sparse import csr_matrix, issparse

    n_comps = min(50, adata.shape[1] - 1, adata.shape[0] - 1)
    sc.pp.pca(adata, n_comps=n_comps)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_comps)
    conn = adata.obsp["connectivities"]
    raw = adata.layers["raw_norm"]
    if not issparse(raw):
        raw = csr_matrix(raw)
    if not issparse(conn):
        conn = csr_matrix(conn)
    adata.layers["M"] = csr_matrix.dot(conn, raw).astype(np.float32).toarray()


def main() -> int:
    ap = argparse.ArgumentParser(description="CNV-based malignant calling via infercnvpy (isolated env).")
    ap.add_argument("--h5ad", required=True, help="section h5ad; X (or layer 'counts') = RAW counts")
    ap.add_argument("--out", required=True, help="output {cell_id, cnv_score, is_malignant} parquet")
    ap.add_argument("--gtf", required=True, help="GTF for gene genomic positions")
    ap.add_argument("--cell-type-key", default="cell_type")
    ap.add_argument("--reference", required=True, help="comma-separated non-tumor reference cell types")
    ap.add_argument("--n-neighbors", type=int, default=100)
    ap.add_argument("--max-cells", type=int, default=0, help="uniform subsample cap (0 = all cells)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    import anndata as ad
    import infercnvpy as cnv
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy.sparse import issparse

    adata = ad.read_h5ad(args.h5ad)
    if "counts" in adata.layers:                       # ensure X is raw counts
        adata.X = adata.layers["counts"].copy()
    if args.max_cells and adata.n_obs > args.max_cells:
        rng = np.random.default_rng(args.seed)
        keep = np.sort(rng.choice(adata.n_obs, args.max_cells, replace=False))
        adata = adata[keep].copy()

    # 1. genomic coordinates
    pos = gene_positions_from_gtf(args.gtf, list(adata.var_names))
    have = [g for g in adata.var_names if g in pos]
    if len(have) < 100:
        print(f"ERROR: only {len(have)} genes got genomic positions - too few for CNV", file=sys.stderr)
        return 2
    adata = adata[:, have].copy()
    adata.var["chromosome"] = pd.Categorical([pos[g][0] for g in adata.var_names])
    adata.var["start"] = [pos[g][1] for g in adata.var_names]
    adata.var["end"] = [pos[g][2] for g in adata.var_names]

    # 2. Jensen neighbor-smoothing (essential on sparse spatial data)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    adata.layers["raw_norm"] = adata.X.copy()
    sc.pp.log1p(adata)                                  # log space for PCA/neighbors only
    _smooth_for_cnv(adata, n_neighbors=args.n_neighbors)
    adata.X = adata.layers["M"].copy()                  # smoothed normalized counts
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    # 3. infercnv against the non-tumor reference (only categories actually present)
    ct = adata.obs[args.cell_type_key].astype(str)
    ref_cats = [c for c in args.reference.split(",") if c and (ct == c).sum() >= 20]
    if not ref_cats:
        print(f"ERROR: no reference category from {args.reference!r} has >=20 cells", file=sys.stderr)
        return 3
    adata.obs[args.cell_type_key] = ct.astype("category")
    cnv.tl.infercnv(adata, reference_key=args.cell_type_key, reference_cat=ref_cats,
                    window_size=100, step=10, dynamic_threshold=1.5)

    # 4. per-cell CNV burden = RMS of the X_cnv row (published-figure convention)
    X_cnv = adata.obsm["X_cnv"]
    if issparse(X_cnv):
        X_cnv = X_cnv.toarray()
    burden = np.sqrt(np.mean(np.square(X_cnv), axis=1)).ravel()

    ref_mask = ct.isin(ref_cats).to_numpy()
    thr = float(np.nanpercentile(burden[ref_mask], 95))   # malignant = above the normal 95th pct
    is_mal = burden > thr

    out = pd.DataFrame({"cell_id": adata.obs_names.astype(str),
                        "cnv_score": burden.astype(float),
                        "is_malignant": is_mal.astype(bool)})
    import pyarrow as pa
    import pyarrow.parquet as pq
    pq.write_table(pa.Table.from_pandas(out, preserve_index=False), args.out)

    print(f"cnv: {adata.n_obs} cells, {len(have)} positioned genes, ref {ref_cats} "
          f"(thr={thr:.4f}), malignant {is_mal.mean():.1%}", file=sys.stderr)
    print(args.out)   # convention: parquet path is the last stdout line
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
