"""SingleR (BiocPy) reference annotation -> per-cell parquet. CPU; runs in the singler env.
Usage: python run_singler.py <spatial.h5ad> <reference.h5ad> <ref_label_key> <out.parquet>
SingleR expects genes x cells (transpose of AnnData) and returns NO cell_id (positional).
"""
import sys

import anndata as ad
import numpy as np
import pandas as pd


def main(spatial_h5, ref_h5, ref_key, out_path):
    adata = ad.read_h5ad(spatial_h5)
    ref = ad.read_h5ad(ref_h5)
    import singler

    result = singler.annotate_single(
        test_data=adata.X.T,                 # genes x cells
        test_features=list(adata.var_names),
        ref_data=ref.X.T,                    # genes x cells
        ref_features=list(ref.var_names),
        ref_labels=ref.obs[ref_key].astype(str).tolist(),
    )
    df = result.to_pandas()
    out = pd.DataFrame({
        "cell_id": np.asarray(adata.obs_names).astype(str),
        "label": np.asarray(df["best"]).astype(str),
        "delta": np.asarray(df["delta"], dtype=float),
    })
    import pyarrow as pa
    import pyarrow.parquet as pq
    pq.write_table(pa.Table.from_pandas(out, preserve_index=False), out_path)  # avoid pandas.to_parquet pyarrow-24 patch bug
    print(out_path)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
