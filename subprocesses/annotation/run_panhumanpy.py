"""Run panhumanpy (Azimuth Neural Network) on a section's counts; write a per-cell parquet.
Reference-free. Runs in its own env (panhumanpy + TensorFlow); never imported by the main app.
Usage: python run_panhumanpy.py <counts.h5ad> <out.parquet>
"""
import sys

import anndata as ad
import pandas as pd


def main(h5ad_path, out_path):
    adata = ad.read_h5ad(h5ad_path)
    import panhumanpy as ph

    az = ph.AzimuthNN(adata)               # high-level interface: runs the annotation pipeline
    annotated = az.pack_adata()            # returns the AnnData with azimuth_* columns in obs
    meta = annotated.obs
    conf_col = "final_level_confidence" if "final_level_confidence" in meta else "annotation_confidence"
    out = pd.DataFrame({
        "cell_id": adata.obs_names.astype(str),
        "broad": meta["azimuth_broad"].to_numpy(),
        "medium": meta["azimuth_medium"].to_numpy(),
        "fine": meta["azimuth_fine"].to_numpy(),
        "confidence": meta[conf_col].to_numpy(),
    })
    import pyarrow as pa
    import pyarrow.parquet as pq
    pq.write_table(pa.Table.from_pandas(out, preserve_index=False), out_path)  # avoid pandas.to_parquet pyarrow-24 patch bug
    print(out_path)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
