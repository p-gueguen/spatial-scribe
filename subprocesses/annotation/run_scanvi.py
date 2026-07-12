"""scANVI (scArches surgery) reference annotation -> per-cell parquet. GPU-preferred; scanvi env.
Usage: python run_scanvi.py <spatial.h5ad> <reference.h5ad> <ref_label_key> <out.parquet>
predict(soft=True) is a DataFrame indexed by obs_names (cell_id ready). label=argmax, conf=max,
entropy=-sum p log p / log K.
"""
import sys

import anndata as ad
import numpy as np
import pandas as pd


def main(spatial_h5, ref_h5, ref_key, out_path):
    query = ad.read_h5ad(spatial_h5)
    ref = ad.read_h5ad(ref_h5)
    import scvi

    if "counts" not in ref.layers:
        ref.layers["counts"] = ref.X.copy()
    scvi.model.SCVI.setup_anndata(ref, layer="counts", labels_key=ref_key)
    scvi_ref = scvi.model.SCVI(ref, n_latent=30)
    scvi_ref.train(max_epochs=200)
    scanvi_ref = scvi.model.SCANVI.from_scvi_model(scvi_ref, labels_key=ref_key, unlabeled_category="Unknown")
    scanvi_ref.train(max_epochs=50)

    query = query[:, [g for g in ref.var_names if g in set(query.var_names)]].copy()
    query.layers["counts"] = query.X.copy()
    scvi.model.SCANVI.prepare_query_anndata(query, scanvi_ref)
    scanvi_q = scvi.model.SCANVI.load_query_data(query, scanvi_ref)
    scanvi_q.train(max_epochs=100, plan_kwargs={"weight_decay": 0.0})

    probs = scanvi_q.predict(soft=True)            # DataFrame indexed by obs_names
    p = probs.to_numpy(dtype=float)
    K = p.shape[1]
    with np.errstate(divide="ignore", invalid="ignore"):
        ent = -(p * np.where(p > 0, np.log(p), 0.0)).sum(1) / np.log(max(K, 2))
    out = pd.DataFrame({
        "cell_id": np.asarray(probs.index).astype(str),
        "label": np.asarray(probs.columns)[p.argmax(1)].astype(str),
        "confidence": p.max(1),
        "entropy": ent,
    })
    import pyarrow as pa
    import pyarrow.parquet as pq
    pq.write_table(pa.Table.from_pandas(out, preserve_index=False), out_path)  # avoid pandas.to_parquet pyarrow-24 patch bug
    print(out_path)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
