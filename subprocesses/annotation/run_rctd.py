"""RCTD doublet-mode annotation -> per-cell parquet. Runs in the rctd-py env; GPU-optional.
Usage: python run_rctd.py <spatial.h5ad> <reference.h5ad> <ref_label_key> <out.parquet>
The DoubletResult is FILTERED-length; map back to cell_id via spatial.obs_names[result.pixel_mask].
first_type/second_type are integer indices into result.cell_type_names -> map to names.
"""
import sys

import anndata as ad
import numpy as np
import pandas as pd


def main(spatial_h5, ref_h5, ref_key, out_path):
    spatial = ad.read_h5ad(spatial_h5)
    import rctd
    from rctd import Reference, RCTDConfig, run_rctd

    ref = Reference(ad.read_h5ad(ref_h5), cell_type_col=ref_key)
    result = run_rctd(spatial, ref, mode="doublet", config=RCTDConfig(device="auto", UMI_min=20),
                      batch_size=5000)
    names = list(result.cell_type_names)
    ids = np.asarray(spatial.obs_names)[np.asarray(result.pixel_mask)].astype(str)

    def _name(idx_arr):
        return [names[int(i)] if i is not None and int(i) >= 0 else None for i in idx_arr]

    weights_doublet = np.asarray(getattr(result, "weights_doublet", np.empty((len(ids), 0))))
    weight = weights_doublet.max(1) if weights_doublet.size else np.full(len(ids), np.nan)
    out = pd.DataFrame({
        "cell_id": ids,
        "first_type": _name(np.asarray(result.first_type)),
        "second_type": _name(np.asarray(result.second_type)),
        "spot_class": [rctd.SPOT_CLASS_NAMES[int(c)] for c in np.asarray(result.spot_class)],
        "weight": np.asarray(weight, dtype=float),
        "singlet_score": np.asarray(getattr(result, "singlet_score", np.full(len(ids), np.nan)), dtype=float),
    })
    import pyarrow as pa
    import pyarrow.parquet as pq
    pq.write_table(pa.Table.from_pandas(out, preserve_index=False), out_path)  # avoid pandas.to_parquet pyarrow-24 patch bug
    print(out_path)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
