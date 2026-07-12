"""Fetch a small CELLxGENE reference (skin / melanoma) via gget - reproducible, CC-licensed.

Run ONCE in a dedicated env (gget + cellxgene-census + tiledbsoma), output cached on /data
(NOT committed - it is data). The committed artifact is this script (the reproducible query).

    uv pip install gget cellxgene-census tiledbsoma      # avoids gget.setup() shelling bare pip
    python scripts/fetch_reference.py --out /data/spatial-scribe/reference_skin.h5ad
"""
from __future__ import annotations

import argparse


def build_query(tissue: str = "skin of body", disease: str | None = None,
                census_version: str = "2025-01-30") -> dict:
    """Reproducible gget.cellxgene query. Pin census_version; keep the CL ontology id."""
    q = {
        "species": "homo_sapiens",
        "tissue": tissue,
        "census_version": census_version,          # pin -> reproducible
        "is_primary_data": True,
        "column_names": ["cell_type", "cell_type_ontology_term_id", "assay", "suspension_type"],
    }
    if disease:
        q["disease"] = disease
    return q


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--tissue", default="skin of body")
    ap.add_argument("--disease", default=None)
    ap.add_argument("--census-version", default="2025-01-30")
    ap.add_argument("--max-cells", type=int, default=30000)
    a = ap.parse_args()

    import gget

    adata = gget.cellxgene(**build_query(a.tissue, a.disease, a.census_version))
    if adata.n_obs > a.max_cells:                    # downsample to a small demo reference
        import numpy as np
        idx = np.random.default_rng(0).choice(adata.n_obs, a.max_cells, replace=False)
        adata = adata[idx].copy()
    adata.write_h5ad(a.out)
    print(a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
