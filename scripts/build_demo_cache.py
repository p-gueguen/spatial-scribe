"""Precompute a processed demo object so the app + demo are instant.

Loads the real public Prime 5K breast Xenium bundle, subsamples for CPU responsiveness,
runs the full pipeline (QC -> cluster -> marker annotation -> confidence -> spatial graph),
and writes a processed ``.h5ad`` the app loads directly.

Run:  pixi run -e main python scripts/build_demo_cache.py
"""

from __future__ import annotations

import os
import sys

os.environ.setdefault("SPATIALSCRIBE_FORCE_CPU", "1")

RAW = os.environ.get(
    "SPATIALSCRIBE_DEMO_RAW",
    "/data/user/spatial_platform_comparison/data/prime5k_breast",
)
OUT = os.environ.get(
    "SPATIALSCRIBE_DEMO_CACHE",
    "/data/spatial-scribe/data/demo_breast_processed.h5ad",
)
N_SUB = int(os.environ.get("SPATIALSCRIBE_DEMO_NSUB", "100000"))
TISSUE = "breast"


def main() -> None:
    import numpy as np

    from spatialscribe.analysis import annotate, cluster, io, markers, panel_check, qc, spatial

    print(f"[cache] loading {RAW}", flush=True)
    sample = io.load(RAW)
    a = sample.adata
    a = a[:, ~sample.control_mask].copy()          # drop control probes from the matrix
    print(f"[cache] {a.n_obs:,} cells x {a.n_vars:,} genes", flush=True)

    if a.n_obs > N_SUB:                             # subsample for CPU responsiveness
        rng = np.random.default_rng(0)
        idx = rng.choice(a.n_obs, N_SUB, replace=False)
        a = a[idx].copy()
        print(f"[cache] subsampled to {a.n_obs:,} cells", flush=True)

    mset = markers.for_tissue(TISSUE)
    print("[cache] QC", flush=True); qc.compute_qc(a)
    print("[cache] panel check", flush=True)
    a.uns["panel_check"] = panel_check.check_panel(a.var_names.tolist(), marker_sets=mset)
    print("[cache] cluster (this is the slow step on CPU)", flush=True)
    cluster.cluster(a, resolution=0.5)   # match the app default so the demo opens on the same view
    print("[cache] annotate (marker consensus, no LLM for the cache)", flush=True)
    annotate.consensus_annotate(a, cluster_key="leiden", use_llm=False,
                                context="human breast carcinoma")
    # re-run marker labels with the tissue set, then confidence. Low-signal cells (leiden NaN) map to
    # NaN -> label "Unassigned" so cell_type has no NaN category (else a "nan" type in the legend +
    # a squidpy nhood crash); apply_confidence still abstains them as "Unassigned: low quality".
    lab = annotate.marker_labels(a, cluster_key="leiden", marker_sets=mset)
    a.obs["cell_type"] = a.obs["leiden"].astype(str).map(lab).fillna("Unassigned").astype("category")
    annotate.apply_confidence(a, cluster_key="cell_type", marker_sets=mset,
                              panel_check_result=a.uns["panel_check"])
    print("[cache] spatial neighbors", flush=True)
    spatial.spatial_neighbors(a)

    a.uns["spatialscribe_demo"] = {"tissue": TISSUE, "source": "10x Prime5K breast (public)",
                                   "n_subsampled": int(a.n_obs)}
    # Persist the platform + panel-size stamps (io.load set them; controls were dropped so
    # n_vars == panel size) so the loaded cache drives the panel-indexed floor / confidence gate.
    a.uns["platform"] = "xenium"
    a.uns["n_panel_genes"] = int(a.n_vars)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    # Use export_h5ad, not raw write_h5ad: it sanitizes uns dict keys (panel_check.coverage is keyed
    # by cell types like 'Epithelial/Tumor', and h5py forbids '/' in group keys).
    from spatialscribe.analysis import export

    export.export_h5ad(a, OUT)
    print(f"[cache] wrote {OUT}  ({a.n_obs:,} cells)", flush=True)


if __name__ == "__main__":
    sys.exit(main())
