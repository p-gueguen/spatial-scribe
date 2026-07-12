# Architecture

SpatialScribe is a React + FastAPI + deck.gl app over a pure-function analysis engine. The engine
is the single source of truth shared by the wizard UI and the Claude copilot, so what the
scientist clicks and what they ask always run the same code.

```
src/spatialscribe/
  analysis/            # <- single source of truth (pure functions)
    io.py              # load() -> SpatialSample  (platform-agnostic ingestion)
    backend.py         # get_backend() -> GPU (rapids_singlecell) | CPU (scanpy)
    qc.py              # per-cell QC, thresholds, region-scoped QC (H1)
    panel_check.py     # panel-adequacy / resolvability (H3)
    cluster.py         # normalize/PCA/neighbors/Leiden/UMAP via backend
    annotate.py        # marker + Claude consensus -> obs['cell_type']
    markers.py         # curated lineage / TAM / cell-state panels
    spatial.py         # squidpy neighborhoods, immune-exclusion, Moran's I
    subcluster.py      # click-to-subcluster + subtype naming (H2)
    states.py          # cell-type x cell-state scoring (H6)
    llm.py             # grounded Claude helpers
    export.py          # annotated .h5ad + HTML report + re-runnable notebook
  agent/               # Claude Agent SDK tools wrapping analysis/ (the copilot, H4)
backend/               # FastAPI over the analysis/ capability registry (cap.run); server-side session
webapp/                # React + Vite + deck.gl SPA (thin): rails | WebGL canvas | streaming copilot
```

## The two contracts

**`SpatialSample`** (from `io.load`) is the one ingestion contract:
`platform, adata, control_mask, panel_genes, has_z, transcripts_path, sdata`. Detection is
by directory signature; the Xenium demo path reads `cell_feature_matrix.h5` + `cells.parquet`
directly (no `squidpy.read.xenium`; avoids the `spatialdata_io` centroid gap). Everything
downstream is platform-agnostic because it only sees a `SpatialSample`.

**`get_backend()`** is the compute contract: every heavy step (normalize/hvg/pca/neighbors/
leiden/umap) goes through it, so GPU (rapids_singlecell) vs CPU (scanpy) is one switch.
`SPATIALSCRIBE_FORCE_CPU=1` pins CPU for the judge-reproducibility gate.

## Data & control flow

1. `backend/` holds a `SpatialSample` + an **action log** server-side, one per session.
2. Each wizard step calls an `analysis/` function and appends `{label, code, params}` to the log.
3. The copilot (`agent/`) calls the *same* `analysis/` functions through a whitelisted tool
   schema bound to the session AnnData - no arbitrary code execution.
4. `export.py` turns the action log into a re-runnable notebook, plus an annotated `.h5ad`
   and a self-contained HTML report.

## Trust model

- Whitelisted tools only; optional guarded `run_python` is shown before it runs.
- Claude reasons only over computed numbers passed into the prompt (annotation markers,
  panel-coverage tables) - it never invents expression. Marker *presence* is framed as
  necessary-not-sufficient (imaging probes drop out).
- Data stays local; only cluster-level summaries go to the Claude API.

## The ovrlpy subprocess (H5)

`ovrlpy` pins clash with the main stack, so it runs in a separate env as a subprocess
(`subprocesses/ovrlpy/run_ovrlpy.py`) that consumes `transcripts.parquet` and writes a
per-cell VSI parquet keyed by `cell_id`; `qc.py` joins it back onto `obs`. Never imported
in the main interpreter.
