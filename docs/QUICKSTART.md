# Quickstart

## 1. Install

SpatialScribe uses [pixi](https://pixi.sh) for the environment.

```bash
git clone <repo> && cd spatial-scribe
pixi install -e main            # CPU stack (scanpy/squidpy/decoupler/celltypist)
```

No GPU is required - the app runs CPU-only out of the box (`backend.py` auto-detects a GPU
and falls back to CPU). To force CPU explicitly (the judge-reproducibility gate):

```bash
export SPATIALSCRIBE_FORCE_CPU=1
```

**GPU (optional, for large sections):** on a CUDA node, install `rapids-singlecell`
(`conda install -c rapidsai -c conda-forge rapids-singlecell`) into the environment; the
backend picks it up automatically.

## 2. Set your Claude key

The copilot and Claude-based annotation need an API key (never commit it):

```bash
export ANTHROPIC_API_KEY=sk-...
export ANTHROPIC_MODEL=claude-sonnet-5   # optional override
```

## 3. Run the app (React SPA + FastAPI)

```bash
# backend - serves /api (PYTHONPATH so the src-layout package resolves)
PYTHONPATH=.:src pixi run python -m uvicorn backend.app:app --port 8000
# frontend - vite dev server on :5173, proxies /api -> backend
cd webapp && npm install && npm run dev
```

Open http://localhost:5173 and click **Load synthetic melanoma (instant demo)** (always works,
no data files) → walk the rails (Panel check → QC → Cluster → Annotate → Spatial exploration →
Report), or ask the copilot on the right in plain English. Production serves the built
`webapp/dist` single-origin from FastAPI on one port (see the internal docs).

## 4. Verify the install

```bash
pixi run smoke                  # imports OK
pixi run test                   # unit + synthetic end-to-end pipeline
```

## Demo dataset

The example is the public, CC BY 4.0
[10x FFPE Human Breast (Xenium Prime 5K)](https://www.10xgenomics.com/datasets).
Point `SPATIALSCRIBE_DEMO_CACHE` at a processed `.h5ad` (build one from a raw download with
`scripts/build_demo_cache.py`), or load any Xenium output folder / `.h5ad` from the Load step.
