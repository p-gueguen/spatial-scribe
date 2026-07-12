# SpatialScribe - agent guide

Self-serve spatial-transcriptomics analysis copilot for wet-lab scientists. A three-pane
React + FastAPI + deck.gl app (guided rails + interactive spatial canvas + Claude copilot) over
one shared analysis engine. Takes a bench scientist from a raw Xenium / CosMx / MERSCOPE / `.h5ad`
section to QC'd, clustered, annotated, spatially-explored data and a shareable report - in plain
language, no terminal or R. GPU-accelerated with a CPU fallback.

This file is the entry point for AI coding agents. `CLAUDE.md` is a thin pointer to it.

## Architecture (respect these boundaries)

- **`src/spatialscribe/analysis/` is the single source of truth.** Pure functions used by BOTH the
  wizard buttons and the copilot, so rails and chat never diverge. Add logic here, not in the UI.
  - `io.py` - `load(path) -> SpatialSample`. The one ingestion contract (platform-agnostic).
  - `backend.py` - `get_backend()` GPU/CPU switch; route normalize/hvg/pca/neighbors/leiden/umap
    through it, never call rapids/scanpy directly in feature code.
  - `capabilities.py` - **the registry**: one declarative `Capability` per operation (pure fn +
    tool schema + `requires`/`produces` keys), consumed by BOTH the wizard and the copilot via
    `cap.run(...)`. Add copilot tools here, not in `agent/`. Conformance is enforced by
    `tests/test_capabilities_conformance.py`.
  - `pipeline.py` - the full spine as ONE degrade-graceful function (`run_pipeline`), driven by
    both `scripts/run.py` and the app's "Run full analysis" job.
- **`src/spatialscribe/agent/`** - `run_copilot` is a thin, endpoint-agnostic, whitelisted tool-use
  loop over `capabilities.copilot_tools()` (not arbitrary exec). It answers only from computed
  results and can render real plots onto the canvas.
- **LLM backend is endpoint-agnostic** (`llm.provider()`, chosen by env): Anthropic by default, or
  any OpenAI-compatible `/v1` server. Every LLM feature gates on `llm.available()`.
- **`backend/` + `webapp/` are the app.** `backend/app.py` is a thin FastAPI over the same registry
  (`cap.run`); a session holds the `SpatialSample` + an action log server-side (the log is the
  re-runnable export). `webapp/` is a React + Vite + deck.gl SPA. Dev = backend `:8000` + vite
  `:5173` (proxies `/api`); prod = FastAPI serves `webapp/dist` single-origin on `:8000`.
- **Out-of-env runners live under `subprocesses/<tool>/`** (annotation = RCTD / SingleR / scANVI /
  panhumanpy, plus CNV / Cancer-Finder / ovrlpy / reference fetch). Each is a standalone script run
  in its own env via `subprocess`, writing a `cell_id`-keyed parquet the main app joins back. They
  **degrade gracefully** (return `{"status": "skipped: ..."}`) when their env is absent - the app
  works fully without them.
- Heavy imports (scanpy, squidpy, rapids) go **inside** functions so the package imports cheaply.

## Run

```bash
pip install "spatial-anno-metrics @ git+https://github.com/p-gueguen/spatial-anno-metrics"
pip install -e ".[dev]"                 # or: pixi install -e main
export SPATIALSCRIBE_FORCE_CPU=1        # force the CPU path (reproducibility gate)
PYTHONPATH=.:src python -m uvicorn backend.app:app --port 8000   # backend
cd webapp && npm install && npm run dev                          # frontend (vite :5173)
python -m pytest -q                     # tests
```

The full spine, headless: `python scripts/run.py --demo --out results/`.

## Grounding and honesty (do not violate)

- **Never invent numbers or labels.** Report only figures a capability returned. Gene *presence* on
  a panel is not *detectability*.
- **Per-cell confidence is a heuristic** driving abstention; it is not calibrated. Abstention is
  honest, never a confident wrong label.
- **Reference choice matters more than method choice**, and a wrong-tissue reference can pass an
  overlap check yet annotate nonsense - match the tissue first. **Panel size gates label
  granularity** - coarsen for small panels.

## Conventions

- No emojis in product surfaces; state is shown with type / color / layout.
- `npm run build` (esbuild) does NOT typecheck - run `npx tsc --noEmit` for real type errors.

## The bundled skill

`.claude/skills/spatialscribe/` teaches a Claude Code agent to drive the engine end to end (the
headless run, the copilot, the panel/QC/reference checks). See its `SKILL.md`. The HTTP surface is
documented in `docs/API.md`.
