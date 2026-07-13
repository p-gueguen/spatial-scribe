# SpatialScribe HTTP API

The backend (`backend/app.py`) is a thin FastAPI layer over the UI-agnostic `spatialscribe.analysis` capability registry. A **session** holds one spatial sample (a `SpatialSample` / AnnData) plus an append-only action log entirely server-side; the browser (or your agent) never touches the AnnData - it asks for viewport points, capability results, verdicts, and copilot answers. Every endpoint runs a capability via `cap.run` (or the copilot's Claude tool-use loop) and returns JSON. Data stays on the server: you drive the analysis by passing a session id around and pulling back small payloads.

## Base URL and prefix

- Base URL: `http://localhost:8000` (single origin - the built SPA is served at `/`, and `/api/*` routes take precedence).
- All API routes are under the `/api` prefix.
- Run it with: `uvicorn backend.app:app --host 0.0.0.0 --port 8000`.
- No auth. The copilot endpoints additionally need an LLM configured in the server's environment (an Anthropic API key or an OpenAI-compatible base URL); check `llm.available` in any state response.
- Sessions are in-memory. A server restart or timeout drops them, and later calls return `404 "unknown session"` - just reload the section.

## The session id

Every session-creating endpoint (`POST /api/load_demo`, `POST /api/load_synthetic`, `POST /api/load_section`) returns a JSON object whose session id lives under the key **`session_id`**. Everything else is addressed as `/api/{sid}/...`, where `sid` is that value. `region_filter` and `rename_celltype` also echo `session_id` back (they replace the session's AnnData in place). Session-create responses spread the full `state` object (see below) alongside `session_id`.

## End-to-end curl walkthrough

```bash
BASE=http://localhost:8000

# 1. Create a session from the bundled demo. The response carries session_id + initial state.
SID=$(curl -s -X POST "$BASE/api/load_demo?name=breast" | jq -r .session_id)
echo "session: $SID"

# 2. Read the per-step summary (only what has been computed so far).
curl -s "$BASE/api/$SID/summary" | jq

# 3. Run a rail-step capability (panel | qc | cluster | annotate | spatial).
curl -s -X POST "$BASE/api/$SID/run/qc" \
     -H 'Content-Type: application/json' -d '{}' | jq '{ok, error}'

# ...or run ANY registered capability by name, getting its value + artifacts back.
curl -s -X POST "$BASE/api/$SID/run_cap/qc_funnel" \
     -H 'Content-Type: application/json' -d '{"params": {}}' | jq '.value'

# 4. Ask the copilot (needs an LLM configured on the server).
curl -s -X POST "$BASE/api/$SID/copilot" \
     -H 'Content-Type: application/json' \
     -d '{"prompt": "colour the map by cell type and tell me what you see"}' \
   | jq '{reply, map_view, mutated}'

# 5. Export. Each is a file download; -OJ keeps the server-supplied filename.
curl -s -OJ "$BASE/api/$SID/export/h5ad"     # annotated.h5ad
curl -s -OJ "$BASE/api/$SID/export/script"   # analysis.py (re-runnable)
curl -s -OJ "$BASE/api/$SID/export/report"   # spatialscribe_report.html
```

## The `state` object

Returned by `GET /api/{sid}/state`, spread into every session-create response, and nested as `state` in run/copilot/pipeline responses. Key fields:

| Field | Meaning |
|-------|---------|
| `n_obs`, `n_vars` | cells and genes in the current section |
| `tissue` | free-text tissue label driving marker/reference choice |
| `is_tumour` | `true`/`false`/`null` (null = decide from the tissue keyword) |
| `obs_fields` | colourable `obs` columns (feed these to `points?color_by=`) |
| `done` | data-derived: which step outputs exist (`load, panel, qc, cluster, annotate, spatial, report`) |
| `ran` | steps the user actually ran this session (the left rail greens off this) |
| `device` | `"GPU"` or `"CPU"` |
| `has_key`, `has_stt` | LLM configured; speech-to-text installed |
| `llm` | `{available, provider, model, providers}` |
| `panel_name` | e.g. `"hMulti_100g (Human, Multi, 480 targets)"` or null |
| `demo` | `{source, role}` for bundled sections, else null |
| `reference` | `{name, n_cells, label_key, n_types}` or null |

## Endpoints

Base path `/api`. `{sid}` is the `session_id`. Request bodies are JSON unless noted.

### Session / Load

| Method | Path | Purpose | Key params / response |
|--------|------|---------|-----------------------|
| GET | `/demos` | Bundled demo sections openable on this server | -> `{demos: [{name, available}]}` |
| POST | `/load_demo` | Open a bundled demo section, create a session | query `name` (default `breast`) -> `{session_id, ...state}` |
| POST | `/load_synthetic` | Open the offline synthetic melanoma demo (no data file) | -> `{session_id, ...state}` |
| POST | `/load_section` | Ingest a real section from a server-side path (Xenium/CosMx/MERSCOPE folder or `.h5ad`) | body `{path, tissue?, is_tumour?}` -> `{session_id, ...state}` |
| GET | `/{sid}/state` | Current session state | -> `state` |
| POST | `/{sid}/is_tumour` | Set the "contains malignant cells?" flag (gates Cancer-Finder / InSituCNV at annotate) | body `{is_tumour}` (bool or null) -> `{ok, ...state}` |
| POST | `/{sid}/reference` | Register a single-cell reference by server-side `.h5ad` path | body `{path, label_key?, gene_name_col?}` -> `{ok, n_ref_cells, label_key, n_labels, match, ...state}` |
| POST | `/{sid}/reference/upload` | Register a reference by uploading a `.h5ad` (multipart) | form field `file`; query `label_key?` -> same as `/reference` |
| POST | `/{sid}/reference/auto` | Free-text tissue -> auto-choose + load the best pre-computed reference (or CELLxGENE fetch) | body `{tissue?, allow_fetch?}` -> `{ok, auto, match, recommended_mode, route, ...state}`; `ok:false` + `recommended_mode:"cluster"` when none fits |

### Analysis / Capabilities

| Method | Path | Purpose | Key params / response |
|--------|------|---------|-----------------------|
| POST | `/{sid}/run/{step}` | Run a rail step: `panel`, `qc`, `cluster`, `annotate`, `spatial` | body `{params?}` -> `{ok, error, state}` |
| POST | `/{sid}/run_cap/{name}` | Run ANY registered capability by name (e.g. `qc_funnel`, `immune_exclusion`, `neighborhood_enrichment`, `state_by_celltype`, `malignant_score`, `discover_programs`, `subcluster`, `rejection_reasons`) | body `{params?}` -> `{ok, error, value, artifacts, state}` |
| POST | `/{sid}/ran/{step}` | Mark a rail step ran (greens the rail after a panel's auto-compute); idempotent | -> `state` |
| GET | `/{sid}/progress` | Coarse progress of the step running now (pollable mid-run) | -> `{running, step, frac, label}` |
| POST | `/{sid}/run_pipeline` | Start the FULL analysis spine as a background job (202) | body `{tumour?, rctd?, split?, resolution?}` -> `{started: true}` |
| GET | `/{sid}/pipeline_status` | Progress + outcome of the background full run | -> `{running, frac, label, stage, done, error, route, stages, summary, state}` |
| POST | `/{sid}/region_qc` | QC summary of a box-selected region vs the whole section | body `{indices: int[]}` -> `{n, region, section}` |
| POST | `/{sid}/region_filter` | Drop (`exclude`) or crop-to (`keep`) selected cells, in place | body `{indices: int[], mode}` -> `{session_id, ...state}` |
| POST | `/{sid}/rename_celltype` | Rename a cell type (reuse a name to merge), then recompute confidence | body `{old, new}` -> `{session_id, ...state}` |

### Views

| Method | Path | Purpose | Key params / response |
|--------|------|---------|-----------------------|
| GET | `/{sid}/points` | Per-cell map coordinates + colours for the canvas | query `color_by` (default `cell_type`), `max_points` (200000), `basis` (`spatial`\|`umap`), `colors_only`, `only_type`, `pos_min` -> `{x, y, ids, rgb, n, color_by, legend, basis, ramp}` (`colors_only=true` omits `x/y/ids`) |
| GET | `/{sid}/summary` | Per-step panel data computed from the current state (only what is present) | -> `{panel?, qc?, cluster?, annotate?, spatial?, report}` |
| GET | `/{sid}/panel_report` | Panel adequacy + reference concordance: coverage lights, identifiability AUC, confusable pairs, merge nudge | -> `{global, per_type, confusable, merge_nudge}` |
| GET | `/{sid}/panel_verdict` | Plain-language LLM panel verdict (cached per section/inputs) | -> `{verdict, typability_basis}` or `{verdict:null, note}` |
| GET | `/{sid}/qc_verdict` | Plain-language LLM QC verdict grounded in the QC funnel + rejections + reference match | -> `{verdict, reference_match}` or `{verdict:null, note, ...}` |
| GET | `/{sid}/verify_report` | Self-verify audit of the labels: marker-argmax agreement + one-vs-rest AUC + failing types | query `neighborhood` (bool) -> `{status, section_agreement, mean_auc, per_type, failed, suggestions}` |

### Copilot

| Method | Path | Purpose | Key params / response |
|--------|------|---------|-----------------------|
| POST | `/{sid}/copilot` | Ask the copilot; runs the Claude tool-use loop and may recolour the map or mutate labels | body `{prompt}` -> `{reply, map_view, figures, mutated, state}` |
| POST | `/{sid}/copilot/stream` | Same, streamed as Server-Sent Events over POST | body `{prompt}` -> `text/event-stream` of `data: {...}` frames |
| POST | `/stt` | Transcribe a dictated question (session-free; audio stays on the server) | multipart form field `audio` -> `{text}` (empty string = silence) |

**SSE frame shape** (`copilot/stream`): each event is `data: <json>\n\n`, where `<json>` has a `type` of:
- `status` -> `{type, text}` (per-tool progress ticks)
- `map_view` -> `{type, color_by, highlight}` (map-recolour directive)
- `figures` -> `{type, figures: [...]}` (rendered PNG artifacts)
- `token` -> `{type, text}` (the answer, accumulated word by word)
- `done` -> `{type, state, loaded, mutated}` (final state; `loaded` = a new section was swapped in, `mutated` = labels changed in place)

### Export

| Method | Path | Purpose | Response |
|--------|------|---------|----------|
| GET | `/{sid}/export/script` | Re-runnable Python script of the session's action log | `text/plain` attachment `analysis.py` |
| GET | `/{sid}/export/h5ad` | Annotated AnnData | `.h5ad` file `annotated.h5ad` |
| GET | `/{sid}/export/report` | Self-contained HTML report | `text/html` attachment `spatialscribe_report.html` |

### Misc

| Method | Path | Purpose | Response |
|--------|------|---------|----------|
| GET | `/health` | Liveness + live session count | `{ok, sessions}` |
| POST | `/llm/provider` | Switch the copilot LLM backend at runtime (process-wide) | body `{provider}` -> `{ok, llm}` |

## Drive it from an agent

A bundled skill teaches an agent to operate this API end to end (create a session, run the pipeline, inspect verdicts, ask the copilot, export): see **`.claude/skills/spatialscribe/`** in this repo (`SKILL.md` plus `scripts/` and `references/`). The pattern is always the same: `POST` a load endpoint, grab `session_id`, then address every later call as `/api/{session_id}/...` and poll `progress` / `pipeline_status` for long-running steps.
