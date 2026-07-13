# Frontend stack: decision + roadmap

> **SUPERSEDED (2026-07-09):** the full React rewrite this doc treats as a "maybe later" already
> happened. The shipped app is a React + Vite + deck.gl SPA over a FastAPI backend; Streamlit was
> removed entirely. Read this as the decision journal that led there (see the 2026-07-09 note under
> "Phase 1"), not the current recommendation.

## TL;DR

Keep **Streamlit** for the app shell, wizard, forms, tables and copilot chat. Move only the one
layer that fights Streamlit - the **interactive specimen viewport** - to a real WebGL component
(**deck.gl**), embedded as a Streamlit custom component. Do a full React rewrite *only* if
SpatialScribe becomes a maintained product, and when you do, keep the Python analysis engine as the
API and use **Vitessce** for the viewport. The clean `analysis/` <-> UI split already in the repo is
what makes any of these a port, not a rewrite.

## Why not "just switch to React"

The valuable half of SpatialScribe is Python and stays Python regardless of the frontend:

- `spatialscribe.analysis` - the capability registry (scanpy / squidpy / rapids, GPU), the single
  source of truth for both the wizard and the copilot.
- the AnnData in-memory data model.
- the Claude copilot tool-use loop.

So "switch to React" really means "stand up a FastAPI service around the registry + a stateful
session that holds the 100k-1M-cell AnnData + a React SPA." That re-implements everything Streamlit
gives for free (session-held data, deploy, state wiring) **and** adds a frontend - a bad trade for a
solo dev on a deadline, and a heavier long-term maintenance/skills burden for a Python-first lab.

## Where Streamlit actually hurts (and where it doesn't)

| Pain | Cause | A WebGL component fixes it |
|---|---|---|
| Floating HUD over the canvas | Streamlit can't overlay widgets on a chart | yes (absolute positioning in-component) |
| Rail tracker / bespoke layout | limited widget composition | (solved here with CSS on buttons) |
| Copilot-drives-map rerun juggling | whole-script reruns, session lag | yes (client-held view state) |
| 1M-cell rendering + snappy box-select | Plotly WebGL + rerun round-trips | yes (deck.gl, client-side interaction) |

The pain is concentrated in the **viewport**. Forms, the wizard, tables and the chat are fine in
Streamlit. So the surgical fix is a viewport component, not a stack change.

## Phase 1 (now): deck.gl specimen viewport as a Streamlit custom component

- `frontend/` - React + Vite + deck.gl + `streamlit-component-lib`. Renders the cells as a deck.gl
  `ScatterplotLayer` (WebGL, scales past 10^6 points), with a floating glass HUD *in the component*
  (trivial in React) and native zoom/pan. box selection posts the selected cell indices back
  to Python via `Streamlit.setComponentValue` - preserving the H1 region-QC round-trip.
- `app/specimen_deck.py` - `declare_component(...)` wrapper: `render(adata, color_by) -> selection`.
- NOTE (2026-07-09): superseded. Streamlit was removed; the deck.gl WebGL canvas is now the ONLY
  canvas, as a native React component in `webapp/src/SpecimenCanvas.tsx` (no `SPATIALSCRIBE_CANVAS`
  flag, no Plotly fallback, no Streamlit component bridge). The notes above are historical.
- Build hygiene: `node_modules` lives on `/data` (home has a 100 GB quota) and is git-ignored;
  only the small source + the built `dist/` bundle are committed, so the deployed app needs no Node.

This is the ~80/20: it removes the viewport pain (overlay HUD, fast rendering, client-side
interaction, copilot-driven recolour without reruns) at ~20% of a rewrite, and everything else stays.

## Phase 2 (only if it becomes a product): React + FastAPI + Vitessce

- **Frontend:** Next.js (React).
- **Backend:** FastAPI wrapping the existing capability registry + a session/task layer, with
  SSE/WebSocket for streaming copilot tokens and SLURM/GPU job progress.
- **Viewport:** the WebGL scatter primitive - **regl-scatterplot** (purpose-built for millions of 2D
  points with native box-select + categorical/continuous colouring; basically made for single-cell
  embeddings) or **deck.gl** `ScatterplotLayer` (more general, GPU picking, easy to add boundary /
  niche layers). Phase 1 already uses deck.gl; benchmark regl-scatterplot for the pure-embedding case
  if the viewport becomes the sole focus. For a whole multi-view framework (linked views + Zarr
  streaming) rather than a single scatter, **Vitessce** is the domain-correct choice and connects
  directly to the DuckDB / PyVips / Zarr out-of-core idea - don't hand-roll what the field solved.
- The copilot stays Python (tool-use loop over `capabilities.copilot_tools()`), exposed as a
  streaming endpoint.

Migrating is a port because `analysis/` is already UI-agnostic: the wizard and the copilot both call
`cap.run(...)`; a React SPA would call the same capabilities over HTTP.

## Decision criteria

- **Who maintains it?** Python-first lab / the cluster a reverse proxy shop -> stay Python-native unless there is
  a frontend owner. A React + FastAPI stack is real ongoing ops.
- **Demo or product?** A demo does not justify the rewrite; a lab's daily tool does.
- **Is the viewport the product?** For spatial-omics, yes - so that is the one place worth real WebGL
  (deck.gl now, Vitessce later), whether inside Streamlit or a React app.

## Also-rans (more flexible than Streamlit, still Python)

- **Panel + Bokeh + Datashader** - if a JS frontend is off the table, this is the strongest
  all-Python path: Datashader is the standard for rendering millions of points server-side with
  pan/zoom, and Panel gives far more layout freedom than Streamlit. You still won't get the bespoke
  HUD/drawer feel as cleanly as React, but you stay in one language.
- **Dash** (React under the hood) - more layout control, more boilerplate; still need a custom
  component for the WebGL canvas.
- **Reflex** - pure-Python compiled to Next.js; promising but younger, and migrating is still a
  rewrite. None clearly beats "Streamlit now, deck.gl (or regl-scatterplot) component next, full
  React only as a product."
