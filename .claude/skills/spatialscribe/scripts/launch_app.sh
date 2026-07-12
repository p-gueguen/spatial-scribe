#!/usr/bin/env bash
# Launch the SpatialScribe app (React SPA + FastAPI, single-origin) in the user's own environment.
# Runs on the user's own data + Anthropic key; nothing phones home.
set -euo pipefail

if [ -z "${ANTHROPIC_API_KEY:-}" ]; then
  echo "note: ANTHROPIC_API_KEY not set - the app runs, but the Claude copilot is disabled." >&2
fi

PORT="${SPATIALSCRIBE_PORT:-8000}"

# FastAPI serves the built SPA (webapp/dist) + /api on one port. Build it once on first launch.
if [ ! -f webapp/dist/index.html ]; then
  echo ">> building the web UI (webapp/dist) - first launch only" >&2
  ( cd webapp && npm install && npm run build )
fi

echo ">> SpatialScribe -> http://localhost:${PORT}" >&2
if command -v pixi >/dev/null 2>&1 && [ -f pixi.toml ]; then
  exec env SPATIALSCRIBE_PORT="$PORT" pixi run python -m uvicorn backend.app:app --host 0.0.0.0 --port "$PORT"
else
  exec env PYTHONPATH="${PYTHONPATH:-}:.:src" python -m uvicorn backend.app:app --host 0.0.0.0 --port "$PORT"
fi
