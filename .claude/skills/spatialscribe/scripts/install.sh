#!/usr/bin/env bash
# Install the SpatialScribe engine into the user's OWN environment and verify it.
# Runs on the user's compute + data + Anthropic key; data never leaves the machine.
set -euo pipefail

REPO="${1:-https://github.com/p-gueguen/spatial-scribe}"
DEST="${2:-spatial-scribe}"

if [ ! -d "$DEST/.git" ]; then
  echo ">> cloning $REPO -> $DEST"
  git clone "$REPO" "$DEST"
fi
cd "$DEST"

if command -v pixi >/dev/null 2>&1; then
  echo ">> pixi install -e main (manages the scanpy/squidpy/streamlit/celltypist stack)"
  pixi install -e main
  echo ">> verify"
  pixi run smoke && pixi run test
else
  echo ">> pixi not found; installing into the active Python env with pip"
  # spatial-anno-metrics is declared by NAME in pyproject (not on PyPI), so pull it from git first;
  # `pip install -e .` then sees it already satisfied.
  python -m pip install "spatial-anno-metrics @ git+https://github.com/p-gueguen/spatial-anno-metrics"
  python -m pip install -e .
  echo ">> verify"
  python -c "import scanpy, squidpy, anndata, spatial_anno_metrics; print('main env OK')"
  python -m pytest -q
fi

echo
echo ">> engine installed in $(pwd)"
echo ">> set your key:  export ANTHROPIC_API_KEY=sk-...   (optional: ANTHROPIC_MODEL=claude-sonnet-5)"
echo ">> GPU is optional (rapids-singlecell on a CUDA node; else CPU). SPATIALSCRIBE_FORCE_CPU=1 forces CPU."
