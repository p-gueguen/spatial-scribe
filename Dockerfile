# SpatialScribe - self-contained CPU image. Builds the React SPA, installs the CPU engine, and
# serves the SPA + /api single-origin on :8000. GPU is optional and not required.

# ---- stage 1: build the web UI ----
FROM node:22-slim AS webbuild
WORKDIR /app/webapp
COPY webapp/package.json webapp/package-lock.json ./
RUN npm ci
COPY webapp/ ./
RUN npm run build

# ---- stage 2: python runtime ----
FROM python:3.11-slim
WORKDIR /app

RUN apt-get update \
 && apt-get install -y --no-install-recommends git build-essential \
 && rm -rf /var/lib/apt/lists/*

# spatial-anno-metrics is a git dependency (declared by name in pyproject); install it first.
RUN pip install --no-cache-dir \
    "spatial-anno-metrics @ git+https://github.com/p-gueguen/spatial-anno-metrics"

COPY . /app
COPY --from=webbuild /app/webapp/dist /app/webapp/dist

# Editable install so the code's source-relative paths (docs/research/*.yaml, subprocesses/) resolve.
RUN pip install --no-cache-dir -e .

ENV PYTHONPATH=/app:/app/src \
    SPATIALSCRIBE_FORCE_CPU=1

EXPOSE 8000
CMD ["python", "-m", "uvicorn", "backend.app:app", "--host", "0.0.0.0", "--port", "8000"]
