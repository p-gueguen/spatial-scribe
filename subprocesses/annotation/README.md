# Annotation subprocess envs

Each cell-type annotation method runs in its **own** env via a subprocess that writes a
`{cell_id, ...}` parquet the main app joins back onto `adata.obs`. Never import these
packages into the main interpreter - torch/TF/R (and their pinned dependency trees) clash
with the main scanpy/spatialdata stack, the same reason `subprocesses/ovrlpy/` is isolated.

All env content lives on `/data/spatial-scribe/` (home has a 100GB
quota) - create the venvs there, not under `~`. Each `run_<m>` in `methods.py` takes an
`env_python` path and no-ops (returns `None`) when it is not supplied, so the main pipeline
degrades gracefully with any subset of envs installed.

Create one `uv venv` per method on `a compute node`; GPU envs (rctd, scanvi) should additionally
be smoke-tested on a GPU node (a compute node / a compute node) before relying on them.

**Every** env additionally needs `anndata pyarrow scipy` (the subprocess reads the section
`.h5ad` and writes the `{cell_id, ...}` parquet the main app joins back) - add them to each
`uv pip install` below. This was verified by the live end-to-end run: SingleR annotated fine
but the write failed with `Unable to find a usable engine ... pyarrow` until pyarrow was added.

## rctd

GPU-optional. Runs `run_rctd.py`.

```bash
uv venv /data/spatial-scribe/rctd_env
uv pip install --python .../rctd_env/bin/python torch --index-url https://download.pytorch.org/whl/cu124
uv pip install --python .../rctd_env/bin/python rctd-py
# or, to run from a live checkout: uv pip install --python .../rctd_env/bin/python -e ~/git/rctd-py
```

## singler

CPU (BiocPy port of SingleR). Runs `run_singler.py`.

```bash
uv venv /data/spatial-scribe/singler_env
uv pip install --python .../singler_env/bin/python singler celldex
```

## scanvi

GPU (torch-backed scvi-tools). Runs `run_scanvi.py`.

```bash
uv venv /data/spatial-scribe/scanvi_env
uv pip install --python .../scanvi_env/bin/python scvi-tools
```

## panhuman

CPU (TensorFlow). Reference-free. Runs `run_panhumanpy.py`.

```bash
uv venv /data/spatial-scribe/panhuman_env
uv pip install --python .../panhuman_env/bin/python panhumanpy
```

## ref-fetch

CPU. Supports `scripts/fetch_reference.py` (reproducible CELLxGENE query for the
reference-based annotators above).

```bash
uv venv /data/spatial-scribe/ref_fetch_env
uv pip install --python .../ref_fetch_env/bin/python gget cellxgene-census tiledbsoma
```
