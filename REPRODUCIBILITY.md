# Reproducibility

SpatialScribe is built so that a given section produces the same analysis on any
machine, with or without a GPU, and with or without an Anthropic API key. This
document explains the guarantees and how to reproduce a run end to end.

## 1. Environment pinning

The core environment is managed by pixi and pinned exactly by `pixi.lock`. That
lockfile records the resolved versions of every conda and PyPI dependency, so:

```bash
pixi install -e main      # reproduces the exact core (CPU/dev) environment - the only pixi env
# GPU is NOT a pixi env: conda install -c rapidsai -c conda-forge rapids-singlecell on a GPU node
# (backend.get_backend() auto-detects it; the result is identical, CPU is the fallback)
```

`pixi.toml` declares the dependency ranges and the task shortcuts; `pixi.lock` is
the source of truth for the exact versions everyone gets. The metrics battery is a
separate package, `spatial-anno-metrics`, wired in as a git-URL dependency;
`analysis/eval_metrics.py` and `analysis/signal_qc.py` are thin re-export shims over
it, so every annotation-quality number flows through one pinned implementation.

Some steps run in **isolated subprocess environments** that are intentionally kept
out of the main resolver because their dependencies clash with the scanpy 1.12 /
spatialdata stack. These are pinned and configured separately:

- **ovrlpy** (VSI contamination signal): its own throwaway env; `subprocesses/ovrlpy/`
  runs it as a subprocess that writes a per-cell parquet the main app joins back.
- **Annotation / CNV methods** (RCTD, SingleR, scANVI, panhumanpy, infercnvpy,
  Cancer-Finder): each runs as an isolated subprocess whose interpreter and model
  paths are supplied via environment variables (for example `CANCERFINDER_PYTHON`,
  `CANCERFINDER_REPO`, `CANCERFINDER_CKPT`, and the `SPATIALSCRIBE_CNV_*` /
  `SPATIALSCRIBE_REF_*` variables). All of these degrade gracefully to a `skipped`
  status when their env is absent, so the core pipeline stays reproducible without them.

## 2. Determinism knobs

- **`SPATIALSCRIBE_FORCE_CPU=1` forces the CPU path** regardless of hardware. This is
  the judge-reproducibility gate: set it and the same section yields the same result
  on any machine, because the GPU backend is skipped entirely.

  ```bash
  export SPATIALSCRIBE_FORCE_CPU=1
  ```

- **The GPU path is an accelerator with a CPU fallback, not a separate result.**
  `analysis/backend.get_backend()` auto-selects `rapids_singlecell` only when a CUDA
  device and the GPU packages are present, and otherwise returns the plain-scanpy CPU
  backend. `SPATIALSCRIBE_FORCE_CPU=1` short-circuits the check. Leiden clustering is
  deliberately run on CPU (scanpy igraph) even on the GPU path, so cluster labels are
  identical to the reproducibility path; the GPU is used only for the behaviour-neutral
  wins (neighbors, UMAP). See the internal docs.

- **Random seeds.** Every stochastic step is seeded (`seed=0` / `random_state=0` /
  `np.random.default_rng(0)`): the stability bootstrap, niche k-means, NMF programs,
  reference subsampling, and the synthetic demo generator all
  fix their seed, so repeated runs match.

## 3. The demo dataset

The demo is a **public, CC BY 4.0** 10x Xenium Prime 5K section (FFPE Human Skin
Primary Dermal Melanoma, with a breast Xenium Prime 5K variant used for the processed
cache). No proprietary data is required.

Two demo entry points exist:

- **Synthetic melanoma demo** (`analysis/demo.py`, always available): a small,
  seeded, Xenium-like AnnData generated in-process. Instant and fully deterministic;
  no download or cache needed.
- **Processed real-data cache** (optional): a subset of the public Prime 5K bundle,
  run once through the full pipeline (QC, cluster, marker annotation, spatial graph)
  by `scripts/build_demo_cache.py` and written to an `.h5ad`. The app reads this cache
  from the `SPATIALSCRIBE_DEMO_CACHE` environment variable only - no path is hardcoded.
  When the variable is unset the app degrades to a clear message and the synthetic demo
  plus "load your own section" remain available.

```bash
export SPATIALSCRIBE_DEMO_CACHE=/path/to/your/processed_demo.h5ad
```

To regenerate the reference atlas used for reference-based annotation, `scripts/fetch_reference.py`
is the deterministic, CC-licensed fetch: it runs a pinned CELLxGENE Census query
(`census_version` is pinned) via `gget` and downsamples with a fixed seed. Only the
script (the reproducible query) is committed; the fetched atlas is data, not committed.

## 4. LLM reproducibility

- The copilot and the Claude-as-annotator both default to **`claude-haiku-4-5`**
  (`llm.DEFAULT_MODEL`). Override with `ANTHROPIC_MODEL` (for example
  `ANTHROPIC_MODEL=claude-sonnet-5`).
- The API key is read from **`ANTHROPIC_API_KEY` (env var only)** - never a literal in
  code or config.
- **The pipeline runs fully without an API key.** Annotation falls back to marker-only
  scoring (`consensus_annotate(..., use_llm=False)` when no key is set), and the copilot
  is simply disabled. So the entire core analysis - load, QC, cluster, annotate, spatial,
  export - is reproducible offline with no LLM in the loop. The LLM only adds a
  plain-language layer on top of numbers that are computed deterministically; it is never
  allowed to invent values (unknown keys return an honest note).

## 5. Reproducing a run

1. Point at a section: a Xenium / CosMx / MERSCOPE output folder or an `.h5ad`
   (`io.load()` returns a platform-agnostic `SpatialSample`).
2. Either walk the guided wizard in the app (`pixi run serve`, or the vite dev flow in the README), or run the exported
   re-runnable script. Every action is recorded in the action log, and
   `export.build_runnable_script(action_log, adata=..., source_path=...)` emits a
   complete, standalone Python script that reproduces the shown pipeline from the
   section loader onward.
3. Durable artifacts of any run:
   - the annotated `.h5ad` (all obs columns: `cell_type`, confidence, verdict, niche);
   - a self-contained HTML report (figures embedded, plus the re-runnable script appended);
   - the re-runnable script itself.

Because the script is derived from the action log and reads the section from a path you
supply, anyone with the same input and the pinned environment reproduces the same
outputs.

## 6. Tests and cross-dataset consistency as reproducibility guards

```bash
pixi run test     # full pytest suite
pixi run check    # cross-dataset invariants (scripts/check_consistency.py)
```

- `pixi run test` runs the whole suite, including registry conformance (every capability
  is covered by construction) and the synthetic-fixture behaviour tests.
- `pixi run check` enforces invariants across the dataset roster (Xenium 480 / 5K, CosMx,
  MERSCOPE, Atera-WTA): panel-size resolution, platform-aware QC flagging, and the
  panel-indexed count floor must behave consistently on all of them. Adding a dataset
  means adding a profile row, so drift is caught.

**Non-interactive shell gotcha.** A background or non-login shell may not have `pixi` on
`PATH`, which makes the wrapper exit non-zero (or a trailing `grep`/`tail` can become the
reported exit code) and report a false green or false red. Confirm a run by reading
pytest's own `N passed` summary line, not the shell / background exit code. Put pixi on
`PATH` first (or `module load` it) if the environment does not do so automatically.

## 7. Performance numbers

Performance figures are measured and maintained separately - see
the internal docs (GPU vs CPU pipeline timings and scalability)
and the internal docs (per-stage wall-clock profile). They are not
duplicated here so there is a single source of truth for each number.
