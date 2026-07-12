# Contributing to SpatialScribe

Thanks for your interest in improving SpatialScribe, a self-serve spatial-transcriptomics
analysis copilot. This guide covers the dev setup, the one architectural contract that keeps
the project coherent, and how to add new functionality so that every frontend gets it for free.

## Dev setup

The project uses [pixi](https://pixi.sh) to manage its environments.

```bash
pixi install -e main     # install the main (CPU/dev) environment
pixi run smoke           # fast import check - confirms the stack loads
pixi run test            # pytest -q - the full test suite
pixi run check           # cross-dataset consistency invariants
export ANTHROPIC_API_KEY=sk-...
pixi run serve           # launch the app (FastAPI serves the built SPA + /api on :8000; see README for the vite dev flow)
```

### Two environments, and why you must never merge them

`ovrlpy` pins `polars` and `umap-learn` versions that clash with the scanpy / spatialdata
stack. They cannot coexist in one environment, so the project keeps them apart:

- **main env** (`pixi -e main`, the only pixi env): everything except ovrlpy. All application
  code, tests, and the app run here. GPU acceleration is NOT a pixi env - `conda install -c
  rapidsai -c conda-forge rapids-singlecell` on a GPU node, which `backend.get_backend()`
  auto-detects (CPU fallback otherwise).
- **isolated subprocess envs**: `ovrlpy` (and the optional learned malignant callers) run as
  **subprocesses in their own throwaway environments**. Each writes a per-cell parquet keyed by
  `cell_id` that the main app joins back. Never `import ovrlpy` (or the CNV/cancer-caller stacks)
  in the main interpreter, and never add their dependencies to the main env.

If you touch anything that shells out to a subprocess env, keep the boundary a parquet hand-off,
not a shared import.

## Architecture boundary contract (read this before writing code)

There is one rule that everything else follows:

> **All analysis logic lives in `src/spatialscribe/analysis/` as UI-agnostic pure functions.
> The frontends are thin and both drive the same capability registry.**

- `src/spatialscribe/analysis/` is the single source of truth. Functions here take an `adata`
  (and explicit params) and return plain, JSON-able values. They know nothing about React,
  HTTP, or the copilot.
- The frontend consumes that engine, and **both its guided rails and its copilot chat go
  through the capability registry via `cap.run(...)`**, so they can never diverge:
  - the React SPA in `webapp/` with its FastAPI service in `backend/` (thin: HTTP over the
    same registry; the session's `SpatialSample` + action log are held server-side). The
    legacy Streamlit app (`src/spatialscribe/app/`) was removed 2026-07-09.
- **Never put analysis logic in a UI.** If a button or an endpoint needs to compute something,
  the computation belongs in `analysis/`, exposed as a capability.
- **Heavy imports go inside functions.** Import scanpy, squidpy, and similar only inside the
  function that needs them, so the package itself imports cheaply.

See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for the full picture; this section is only the
contract you must not break.

## How to add a capability

A capability is one declarative entry that binds a pure function to a wizard step and a copilot
tool at the same time. Adding one is four small, local edits, and both frontends plus the copilot
pick it up automatically.

1. **Write the pure function** in the relevant `analysis/*.py` module (or reuse an existing one).
   Keep it UI-agnostic and put heavy imports inside it.
2. **Write a thin adapter** in `analysis/capabilities.py`:

   ```python
   def _cap_<name>(adata, ctx, **params):
       from . import <module> as _m       # heavy import stays inside
       return _m.<pure_function>(adata, ...)   # return a JSON-able value
   ```

   The adapter is the one place that translates `(adata, ctx, **params)` into your real call;
   the pure function keeps its own signature.
3. **Append a `Capability(...)` entry** to `_CAPABILITIES` in the same file, wiring up:
   - `name`, `label`, `description` (the description doubles as the copilot tool schema, so make
     it clear and grounded),
   - a JSON-schema `params` block and any `required_params`,
   - `requires` / `produces` keys (the explicit inter-stage data contract),
   - `copilot_exposed=True` if the copilot should be able to call it.
4. **Declare any new obs/uns/obsm keys** in `analysis/keys.py`. Every key referenced by a
   capability's `requires` / `produces` must be defined there - the conformance test enforces it.

That is it. `tests/test_capabilities_conformance.py` covers your capability by construction (it
walks the registry), and both frontends and the copilot get it through `cap.run(...)`. Do not
register copilot tools anywhere else.

**A new or reworded copilot tool can steal a neighbouring intent, and `pytest` cannot catch it**
(the suite mocks the LLM; routing lives in the model and is provider-specific). After you add or
change a `copilot_exposed` tool's `name` / `description` / `params`, run the per-provider routing
regression and confirm both the new intent routes AND the neighbours still do:
`pixi run python an internal LLM smoke-check --provider qwen` (and `--provider anthropic`). See
docs/LLM_EVAL.md.

**Annotation-quality scores are not yours to reimplement.** All spatial cell-type
annotation-quality metrics live in the separate `spatial-anno-metrics` package (a git-URL
dependency, see `pixi.toml`). `analysis/eval_metrics.py` and `analysis/signal_qc.py` are thin
re-export shims over it - import through them, and add metric logic to that package, never to the
shims.

## Testing conventions

- **TDD for bugs and features.** When fixing a bug, first write a test that reproduces it and
  fails, then make it pass. When adding a feature, add tests alongside it.
- **Keep the suite green.** Run `pixi run test` before you push; `pixi run check` guards the
  cross-dataset consistency invariants (add a dataset -> add its profile row, or the consistency
  test fails).
- **Confirm the result, not the exit code.** A non-interactive or background shell may not have
  `pixi` on `PATH`; the login profile that adds it is not always sourced. When that happens
  `pixi run test` exits **127** ("pixi: No such file or directory") and a wrapper can report a
  **false green**. Put pixi on `PATH` first (via your module system or the pixi install's `bin`
  directory), and verify a run by reading pytest's **`N passed`** summary line - never trust the
  shell or background exit code (a trailing `grep`/`tail`/`&& echo` becomes the reported exit and
  masks the real pytest result).

## Style

- **Match the surrounding code.** Follow the conventions already in the file you are editing.
- **Every function documents what it does, how to use it, and what it depends on** (the
  isolation principle used throughout `analysis/`).
- **Hyphens, not em-dashes** in user-facing prose (app text, reports, docstrings surfaced to
  users). Use `-` or `--`.
- **No emojis in product surfaces.** State is communicated with type, color, and layout, not
  emoji.

## Open-source hygiene (important)

This is a public repository. Committed code must contain **no site-specific or internal paths and
no credentials**.

- Any concrete data path must come from an **environment variable with an EMPTY public default**,
  e.g. `os.environ.get("SPATIALSCRIBE_REF_*", "")`, `SPATIALSCRIBE_CNV_*`,
  `SPATIALSCRIBE_DEMO_CACHE`, with a pointer in `docs/DATASETS.md`. Never hardcode an internal
  path as a string literal - and keep it out of comments too (leak scans are case-sensitive and
  match lowercase tokens).
- `ANTHROPIC_API_KEY` is read from the environment only. Never commit a key.
- Re-run the leak scan **after every new file lands**, not once: a fresh module can silently
  reintroduce an internal path long after an earlier scan was clean.
- Annotation-quality metrics live in a separate package, `spatial-anno-metrics`
  ([github.com/p-gueguen/spatial-anno-metrics](https://github.com/p-gueguen/spatial-anno-metrics));
  contribute metric changes there, not to the re-export shims here.

## Progress tracking

`PROGRESS.html` is a **generated** dashboard - do not hand-edit it. Edit the source of truth
`progress.yaml` and regenerate:

```bash
python the progress renderer            # progress.yaml -> PROGRESS.html
python the progress renderer --check    # exit 1 if PROGRESS.html is stale (CI gate)
```

Section percentages and progress bars are computed from each card's one-word `status`
(`done` / `partial` / `todo` / `roadmap`); never write percentages by hand.

The architecture docs are generated the same way, but from the **live code** rather than a YAML
file: `scripts/render_architecture.py` reads `capabilities.REGISTRY` and the `analysis/*.py`
docstrings, then writes `docs/CODE_ARCHITECTURE.html`, the module fence in `docs/ARCHITECTURE.md`,
the `<!-- generated:modules -->` block in `CLAUDE.md`, and the full capability catalog at
`.claude/skills/spatialscribe/references/capabilities.md`.

```bash
python scripts/render_architecture.py          # regenerate all three
python scripts/render_architecture.py --check  # exit 1 if any is stale (CI gate)
```

Both generators run automatically via a `PostToolUse` hook in `.claude/settings.json` when an agent
saves `progress.yaml` or `src/spatialscribe/analysis/*.py`, and `tests/test_render_architecture.py`
gates the committed outputs. Add a capability and the docs follow; you never update a count by hand.

## PR checklist

Before opening a pull request, confirm:

- [ ] `pixi run test` is green (checked via pytest's `N passed` line, not the shell exit code).
- [ ] `pixi run check` passes if you touched datasets, thresholds, or QC invariants.
- [ ] New analysis logic lives in `analysis/` (not in a UI), with heavy imports inside functions.
- [ ] New capabilities are registered in `capabilities.py` and any new keys are declared in `keys.py`.
- [ ] If you added/reworded a `copilot_exposed` tool: `an internal LLM smoke-check --provider qwen` (and
      `--provider anthropic`) still green - the new intent routes and no neighbour regressed.
- [ ] New behavior has tests; a bug fix has a regression test that failed before the fix.
- [ ] No internal paths or credentials committed; concrete paths come from env vars with empty
      public defaults. Leak scan is clean.
- [ ] User-facing text uses hyphens (no em-dashes) and no emojis.
- [ ] Docs updated if behavior changed; `progress.yaml` updated (never `PROGRESS.html` by hand).
