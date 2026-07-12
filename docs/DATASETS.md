# Datasets

## Bundled demos

| Demo | What | How to load |
|------|------|-------------|
| **Synthetic melanoma** | An instant, reproducible synthetic section built in code (`analysis/demo.py`). No data file, no GPU, no key. | "Load synthetic demo" in the app, or `python scripts/run.py --demo` |
| **breast** | A public 10x **Xenium Prime 5K** FFPE human breast section (CC BY 4.0), shipped downsampled at `data/demo_public.h5ad`. ~5,000-gene panel at shallow depth - the honest hard case for panel adequacy. | "Load breast example" in the app, or `POST /api/load_demo?name=breast` |

The synthetic demo always works and is the right target for a quick look or for reproducing a bug.
The breast demo is a real section, so the panel-adequacy and confidence stories are genuine.

Point the breast demo at a different file with `SPATIALSCRIBE_DEMO_CACHE=/path/to.h5ad`.

## Bring your own section

`io.load(path)` ingests, platform-agnostically:

- a **Xenium** output folder (`cell_feature_matrix.h5` + `cells.parquet`), including Prime 5K and
  the Atera WTA panel,
- a **CosMx** AtoMx export,
- a **MERSCOPE** output,
- a Flex scRNA matrix or VisiumHD `binned_outputs` / `segmented_outputs`,
- or any `.h5ad` with cell x gene counts and `obsm['spatial']`.

A transcripts-only export (just `transcripts.parquet`, no segmentation) is not loadable by design -
a molecule table is not a cell x gene matrix.

In the app, use **Load your own section** and paste a server-side path; from the copilot, ask it to
`load the section at <path>`; headless, pass `--path <dir_or_.h5ad> --tissue "<tissue>"` to
`scripts/run.py`.

## References (optional)

Reference-based annotation (RCTD / SingleR / scANVI / TACCO) and reference transfer are **optional**
and **bring-your-own**: pass a single-cell/-nucleus reference `.h5ad` (or a spacexr `Reference.rds`)
via the app's reference upload, or `--reference ref.h5ad` on the CLI. Without a reference, the app
annotates from curated markers + the LLM inside a closed vocabulary (see
[ANNOTATION_SCHEME.md](ANNOTATION_SCHEME.md)). **Match the tissue first** - a wrong-tissue reference
can pass an overlap check yet annotate nonsense.
