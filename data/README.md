# Demo data

## `demo_public.h5ad`

A downsampled, pre-processed **10x Genomics Xenium Prime 5K** FFPE human breast section
(20,051 cells x 5,101 genes), shipped so the app opens on a real section out of the box.

- **Source:** 10x Genomics public datasets (Xenium Prime 5K human breast), <https://www.10xgenomics.com/datasets>
- **License:** [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) (10x Genomics). Attribution: 10x Genomics.
- **How it was made:** a 100k-cell processed section stratified-subsampled to ~20k cells (all cell
  types kept, >=40 per type), preserving the annotation, spatial coordinates, neighbor graph, and
  per-cluster markers. Reproduce with `python scripts/make_public_demo.py <source.h5ad> demo_public.h5ad`.

The default "breast" demo points here (`SPATIALSCRIBE_DEMO_CACHE` overrides). An instant synthetic
melanoma section (no data file) is also built in - see `docs/DATASETS.md`.
