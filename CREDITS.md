# Credits

SpatialScribe stands on a lot of open-source science. It orchestrates these tools rather than
re-implementing them; each keeps its own license and should be cited when you publish.

## Bundled demo data

- **10x Genomics** - Xenium Prime 5K FFPE human breast section (`data/demo_public.h5ad`, downsampled),
  used under **CC BY 4.0**. See [data/README.md](data/README.md).

## Methods & libraries SpatialScribe drives

| Area | Project |
|------|---------|
| Single-cell / spatial core | [scanpy](https://scanpy.readthedocs.io), [anndata](https://anndata.readthedocs.io), [squidpy](https://squidpy.readthedocs.io) |
| GPU acceleration (optional) | [rapids-singlecell](https://rapids-singlecell.readthedocs.io) |
| Annotation-quality metrics | [spatial-anno-metrics](https://github.com/p-gueguen/spatial-anno-metrics) |
| Reference-based annotation | [CellTypist](https://www.celltypist.org), [TACCO](https://github.com/simonwm/tacco), RCTD ([spacexr](https://github.com/dmcable/spacexr)), SingleR, [scvi-tools / scANVI](https://scvi-tools.org) |
| Contamination / purity | [SPLIT](https://github.com/buenrostrolab/SPLIT), [ovrlpy](https://github.com/HiDiHlabs/ovrlpy) |
| Malignant / CNV calling | [infercnvpy](https://github.com/icbi-lab/infercnvpy), Cancer-Finder |
| Marker grounding | [CZI CellGuide](https://cellxgene.cziscience.com/cellguide) |
| App | [React](https://react.dev), [Vite](https://vitejs.dev), [deck.gl](https://deck.gl), [FastAPI](https://fastapi.tiangolo.com) |
| Dictation (optional) | [faster-whisper](https://github.com/SYSTRAN/faster-whisper) |

Built with [Claude Code](https://claude.com/claude-code). Licensed MIT (see [LICENSE](LICENSE)).
