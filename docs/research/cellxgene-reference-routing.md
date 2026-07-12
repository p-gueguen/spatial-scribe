# CELLxGENE reference routing: streamlined any-tissue fetch + chemistry-aware selection

*How SpatialScribe pulls a single-cell reference from the CZI CELLxGENE census, why the old path
was weak, and what the metadata-first routing now does. Grounds `reference.fetch_cellxgene_reference`
+ `subprocesses/reference_fetch/run_fetch_cellxgene.py`.*

## The goal

Load, in a streamlined way, a good single-cell reference for **any tissue** - so a Xenium/MERSCOPE
section of a tissue we have no local atlas for can still route to supervised reference transfer
instead of falling back to de-novo clustering. "Good" is not "present": a wrong-tissue or
low-quality reference passes gene-overlap yet annotates nonsense (the "match the tissue first"
rule), and the reference chosen matters more than the transfer method.

## What was wrong with the old fetch

The previous path was a one-liner over `gget.cellxgene(species, tissue)` that **materialised the
whole tissue** into memory and then `np.random.choice`-downsampled it:

1. **Cannot scale.** A broad tissue does not fit. Human `embryo` alone is ~35M cells in the
   2025-01-30 census; densifying it to downsample is impossible on a normal node. The mouse-embryo
   annotation in this session had to bypass `gget` entirely for exactly this reason.
2. **Random downsample drowns rare types.** A uniform sample of a tissue whose composition is 40%
   T cells returns ~40% T cells and a handful of the rare lineages a reference most needs to teach.
3. **No chemistry preference.** A 2015 inDrop dataset and a 2023 10x 3' v3 dataset were equally
   likely to be sampled, even though the newer chemistry captures far more genes per cell and makes
   a better reference for a targeted panel.
4. **No recency preference** for the *paper* / dataset behind the cells.
5. **Includes spatial assays.** Slide-seq / Visium / MERFISH cells are in the census and were
   eligible - but a spatial spot/bin is not a dissociated single-cell reference.
6. **Free-text tissue only, no cache** - every call re-queried, and a query that did not exactly
   match a census `tissue_general` label returned nothing.

## What the census actually contains (measured, census 2025-01-30)

Assay inventory by cells (top entries), which sets what "recent chemistry" can mean today:

| Human (top assays) | cells | Mouse (top assays) | cells |
|---|---:|---|---:|
| **10x 3' v3** | 59.7M | **sci-RNA-seq3** | 34.3M |
| 10x 3' v2 | 22.7M | 10x 3' v2 | 3.7M |
| 10x 5' v1 | 7.4M | 10x 3' v3 | 3.6M |
| sci-RNA-seq3 | 5.0M | Slide-seqV2 *(spatial)* | ... |
| Slide-seqV2 / Visium *(spatial)* | ... | Smart-seq2 / inDrop / microwell-seq | ... |

Two consequences that the routing has to respect:

- **Flex and GEM-X are not in the 2025-01-30 census yet.** The newest 10x chemistry present at
  scale is `10x 3' v3`. So a preference that hard-required Flex would select nothing. The ranking
  puts Flex/GEM-X at the top so they *win once they appear*, but the operative preference today is
  `10x 3' v3 > v2 > v1 > legacy`.
- **Chemistry preference must be per-tissue, not global.** The mouse embryo has essentially only
  `sci-RNA-seq3` (the MOCA atlas); there is no 10x v3 embryo. A preference that excluded non-10x
  would leave the embryo with no reference. The routing therefore *biases within a tissue's
  available assays* and never drops the only chemistry a type has.

## What the routing does now (implemented)

`run_fetch_cellxgene.py` is **metadata-first**. It reads only the census OBS table (cheap even at
tens of millions of cells), decides which cells to pull, and materialises **only those**:

1. **Resolve the tissue** (`reference.resolve_census_tissue`): free-text -> a real census
   `tissue_general` label (exact -> two-way substring, "mouse embryo" -> `embryo` -> token overlap,
   "skin" -> `skin of body`). An honest miss returns nothing rather than a wrong tissue.
2. **Pull metadata only** for `is_primary_data == True and tissue_general == <resolved>`
   (+ optional `disease` / `development_stage` filters, e.g. an embryo Theiler/Carnegie stage).
3. **Route the cells** (`reference.select_reference_cells`, a pure pandas function so it is unit
   tested offline and reused by both the app and the subprocess):
   - **Exclude spatial** assays (`_is_spatial_assay`), but never strip a spatial-only tissue to nothing.
   - **Stratify by `cell_type`**, capping each type at `max(min_cells_per_type, target // n_types)`
     so rare lineages survive; a type with fewer cells keeps all of them.
   - **Prefer recent chemistry** (`_assay_rank`): within each type take the highest-ranked assay
     first (`10x 3' v3 > v2 > v1 > sci-RNA-seq3 > BD/Drop-seq > Smart-seq > microwell/inDrop`),
     ties broken reproducibly.
4. **Materialise only the chosen soma joinids** via `get_anndata`, de-Ensembl to symbols, stamp
   `obs['tissue']` so the tissue-consistency guard confirms, and **cache** the `.h5ad` (a repeat
   call for the same tissue/species/size skips the whole census round-trip).

**Verified live** on human `pancreas` (target 5,000): 66 cell types stratified, **top assay
`10x 3' v3` at 57.6%**, zero spatial assays. The residual `microwell-seq` (17%) and `sci-RNA-seq3`
(16%) are precisely the types that *only* exist in those older-chemistry datasets - kept alive by
per-type stratification instead of dropped, which is the intended behaviour.

`fetch_cellxgene_reference(...)` exposes the knobs: `development_stage`, `prefer_recent_chemistry`,
`exclude_spatial`, `min_cells_per_type`, `use_cache`.

## What is proposed but not yet built

- **Dataset / paper recency as an explicit signal.** The census `census_info/datasets` table
  carries `collection_id` and dataset metadata; joining `dataset_id` -> collection publication date
  would let the routing prefer *recent papers* directly (a second sort key after chemistry), rather
  than using chemistry as a recency proxy. This is the most direct answer to "preferring recent
  papers" and is a clean follow-up - the `dataset_id` is already pulled into the obs metadata.
- **Quality-weighted dataset choice.** Prefer larger, more-cited, or curated collections; down-rank
  datasets with few genes/cell or a suspiciously flat composition.
- **Registry-first, fetch-second, automatically.** `auto_select_reference` already prefers a local
  registry ref and only fetches when the keyword match is weak; the fetched reference could be
  *written back into the registry* (with its tissue + provenance) so the second section of the same
  tissue is instant.
- **Embedding-based tissue routing.** Resolve free-text tissue against the census `tissue_general`
  vocabulary by embedding similarity rather than substring/token overlap, for queries that share no
  literal token with the census label.
- **Disease-aware defaults for tumour contexts** (pull the matched normal + the disease tissue and
  let the malignant caller separate them).

## Files

- `src/spatialscribe/analysis/reference.py` - `select_reference_cells`, `_assay_rank`,
  `_is_spatial_assay`, `resolve_census_tissue`, `fetch_cellxgene_reference` (the wrapper).
- `subprocesses/reference_fetch/run_fetch_cellxgene.py` - the isolated metadata-first fetch.
- `tests/test_reference.py` - offline routing tests (stratification, chemistry preference, spatial
  exclusion, spatial-only fallback, tissue resolution).
