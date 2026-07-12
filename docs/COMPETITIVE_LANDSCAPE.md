# SpatialScribe - Competitive Landscape

**Compiled:** 2026-07-08 - from two structured, adversarially-verified research passes (a named-competitor
deep dive + a 6-category market sweep of ~120 tools). Preprints/self-benchmarks weighted accordingly.
This doc positions SpatialScribe against the field and states, plainly, where its moat is and is not.

> **One-line thesis.** Across ~120 surveyed tools, **none combines self-serve + spatial + *decomposed,
> plain-language* annotation uncertainty + data-driven merge/split + panel-adequacy + an LLM copilot that
> narrates all of it grounded in real numbers.** Competitors either force a label or expose a single opaque
> score. That specific *stack* is the differentiator - not any one piece, and not "it has a copilot".

---

## 1. The named competitors

| Tool | What it is | The gap that matters |
|------|-----------|----------------------|
| **BioTuring SpatialX** (BBrowserX / Talk2Data) | No-code multi-platform viewer; huge public reference DB; protein-guided annotation | Annotation is **cluster-level majority vote**; no per-cell confidence/abstention, no panel-adequacy, no cell-state library. Talk2Data's copilot queries *their* public DB, not your uploaded run. Sales-gated. |
| **Aspect Analytics Weave** | Cross-modality (MSI/IMC/spatial) fusion, enterprise governance | "Flexible cell annotation" with no method detail. No uncertainty, no panel check, no cell-state, **no shipped copilot** (only "AI-ready data"). Enterprise/consortium. |
| **10x Xenium Explorer** | Free viewer everyone already has; excellent raw-signal visual QC | **Zero automated annotation** (only unsupervised numbered clusters); no confidence, no abstention, no wizard, no copilot. A viewer, not a pipeline. |
| **Jetomics** | Windows desktop no-code GUI; RCTD-or-manual typing; offline pitch | No uncertainty, no panel check, no cell-state, no copilot. Xenium-only, install + account required. |
| **Hevelion** | Marco Varrone (CellCharter author, UNIL/SIB); web platform, "automated annotation" + niches | Behind login/waitlist. No public evidence of a copilot, per-cell uncertainty, cell-state, or panel-adequacy. Closest Swiss-academic competitor; niche/domain-focused. |

---

## 2. The finding that reframes the pitch: a 2026 wave of agentic spatial copilots

**"LLM copilot for spatial" is no longer a differentiator by itself.** A whole 2026 cohort ships the
"LLM runs grounded real spatial analysis" pattern:

- **Latch Spatial Agents** (LatchBio) - per-platform NL agents (Xenium/MERSCOPE/Trekker), commercial,
  assay-manufacturer partnerships. **The real product threat** and the closest architecture to SpatialScribe.
- **ChatSpatial** (bioRxiv Feb 2026) - MCP server for Claude Desktop, 60+ spatial methods, schema-enforced
  (non-hallucinating) tool calls. Nearly identical core mechanic; same MCP family as **stMCP**.
- **spatiAlytica** (bioRxiv May 2026) - napari-embedded agent with a **Spatial-VQA sub-agent that reads the
  live canvas** - directly analogous to SpatialScribe's canvas+copilot loop.
- **CellAgent** (spatial-domain ID + imputation via NL), **STAT** ("shared interactive tissue viewer" +
  staged pipeline), **STAnalyzer**, **PhenoGraph**, **Omega / napari-chatgpt + Napari MCP Server** (LLM
  drives napari + real Cellpose/StarDist) - more of the same wave.

**Implication:** do not lead with the copilot. Lead with the *honest-uncertainty stack* the copilot narrates.

---

## 3. Honest-uncertainty precedents (benchmark against these; none is spatial + self-serve + shipped)

| Tool | Uncertainty mechanism | Why it matters to us |
|------|----------------------|----------------------|
| **popV** (Nat Genet 2024) | **Consensus score** = # of an ~8-method ensemble that agree (ontology-aware); tracks accuracy tightly (>=6 -> >90% correct). Explicitly *refuses* to trust single-method self-scores. | The peer for "trust cross-method agreement". **Being ported** (see `docs/POPV_UNCERTAINTY.md`) - it also corrects a real flaw: SpatialScribe weighted the wrong signal. No spatial, no self-serve UI. |
| **scConform** (Harvard CCB) | Conformal prediction **sets** with distribution-free coverage guarantees, ontology-consistent. | The **statistical-rigor benchmark** our heuristic confidence will be compared to. Concrete upgrade path (conformal set-size -> plain sentence). No product surface. |
| **mLLMCelltype** | Multi-LLM debate -> Consensus Proportion + Shannon Entropy as uncertainty. | Same "LLM consensus as uncertainty" instinct; scRNA-only. |
| **GPTAnno** (bioRxiv Nov 2025) | Ontology-tree "uncertainty-aware" annotator; auto-selects resolution; flags ambiguous clusters. | Closest scRNA match to SpatialScribe's per-cluster honest-uncertainty philosophy. |
| **InSituType** (CosMx engine) | Per-cell posterior; docs recommend an 80% "unclassified" threshold. | Notably documents its *own overconfidence* - but a bare number, not reason-codes. |
| **STELLAR** (Stanford) | Spatial GNN with a built-in novel-cell **reject** path. | The closest *spatial-native* abstention precedent; 2022 research code, no product. |

**None ships per-cell, plain-language abstention *reason codes* in a self-serve spatial product.** That
specific instantiation (`rejection.py`, 11 codes) is still unoccupied.

---

## 4. Merge/split prior art (validates the `cluster_confidence` feature)

- **sc-SHC** (Grabski & Purdom, Genome Biol 2023) - model-based hypothesis test for whether a cluster split
  is a genuinely distinct population; usable post-hoc on any clustering.
- **CHOIR** (Corces Lab, Nat Genet 2025) - iterative **random-forest permutation testing** across a cluster
  hierarchy for cell types *and* states.

Both are R-only, bioinformatician-facing. SpatialScribe's `cluster_confidence` reimplements their core in
scikit-learn and wraps it in a no-code copilot + wizard - **no self-serve tool surfaces merge/split nudges.**

---

## 5. Self-serve spatial wizards (structurally close; no uncertainty / no copilot)

- **Recognize** (Resolve Biosciences) - browser segment -> matrix -> cluster -> type -> explore; the closest
  end-to-end self-serve wizard, but no uncertainty/copilot and locked to Molecular Cartography.
- **AtoMx SIP** (NanoString/Bruker) - genuinely no-code cloud wizard with bundled **InSituType** confidence,
  but a bare number, platform-locked to CosMx/GeoMx, no LLM narration.
- **i-stLearn** - no-code web wizard (load -> QC -> cluster -> trajectory -> CCI).
- **Partek Flow**, **Trailmaker** (Parse Biosciences), **Pluto Bio** - mature/funded self-serve multi-omics
  incumbents adding spatial breadth; compete on maturity, not uncertainty or narration.

---

## 6. Pathology-AI with real confidence UX (adjacent; protein/H&E, not transcriptomics, no narration)

**Visiopharm Phenoplex** (tunable per-object confidence 50-100%), **HALO / HALO AI** (per-object probability
markup), **Aiforia** (no-code "train-by-drawing"), **Nucleai**, **Aignostics**, **QuPath + QuST**. These have
confidence + (some) self-serve, but for protein/H&E images, no transcriptomics, and no conversational layer.
The lesson to borrow: a *tunable confidence threshold* the user controls is an expected UX.

---

## 7. Close to home (Swiss / UNIL cluster)

- **CellCharter** (Ciriello Lab, UNIL) - spatial niche/domain identification; the engine behind **Hevelion**.
- **ProjecTILs / SPICA** and **UCell** (Carmona Lab, UNIL) - continuum cell-**state** modeling + program
  scoring (UCell already used at a genomics facility).
- **Scailyte** (Swiss) - interpretability via integrated gradients (same "which genes drove the call" instinct).

Relevant because SpatialScribe is built at a university/the cluster; these are the nearest local precedents and potential
collaborators-or-competitors.

---

## 8. Strategic takeaways

1. **Sell the stack, not the copilot.** The copilot is table stakes as of 2026. The moat is
   cluster -> cell -> panel decomposed honest uncertainty + cell-state, *narrated* by the copilot.
2. **Pre-empt the calibration question.** Be ready for "how is your confidence calibrated vs popV / scConform?"
   The spec already flags calibration as unsolved; porting popV's consensus score (empirically accuracy-correlated)
   and, later, conformal set-size, is the credible answer.
3. **The strategic threats are commercial, not academic.** Latch (partnerships), Noetik (capital + a 40M-cell
   spatial FM), and the instrument vendors' bundled tools have distribution SpatialScribe cannot match; the
   defensible edge is depth of *honest-uncertainty UX*, not breadth or scale.
4. **Single-axis competitors abound; the combination does not.** popV = uncertainty w/o spatial/UI; CHOIR =
   merge/split w/o UI; Recognize = wizard w/o uncertainty; InSituType = a confidence number w/o reasons. The
   demo should show all rungs of the ladder in one flow, because that is what nobody else has.

---

## Appendix - lower-relevance tools surfaced (names only)

*Instrument/vendor & IMC:* GeoMx DSP Data Analysis Suite, nCounter nSolver, MERSCOPE Vizualizer, Vizgen Gene
Panel Design Portal, MIBIsight/MIBI-O, Polylux, MCD Viewer, steinbock + imcRtools + cytomapper, CODEX MAV,
Curio Seeker/Trekker pipelines, napari-CosMx, MOSAIK, PowerOMX/CellScape Navigator.
*Cloud bioinformatics:* Velsera/Seven Bridges, DNAnexus, Lifebit, Code Ocean, Genestack, Qlucore, Basepair,
Watershed Bio, Genular, Biocartesian.
*Open-source GUI/viewers:* Samui, Cirrocumulus, Vitessce, scimap, napari-spatialdata, napari-clusters-plotter,
cellxgene VIP, ShinyCell2, Giotto Viewer/Suite, SPIAT, UCSC Cell Browser, iSEE, SpatialData, EasyCellType,
Clustifyr, ACT web server.
*Pathology-AI:* Indica HALO, Proscia Concentriq, Mindpeak, PathChat, Paige.AI, PathAI, Owkin MOSAIC, Aiosyn,
Deepcell, Lunaphore COMET/HORIZON, Ultivue, NeoGenomics, OracleBio.
*Annotation methods / agents:* Azimuth (retreating from its web UI), ScType (now spatial), CellTypist (web
portal + conf_score), Symphony, SingleR, Tangram/cell2location, STELLAR, Besca, UCell, CellCharter, scChat,
Biomni, CytoAnalyst, Noetik OCTO-VirtualCell, "Analyst-Consensus-Reviewer" ST framework, GPTCelltype,
Nicheformer, scGPT-spatial, scDblFinder, InstructCell, CellWhisperer.

## Sources (primary)

BioTuring (bioturing.com/spatialx, /talk2data), Aspect Analytics (aspect-analytics.com/weave), 10x Xenium
Explorer, jetomics.com, Hevelion; popV [Nat Genet 2024, 10.1038/s41588-024-01993-3]; scConform
(ccb-hms.github.io/scConform); sc-SHC (github.com/igrabski/sc-SHC); CHOIR (choirclustering.com, Nat Genet 2025);
LatchBio "Agents for Spatial Biology" (blog.latch.bio); ChatSpatial / spatiAlytica / CellAgent / STAT
(bioRxiv 2026); NanoString AtoMx + InSituType (github.com/Nanostring-Biostats/InSituType); CyteType (Nygen,
bioRxiv 2025); Visiopharm Phenoplex; Indica HALO; Aiforia; CellCharter (Ciriello Lab); ProjecTILs/UCell
(Carmona Lab). Full per-tool profiles + threat/net-new flags in the two research-pass artifacts.
