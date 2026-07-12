# SpatialScribe capability catalog

Auto-generated from `capabilities.REGISTRY` by `scripts/render_architecture.py` (do not hand-edit;
run the generator) - the single source of truth both the wizard and the copilot use. Call
`cap.run(adata, name, params, ctx)` or `cap.ensure(adata, name, ctx, params=...)`.
`requires`/`produces` are data-contract keys; a missing prerequisite returns a structured
`prerequisite_missing` error whose hint names the step to run first.

## `load_section` - Load a section from a path
- Load a spatial section from a SERVER-SIDE path and make it the active section. Use this when the user asks to open/load/read a dataset at a path (a Xenium / CosMx / MERSCOPE output folder, or a .h5ad). Auto-detects the platform, INFERS the tissue from the panel metadata (e.g. a Xenium Mouse Brain panel -> 'mouse brain'), and - unless auto_reference is false - automatically SELECTS the best-matched single-cell reference for that tissue and recommends supervised label transfer vs de-novo clustering. Returns what it loaded and which reference it chose. After this, the section is swapped in - subsequent steps run on it.
- requires: -  |  produces: -
- params: allow_fetch, auto_reference, path, tissue  |  copilot-exposed: yes

## `describe_sample` - Describe sample
- Summary of the loaded section: n cells, cell-type composition, QC medians.
- requires: -  |  produces: -
- params: -  |  copilot-exposed: yes

## `panel_check` - Panel check
- Which cell types the panel can/cannot resolve (panel-adequacy check), plus pairs it cannot separate. Marker presence is necessary, not sufficient.
- requires: -  |  produces: uns:panel_check
- params: -  |  copilot-exposed: yes

## `reference_match` - Reference match
- How well does a single-cell reference match this panel? Global panel-gene overlap plus depth-matched per-type resolvability (which cell types the panel can confidently transfer from this reference and which it cannot), and a clustering nudge when the reference is a poor fit. Grounded in eval_metrics.panel_resolvability (supersedes identifiability AUC). Pass reference_path for a .h5ad, else it uses the reference chosen at the Panel-check step. Needs no prior clustering or annotation.
- requires: -  |  produces: uns:reference_match
- params: label_key, reference_path  |  copilot-exposed: yes

## `auto_select_reference` - Auto-select reference (free text)
- Given a FREE-TEXT tissue/tumour context (e.g. 'uveal melanoma', 'mouse brain', 'lung adenocarcinoma'), automatically CHOOSE and LOAD the best-matched pre-computed single-cell reference for this panel - not just rank them. Picks the top registry reference whose keywords + panel-gene overlap fit the tissue and loads it (skips a wrong-tissue atlas even if it is the only one available), or when allow_fetch is set and nothing local fits, fetches a small CELLxGENE reference live via gget. Then scores the reference<->panel match. Use this when the user names a tissue instead of uploading a reference. Degrades to a 'no_reference' note (cluster instead) when nothing suitable is found.
- requires: -  |  produces: -
- params: allow_fetch, tissue  |  copilot-exposed: yes

## `annotation_strategy` - Annotation strategy (supervised vs clustering)
- A REFERENCE tool: answers 'is my single-cell REFERENCE good enough to annotate from, or should I cluster instead?'. Self-verifies the reference<->panel match and reruns a ladder (coarsen the labels the panel cannot separate, then try the best-overlap tissue-matched reference), then RECOMMENDS supervised transfer vs de-novo clustering. NEEDS a reference (from Load or reference_path); with none it returns status 'no_reference' (an empty ladder) and recommends marker-based annotation - 'no_reference' means NO ATLAS IS ATTACHED, it does NOT mean no section is loaded, so never tell the user to load data because of it. To IMPROVE annotation quality by MERGING confusable cell types with no reference (e.g. on the demo), use `merge_types` or `self_heal` (merge_confusable=true) instead. Writes uns['annotation_route'].
- requires: -  |  produces: uns:annotation_route
- params: label_key, max_rounds, reference_path  |  copilot-exposed: yes

## `compute_qc` - Compute QC
- Per-cell QC metrics (counts, genes, % control) and the section summary.
- requires: -  |  produces: obs:total_counts, obs:n_genes_by_counts, obs:pct_counts_control
- params: -  |  copilot-exposed: no

## `qc_funnel` - QC funnel
- Full six-layer QC funnel headline: segmentation, panel-indexed count floor, purity, spatial coherence, and the annotatability breakdown (pct pass/warn/abstain + top abstention reasons).
- requires: -  |  produces: -
- params: -  |  copilot-exposed: yes

## `cluster` - Cluster
- Normalize, embed and Leiden-cluster the section at a given resolution.
- requires: -  |  produces: obs:leiden, obsm:X_umap, uns:rank_genes_groups
- params: resolution  |  copilot-exposed: no

## `cluster_markers` - Top DEGs per cluster
- Top-N differentially-expressed genes for every Leiden cluster, drawn as a dot-plot (colour = mean expression, size = % of cells expressing). Use for 'what are the marker genes of each cluster' or 'top 2 DEGs per cluster'. Reuses the ranking computed during clustering, so it is near-instant; recomputes only when the clusters changed.
- requires: obs:leiden  |  produces: -
- params: n_genes  |  copilot-exposed: yes

## `cluster_confidence` - Cluster confidence
- Data-driven over/under-clustering check: which cluster PAIRS are statistically indistinguishable (merge) and which single clusters hide substructure (split). Grounded in a RandomForest-vs-permutation test (sc-SHC/CHOIR-style p-value + accuracy) and a bimodality coefficient - advisory nudges, never auto-applied. The cluster rung of the uncertainty ladder (cluster then cell then panel).
- requires: obs:leiden  |  produces: -
- params: -  |  copilot-exposed: yes

## `annotate` - Annotate
- Marker + consensus cell-type annotation with per-cell confidence and abstention.
- requires: obs:leiden  |  produces: obs:cell_type, obs:cell_type_final, obs:annotation_confidence, obs:annotation_verdict, obs:annotation_reason
- params: reliability_weighted, truth_key  |  copilot-exposed: no

## `immune_exclusion` - Immune exclusion
- Neighborhood-enrichment z-score between two cell types; tells whether one is spatially excluded from (negative) or infiltrating (positive) the other.
- requires: obs:cell_type  |  produces: -
- params: type_a, type_b  |  copilot-exposed: yes

## `neighborhood_enrichment` - Neighborhood enrichment
- Full cell-type x cell-type neighborhood-enrichment z-score matrix.
- requires: obs:cell_type  |  produces: -
- params: -  |  copilot-exposed: yes

## `niches` - Niches
- Call TME spatial niches (neighborhood composition) and list them.
- requires: obs:cell_type  |  produces: obs:niche
- params: -  |  copilot-exposed: yes

## `co_occurrence` - Co-occurrence
- Cell-type co-occurrence probability vs. spatial distance - which cell types are found together, and at what radius (squidpy co_occurrence).
- requires: obs:cell_type  |  produces: -
- params: -  |  copilot-exposed: yes

## `spatial_genes` - Spatially variable genes
- Top spatially variable genes by Moran's I - genes whose expression is spatially structured across the section rather than randomly distributed.
- requires: -  |  produces: -
- params: n_top  |  copilot-exposed: no

## `state_by_celltype` - Cell states
- Cell-type x cell-state (cycling / IFN / hypoxia / exhaustion ...) mean-score matrix.
- requires: obs:cell_type  |  produces: uns:state_columns
- params: -  |  copilot-exposed: no

## `assign_cell_states` - Type cell states
- Type each cell with its DOMINANT cell-state program (cycling / interferon / hypoxia / stress / EMT / T-exhaustion / T-cytotoxicity / ECM-remodeling / antigen-presentation), CyteType-style - a colourable obs['cell_state'] label, plus the state distribution and each cell type's dominant state. Names the states in plain language via the LLM when a key is set. States are orthogonal to lineage identity (who a cell IS vs what it is DOING).
- requires: obs:cell_type  |  produces: obs:cell_state
- params: min_z  |  copilot-exposed: yes

## `malignant_score` - Malignant score
- Where is the tumor? Marker-based malignant score per cell; returns the fraction of high-malignant cells.
- requires: -  |  produces: obs:malignant_score
- params: -  |  copilot-exposed: yes

## `discover_programs` - De-novo programs
- Discover data-driven gene programs (NMF) beyond the fixed marker/state lists; returns the top genes per program plus a PLAID enrichment score (how well each program's cells express its own signature: mean in vs out + AUROC specificity).
- requires: -  |  produces: obsm:programs, obsm:program_scores, obs:program, varm:program_loadings
- params: n_programs  |  copilot-exposed: yes

## `name_programs` - Label programs
- Name each de-novo NMF program in plain language from its top loading genes (grounded; degrades to 'Program k' with no API key). Relabels obs['program'] so the map legend reads by biological program, and returns the program table with a stable program_id (the program_score_<k> colour field is preserved).
- requires: -  |  produces: obs:program
- params: -  |  copilot-exposed: yes

## `subcluster` - Subcluster
- Subcluster ONE named cell type into finer subtypes and name them - on demand, for the specific type the user asks about. Use this whenever the user says to subcluster / break down / re-cluster / find subpopulations (or subtypes / substructure) WITHIN a named cell type. (To instead auto-subcluster only the heterogeneous or marker-failing types across the WHOLE annotation, use self_heal.)
- requires: obs:cell_type  |  produces: obs:subtype
- params: cell_type, resolution  |  copilot-exposed: yes

## `annotation_methods` - Annotation methods
- Which reference-based annotation methods ran (RCTD / SingleR / scANVI / panhumanpy), their per-cell coverage, and the multi-method consensus-agreement distribution. Reports only computed coverage/agreement - it does not assign labels.
- requires: -  |  produces: -
- params: -  |  copilot-exposed: yes

## `rejection_reasons` - Rejection reasons
- Why weren't some cells confidently typed? Granular, plain-language reasons (too few transcripts, spatial doublet, panel lacks markers, mixed lineages, no clear winner, ...) with the count and % of untyped cells for each.
- requires: obs:annotation_verdict  |  produces: obs:rejection_reason, obs:rejection_detail
- params: -  |  copilot-exposed: yes

## `self_verify` - Self-verify annotation
- Check whether each cell type's cells actually score highest on that type's canonical markers (marker-argmax agreement + one-vs-rest AUC), flag the types that fail, and suggest concrete grounded fixes (subcluster / abstain / merge a confusable pair). Advisory - it never changes labels.
- requires: obs:cell_type  |  produces: uns:annotation_verification
- params: neighborhood  |  copilot-exposed: yes

## `trust_ledger` - Trust ledger
- Three INDEPENDENT per-cell-type verdicts and their disagreements: resolvable? (panel adequacy) x coherent? (cells express the type's markers) x agreed? (the reference methods back the label, when >=3 voted). The informative rows are the CONTRADICTIONS - especially 'coherent but DISPUTED', the coherent-whole-cluster-mislabel the AQI index alone cannot see. Advisory - never changes labels.
- requires: obs:cell_type  |  produces: -
- params: -  |  copilot-exposed: yes

## `merge_types` - Merge confusable types
- IMPROVE / optimize annotation QUALITY by MERGING the cell types the panel cannot separate (no private on-panel marker) into ONE coarser label, and report the quality delta (marker-agreement + number of failing types, before vs after). Reference-FREE: works on ANY already-annotated section (the loaded demo included) - no atlas needed and it does NOT reload the data. This is the tool for 'merge / coarsen the confusable cell types to raise annotation quality'. APPLIES the merges (self_verify only advises). Defaults to the groups panel_check flagged; pass `groups` to merge specific sets. A coin-flip between indistinguishable types becomes one defensible call.
- requires: obs:cell_type  |  produces: -
- params: dry_run, groups  |  copilot-exposed: yes

## `calibrate_confidence` - Calibrate confidence
- Fit a post-hoc ISOTONIC calibration of the per-cell annotation confidence against a ground-truth label column, so the score finally means P(correct). Writes annotation_confidence_calibrated plus a reliability report (ECE and Brier before/after, reliability curves, and an AUC discrimination check); the raw heuristic is never overwritten. Benchmarked on a 12k CosMx section the raw score had ECE 0.33 with accuracy NON-monotonic in confidence. Skips honestly with no labels - calibration is claimed only when it is fit against real labels. Note the AUC: calibration fixes what the number MEANS but cannot change how it RANKS cells, and abstention depends on the ranking.
- requires: obs:annotation_confidence  |  produces: -
- params: cal_frac, n_bins, pred_key, truth_key  |  copilot-exposed: yes

## `highlight_cells` - Highlight cells
- Render a NEW spatial view on the canvas with cells matching a criterion lit up. Call this whenever the user asks to SEE / show / highlight / where-are a population (e.g. 'highlight low-quality cells', 'show the T cells', 'where are the malignant cells'). `criterion` accepts: a cell-type name, 'low quality', 'low confidence' / 'uncertain' / 'abstained', or 'malignant'. Returns the matched count + an honest note (it states any fallback, e.g. when no cells were flagged) and draws the plot - then describe what it shows.
- requires: -  |  produces: -
- params: criterion  |  copilot-exposed: yes

## `show_spatial` - Show spatial
- Draw the spatial map colored by a GENE's expression or a field. Use for 'color/show the tissue by <gene>', 'show <gene> expression', 'color by total counts / malignant / niche / cell type'. `color_by` is a gene symbol or a field name. Only on-panel genes and real fields are accepted (it says so if the key is unknown).
- requires: -  |  produces: -
- params: color_by  |  copilot-exposed: yes

## `marker_dotplot` - Marker dot-plot
- Draw a dot-plot of marker genes across cell types (dot color = expression scaled 0-1 PER GENE i.e. relative across types, size = % of cells expressing). Use for 'dotplot of <cell type> markers' or 'dotplot of GENE1, GENE2 ...'. Off-panel genes are dropped (reported). Once spillover purification has run (a split_corrected layer exists) it automatically draws the raw + SPLIT-purified pair.
- requires: obs:cell_type  |  produces: -
- params: cell_type, corrected, genes  |  copilot-exposed: yes

## `annotate_rctd` - Annotate with RCTD
- Run RCTD doublet-mode deconvolution (rctd-py subprocess, GPU-optional) against the uploaded single-cell reference and write obs['rctd_first_type'] (+ second_type, spot_class, weight, singlet_score). RCTD was the most accurate annotator on the non-circular benchmark, so its per-cell labels feed the consensus vote; set_as_primary adopts them as the primary cell_type. Needs the rctd-py env (SPATIALSCRIBE_RCTD_PYTHON) and a reference; skips honestly otherwise. Slower than the default TACCO annotation - use when you want RCTD's doublet / contamination model.
- requires: -  |  produces: obs:cell_type
- params: max_cells, set_as_primary  |  copilot-exposed: yes

## `split_purify` - SPLIT spillover purify
- Decontaminate transcript spillover between neighbouring cells (SPLIT reference-path purification with residual-contamination removal): deconvolve the section against the reference, then SPLIT::rctd_free_purify, writing the purified counts to layers['split_corrected']. Deconvolution weights come from TACCO by default (fast, in-env); weights_engine='rctd' uses the rctd-py subprocess instead. Reports the median library-size reduction (spillover removed). When the SPLIT-R env or a reference are absent it falls back to an in-app neighbour+marker decontamination (method='marker_neighbour'). Then compare markers with marker_dotplot(corrected=true).
- requires: obs:cell_type  |  produces: -
- params: max_cells, weights_engine  |  copilot-exposed: yes

## `expression_violin` - Expression violin
- Draw a violin of one gene's expression across cell types. Use for 'how is <gene> expressed across cell types' / 'violin of <gene>'. Off-panel gene -> reported, no plot.
- requires: obs:cell_type  |  produces: -
- params: gene  |  copilot-exposed: yes

## `composition_chart` - Composition chart
- Draw a bar chart of the cell-type composition (counts per type). Use for 'show the composition' / 'what are the cell-type proportions'.
- requires: obs:cell_type  |  produces: -
- params: -  |  copilot-exposed: yes

## `show_segmentation` - Show segmented cells
- Draw the actual cell-segmentation POLYGONS (not centroids) colored by a field. Use for 'show the segmented cells' / 'show cell boundaries'. Needs a Xenium run directory with cell_boundaries.parquet - says so honestly if only centroids are available.
- requires: -  |  produces: -
- params: color_by  |  copilot-exposed: yes

## `reference_transfer` - Reference transfer
- Transfer cell-type labels from a user single-cell reference onto the section (reference-anchored annotation). Trains a CellTypist model on the reference (the reliable in-env arm) and, when the optional `tacco` package is installed, also runs TACCO optimal-transport transfer. Writes per-cell reference labels and reports coverage + agreement with the current cell_type. Uses the reference chosen at the Panel-check step (or pass reference_path); degrades to 'no_reference'. Set set_as_primary=true to adopt the CellTypist labels as cell_type.
- requires: -  |  produces: -
- params: label_key, reference_path, set_as_primary  |  copilot-exposed: yes

## `malignant_concordance` - Malignant concordance
- Where is the tumour, cross-checked? Runs the malignant callers - the always-on marker score plus the LEARNED Cancer-Finder caller (isolated env) - and reports each caller's %-malignant. Needs a tumour tissue context; Cancer-Finder degrades to a 'skipped' status when its env is not configured. (InSituCNV is not run interactively - its infercnvpy pass is minutes-long; use it offline/in batch.)
- requires: -  |  produces: obs:malignant_score
- params: cf_threshold, max_cells  |  copilot-exposed: yes

## `self_heal` - Self-verify + re-run
- Self-verify the annotation, then RE-RUN what fails: flags cell types whose cells disagree with their canonical markers and automatically subclusters the heterogeneous ones, re-verifying each round (up to max_rounds). Set merge_confusable to ALSO apply the panel-gap coarsening first (collapse the types the panel cannot separate into one label) and report the quality delta; relabel stays advisory. Returns the per-round action log, the merge delta, and the before/after failing types.
- requires: obs:cell_type  |  produces: -
- params: max_rounds, merge_confusable  |  copilot-exposed: yes
