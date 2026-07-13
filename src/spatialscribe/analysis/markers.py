"""Curated marker-gene panels and cell-state signatures.

Curated marker panels for melanoma / skin lineages plus macrophage/TAM cell states.
Used for marker-based annotation (score_genes) and as the
canonical-marker reference for the panel-adequacy check (H3).

All scoring on Xenium must use ``maxRank=150`` (UCell) or the score_genes equivalent -
the default (1500) dilutes scores into zero-count noise because most Xenium cells have
<500 detected genes.

Known panel gap (H3 proof case): on the Xenium 5K panel only MLANA/MITF/DCT/SOX10 are
present; PMEL/TYR/TYRP1 are absent - so a naive melanocyte score is under-powered.
"""

from __future__ import annotations

import hashlib

# 8 canonical lineages for melanoma / skin / uveal-melanoma tissue (+ keratinocyte).
LINEAGE_MARKERS: dict[str, list[str]] = {
    "Malignant/Melanocyte": ["MLANA", "PMEL", "TYR", "TYRP1", "DCT", "MITF", "SOX10", "PAX3"],
    "T cell": ["CD3D", "CD3E", "CD2", "CD8A", "CD4", "IL7R", "TRAC"],
    "NK cell": ["NKG7", "GNLY", "KLRD1", "NCAM1"],
    "B/Plasma": ["MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN"],
    "Myeloid": ["CD68", "LYZ", "C1QA", "C1QB", "AIF1", "ITGAX", "CD14"],
    # Mast cells are abundant in dermis and were MISSING here while the sibling EPITHELIAL_LINEAGES
    # had them. Survivable while the LLM could invent a label; once the vocabulary is CLOSED
    # (`annotate.consensus_annotate`) a missing lineage has nowhere to go but "Novel / unknown" - or
    # worse, gets misfiled into Myeloid. A lineage whose markers are all off-panel is dropped by
    # `present()`, so adding it costs nothing on panels that cannot see it.
    "Mast": ["TPSAB1", "TPSB2", "CPA3", "MS4A2"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5"],
    "Stromal/CAF": ["COL1A1", "COL1A2", "DCN", "PDGFRB", "ACTA2", "RGS5"],
    "Keratinocyte": ["KRT14", "KRT5", "KRT1", "KRT10"],
    "Glial/Neural": ["PLP1", "S100B", "MBP", "GFAP"],
}

# Macrophage / TAM states (internal) - orthogonal states within the Myeloid lineage.
# Kept for `panel_check._load_marker_db`'s curated fallback and `demo.py`, but NO LONGER part of the
# default cell-state set (`states.score_states` now defaults to CELL_STATES only): the three TAM
# programs are Myeloid-only, so they were dead weight on non-myeloid sections. CELL_STATES below
# replaces them with programs that carry across tumour contexts (T-cytotoxicity, ECM, antigen
# presentation).
TAM_STATES: dict[str, list[str]] = {
    "IFN-TAM": ["CXCL9", "CXCL10", "CXCL11", "IDO1", "GBP1", "GBP2", "ISG15", "STAT1", "IRF1", "IFI44L"],
    "LA-TAM": ["SPP1", "APOE", "APOC1", "TREM2", "LIPA", "FABP5", "GPNMB"],
    "Reg-TAM": ["CD163", "MRC1", "MAFB", "MSR1", "CCL18", "STAB1", "SLC40A1"],
}

# Generic cell-state / program signatures (H6) - orthogonal to lineage identity; the single default
# set for `states.score_states`. All symbols in HUMAN UPPERCASE: `_resolve_one` folds case
# unambiguously onto a MOUSE title-case panel (GZMB -> Gzmb), so one set serves human AND mouse - do
# NOT add a parallel mouse dict. Two SignatuR-sourced programs (verbatim gene lists, see comments);
# the rest curated. A program needs >=2 on-panel genes to score (see `score_states`).
CELL_STATES: dict[str, list[str]] = {
    "Cycling": ["MKI67", "TOP2A", "PCNA", "CCNB1", "CDK1", "BIRC5", "UBE2C"],
    "IFN/ISG": ["ISG15", "IFI6", "MX1", "OAS1", "STAT1", "IRF7", "IFIT1", "IFIT3"],
    "Hypoxia": ["VEGFA", "CA9", "HIF1A", "SLC2A1", "NDRG1", "LDHA"],
    # SignatuR Hs/Programs/HeatShock ("Curated HSP") - the full chaperone complement, verbatim.
    # https://raw.githubusercontent.com/carmonalab/SignatuR/master/data-raw/mySignatuR.csv
    # The SignatuR *Mm* "HeatShock" row is a DIFFERENT signature (nucleolar / ribosome-biogenesis
    # genes), not the mouse ortholog set of these chaperones, so a per-gene Mm case-variant check does
    # not apply here; the standard mouse orthologs of these symbols ARE case-variants (HSPA1A ->
    # Hspa1a, DNAJB1 -> Dnajb1, HSP90AA1 -> Hsp90aa1) and resolve via `_resolve_one` on a mouse panel.
    # Human-specific pseudogenes in the list (HSPA7, HSP90AA3P) simply do not resolve on any panel -
    # harmless, and kept only to preserve the curated list verbatim ("do not invent genes").
    "Stress/HSP": [
        "BBS10", "BBS12", "TCP1", "CCT2", "CCT3", "CCT4", "CCT5", "CCT6A", "CCT6B", "CCT7", "CCT8",
        "CLPB", "HSPD1", "HSPE1", "MKKS", "DNAJA1", "DNAJA2", "DNAJA3", "DNAJA4", "DNAJB1", "DNAJB11",
        "DNAJB12", "DNAJB13", "DNAJB14", "DNAJB2", "DNAJB3", "DNAJB4", "DNAJB5", "DNAJB6", "DNAJB7",
        "DNAJB8", "DNAJB9", "DNAJC1", "DNAJC10", "DNAJC11", "DNAJC12", "DNAJC13", "DNAJC14", "DNAJC15",
        "DNAJC16", "DNAJC17", "DNAJC18", "DNAJC19", "DNAJC2", "HSCB", "DNAJC21", "DNAJC22", "SEC63",
        "DNAJC24", "DNAJC25", "GAK", "DNAJC27", "DNAJC28", "SACS", "DNAJC3", "DNAJC30", "DNAJC4",
        "DNAJC5", "DNAJC5B", "DNAJC5G", "DNAJC6", "DNAJC7", "DNAJC8", "DNAJC9", "HSPA12A", "HSPA12B",
        "HSPA13", "HSPA14", "HSPA1A", "HSPA1B", "HSPA1L", "HSPA2", "HSPA4", "HSPA4L", "HSPA5", "HSPA6",
        "HSPA7", "HSPA8", "HSPA9", "HSPH1", "HYOU1", "HSP90AA1", "HSP90AA3P", "HSP90AB1", "HSP90B1",
        "TRAP1", "HSPB1", "ODF1", "HSPB11", "HSPB2", "HSPB3", "CRYAA", "CRYAB", "HSPB6", "HSPB7",
        "HSPB8", "HSPB9",
    ],
    "EMT": ["VIM", "ZEB1", "SNAI2", "FN1", "TWIST1"],
    "T-exhaustion": ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "TOX", "ENTPD1"],
    # SignatuR Hs/Programs/Tcell.cytotoxicity (Carmona et al. 2020 Oncoimmunology): [GZMB, PRF1, FASLG].
    # FASLG dropped: its SignatuR Mm ortholog is `Fasl` (not `Faslg`), i.e. NOT a case-variant, so it
    # would fail `_resolve_one` on a mouse title-case panel. GZMB/PRF1 map cleanly (Gzmb/Prf1) and are
    # the canonical cytotoxic-effector pair, so the human-uppercase set stays valid for both species.
    "T-cytotoxicity": ["GZMB", "PRF1"],
    # Curated (NOT SignatuR - it has no ECM signature): canonical matrix-remodeling genes (MMPs +
    # cross-linkers + matricellular). Mouse orthologs are all case-variants (Mmp9, Lox, Postn, ...).
    "ECM remodeling": ["MMP2", "MMP9", "MMP14", "TIMP1", "LOX", "LOXL2", "FN1", "SPARC", "POSTN",
                       "COL1A1", "COL3A1", "TNC"],
    # Curated (NOT SignatuR): antigen processing AND presentation (GO:0019882) - both MHC arms plus
    # the shared machinery. Class II (HLA-D*, CD74) with its transactivator CIITA; class I (HLA-A/B/C,
    # B2M) with its transactivator NLRC5; the peptide-loading / immunoproteasome path (TAP1, TAP2,
    # TAPBP, PSMB8, PSMB9). Scoped to class II alone only CIITA is on a Xenium Prime 5K panel - 1 gene,
    # under the >=2 floor - so the program silently vanished on BOTH bundled demos. The machinery genes
    # are canonical members of the axis, not padding to force a score. Consequence to state when
    # interpreting it: on a targeted panel the surviving genes are mostly the processing arm, so the
    # score reads as processing capacity rather than MHC-II surface load.
    # Human symbols only. The mouse class-I/II orthologs (H2-K1/D1, H2-Aa/Ab1) are NOT case-variants of
    # the HLA symbols, so `_resolve_one` drops them on a mouse panel; the machinery genes (Ciita, Cd74,
    # Nlrc5, Tap1, Tap2, Tapbp, Psmb8, Psmb9, B2m) ARE case-variants, and carry the program there.
    "Antigen presentation": ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
                             "HLA-DMA", "HLA-DMB", "CD74", "CIITA",
                             "HLA-A", "HLA-B", "HLA-C", "B2M", "NLRC5",
                             "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9"],
}

# Epithelial / carcinoma lineages (breast + general solid-tumor panels, e.g. Xenium Prime).
# Universal immune / stromal / vascular lineages are shared with LINEAGE_MARKERS.
EPITHELIAL_LINEAGES: dict[str, list[str]] = {
    "Epithelial/Tumor": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "ELF3", "CDH1"],
    "Basal/Myoepithelial": ["KRT14", "KRT5", "ACTA2", "OXTR", "TP63"],
    "T cell": ["CD3D", "CD3E", "CD2", "CD8A", "CD4", "IL7R", "TRAC"],
    "NK cell": ["NKG7", "GNLY", "KLRD1", "NCAM1"],
    "B/Plasma": ["MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN"],
    "Myeloid": ["CD68", "LYZ", "C1QA", "C1QB", "AIF1", "ITGAX", "CD14"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3", "MS4A2"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5"],
    "Stromal/CAF": ["COL1A1", "COL1A2", "DCN", "PDGFRB", "ACTA2", "RGS5"],
    "Adipocyte": ["ADIPOQ", "LEP", "FABP4", "PLIN1"],
}


# CNS / brain-parenchyma lineages. Written in HUMAN uppercase like every other set; the
# species-adaptive matcher in `present`/`on_panel` resolves them onto a MOUSE title-case panel
# (SLC17A7 -> Slc17a7), so this one set serves both human and mouse brain sections. Marker choice
# is deliberately over-inclusive per lineage (more synonyms than any single panel carries) so
# coverage degrades gracefully across panels - the panel-adequacy check reports what actually lands.
BRAIN_MARKERS: dict[str, list[str]] = {
    "Excitatory neuron": ["SLC17A7", "SLC17A6", "SATB2", "RORB", "FEZF2", "NEUROD6", "TBR1", "CUX2"],
    "Inhibitory neuron": ["GAD1", "GAD2", "SLC32A1", "PVALB", "SST", "VIP", "LAMP5", "CALB1", "CALB2", "NPY", "CCK"],
    "Astrocyte": ["AQP4", "GJA1", "SLC1A3", "SLC1A2", "GFAP", "ALDH1L1", "ACSBG1", "AGT"],
    "Oligodendrocyte": ["PLP1", "MBP", "MOG", "MOBP", "MAG", "CNP", "CLDN11", "ASPA", "ERMN", "OPALIN", "ST18"],
    "OPC": ["PDGFRA", "CSPG4", "OLIG1", "OLIG2", "SOX10"],
    "Microglia": ["CX3CR1", "P2RY12", "TMEM119", "HEXB", "CTSS", "C1QA", "C1QB", "CSF1R", "TREM2", "TYROBP", "LAPTM5", "SIGLECH", "CD68", "AIF1"],
    "Endothelial": ["CLDN5", "PECAM1", "FLT1", "SLCO1A4", "SLC2A1", "MECOM", "LEF1", "LY6C1"],
    "Mural/Pericyte": ["PDGFRB", "RGS5", "KCNJ8", "ACTA2", "MYH11", "TAGLN", "NOTCH3", "VTN", "CARMN", "ATP13A5"],
    "Ependymal": ["FOXJ1", "TMEM212", "CCDC153", "PIFO", "SPAG16", "HDC", "TTR"],
}


def _match_curated(tissue: str) -> dict[str, list[str]] | None:
    """Return the curated marker set for a KNOWN tissue context, or None if unrecognised.

    Species is handled downstream by the case-adaptive matcher, so the same set serves human and
    mouse (a mouse brain section still routes to BRAIN_MARKERS). None means "no curated set" - the
    caller decides whether to ask Claude for one (:func:`resolve_markers`) or fall back.
    """
    t = (tissue or "").lower()
    if any(k in t for k in ("brain", "cortex", "cerebr", "cerebell", "hippocamp", "striatum",
                            "thalamus", "cns", "neural", "neuron")):
        return BRAIN_MARKERS
    if any(k in t for k in ("melanoma", "skin", "uveal", "melanocyt", "dermal", "cutaneous")):
        return LINEAGE_MARKERS
    if any(k in t for k in ("breast", "mammary", "ductal", "luminal", "brca", "epithelial", "carcinoma")):
        return EPITHELIAL_LINEAGES
    return None


def tissue_has_curated_set(tissue: str) -> bool:
    """True if a CURATED lineage marker set exists for this tissue (brain / skin-melanoma / breast).

    False means :func:`resolve_markers` (offline / no LLM) falls back to the generic
    :data:`EPITHELIAL_LINEAGES` breast-carcinoma set - which mislabels a non-breast organ (e.g. a
    kidney tubule as 'Epithelial/Tumor'). Callers use this to surface a wrong-tissue caveat.
    """
    return _match_curated(tissue) is not None


def for_tissue(tissue: str = "melanoma") -> dict[str, list[str]]:
    """Curated lineage marker set for a tissue context (offline; no LLM).

    Known contexts route to their curated set (brain/skin-melanoma/breast); anything unrecognised
    defaults to the epithelial (general) set. This is the deterministic fast path; for arbitrary
    tissues with a loaded panel and an API key, prefer :func:`resolve_markers`.
    """
    return _match_curated(tissue) or EPITHELIAL_LINEAGES


# Cache generated panels per (tissue, panel-signature) so the LLM is called once, and every
# downstream step (panel-check / annotate / verify) sees an identical, stable marker set.
_TISSUE_PANEL_CACHE: dict[tuple[str, str], dict[str, list[str]]] = {}


def _panel_signature(panel_genes) -> str:
    return hashlib.sha1(",".join(sorted(str(g) for g in panel_genes)).encode()).hexdigest()[:12]


def _generate_llm_panel(tissue: str, panel_genes, cell_types=None) -> dict[str, list[str]]:
    """Claude-generated ``{cell_type: [markers]}`` for ``tissue``, grounded to ``panel_genes``.

    With ``cell_types`` (the section's current categories) it asks for markers for exactly those
    types. Keeps only lineages with >=1 gene on the panel (original gene lists preserved, so
    panel_check can still report present-vs-missing). Returns ``{}`` on any failure (no key / parse
    error / nothing on-panel) so callers fall back to a curated set.
    """
    try:
        from . import llm
        raw = llm.generate_marker_panel(tissue, list(panel_genes), cell_types=cell_types)
    except Exception:
        return {}
    if not isinstance(raw, dict):
        return {}
    grounded: dict[str, list[str]] = {}
    for ct, genes in raw.items():
        if isinstance(genes, list) and present(panel_genes, {ct: [str(g) for g in genes]})[ct]:
            grounded[str(ct)] = [str(g) for g in genes]
    return grounded


def resolve_markers(tissue: str, panel_genes=None, use_llm: bool = False) -> dict[str, list[str]]:
    """Lineage markers for ANY tissue: curated when known, else Claude-generated (panel-grounded).

    Resolution order: (1) a curated set if ``tissue`` matches a known context; (2) for an
    unrecognised tissue with a panel and ``use_llm``, a Claude-proposed panel grounded to the
    section's genes, cached per (tissue, panel); (3) otherwise the epithelial fallback (so the
    offline / no-key path behaves exactly as :func:`for_tissue`). Never raises.
    """
    curated = _match_curated(tissue)
    if curated is not None:
        return curated
    if use_llm and panel_genes is not None and len(panel_genes):
        key = ((tissue or "").strip().lower(), _panel_signature(panel_genes))
        if key not in _TISSUE_PANEL_CACHE:
            gen = _generate_llm_panel(tissue, panel_genes)
            if len(gen) >= 2:
                _TISSUE_PANEL_CACHE[key] = gen
        if key in _TISSUE_PANEL_CACHE:
            return _TISSUE_PANEL_CACHE[key]
    return EPITHELIAL_LINEAGES


_TYPES_PANEL_CACHE: dict[tuple, dict[str, list[str]]] = {}


def infer_organism(panel_genes) -> str:
    """"mouse" if the panel symbols read as mouse Title-case (``Cd3d``), else "human" (``CD3D``).

    A cheap heuristic on symbol case - the same signal the marker resolver already relies on - so
    CellGuide is queried for the right organism without threading an explicit flag through every call.
    Judges only the LETTERS (gene symbols carry digits: ``CD3D`` / ``Cd3d`` / ``MS4A1``), so
    ``str.isalpha`` cannot be used to gate them.
    """
    title = upper = 0
    # `panel_genes or []` raises on a pandas Index / numpy array ("truth value is ambiguous"), and
    # a.var_names (an Index) is the natural argument, so guard on None explicitly.
    for g in ([] if panel_genes is None else panel_genes):
        s = str(g)
        if s.upper().startswith("ENS"):
            continue                                            # Ensembl id, not a symbol - no case signal
        letters = [c for c in s if c.isalpha()]
        if len(letters) < 2:
            continue
        if all(c.isupper() for c in letters):
            upper += 1                                          # CD3D  -> human
        elif letters[0].isupper() and all(c.islower() for c in letters[1:]):
            title += 1                                          # Cd3d  -> mouse
    return "mouse" if title > upper else "human"


def markers_for_types(cell_types, panel_genes=None, tissue: str = "", use_llm: bool = False,
                      organism: str | None = None) -> dict[str, list[str]]:
    """``{cell_type: [markers]}`` for EXACTLY the section's current categories, cached per
    (organism, tissue, panel, types).

    Lets panel-check report marker coverage / identifiability against the labels actually on the
    section (renamed, subclustered, reference-transferred types the curated dictionaries never
    anticipated). Resolution order per type:

      1. **CellGuide** (``cellguide.markers_for_label``): Cell-Ontology-grounded CANONICAL markers
         (literature-curated) or, failing that, COMPUTATIONAL (census effect-size ranked). Deterministic,
         needs no API key, and cannot hallucinate a gene that is not a real marker.
      2. **LLM** for any type CellGuide could not resolve, and only when ``use_llm``.

    Full marker lists are returned (NOT pre-filtered to the panel) so panel_check computes
    present-vs-missing coverage itself. Keys are the requested category names verbatim. Returns ``{}``
    only when nothing grounds at all, so callers fall back to the curated resolver.
    """
    cts = [str(c) for c in (cell_types or []) if str(c).strip()]
    if not (panel_genes is not None and len(panel_genes) and cts):
        return {}
    org = organism or infer_organism(panel_genes)
    # use_llm is in the key: the LLM-fallback portion is use_llm-dependent, so a use_llm=False call
    # must NOT receive LLM-derived markers cached by an earlier use_llm=True call.
    key = (org, bool(use_llm), (tissue or "").strip().lower(),
           _panel_signature(panel_genes), tuple(sorted(cts)))
    if key in _TYPES_PANEL_CACHE:
        return _TYPES_PANEL_CACHE[key]

    out: dict[str, list[str]] = {}
    try:
        from . import cellguide as _cg
        grounded = _cg.markers_for_labels(cts, organism=org, tissue=(tissue or None))
    except Exception:
        grounded = {}
    unresolved = []
    for c in cts:
        genes = (grounded.get(c) or {}).get("markers") or []
        if genes:
            out[c] = genes
        else:
            unresolved.append(c)

    if unresolved and use_llm:                       # LLM only for what CellGuide could not ground
        gen = _generate_llm_panel(tissue, panel_genes, cell_types=unresolved)
        want = {c.lower(): c for c in unresolved}
        for k, v in gen.items():
            if k.lower() in want:
                out[want[k.lower()]] = v

    _TYPES_PANEL_CACHE[key] = out
    return out


# T-cell subtypes for the click-to-subcluster demo (H2).
TCELL_SUBTYPES: dict[str, list[str]] = {
    "CD8 T": ["CD8A", "CD8B", "GZMK", "GZMB"],
    "CD4 T": ["CD4", "IL7R", "CD40LG"],
    "Treg": ["FOXP3", "IL2RA", "CTLA4", "IKZF2"],
    "Exhausted T": ["PDCD1", "HAVCR2", "LAG3", "TOX"],
    "gd/NKT": ["TRDC", "TRGC1", "KLRD1"],
}


def _panel_index(panel_genes):
    """(exact set, lowercase->[symbols]) index over a panel, built once per resolve call."""
    exact = set(panel_genes)
    lower: dict[str, list[str]] = {}
    for g in panel_genes:
        lower.setdefault(str(g).lower(), []).append(str(g))
    return exact, lower


def _resolve_one(gene: str, exact: set, lower: dict):
    """Resolve one marker symbol to the panel's ACTUAL symbol, or None if absent.

    Exact (case-sensitive) match first - so a human panel is unchanged and two symbols that
    differ only in case are never conflated. On a miss, fall back to a case-insensitive match
    and return the panel's real symbol, but ONLY when it is unambiguous (exactly one panel gene
    shares that lowercase form). This is what lets the human-uppercase marker sets resolve onto a
    mouse title-case panel (ACTA2 -> Acta2) without ever mutating the panel or emitting a symbol
    that is not in ``adata.var_names``.
    """
    if gene in exact:
        return gene
    hits = lower.get(str(gene).lower())
    return hits[0] if hits and len(hits) == 1 else None


def on_panel(panel_genes, genes) -> list[str]:
    """Resolve a flat gene list to the panel's actual symbols (species/case-adaptive)."""
    exact, lower = _panel_index(panel_genes)
    return [r for r in (_resolve_one(g, exact, lower) for g in genes) if r is not None]


def present(panel_genes, markers: dict[str, list[str]]) -> dict[str, list[str]]:
    """Restrict each signature to genes on the panel, order preserved, species/case-adaptive.

    Returns the PANEL's actual symbols (so downstream ``adata[:, genes]`` indexing is valid on a
    mouse title-case panel). See :func:`_resolve_one` for the exact-then-unambiguous-fold rule.
    """
    exact, lower = _panel_index(panel_genes)
    return {k: [r for r in (_resolve_one(g, exact, lower) for g in v) if r is not None]
            for k, v in markers.items()}
