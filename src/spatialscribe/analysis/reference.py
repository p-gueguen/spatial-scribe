"""Reference <-> panel matching: pick a single-cell reference and score how well it fits.

What it does
------------
Answers, in one module, the reference-side questions a wet-lab scientist has before they
transfer labels onto a spatial section:

- ``load_reference(path)``      - read a custom scRNA/snRNA reference ``.h5ad``, auto-detecting
                                  its cell-type label column and de-Ensembling gene names.
- ``choose_reference(tissue)``  - given a free-text tissue/tumour context (and optionally the
                                  panel genes), rank the pre-computed references in
                                  :data:`REFERENCE_REGISTRY` by keyword match and panel-gene overlap.
- ``reference_panel_match(ref, panel_genes, label_key)`` - THE metric: a single global
                                  reference<->panel match score plus per-cell-type resolvability
                                  (which types this panel can confidently transfer and which it
                                  cannot), built on the depth-matched per-class-F1
                                  :func:`eval_metrics.panel_resolvability` (which supersedes the
                                  one-vs-rest identifiability AUC), with a clustering nudge when
                                  the reference is a poor fit.

How to use it
-------------
>>> ref, key = load_reference("atlas.h5ad")               # (AnnData, "cell_type")
>>> m = reference_panel_match(ref, panel_genes, key)       # {global, per_type, clustering_nudge}
>>> m["global"]["verdict"], m["clustering_nudge"]

Depends on
----------
anndata (inside functions), numpy, difflib, re; and
:func:`spatialscribe.analysis.eval_metrics.panel_resolvability` for the resolvability metric.
Heavy imports stay inside functions so importing the registry stays cheap.

Notes
-----
* ``panel_resolvability`` matching is CASE-SENSITIVE (``str(g) in panel``, no ``upper()``); this
  module passes ``panel_genes`` through unchanged so the overlap count matches the metric's.
* Every registry path is OPTIONAL and existence-checked lazily, so a public checkout with no
  local reference data still ranks references by keyword (the file just shows as unavailable).
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import anndata

# Ensembl gene-id shape (human ENSG..., mouse ENSMUSG..., etc.).
_ENSEMBL = re.compile(r"^ENS[A-Z]*G\d{6,}")
# var columns that commonly hold gene SYMBOLS when var_names are Ensembl ids (cellxgene: feature_name).
_SYMBOL_COLS = ("feature_name", "gene_name", "gene_symbols", "gene_symbol",
                "feature_symbol", "symbol", "Symbol")
# obs columns that commonly hold the cell-type label, best-first.
_LABEL_CANDIDATES = ("cell_type", "author_cell_type", "celltype", "cell_types", "CellType",
                     "cell_type_fine", "annotation", "annotations", "subclass", "class",
                     "labels", "label")


# --------------------------------------------------------------------------- #
# Reference ingestion
# --------------------------------------------------------------------------- #
def _looks_like_label(adata, col: str) -> bool:
    """A column looks like a cell-type label: categorical-ish, not near-unique, not an ontology id."""
    if col.endswith("_ontology_term_id"):
        return False
    n = int(adata.n_obs)
    try:
        nu = int(adata.obs[col].nunique(dropna=True))
    except Exception:
        return False
    # >= 2 types and not near-unique (rejects cell-id / barcode columns). No 0.5*n cap: a small or
    # heavily-subsampled reference can legitimately have more fine-grained types than half its cells.
    return 2 <= nu and (nu / max(1, n)) < 0.9


def detect_label_key(adata) -> str | None:
    """First obs column that looks like a cell-type label (see :func:`_looks_like_label`).

    Priority: the curated :data:`_LABEL_CANDIDATES` that pass the heuristic, then any other
    object/categorical column that does. Returns ``None`` if nothing qualifies.
    """
    for c in _LABEL_CANDIDATES:
        if c in adata.obs.columns and _looks_like_label(adata, c):
            return c
    for c in adata.obs.columns:
        s = adata.obs[c]
        if (s.dtype == object or str(s.dtype) == "category") and _looks_like_label(adata, c):
            return c
    return None


def _symbolize_var_names(adata, gene_name_col: str | None = None) -> str | None:
    """If ``var_names`` look Ensembl, swap them for a gene-SYMBOL column (cellxgene gotcha).

    Fires only when >50% of the first 50 var_names match :data:`_ENSEMBL` and a symbol column
    exists (``gene_name_col`` first, else :data:`_SYMBOL_COLS`). Sets ``var_names`` to that column
    (as str) + ``var_names_make_unique()``. Returns the column used, or ``None`` (no change).
    """
    import numpy as np

    vn = [str(v) for v in list(adata.var_names[:50])]
    if not vn or float(np.mean([bool(_ENSEMBL.match(v)) for v in vn])) <= 0.5:
        return None
    col = gene_name_col if (gene_name_col and gene_name_col in adata.var.columns) else None
    if col is None:
        col = next((c for c in _SYMBOL_COLS if c in adata.var.columns), None)
    if col is None:
        return None
    adata.var_names = [str(x) for x in adata.var[col]]
    adata.var_names_make_unique()
    return col


def _rctd_rds_to_adata(path: str | Path) -> "anndata.AnnData":
    """Convert an R RCTD / spacexr ``Reference`` ``.rds`` to AnnData (cells x genes, obs['cell_type']).

    RCTD references are R objects (not HDF5), so ``anndata.read_h5ad`` fails on them with
    "file signature not found". This reads the object with rds2py and rebuilds the reference from the
    spacexr ``Reference`` slots: ``counts`` (a dgCMatrix, genes x cells), ``cell_types`` (an R factor).
    Accepts either a bare ``Reference`` or a bundle list that holds one under a ``reference`` slot.
    """
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    try:
        from rds2py import read_rds
    except ImportError as e:  # pragma: no cover - optional dep, present in the deploy env
        raise ImportError("loading an RCTD .rds reference needs the 'rds2py' package "
                          "(uv pip install rds2py); or convert the reference to .h5ad first.") from e

    def _vec(node):
        return node.get("data") if isinstance(node, dict) else node

    obj = read_rds(str(path))
    ref = obj if (isinstance(obj, dict) and obj.get("class_name") == "Reference") else (
        obj.get("reference") if isinstance(obj, dict) else None)
    if not (isinstance(ref, dict) and ref.get("class_name") == "Reference"):
        raise ValueError(f"{path} is not a spacexr/RCTD 'Reference' .rds (expected an S4 Reference "
                         "object, or a bundle list with a 'reference' slot).")
    attrs = ref["attributes"]
    c = attrs["counts"]["attributes"]                                # dgCMatrix (CSC) slots
    i = np.asarray(_vec(c["i"]), dtype=np.int64)
    p_ = np.asarray(_vec(c["p"]), dtype=np.int64)
    x = np.asarray(_vec(c["x"]), dtype=np.float32)
    n_genes, n_cells = (int(v) for v in np.asarray(_vec(c["Dim"]), dtype=np.int64))
    counts = sp.csc_matrix((x, i, p_), shape=(n_genes, n_cells))     # genes x cells
    dn = c.get("Dimnames")
    dn_items = dn.get("data") if isinstance(dn, dict) else dn
    genes = _vec(dn_items[0]) if dn_items is not None and len(dn_items) > 0 else None
    cells = _vec(dn_items[1]) if dn_items is not None and len(dn_items) > 1 else None
    ct = attrs["cell_types"]                                          # R factor: 1-based codes + levels
    codes = np.asarray(_vec(ct), dtype=np.int64)
    lv = ct["attributes"]["levels"]
    levels = np.asarray(_vec(lv) if isinstance(lv, dict) else lv, dtype=object)
    labels = levels[codes - 1]

    adata = ad.AnnData(X=counts.T.tocsr())                            # cells x genes
    adata.var_names = ([str(g) for g in genes] if genes is not None
                       else [f"g{j}" for j in range(n_genes)])
    adata.obs_names = ([str(v) for v in cells] if cells is not None
                       else [f"c{j}" for j in range(n_cells)])
    adata.var_names_make_unique()
    adata.obs["cell_type"] = pd.Categorical([str(v) for v in labels])
    return adata


def load_reference(path: str | Path, *, label_key: str | None = None,
                   gene_name_col: str | None = None) -> tuple["anndata.AnnData", str]:
    """Read a scRNA/snRNA reference -> ``(AnnData, label_key)``.

    Accepts an ``.h5ad`` (AnnData) or an ``.rds`` (an R RCTD / spacexr ``Reference``, converted via
    :func:`_rctd_rds_to_adata`). De-Ensembles gene names when needed and auto-detects the cell-type
    column unless ``label_key`` is given. Raises ``FileNotFoundError`` if ``path`` is absent and
    ``ValueError`` (listing obs columns) if no label column can be found - callers wrap and degrade.
    """
    import anndata as ad

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"reference not found: {path}")
    adata = _rctd_rds_to_adata(p) if p.suffix.lower() == ".rds" else ad.read_h5ad(p)
    _symbolize_var_names(adata, gene_name_col)
    key = label_key or detect_label_key(adata)
    if key is None or key not in adata.obs.columns:
        raise ValueError(
            f"no cell-type label column detected in reference {path}; "
            f"obs columns: {list(adata.obs.columns)}. Pass label_key=... .")
    return adata, key


# --------------------------------------------------------------------------- #
# Tissue-consistency guard: does the reference's declared tissue match the section's tissue?
#
# The reference<->panel resolvability metric is reference-INTERNAL cross-validation, so a WRONG-tissue
# reference (a kidney atlas on a breast section) passes it with a high verdict and can outrank the
# correct atlas - it silently annotates nonsense. Most standard references (CellTypist / cellxgene /
# DISCO / Azimuth) DECLARE their tissue in obs; comparing it to the section tissue via an organ-synonym
# map catches a clear mismatch definitively and cheaply, without the false alarms of a raw token match
# (uveal melanoma <-> eye stays consistent). No declared tissue -> undecidable -> an honest caveat.
# --------------------------------------------------------------------------- #
_ORGAN_KEYWORDS: dict[str, tuple[str, ...]] = {
    "breast": ("breast", "mammary", "brca", "ductal", "luminal"),
    "kidney": ("kidney", "renal", "nephron", "glomerul"),
    "lung": ("lung", "pulmonary", "airway", "nsclc", "luad", "lusc", "alveolar", "bronch"),
    "liver": ("liver", "hepatic", "hepato"),
    "brain": ("brain", "cortex", "cerebr", "cerebell", "hippocamp", "striatum", "thalamus", "glioma", "glioblastoma"),
    "skin": ("skin", "dermal", "cutaneous", "epiderm", "melanoma", "melanocyt"),
    "eye": ("eye", "uveal", "retina", "ocular", "choroid"),
    "pancreas": ("pancrea", "islet"),
    "intestine": ("intestin", "colon", "duoden", "ileum", "jejunum", "bowel", "rectum", "colorect"),
    "stomach": ("stomach", "gastric"),
    "heart": ("heart", "cardiac", "myocard"),
    "spleen": ("spleen", "splenic"),
    "lymph_node": ("lymph node", "lymph_node", "tonsil"),
    "bone_marrow": ("bone marrow", "bone_marrow", "marrow"),
    "blood": ("pbmc", "peripheral blood", "leukocyte"),
    "prostate": ("prostate", "prostatic"),
    "bladder": ("bladder", "urothel"),
    "thyroid": ("thyroid",),
}
_REF_TISSUE_COLS = ("tissue", "tissue_general", "Tissue", "organ", "Organ", "tissue_label")


def _organ_of(text: str) -> set[str]:
    """The coarse organ(s) a free-text tissue string names (via :data:`_ORGAN_KEYWORDS`). Organism
    words (human/mouse) carry no organ, so they never match - only the ORGAN drives consistency."""
    t = " " + re.sub(r"[^a-z0-9]+", " ", str(text or "").lower()).strip() + " "
    return {organ for organ, kws in _ORGAN_KEYWORDS.items() if any(k in t for k in kws)}


def detect_reference_tissue(reference) -> str | None:
    """The tissue a reference DECLARES (most common value of an obs tissue column, else None).

    Standard atlases (CellTypist / cellxgene / DISCO) carry ``obs['tissue']``; a bare custom reference
    may not. Returns the label (e.g. ``"kidney"``) or ``None`` when nothing declares it."""
    try:
        obs = reference.obs
    except Exception:
        return None
    for c in _REF_TISSUE_COLS:
        if c in getattr(obs, "columns", []):
            try:
                vc = obs[c].astype(str).value_counts()
                if len(vc):
                    top = str(vc.index[0]).strip()
                    if top and top.lower() not in ("nan", "none", "", "unknown", "na"):
                        return top
            except Exception:
                continue
    return None


def tissue_consistency(query_tissue: str, reference) -> dict:
    """Does the reference's DECLARED tissue match the section's ``query_tissue``?

    ``consistent`` is ``True`` (organs overlap), ``False`` (both name a KNOWN, different organ - a real
    wrong-tissue reference), or ``None`` (undecidable: the reference declares no tissue, or an organ
    cannot be mapped). Only a confident ``False`` should downgrade a verdict; ``None`` is a caveat."""
    ref_tissue = detect_reference_tissue(reference)
    if not ref_tissue:
        return {"reference_tissue": None, "consistent": None,
                "note": "The reference declares no tissue in its metadata, so tissue match cannot be "
                        "confirmed - the verdict measures separability only. Verify it matches your section."}
    q, r = _organ_of(query_tissue), _organ_of(ref_tissue)
    if not q or not r:
        return {"reference_tissue": ref_tissue, "consistent": None, "query_organs": sorted(q),
                "ref_organs": sorted(r),
                "note": f"The reference tissue is '{ref_tissue}'; could not map it and '{query_tissue}' to "
                        "comparable organs, so tissue match is unconfirmed."}
    consistent = bool(q & r)
    note = (f"Reference tissue '{ref_tissue}' matches the section tissue '{query_tissue}'." if consistent
            else f"TISSUE MISMATCH: the reference is '{ref_tissue}' but the section is '{query_tissue}'. "
                 "A wrong-tissue reference passes the internal separability check yet annotates nonsense - "
                 "prefer unsupervised clustering or a tissue-matched reference.")
    return {"reference_tissue": ref_tissue, "consistent": consistent,
            "query_organs": sorted(q), "ref_organs": sorted(r), "note": note}


# --------------------------------------------------------------------------- #
# The pre-computed reference registry (populated from docs/DATASETS.md)
# --------------------------------------------------------------------------- #
# Public-safe: the registry carries only tissue METADATA (keywords / label key / description). The
# concrete on-disk reference PATH for each entry is read from an env var with an empty default, so a
# committed public checkout has NO site-specific paths (an unset entry simply ranks by keyword and
# shows as unavailable). Set SPATIALSCRIBE_REF_<KEY> for a single reference, or point
# $SPATIALSCRIBE_REFERENCE_REGISTRY at a YAML (merged over these) to repoint many at once. The the cluster
# locations for these atlases are inventoried in docs/DATASETS.md, not hardcoded here.
REFERENCE_REGISTRY: dict[str, dict] = {
    "breast": {
        "path": os.environ.get("SPATIALSCRIBE_REF_BREAST", ""),
        "label_key": "cell_type", "gene_name_col": "feature_name",   # Ensembl var_names
        "description": "Global breast atlas (621k cells, cell_type / author_cell_type).",
        "panel_hint": "coarse labels for 353-405-gene panels; fine only on 5K/WTA.",
        "keywords": ["breast", "mammary", "brca", "ductal", "luminal", "carcinoma"]},
    "skin": {
        "path": os.environ.get("SPATIALSCRIBE_REF_SKIN", ""),
        "label_key": "cell_type", "gene_name_col": "feature_name",
        "description": "CELLxGENE skin-of-body reference (gget cellxgene, CC-licensed).",
        "panel_hint": "coarsen to lineage on <500-gene panels.",
        "keywords": ["skin", "melanoma", "dermal", "cutaneous", "epidermis", "uveal", "melanocyte"]},
    "kidney": {
        "path": os.environ.get("SPATIALSCRIBE_REF_KIDNEY", ""),
        "label_key": "cell_type", "gene_name_col": None,
        "description": "DISCO kidney reference.",
        "panel_hint": "collapses on <500-gene panels.",
        "keywords": ["kidney", "renal", "nephron"]},
    "lung": {
        "path": os.environ.get("SPATIALSCRIBE_REF_LUNG", ""),
        "label_key": None, "gene_name_col": None,
        "description": "CellTypist human lung.",
        "panel_hint": "", "keywords": ["lung", "pulmonary", "airway", "nsclc", "lusc", "luad"]},
    "lymph_node": {
        "path": os.environ.get("SPATIALSCRIBE_REF_LYMPH_NODE", ""),
        "label_key": None, "gene_name_col": None,
        "description": "CellTypist lymph node (immune-rich).",
        "panel_hint": "", "keywords": ["lymph", "node", "immune", "tonsil", "lymphoid"]},
    "blood": {
        "path": os.environ.get("SPATIALSCRIBE_REF_BLOOD", ""),
        "label_key": None, "gene_name_col": None,
        "description": "CellTypist blood / PBMC.",
        "panel_hint": "", "keywords": ["blood", "pbmc", "peripheral", "leukocyte"]},
    "glioblastoma": {
        "path": os.environ.get("SPATIALSCRIBE_REF_GLIOBLASTOMA", ""),
        "label_key": None, "gene_name_col": None,
        "description": "ArchMap glioblastoma (disease-specific).",
        "panel_hint": "", "keywords": ["glioblastoma", "gbm", "glioma", "brain tumour", "brain tumor"]},
    "mouse_brain": {
        "path": os.environ.get("SPATIALSCRIBE_REF_MOUSE_BRAIN", ""),
        "label_key": "cell_type", "gene_name_col": None,   # Allen atlas ships mouse symbols in var_names
        "description": "Allen mouse brain atlas (snRNA-seq; mouse gene symbols, title-case).",
        "panel_hint": "mouse title-case symbols; coarsen to class on <500-gene panels.",
        "keywords": ["mouse", "brain", "cortex", "cerebrum", "cerebellum", "hippocampus",
                     "striatum", "thalamus", "cns", "neural", "neuron"]},
    "pan_tissue": {
        "path": os.environ.get("SPATIALSCRIBE_REF_PAN_TISSUE", ""),
        "label_key": "cell_type", "gene_name_col": None,
        "description": "Pan-tissue reference - generic fallback.",
        "panel_hint": "very granular; coarsen for targeted panels.", "keywords": []},
}


def _load_registry() -> dict:
    """The built-in registry, with a ``$SPATIALSCRIBE_REFERENCE_REGISTRY`` YAML merged over it.

    Each entry's on-disk PATH is re-resolved from ``$SPATIALSCRIBE_REF_<KEY>`` at CALL time (not the
    import-time default baked into :data:`REFERENCE_REGISTRY`), so an env var exported after this
    module is imported - e.g. by a deploy wrapper or a test - is honoured.
    """
    reg = {k: dict(v) for k, v in REFERENCE_REGISTRY.items()}
    for k, v in reg.items():
        env_path = os.environ.get(f"SPATIALSCRIBE_REF_{str(k).upper()}", "").strip()
        if env_path:
            v["path"] = env_path
    ypath = os.environ.get("SPATIALSCRIBE_REFERENCE_REGISTRY")
    if ypath and Path(ypath).exists():
        try:
            import yaml
            extra = yaml.safe_load(Path(ypath).read_text()) or {}
            for k, v in extra.items():
                reg[k] = {**reg.get(k, {}), **v}
        except Exception:
            pass
    return reg


def _ref_gene_overlap(path: str, gene_name_col: str | None, panel_set: set[str]) -> float:
    """Fraction of ``panel_set`` present in the reference at ``path`` (cheap backed var read)."""
    import anndata as ad

    r = ad.read_h5ad(path, backed="r")
    try:
        _symbolize_var_names(r, gene_name_col)
        ref_genes = {str(g) for g in r.var_names}
    finally:
        try:
            r.file.close()
        except Exception:
            pass
    return len(ref_genes & panel_set) / max(1, len(panel_set))


def choose_reference(tissue_freetext: str, panel_genes=None, registry: dict | None = None,
                     top_n: int = 3) -> list[dict]:
    """Rank the pre-computed references for a free-text tissue/tumour context (best first).

    Keyword score = best token match (substring -> 1.0, else difflib ratio) of the query against
    each entry's ``keywords`` + key + description. When ``panel_genes`` is given and the file is
    present, blend in the panel-gene overlap: ``score = 0.6*keyword + 0.4*overlap``. Never raises
    on a missing file (a public checkout still ranks by keyword).
    """
    import difflib

    reg = registry if registry is not None else _load_registry()
    q = (tissue_freetext or "").lower().strip()
    q_tokens = [t for t in re.split(r"[^a-z0-9]+", q) if t] or ([q] if q else [])
    panel_set = {str(g) for g in panel_genes} if panel_genes is not None else None   # built once

    out: list[dict] = []
    for key, entry in reg.items():
        path = entry.get("path") or ""
        available = bool(path) and Path(path).exists()
        # Match against the CURATED surface only - the key and keywords - NOT the free-text
        # description. Tokenising the description let a word that is merely descriptive hijack an
        # unrelated query: a "Melanoma + TME" description made every melanoma query rank the breast
        # UM ref, and "Allen mouse cortex" made a human "cortex" query pick the mouse atlas. Keywords
        # are what an entry claims to match; the description is for the human reading the ranked list.
        cand = {str(key).lower()} | {str(k).lower() for k in entry.get("keywords", [])}
        cand = {c for c in cand if c}
        kw = 0.0
        coverage = 0                                     # DISTINCT query tokens with a strong match
        for qt in q_tokens:
            best = 0.0
            for ct in cand:
                # Substring match only for tokens >= 3 chars, so a stray 1-2 char token (e.g. "k"
                # from a "K=29" description, or a numeral) cannot spuriously score 1.0 against a
                # query like "skin" and hijack the ranking. Shorter pairs fall back to fuzzy ratio.
                sub = len(qt) >= 3 and len(ct) >= 3 and (qt in ct or ct in qt)
                best = max(best, 1.0 if sub else difflib.SequenceMatcher(None, qt, ct).ratio())
            kw = max(kw, best)
            if best >= 0.8:
                coverage += 1
        overlap = None
        if panel_set is not None and available:
            try:
                overlap = _ref_gene_overlap(path, entry.get("gene_name_col"), panel_set)
            except Exception:
                overlap = None
        # score IS the keyword match. Gene overlap is a POSITIVE tertiary sort key (below), never
        # blended into score: the old 0.6*kw + 0.4*overlap could only SUBTRACT from a kw=1.0 tissue
        # match, so a gene-verified .h5ad ref (overlap<1) lost to an unreadable .rds ref (overlap=None
        # -> kept full kw), demoting the species+organ-matched atlas. Tissue match first, then overlap.
        score = kw
        out.append({
            "tissue_key": key, "path": path or None, "label_key": entry.get("label_key"),
            "gene_name_col": entry.get("gene_name_col"), "description": entry.get("description", ""),
            "panel_hint": entry.get("panel_hint", ""), "available": available,
            "keyword_score": round(kw, 3), "match_coverage": coverage,
            "gene_overlap_frac": (round(overlap, 3) if overlap is not None else None),
            "score": round(score, 3)})
    # Rank by score (= keyword match), then break ties on COVERAGE: a reference whose keywords match
    # MORE of the query tokens (organ AND disease) beats a generic one that matched a single token to
    # the same 1.0. Without this, "uveal melanoma liver metastasis" cannot pick the liver reference over
    # the eye reference (both hit 1.0 on the melanoma token). Then break remaining ties on gene OVERLAP
    # (a positive key: a measured, higher-overlap ref wins; an unmeasured overlap sorts last, never
    # ahead). Pure secondary/tertiary keys: they never reorder a non-tied primary ranking.
    out.sort(key=lambda r: (r["score"], r["match_coverage"],
                            r["gene_overlap_frac"] if r["gene_overlap_frac"] is not None else -1.0),
             reverse=True)
    return out[:top_n]


# --------------------------------------------------------------------------- #
# CELLxGENE routing: which CELLS to pull for a reference (chemistry + stratify)
# --------------------------------------------------------------------------- #
# Assay -> preference rank (higher = preferred), for picking the BEST chemistry a tissue
# offers rather than a random slice. Matched by lowercased substring so census EFO label
# variants ("10x 3' v3", "10x 3' v3 transcription profiling", ...) all resolve. Newest 10x
# chemistries (Flex / GEM-X) rank highest so they WIN once they enter the census - they are
# NOT in the 2025-01-30 census yet (measured: the newest widely-present 10x label is `10x 3' v3`),
# so today this preference mostly separates 10x v3 > v2 > v1 and demotes plate/legacy methods.
# The rank only biases WITHIN a tissue's available assays - if a tissue only has sci-RNA-seq3
# (e.g. the mouse-embryo MOCA atlas), that assay still wins by being the only option.
_ASSAY_RANK: dict = {}  # populated by _assay_rank via ordered rules; kept as a doc anchor.
_SPATIAL_ASSAY_TOKENS = ("slide-seq", "visium", "merfish", "xenium", "cosmx", "curio", "spatial")


def _assay_rank(assay) -> int:
    """Chemistry-recency preference for a census ``assay`` label (higher = preferred)."""
    a = str(assay).lower()
    if any(t in a for t in ("flex", "fixed rna", "gem-x", "gemx")):
        return 9                                      # newest 10x - not in census 2025-01-30 yet
    if "10x 3' v3" in a or "10x 3' v4" in a or "10x 5' v3" in a:
        return 8
    if "10x 5' v2" in a or "10x multiome" in a:
        return 7
    if "10x 3' v2" in a or "10x 5' v1" in a:
        return 6
    if "10x 3'" in a or "10x 5'" in a or "10x transcription" in a:
        return 5                                      # 10x v1 / generic "transcription profiling"
    if "scalebio" in a or "sci-rna-seq3" in a:
        return 4                                      # recent high-throughput non-10x
    if "bd rhapsody" in a or "seq-well" in a or "drop-seq" in a:
        return 3
    if "smart-seq" in a:
        return 2
    if any(t in a for t in ("indrop", "cel-seq", "mars-seq", "microwell", "quartz-seq", "sci-rna-seq")):
        return 1                                      # older / plate / early combinatorial
    return 3                                          # unknown droplet-ish: middle, still usable


def _is_spatial_assay(assay) -> bool:
    """A spatial assay is not a dissociated single-cell reference - exclude it by default."""
    return any(t in str(assay).lower() for t in _SPATIAL_ASSAY_TOKENS)


def select_reference_cells(obs, *, target_cells: int = 40000, min_cells_per_type: int = 25,
                           prefer_recent_chemistry: bool = True, exclude_spatial: bool = True,
                           cell_type_col: str = "cell_type", assay_col: str = "assay",
                           seed: int = 0):
    """Pick ~``target_cells`` rows from a census ``obs`` frame, stratified by cell type and biased
    toward recent chemistry - the routing that makes a fetched reference GOOD, not just present.

    Metadata-first: operate on the (cheap) obs table so a multi-million-cell tissue never has to be
    materialised. Three rules, in order:

    1. **Exclude spatial** (``exclude_spatial``): Slide-seq / Visium / MERFISH etc. are not a
       dissociated reference. Never drops the last cell - a spatial-only tissue keeps its cells.
    2. **Stratify by cell type**: cap each ``cell_type`` at ``max(min_cells_per_type,
       target_cells // n_types)`` so rare types are represented instead of drowned by a random
       downsample (the old failure). A type with fewer than the cap keeps all its cells.
    3. **Prefer recent chemistry** (``prefer_recent_chemistry``): within each type, take the
       highest-:func:`_assay_rank` cells first (10x v3 > v2 > ... > legacy), ties broken
       reproducibly by ``seed``. If a type only has an older assay, it is used unchanged.

    Pure pandas (no census/anndata import) so it is unit-testable offline and importable in the
    isolated fetch subprocess. Returns the selected sub-frame (a copy), index-sorted.
    """
    import numpy as np
    import pandas as pd

    df = obs if isinstance(obs, pd.DataFrame) else pd.DataFrame(obs)
    df = df.copy()
    if df.empty:
        return df
    if exclude_spatial and assay_col in df.columns:
        keep = ~df[assay_col].map(_is_spatial_assay)
        if keep.any():                                 # never strip the tissue down to nothing
            df = df[keep]
    types = df[cell_type_col].astype(str)
    n_types = max(1, types.nunique())
    cap = max(int(min_cells_per_type), int(target_cells) // n_types)
    rng = np.random.default_rng(seed)
    df = df.assign(
        _ct=types.to_numpy(),
        _arank=(df[assay_col].map(_assay_rank) if (prefer_recent_chemistry and assay_col in df.columns) else 0),
        _r=rng.random(len(df)),
    ).sort_values(["_arank", "_r"], ascending=[False, True])
    picked = df.groupby("_ct", sort=False, observed=True).head(cap)
    if len(picked) < int(target_cells):               # top up toward target from the best leftover chemistry
        need = int(target_cells) - len(picked)
        picked = pd.concat([picked, df.drop(index=picked.index).head(need)])
    return picked.drop(columns=["_ct", "_arank", "_r"], errors="ignore").sort_index()


def resolve_census_tissue(query: str, available) -> str | None:
    """Map a free-text ``query`` to a census ``tissue_general`` label from ``available``.

    Exact (case-insensitive) match wins; else a two-way substring match (query contains a label
    or vice versa), preferring the longest such label ("mouse embryo" -> "embryo"); else the label
    with the most shared word tokens; else ``None`` (caller degrades). Pure - no census import.
    """
    q = (query or "").strip().lower()
    if not q:
        return None
    avail = [str(a) for a in available if str(a).strip()]
    low = {a.lower(): a for a in avail}
    if q in low:
        return low[q]
    subset = [a for a in avail if a.lower() in q or q in a.lower()]
    if subset:
        return sorted(subset, key=len, reverse=True)[0]
    qtok = set(re.findall(r"[a-z0-9]+", q))
    scored = [(len(qtok & set(re.findall(r"[a-z0-9]+", a.lower()))), a) for a in avail]
    scored = [s for s in scored if s[0] > 0]
    if scored:
        return sorted(scored, key=lambda s: (s[0], len(s[1])), reverse=True)[0][1]
    return None


# --------------------------------------------------------------------------- #
# Reference-derived marker sets - the fix for the AQI purity/fidelity COVERAGE artifact
# on a fine reference-transfer vocabulary (MOCA / Allen atlas labels).
# --------------------------------------------------------------------------- #
def _lognorm_for_ranking(adata) -> None:
    """If ``adata.X`` looks like raw counts, ``normalize_total`` + ``log1p`` IN PLACE so a per-label
    t-test ranks on stabilised values. A no-op when the data is already normalised."""
    import numpy as np
    import scanpy as sc

    sample = adata.X[:200]
    sample = sample.toarray() if hasattr(sample, "toarray") else np.asarray(sample)
    if (sample.size and float(sample.min()) >= 0
            and float(sample.max()) > 20 and np.allclose(sample, np.round(sample))):
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)


def reference_marker_sets(reference, label_key: str = "cell_type", panel_genes=None, *,
                          top_n: int = 8, min_cells: int = 20, min_markers: int = 2,
                          cache: bool = True) -> dict:
    """Per-label marker sets from a REFERENCE's own DEGs, restricted to the panel.

    The fix for the AQI purity/fidelity COVERAGE artifact on a fine reference-transfer vocabulary
    (MOCA / Allen atlas labels the curated + CellGuide dicts never anticipated): ground each label's
    markers in the reference's OWN biology - independent of THIS section's label assignment (unlike
    the section's own DEGs, which are circular and inflate purity). Restricts the reference to the
    panel genes FIRST, then ranks (``rank_genes_groups``, t-test), so the top_n are guaranteed
    on-panel discriminators. A label with < ``min_cells`` cells or < ``min_markers`` on-panel markers
    is dropped - it is genuinely unresolvable on this panel, which belongs in the ADEQUACY term (A),
    not a fake purity penalty. Cached on ``reference.uns``. NEVER raises - returns ``{}`` when it
    cannot compute (no reference / no label column / no panel overlap), so the caller keeps its set.
    """
    if reference is None or label_key not in getattr(getattr(reference, "obs", None), "columns", []):
        return {}
    panel_set = {str(g) for g in (panel_genes or [])}
    shared = [str(g) for g in reference.var_names if (str(g) in panel_set or not panel_set)]
    if len(shared) < 2:
        return {}
    ckey = f"_ref_markers::{label_key}::{len(shared)}::{hash(frozenset(shared)) & 0xffffffff}::{top_n}"
    if cache and ckey in getattr(reference, "uns", {}):
        return reference.uns[ckey]
    try:
        import scanpy as sc

        sub = reference[:, shared].copy()
        vc = sub.obs[label_key].astype(str).value_counts()
        keep = vc[vc >= int(min_cells)].index.tolist()
        sub = sub[sub.obs[label_key].astype(str).isin(keep)].copy()
        if sub.obs[label_key].astype(str).nunique() < 2:
            return {}
        _lognorm_for_ranking(sub)
        sub.obs["_k"] = sub.obs[label_key].astype(str).astype("category")
        sc.tl.rank_genes_groups(sub, "_k", method="t-test", n_genes=int(top_n))
        names = sub.uns["rank_genes_groups"]["names"]
        out: dict = {}
        for t in names.dtype.names:
            genes = [str(g) for g in list(names[t])[:int(top_n)] if (str(g) in panel_set or not panel_set)]
            if len(genes) >= int(min_markers):
                out[str(t)] = genes
    except Exception:                                        # rank failure / degenerate reference
        return {}
    if cache:
        try:
            reference.uns[ckey] = out
        except Exception:
            pass
    return out


# --------------------------------------------------------------------------- #
# Live CELLxGENE fetch (gget) - the "download the best-matched reference" path
# --------------------------------------------------------------------------- #
def _resolve_gget_runner() -> list[str] | None:
    """The command PREFIX that runs the committed census fetch subprocess in an env that HAS census.

    cellxgene-census + tiledbsoma clash with the scanpy/spatialdata main env, so the fetch runs OUT
    of process (same isolation rule as ovrlpy / the annotators). Resolution order:
      1. ``$SPATIALSCRIBE_GGET_PYTHON`` - a python interpreter that already has cellxgene-census
         (empty public default; the name is kept for backward compatibility);
      2. ``uv run --with cellxgene-census,tiledbsoma python`` - a throwaway env, if ``uv`` is on PATH.
    Returns the argv prefix (script path appended by the caller), or ``None`` if neither is available.
    """
    import shutil

    py = os.environ.get("SPATIALSCRIBE_GGET_PYTHON", "").strip()
    if py and Path(py).exists():
        return [py]
    if shutil.which("uv"):
        # --no-project: resolve ONLY the --with packages, NOT the surrounding spatialscribe pyproject.
        # Without it, `uv run` inside the checkout tries to resolve spatialscribe[dev] against the
        # ambient (e.g. 3.13) interpreter and fails ("requires spatialscribe[dev] ... unsatisfiable").
        return ["uv", "run", "--no-project", "--with", "cellxgene-census",
                "--with", "tiledbsoma", "python"]
    return None


def fetch_cellxgene_reference(tissue_query: str, *, disease: str | None = None,
                              species: str = "homo_sapiens", max_cells: int = 20000,
                              census_version: str = "2025-01-30", out_path: str | None = None,
                              development_stage: str | None = None,
                              prefer_recent_chemistry: bool = True, exclude_spatial: bool = True,
                              min_cells_per_type: int = 25, use_cache: bool = True,
                              timeout: int = 1800) -> dict:
    """Fetch a small CELLxGENE reference for a free-text ``tissue_query`` (live, metadata-first).

    Runs the committed ``subprocesses/reference_fetch/run_fetch_cellxgene.py`` in an isolated env
    (see :func:`_resolve_gget_runner`) so the heavy cellxgene-census/tiledbsoma stack never enters
    the main interpreter, then returns the written ``.h5ad`` path + a summary. The subprocess reads
    the census **metadata first** and picks cells via :func:`select_reference_cells` (cell-type
    stratified + recent-chemistry preferred), so a multi-million-cell tissue is never materialised
    and rare types survive. Requires network access to the CZI CELLxGENE census. NEVER raises -
    returns ``{"status": "skipped"|"error", "message": ...}`` when the env/network is unavailable,
    so the app degrades to the local registry or clustering.

    ``development_stage`` narrows a broad tissue (e.g. embryo -> a Theiler/Carnegie stage);
    ``prefer_recent_chemistry`` / ``exclude_spatial`` / ``min_cells_per_type`` steer the routing;
    ``use_cache`` reuses a previously written ``.h5ad`` at the same path instead of re-querying.

    Returns on success ``{"status": "ok", "path", "n_obs", "n_vars", "label_key", "tissue", ...}``.
    """
    import json
    import subprocess
    import tempfile

    cache = out_path or os.path.join(
        os.environ.get("SPATIALSCRIBE_REF_CACHE", tempfile.gettempdir()),
        "ss_cellxgene_" + re.sub(r"[^a-z0-9]+", "_",
                                 f"{species}_{tissue_query or 'ref'}_{max_cells}".lower()).strip("_") + ".h5ad")
    if use_cache and Path(cache).exists():             # streamlined re-load: skip the whole census round-trip
        return {"status": "ok", "path": cache, "label_key": "cell_type",
                "tissue": tissue_query, "cached": True}

    runner = _resolve_gget_runner()
    if runner is None:
        return {"status": "skipped",
                "message": ("no census env: set $SPATIALSCRIBE_GGET_PYTHON to a python with "
                            "cellxgene-census+tiledbsoma, or install uv, to enable the live fetch.")}
    script = Path(__file__).resolve().parents[3] / "subprocesses" / "reference_fetch" / "run_fetch_cellxgene.py"
    if not script.exists():
        return {"status": "error", "message": f"fetch subprocess missing: {script}"}
    argv = [*runner, str(script), "--out", cache, "--tissue", tissue_query or "",
            "--species", species, "--census-version", census_version, "--max-cells", str(int(max_cells)),
            "--min-cells-per-type", str(int(min_cells_per_type))]
    if disease:
        argv += ["--disease", disease]
    if development_stage:
        argv += ["--development-stage", development_stage]
    if not prefer_recent_chemistry:
        argv += ["--no-prefer-recent-chemistry"]
    if not exclude_spatial:
        argv += ["--keep-spatial"]
    try:
        proc = subprocess.run(argv, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"status": "error", "message": f"gget fetch timed out after {timeout}s"}
    except Exception as exc:                                    # env/launch failure -> degrade
        return {"status": "error", "message": f"{type(exc).__name__}: {exc}"}
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-4:]
        return {"status": "error", "message": "gget fetch failed: " + " | ".join(tail)}
    # The subprocess prints one JSON summary line last; parse it, else fall back to the cache path.
    line = next((ln for ln in reversed((proc.stdout or "").splitlines()) if ln.strip().startswith("{")), "")
    try:
        summary = json.loads(line)
    except Exception:
        summary = {"status": "ok", "path": cache}
    summary.setdefault("path", cache)
    summary.setdefault("status", "ok")
    return summary


# --------------------------------------------------------------------------- #
# Free-text -> auto-CHOSEN + LOADED reference (the missing "auto" half)
# --------------------------------------------------------------------------- #
def auto_select_reference(tissue_freetext: str, panel_genes=None, *, registry: dict | None = None,
                          allow_fetch: bool = False, min_keyword_score: float = 0.34,
                          fetch_kwargs: dict | None = None) -> dict:
    """Free-text tissue -> the reference is auto-CHOSEN from the registry AND LOADED (not just ranked).

    This is the piece :func:`choose_reference` was missing: it takes the best-ranked AVAILABLE registry
    entry, actually :func:`load_reference`\\ s it, and hands back the in-memory ``AnnData`` ready to
    register on a session. If nothing local is available (or the best keyword match is below
    ``min_keyword_score``) and ``allow_fetch`` is set, it falls back to a live
    :func:`fetch_cellxgene_reference` and loads that. Pure and NEVER raises: ``ref`` is ``None`` when
    nothing could be loaded, so the caller degrades to unsupervised clustering.

    Returns ``{"status", "ref", "label_key", "chosen", "ranked", "source", "message"}`` where
    ``status`` is ``registry`` (loaded a local ref), ``fetched`` (downloaded one), or
    ``no_reference`` (fall back to clustering).
    """
    ranked = choose_reference(tissue_freetext, panel_genes=panel_genes, registry=registry, top_n=5)
    ranked_public = [{k: r[k] for k in ("tissue_key", "score", "keyword_score",
                                        "gene_overlap_frac", "available", "description")}
                     for r in ranked]
    base = {"ref": None, "label_key": None, "chosen": None, "ranked": ranked_public}

    # 1) best AVAILABLE local registry entry - but only if it is a real tissue match. A low
    # keyword score means the only available atlas is unrelated (e.g. a lung ref for a "kidney"
    # query); auto-loading it would silently annotate nonsense, so we skip it and fall through to
    # a live fetch / clustering instead of trusting a wrong-tissue reference (the "match the tissue
    # first" rule). gene-overlap alone does not rescue a wrong tissue.
    for cand in ranked:
        # Require the candidate to cover at least one query token STRONGLY (match_coverage >= 1), not
        # merely clear the fuzzy keyword floor: the difflib ratio inflates on the shared 'human' token,
        # so a lung ref scores ~0.44 (> 0.34) for a 'human kidney' query yet covers 0 tokens. Coverage
        # is the real wrong-tissue guard ("match the tissue first").
        if (cand.get("available") and cand.get("path") and cand.get("match_coverage", 0) >= 1
                and cand.get("keyword_score", 0.0) >= min_keyword_score):
            try:
                ref, key = load_reference(cand["path"], label_key=cand.get("label_key"),
                                          gene_name_col=cand.get("gene_name_col"))
            except Exception as exc:
                base.setdefault("load_errors", []).append(f"{cand['tissue_key']}: {exc}")
                continue
            return {**base, "status": "registry", "ref": ref, "label_key": key,
                    "chosen": cand["tissue_key"], "source": f"registry:{cand['tissue_key']} ({cand['path']})",
                    "message": f"auto-selected the '{cand['tissue_key']}' reference for '{tissue_freetext}'."}

    # 2) live CELLxGENE fetch (opt-in) when nothing local is available.
    if allow_fetch:
        fetched = fetch_cellxgene_reference(tissue_freetext, **(fetch_kwargs or {}))
        if fetched.get("status") == "ok" and fetched.get("path") and Path(fetched["path"]).exists():
            try:
                ref, key = load_reference(fetched["path"], label_key=fetched.get("label_key"))
                return {**base, "status": "fetched", "ref": ref, "label_key": key,
                        "chosen": "cellxgene", "source": f"gget cellxgene: {fetched['path']}",
                        "message": f"fetched a CELLxGENE reference for '{tissue_freetext}'."}
            except Exception as exc:
                return {**base, "status": "no_reference",
                        "message": f"fetched a reference but could not load it: {exc}"}
        return {**base, "status": "no_reference",
                "message": "live fetch unavailable: " + str(fetched.get("message", fetched.get("status")))}

    top = ranked_public[0]["tissue_key"] if ranked_public else "?"
    return {**base, "status": "no_reference",
            "message": (f"no local reference is available for '{tissue_freetext}' (best keyword match: "
                        f"'{top}'). Set its SPATIALSCRIBE_REF_* path, upload a .h5ad, or enable the "
                        f"live CELLxGENE fetch. Falling back to unsupervised clustering.")}


# --------------------------------------------------------------------------- #
# The metric: reference <-> panel match
# --------------------------------------------------------------------------- #
def reference_panel_match(reference, panel_genes, label_key: str,
                          target_depth: float | None = None, tissue: str | None = None) -> dict:
    """Global reference<->panel match score + per-cell-type resolvability.

    Built on :func:`eval_metrics.panel_resolvability` (depth-matched per-class F1, which
    supersedes one-vs-rest identifiability AUC). Returns a global headline (panel-gene overlap,
    mean per-type F1, fraction resolvable, a good/fair/poor verdict), a per-type resolvable flag,
    and a plain-language ``clustering_nudge`` that fires when the reference is a poor fit (prefer
    unsupervised clustering / coarser labels). All strings use hyphens (never em-dashes).
    """
    from . import eval_metrics as _em

    panel_list = [str(g) for g in panel_genes]     # do NOT uppercase - the metric is case-sensitive
    n_panel = len(panel_list)
    # Finite-positive default substitution: bool(nan) is True, so a plain truthiness check would let a
    # NaN depth (empty/corrupt total_counts -> np.median NaN) through and crash panel_resolvability's
    # rng.binomial thinning. Require a real positive number, else fall back to 50.
    td = float(target_depth) if (target_depth is not None and target_depth == target_depth
                                 and target_depth > 0) else 50.0
    pr = _em.panel_resolvability(reference, label_key, panel_list, target_depth=td)

    n_shared = int(pr.get("n_shared_genes", 0))
    overlap = n_shared / max(1, n_panel)
    status = pr.get("status")

    per_type: dict[str, dict] = {}
    if status == "ok":
        per = pr.get("per_type", {})
        f1s = [float(d["f1"]) for d in per.values()]
        mean_f1: float | None = float(sum(f1s) / len(f1s)) if f1s else 0.0
        n_types = int(pr.get("n_types", len(per)))
        frac_resolvable = float(pr.get("frac_resolvable", 0.0))
        n_resolvable = int(sum(1 for d in per.values() if d.get("tier") == "resolvable"))
        for ct, d in per.items():
            per_type[str(ct)] = {
                "f1": float(d.get("f1", 0.0)), "resolvable": d.get("tier") == "resolvable",
                "tier": d.get("tier"), "confused_with": d.get("confused_with"),
                "confused_frac": float(d.get("confused_frac", 0.0))}
        if overlap >= 0.5 and mean_f1 >= 0.5 and frac_resolvable >= 0.5:
            verdict = "good"
        elif overlap >= 0.3 and mean_f1 >= 0.35:
            verdict = "fair"
        else:
            verdict = "poor"
    else:
        mean_f1, n_types, frac_resolvable, n_resolvable, verdict = None, int(pr.get("n_types", 0)), 0.0, 0, "poor"

    poor_fit = (status != "ok") or overlap < 0.3 or (mean_f1 is not None and mean_f1 < 0.35) \
        or frac_resolvable < 0.3
    nudge: str | None = None
    if poor_fit:
        f1_txt = f"{mean_f1:.2f}" if mean_f1 is not None else "n/a"
        nudge = (f"Reference is a poor fit for this panel ({overlap:.0%} panel-gene overlap, "
                 f"mean per-type F1 {f1_txt}, {n_resolvable}/{n_types} types resolvable). "
                 f"Prefer unsupervised clustering or coarser labels over supervised annotation.")
    elif verdict == "fair":
        nudge = (f"Reference is a fair fit ({overlap:.0%} panel-gene overlap, "
                 f"mean per-type F1 {mean_f1:.2f}). Consider coarsening the labels the panel "
                 f"cannot separate before transferring.")

    # A SINGLE 0-1 reference<->panel matching score - the one headline number a scientist can read
    # ("can I transfer this reference onto this panel?"). A weighted mean of the three components,
    # F1 weighted highest because depth-matched per-type F1 is the strongest evidence (gene presence
    # alone over-promises; see the panel_resolvability gotcha). Falls back to overlap-only when the
    # resolvability metric could not run (status != ok). Aligns with the good/fair/poor verdict.
    if status == "ok" and mean_f1 is not None:
        match_score = round(0.5 * float(mean_f1) + 0.3 * float(frac_resolvable) + 0.2 * float(overlap), 3)
    else:
        match_score = round(0.2 * float(overlap), 3)   # no resolvability -> overlap is the only signal

    # Structured recommendation (the "so what do I do?" decision, from the verdict). This is the
    # deliberate answer to "if the reference-to-dataset match is low, what next?": a GOOD match ->
    # transfer as-is; a FAIR match -> coarsen the labels the panel cannot separate, then transfer; a
    # POOR match -> the panel cannot resolve this reference's types at the section's real depth, so
    # prefer UNSUPERVISED clustering (with de-novo markers), or switch to a better-matched / coarser
    # reference, or coarsen the labels before transferring. The UI can act on `action`, not just prose.
    _REC = {
        "good": ("supervised_transfer",
                 "Supervised label transfer from this reference is reliable on this panel."),
        "fair": ("coarsen_then_transfer",
                 "Transfer is workable, but coarsen the labels the panel cannot separate first - "
                 "fine-grained types will be confused."),
        "poor": ("unsupervised_clustering",
                 "This panel cannot resolve this reference's cell types at your section's depth. "
                 "Prefer unsupervised clustering (with de-novo markers) over supervised transfer, "
                 "or switch to a better-matched or coarser reference, or coarsen the labels first."),
    }
    _action, _head = _REC.get(verdict, _REC["poor"])
    _why = (f"match score {match_score:.2f}: {overlap:.0%} panel-gene overlap, "
            f"mean per-type F1 {mean_f1:.2f}, {n_resolvable}/{n_types} types resolvable."
            if (status == "ok" and mean_f1 is not None)
            else f"match score {match_score:.2f}: resolvability could not be computed "
                 f"({overlap:.0%} panel-gene overlap).")
    recommendation = {"action": _action, "headline": _head, "why": _why, "verdict": verdict}

    # Tissue-consistency guard (the resolvability metric is reference-INTERNAL and tissue-blind): when
    # the section tissue is known and the reference DECLARES a clearly different organ, this is a
    # wrong-tissue reference no matter how well its types separate. Override the recommendation toward
    # unsupervised clustering and flag it - never silently endorse a mismatched atlas.
    tissue_check = tissue_consistency(tissue, reference) if tissue else None
    tissue_mismatch = bool(tissue_check and tissue_check.get("consistent") is False)
    if tissue_mismatch:
        recommendation = {"action": "unsupervised_clustering",
                          "headline": tissue_check["note"], "why": recommendation["why"],
                          "verdict": verdict, "tissue_mismatch": True}
        nudge = tissue_check["note"] + (f" ({nudge})" if nudge else "")

    return {
        "status": status,
        "match_score": match_score,        # single 0-1 headline (also mirrored in `global`)
        "recommendation": recommendation,
        "tissue_check": tissue_check,
        "global": {
            "match_score": match_score,
            "gene_overlap_frac": round(overlap, 3), "n_shared_genes": n_shared,
            "n_ref_genes": int(getattr(reference, "n_vars", 0)), "n_panel_genes": n_panel,
            "mean_f1": (round(mean_f1, 3) if mean_f1 is not None else None),
            "macro_f1": pr.get("macro_f1"), "balanced_accuracy": pr.get("balanced_accuracy"),
            "n_resolvable": int(n_resolvable), "n_types": int(n_types),
            "frac_resolvable": round(float(frac_resolvable), 3), "verdict": verdict,
            "tissue_mismatch": tissue_mismatch, "reference_tissue": (tissue_check or {}).get("reference_tissue"),
            "target_depth": td},
        "per_type": per_type,
        "clustering_nudge": nudge,
    }


# --------------------------------------------------------------------------- #
# Self-verify + auto-rerun ladder: coarsen -> reselect reference -> route
# (supervised label transfer vs de-novo clustering)
# --------------------------------------------------------------------------- #
_VERDICT_RANK = {"poor": 0, "fair": 1, "good": 2}


def _verdict_key(match: dict) -> tuple:
    """(verdict_rank, mean_f1) for comparing two match results - higher is a better fit."""
    g = match.get("global", {}) if isinstance(match, dict) else {}
    f = g.get("mean_f1")
    return (_VERDICT_RANK.get(g.get("verdict"), -1), float(f) if f is not None else -1.0)


def merge_groups_from_confusion(per_type: dict, min_frac: float = 0.25) -> dict:
    """Union-find merge map over cell types that confuse each other on the panel.

    Two types join when one lists the other as ``confused_with`` with ``confused_frac >= min_frac``
    (the resolvability metric's dominant-confuser signal). Returns ``{original: merged}`` for EVERY
    input label (identity for unmerged); a merged label is its members sorted and joined with ' + '.
    Only groups of size >= 2 are merged - this is the panel-driven label coarsening that lets a
    reference the panel cannot resolve at full granularity still transfer at a coarser one.
    """
    labels = list(per_type)
    parent = {t: t for t in labels}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for t, d in per_type.items():
        cw = d.get("confused_with")
        if cw in parent and cw != t and float(d.get("confused_frac", 0.0)) >= min_frac:
            parent[find(t)] = find(cw)

    groups: dict[str, list] = {}
    for t in labels:
        groups.setdefault(find(t), []).append(t)
    mapping: dict[str, str] = {}
    for members in groups.values():
        merged = " + ".join(sorted(members)) if len(members) >= 2 else members[0]
        for m in members:
            mapping[m] = merged
    return mapping


def _safe_match(reference, panel_genes, label_key, target_depth, tissue=None) -> dict:
    """:func:`reference_panel_match` that never raises - a scoring failure (e.g. a ``label_key``
    that is not a reference column, or degenerate data) degrades to a well-formed poor/error match
    so the ladder keeps its 'never raises' contract."""
    try:
        return reference_panel_match(reference, panel_genes, label_key, target_depth=target_depth,
                                     tissue=tissue)
    except Exception as exc:
        return {"status": "error", "error": str(exc)[:200], "per_type": {}, "clustering_nudge": None,
                "global": {"verdict": "poor", "gene_overlap_frac": None, "mean_f1": None,
                           "frac_resolvable": 0.0, "n_resolvable": 0, "n_types": 0}}


def _rescore_coarsened(reference, panel_genes, label_key, mapping, target_depth, tissue=None) -> dict:
    """Re-run the match with labels coarsened by ``mapping`` (transient col; never raises)."""
    import pandas as pd

    tmp = "_ss_coarse_label"
    orig = reference.obs[label_key].astype(str)
    reference.obs[tmp] = pd.Categorical(orig.map(lambda x: mapping.get(x, x)))
    try:
        return _safe_match(reference, panel_genes, tmp, target_depth, tissue=tissue)
    finally:
        try:
            del reference.obs[tmp]
        except Exception:
            pass


def plan_annotation_strategy(reference, panel_genes, label_key, *, target_depth=None,
                             tissue=None, registry=None, max_rounds=3,
                             coarsen_min_frac=0.25) -> dict:
    """Self-verify the reference<->panel match, auto-rerun a remediation ladder, and RECOMMEND
    supervised label transfer vs unsupervised de-novo clustering.

    Ladder (each rung re-scores :func:`reference_panel_match`; stops once the verdict reaches
    good/fair; NEVER routes to clustering on a single low metric - frac_resolvable is granularity-
    dependent, so it routes through coarsening first):

      1. score the reference as given;
      2. if fair/poor with confusable types -> COARSEN the confused groups and re-score;
      3. if still poor and a ``registry`` + ``tissue`` are given -> try the next best-overlap
         AVAILABLE tissue-matched reference (:func:`choose_reference` + :func:`load_reference`);
      4. decide: good/fair -> ``reference_transfer`` (supervised, at the coarsened granularity if
         step 2 helped); genuinely uninformative (still poor / overlap below the classifier floor)
         -> ``cluster`` (de-novo + marker naming).

    Pure and never raises; returns the decision plus the full per-rung ladder for the report/UI.
    """
    attempts: list[dict] = []

    def _snapshot(step, action, match, extra=None):
        g = match.get("global", {}) if isinstance(match, dict) else {}
        row = {"step": step, "action": action, "status": match.get("status"),
               "verdict": g.get("verdict"), "gene_overlap_frac": g.get("gene_overlap_frac"),
               "mean_f1": g.get("mean_f1"), "frac_resolvable": g.get("frac_resolvable"),
               "n_resolvable": g.get("n_resolvable"), "n_types": g.get("n_types")}
        if extra:
            row.update(extra)
        attempts.append(row)

    best = _safe_match(reference, panel_genes, label_key, target_depth, tissue=tissue)
    _snapshot("initial", "score", best, {"label_key": str(label_key)})
    best_ref_label = "provided"
    coarsen_map = None

    # Rung 2: coarsen the confusable labels.
    if best.get("status") == "ok" and best.get("global", {}).get("verdict") in ("fair", "poor"):
        mapping = merge_groups_from_confusion(best.get("per_type", {}), coarsen_min_frac)
        if any(k != v for k, v in mapping.items()):
            m2 = _rescore_coarsened(reference, panel_genes, label_key, mapping, target_depth, tissue=tissue)
            _snapshot("coarsen", "merge_confusable", m2,
                      {"merged_groups": sorted({v for k, v in mapping.items() if k != v})})
            if _verdict_key(m2) > _verdict_key(best):
                best, coarsen_map = m2, {k: v for k, v in mapping.items() if k != v}

    # Rung 3: reselect a better-overlap, tissue-matched, AVAILABLE reference.
    if best.get("global", {}).get("verdict") == "poor" and tissue and registry:
        tried: set = set()
        for cand in choose_reference(tissue, panel_genes, registry, top_n=max(3, int(max_rounds))):
            if sum(1 for a in attempts if a["step"] == "swap_reference") >= int(max_rounds):
                break
            key_c, path_c = cand.get("tissue_key"), cand.get("path")
            if not cand.get("available") or not path_c or key_c in tried:
                continue
            tried.add(key_c)
            try:
                alt, alt_key = load_reference(path_c, label_key=cand.get("label_key"),
                                              gene_name_col=cand.get("gene_name_col"))
                m3 = _safe_match(alt, panel_genes, alt_key, target_depth, tissue=tissue)
            except Exception as exc:            # a bad/unreadable candidate must not break the loop
                _snapshot("swap_reference", "load_failed", {"status": "error", "global": {}},
                          {"reference": key_c, "error": str(exc)[:200]})
                continue
            _snapshot("swap_reference", "score", m3, {"reference": key_c})
            if _verdict_key(m3) > _verdict_key(best):
                best, best_ref_label, coarsen_map = m3, key_c, None
            if best.get("global", {}).get("verdict") != "poor":
                break

    # Rung 4: route.
    g = best.get("global", {})
    verdict = g.get("verdict")
    overlap = g.get("gene_overlap_frac")
    overlap_txt = f"{overlap:.0%}" if overlap is not None else "n/a"
    below_floor = best.get("status") != "ok" or (overlap is not None and overlap < 0.2)
    tissue_mismatch = bool(best.get("global", {}).get("tissue_mismatch"))
    if tissue_mismatch:
        # A declared wrong-tissue reference: the internal separability verdict is irrelevant, do not
        # transfer it. This is the guard for the silent-nonsense case (a kidney atlas on a breast panel).
        mode = "cluster"
        reason = (best.get("tissue_check", {}) or {}).get("note", "") + \
            " De-novo clustering with marker-based naming is more honest than transferring a " \
            "wrong-tissue reference (whose internal separability says nothing about tissue match)."
    elif below_floor or verdict == "poor":
        mode = "cluster"
        reason = ("The reference cannot resolve this panel even after coarsening the confusable labels "
                  "and trying the best available tissue-matched reference "
                  f"(final verdict '{verdict}', {overlap_txt} gene overlap). De-novo clustering with "
                  "marker-based naming is more honest than forcing supervised labels onto it.")
    else:
        gran = " at the coarsened granularity the panel can support" if coarsen_map else ""
        mode = "reference_transfer"
        reason = (f"The reference resolves the panel ({verdict}: mean per-type F1 {g.get('mean_f1')}, "
                  f"{g.get('n_resolvable')}/{g.get('n_types')} types resolvable){gran}. Use supervised "
                  "label transfer, keeping per-cell abstention.")

    return {
        "status": best.get("status"),
        "recommended_mode": mode,               # 'reference_transfer' | 'cluster'
        "reason": reason,
        "reference_used": best_ref_label,
        "coarsen_map": coarsen_map,
        "initial_verdict": attempts[0].get("verdict"),
        "final_verdict": verdict,
        "final_global": g,
        "final_per_type": best.get("per_type", {}),
        "ladder": attempts,
        "clustering_nudge": best.get("clustering_nudge"),
        "tissue_mismatch": tissue_mismatch,
        "tissue_check": best.get("tissue_check"),
    }
