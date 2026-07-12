"""CellGuide marker provider - Cell Ontology-grounded markers from the CZI cell-type corpus.

The curated ``markers.LINEAGE_MARKERS`` stay the benchmarked default for known lineages (they are
tuned for targeted imaging panels; CellGuide's raw top computational markers are census-effect-ranked,
not panel markers). CellGuide adds two things the hardcoded dicts cannot:

  1. **Cell Ontology grounding** - :func:`resolve_cl` maps an arbitrary label ("T cell", "Treg",
     "fibroblast") onto a ``CL:xxxx`` id via the CellGuide metadata (name + synonyms), so novel /
     renamed / reference-transferred types get a real ontology anchor.
  2. **Grounded markers for those types** - :func:`markers_for_label` returns CANONICAL markers
     (CellGuide's literature-curated "canonical" set) when CellGuide has them, else COMPUTATIONAL markers
     (census-derived, ranked by effect size = ``marker_score`` on the tissue-agnostic aggregate),
     else nothing. This replaces "ask the LLM to invent marker genes" with a deterministic,
     corpus-grounded source that cannot hallucinate an off-panel gene.

Everything is fetch-with-fallback and disk-cached: the app must degrade to curated markers offline,
so no call here ever raises - failures return ``None`` / ``[]`` and the caller keeps the curated set.

Data source (undocumented static JSON, versioned by snapshot; treat as best-effort, never load-bearing):
  https://cellguide.cellxgene.cziscience.com/<snapshot>/{celltype_metadata,canonical_marker_genes/<CL>,computational_marker_genes/<CL>}.json
Ranking note: CellGuide's UI "Effect Size" column IS ``marker_score``; the canonical-quality T-cell
markers (CD3E/CD3D/IL7R) appear only on the tissue-AGNOSTIC aggregate - mixing tissue-specific records
surfaces census-abundant non-markers (SKAP1/RIPOR2), so :func:`markers_for_label` filters to it.
"""
from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import Optional
from urllib.request import Request, urlopen

_BASE = "https://cellguide.cellxgene.cziscience.com"
_TIMEOUT = float(os.environ.get("CELLGUIDE_TIMEOUT", "15"))
# Cache on shared /data (NEVER home - 100 GB quota). Overridable; the snapshot id namespaces it.
_CACHE = Path(os.environ.get(
    "CELLGUIDE_CACHE",
    os.path.join(os.path.expanduser("~"), ".cache", "spatialscribe", "cellguide"),
))
# A cell must express a computational marker in this fraction of cells AND clear this specificity to
# count - keeps the census tail (ubiquitous or noisy genes) out of the returned set.
_MIN_PC = 0.10
_MIN_SPECIFICITY = 0.5
_MIN_MARKER_SCORE = 0.0
# The raw computational file is ~1.3 MB (every organism x tissue). We never keep it: distil to the top
# markers per organism (tissue-agnostic aggregate) and cache only that - a few hundred bytes per type.
_CACHE_TOP_N = 20
# Canonical must contribute at least this many symbols to be preferred; below it, the richer
# computational set is used instead (a canonical set of 1-2 genes is thinner than computational).
_MIN_CANONICAL = 3


def _is_human(organism: str) -> bool:
    return "sapiens" in str(organism).lower() or str(organism).lower() in ("human", "homo sapiens")


def _canonical_organism(organism: str) -> str:
    """Fold 'human'/'mouse'/etc onto CellGuide's ``organism_ontology_term_label`` values."""
    o = str(organism).lower()
    if "sapiens" in o or o == "human":
        return "Homo sapiens"
    if "musculus" in o or o == "mouse":
        return "Mus musculus"
    return organism                                            # already a canonical label, or unknown

_snapshot: Optional[str] = None
_meta_index: Optional[dict] = None


def _offline() -> bool:
    """Read the env each call (not at import), so the hermetic test suite can pin it after import."""
    return os.environ.get("CELLGUIDE_OFFLINE") == "1"


def _http_get(url: str) -> Optional[str]:
    if _offline():
        return None
    try:
        req = Request(url, headers={"Accept-Encoding": "gzip", "User-Agent": "spatialscribe"})
        with urlopen(req, timeout=_TIMEOUT) as r:               # noqa: S310 (fixed https host)
            raw = r.read()
            if r.headers.get("Content-Encoding") == "gzip":
                import gzip
                raw = gzip.decompress(raw)
            return raw.decode("utf-8") if raw else None
    except Exception:
        return None


def snapshot() -> Optional[str]:
    """The current CellGuide snapshot id (memoised per process).

    Online: fetch it and remember it in ``<cache>/.latest`` so the cache survives an outage. Offline
    / unreachable: reuse that last-known id, so files cached while online are still served (stale is
    better than nothing). ``None`` only when offline AND nothing was ever cached.
    """
    global _snapshot
    if _snapshot is None:
        s = _http_get(f"{_BASE}/latest_snapshot_identifier")
        if s and s.strip().isdigit():
            _snapshot = s.strip()
            try:
                _CACHE.mkdir(parents=True, exist_ok=True)
                (_CACHE / ".latest").write_text(_snapshot)
            except Exception:
                pass
        else:
            try:
                _snapshot = (_CACHE / ".latest").read_text().strip() or None
            except Exception:
                _snapshot = None
    return _snapshot


def _fetch_json(rel: str, cache_raw: bool = True) -> Optional[object]:
    """GET ``<base>/<snapshot>/<rel>`` with a disk cache. Returns the parsed JSON or None.

    The cache is keyed by snapshot, so a new snapshot transparently invalidates it. A test can seed
    the cache dir (``CELLGUIDE_CACHE``) + pin ``CELLGUIDE_OFFLINE=1`` to run fully hermetically.
    ``cache_raw=False`` skips writing the body (the ~1.3 MB computational file is distilled + cached
    small by the caller instead), but still reads a pre-seeded cache file if one exists.
    """
    snap = snapshot() or "offline"                             # offline: pre-seeded cache namespace
    cache_file = _CACHE / snap / rel
    if cache_file.exists():
        try:
            return json.loads(cache_file.read_text())
        except Exception:
            pass
    body = _http_get(f"{_BASE}/{snap}/{rel}")
    if not body:
        return None
    try:
        obj = json.loads(body)
    except Exception:
        return None
    if cache_raw:
        try:
            cache_file.parent.mkdir(parents=True, exist_ok=True)
            cache_file.write_text(body)
        except Exception:
            pass                                                # cache is best-effort
    return obj


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(s).lower()).strip()


def _metadata_index() -> dict:
    """``{normalised name/synonym: CL_id}`` from celltype_metadata.json. Empty when unreachable."""
    global _meta_index
    if _meta_index is not None:
        return _meta_index
    meta = _fetch_json("celltype_metadata.json")
    idx: dict[str, str] = {}
    if isinstance(meta, dict):
        items = [(cid, v) for cid, v in meta.items() if isinstance(v, dict)]
        # TWO passes so a canonical NAME always wins globally: an earlier record's synonym must not
        # shadow a later record's name (e.g. beta cell's synonym "insulin-secreting cell" vs the
        # distinct term whose NAME is "insulin secreting cell").
        for cid, v in items:
            if v.get("name"):
                idx.setdefault(_norm(v["name"]), cid)
        for cid, v in items:
            for syn in (v.get("synonyms") or []):
                idx.setdefault(_norm(syn), cid)
    _meta_index = idx
    return idx


# A composite curated lineage is not a Cell Ontology term ("Stromal/CAF", "Malignant/Melanocyte");
# anchor those explicitly. Keys are NORMALISED (via _norm) so the lookup below matches; values are
# ontology names/synonyms resolved through the metadata index. Anything not here falls through to
# direct name/synonym matching.
_LINEAGE_ALIAS = {
    "malignant melanocyte": "melanocyte",
    "stromal caf": "fibroblast",
    "b plasma": "B cell",
    "glial neural": "glial cell",
    "epithelial tumor": "epithelial cell",
    "basal myoepithelial": "myoepithelial cell",
}


def resolve_cl(label: str) -> Optional[str]:
    """Map a cell-type label onto a ``CL:xxxx`` id via CellGuide metadata. None if no match.

    Exact (case/punctuation-insensitive) matching on the ontology name and its synonyms - never fuzzy
    nearest-string (same rule as annotate._canonicalise: a wrong ontology anchor is worse than none).
    """
    if not label:
        return None
    idx = _metadata_index()
    if not idx:
        return None
    # Try the label, its composite alias, and a "<label> cell" suffix - many section labels are the
    # bare lineage ("Mast", "Myeloid", "Endothelial") while the ontology name is "mast cell" etc.
    # Still EXACT match against real names/synonyms, so the suffix only resolves a genuine CL term.
    nl = _norm(label)
    for cand in (nl, _norm(_LINEAGE_ALIAS.get(nl, "")), f"{nl} cell"):
        if cand and cand in idx:
            return idx[cand]
    return None


def _cl_path(cl_id: str) -> str:
    return cl_id.replace(":", "_")


def _canonical_freq(cl_id: str, tissue: Optional[str]) -> dict:
    """``{symbol: cross-tissue citation count}`` from the canonical (literature-curated) set.

    Returns COUNTS, not an order: cross-tissue frequency alone floods the head with ubiquitous
    genes and buries flagship markers (COL1A1 for fibroblast), so the caller re-ranks by effect size.
    Optionally restricts to a requested ``tissue``.
    """
    recs = _fetch_json(f"canonical_marker_genes/{_cl_path(cl_id)}.json")
    if not isinstance(recs, list):
        return {}
    from collections import Counter
    if tissue:
        recs = [r for r in recs if str(r.get("tissue", "")).lower() == tissue.lower()] or recs
    return dict(Counter(r.get("symbol") for r in recs
                        if isinstance(r, dict) and r.get("symbol")))


def _distil_computational(recs: list) -> dict:
    """Raw computational records -> ``{organism_label: [[symbol, marker_score], ...top _CACHE_TOP_N]}``.

    Tissue-agnostic ONLY: the CellGuide UI's canonical-quality ranking is the "All Tissues" aggregate
    (tissue == None); mixing tissue-specific records surfaces census-abundant non-markers. Ranked by
    ``marker_score`` (the UI "Effect Size"). The SCORE is kept (not just the symbol) so canonical
    markers can be re-ranked by effect size too. Mouse records carry mouse-case symbols, so one distil
    serves both organisms.
    """
    from collections import defaultdict
    by_org: dict[str, list] = defaultdict(list)
    for r in recs:
        if not isinstance(r, dict):
            continue                                           # honour the "never raises" contract
        d = r.get("groupby_dims", {})
        if d.get("tissue_ontology_term_label") is not None:
            continue
        if (r.get("pc", 0) >= _MIN_PC and r.get("specificity", 0) >= _MIN_SPECIFICITY
                and r.get("marker_score", 0) > _MIN_MARKER_SCORE and r.get("symbol")):
            by_org[d.get("organism_ontology_term_label", "?")].append(r)
    out: dict[str, list] = {}
    for org, rs in by_org.items():
        rs.sort(key=lambda r: -r.get("marker_score", 0))
        seen, ranked = set(), []
        for r in rs:
            g = r["symbol"]
            if g not in seen:
                seen.add(g)
                ranked.append([g, round(float(r.get("marker_score", 0)), 4)])
            if len(ranked) >= _CACHE_TOP_N:
                break
        out[org] = ranked
    return out


def _computational_ranked(cl_id: str, organism: str) -> list:
    """``[[symbol, marker_score], ...]`` for ``organism`` (tissue-agnostic), via the distilled cache.

    The distilled ``{organism: [[sym, score], ...]}`` (a few hundred bytes) is cached instead of the
    ~1.3 MB raw file. Both organisms come from one fetch, so mouse and human share the artifact.
    """
    snap = snapshot() or "offline"
    dist_file = _CACHE / snap / "computational_top" / f"{_cl_path(cl_id)}.json"
    dist = None
    if dist_file.exists():
        try:
            dist = json.loads(dist_file.read_text())
        except Exception:
            dist = None
    if dist is None:
        raw = _fetch_json(f"computational_marker_genes/{_cl_path(cl_id)}.json", cache_raw=False)
        if not isinstance(raw, list):
            return []
        dist = _distil_computational(raw)
        try:
            dist_file.parent.mkdir(parents=True, exist_ok=True)
            dist_file.write_text(json.dumps(dist))
        except Exception:
            pass
    return dist.get(organism, []) or []


def markers_for_label(label: str, organism: str = "Homo sapiens",
                      tissue: Optional[str] = None, n: int = 12) -> dict:
    """Ontology-grounded markers for one cell-type label.

    Returns ``{"cl_id", "source", "markers"}``. ``source`` is ``"canonical"`` (literature-curated,
    preferred), ``"computational"`` (census-derived, effect-size ranked), or ``"none"`` when
    CellGuide cannot resolve or has nothing - the caller then keeps its curated markers. Never raises.

    Both sources are ranked by EFFECT SIZE (the computational ``marker_score``, the CellGuide UI
    column). Canonical carries only symbols and cross-tissue counts, so its genes are ORDERED by their
    computational effect size (flagship markers first; a canonical gene absent from the computational
    top falls to the tail by citation count) rather than by the misleading alphabetical/frequency order.

    Organism: canonical markers are HUMAN symbols with no organism dimension, so canonical is used only
    for human; mouse resolves straight to the computational aggregate (mouse-case symbols, ``Cd3g``).
    Canonical must yield >= ``_MIN_CANONICAL`` symbols to be trusted, else the richer computational set
    is used. Both are ontology-grounded on the same CL id.
    """
    cl = resolve_cl(label)
    if cl is None:
        return {"cl_id": None, "source": "none", "markers": []}
    org = _canonical_organism(organism)
    ranked = _computational_ranked(cl, org)                    # [[symbol, score], ...] effect-ordered
    score = {s: sc for s, sc in ranked}

    if _is_human(org):
        freq = _canonical_freq(cl, tissue)
        if len(freq) >= _MIN_CANONICAL:
            # order canonical genes by computational effect size, then by citation count for the
            # canonical-only tail; this keeps flagship markers ahead of the alphabetical head.
            genes = sorted(freq, key=lambda g: (-score.get(g, -1.0), -freq[g], g))
            return {"cl_id": cl, "source": "canonical", "markers": genes[:n]}

    comp = [s for s, _ in ranked][:n]
    if comp:
        return {"cl_id": cl, "source": "computational", "markers": comp}
    return {"cl_id": cl, "source": "none", "markers": []}


def markers_for_labels(labels, organism: str = "Homo sapiens",
                       tissue: Optional[str] = None, n: int = 12) -> dict:
    """``{label: {cl_id, source, markers}}`` for many labels (one metadata fetch, cached per CL)."""
    return {str(x): markers_for_label(str(x), organism, tissue, n) for x in dict.fromkeys(labels)}


# --------------------------------------------------------------------------- #
# Cell Ontology (CL) hierarchy - "which cell types may be merged, and which may not".
#
# CellGuide's celltype_metadata.json carries only name/synonyms/description (no is-a edges), so the
# hierarchy comes from the CL `cl-basic.obo` (fetched once, cached next to the CellGuide cache). Two
# cell types are MERGEABLE only when they descend from a COMMON major-lineage anchor - so CD4 and CD8
# T cells (both leukocytes) may coarsen to "T cell", but a mast cell (leukocyte) and an epithelial
# tumour cell (epithelial) share NO anchor and must NOT be merged. Positive anchors (not a "too-broad"
# blocklist) are robust to cl-basic's flat upper terms (motile cell / nucleate cell).
# --------------------------------------------------------------------------- #
_CL_OBO_URL = "https://github.com/obophenotype/cell-ontology/releases/latest/download/cl-basic.obo"
# Major-lineage anchors (CL id -> readable compartment). Verified against cl-basic 2026-07.
_COMPARTMENT_ANCHORS: dict[str, str] = {
    "CL:0000988": "hematopoietic", "CL:0000738": "immune",       # hematopoietic cell / leukocyte
    "CL:0000066": "epithelial", "CL:0000115": "endothelial",
    "CL:0000499": "stromal", "CL:0008019": "stromal", "CL:0002320": "stromal",  # stromal/mesenchymal/connective
    "CL:0000540": "neural", "CL:0000125": "neural", "CL:0002319": "neural",     # neuron/glia/neural
    "CL:0000187": "muscle",
}
_cl_isa: Optional[dict] = None


def _cl_obo_path() -> Path:
    return Path(os.environ.get("SPATIALSCRIBE_CL_OBO", str(_CACHE / "cl-basic.obo")))


def _cl_graph() -> dict:
    """``{CL_id: set(parent CL_ids)}`` from cl-basic.obo (is_a edges only). Memoised.

    Reads the cached obo; if absent and online, fetches + caches it once. Returns ``{}`` when the obo
    is unavailable (offline, no cache) so callers degrade to the hardcoded compartment map. Never raises.
    """
    global _cl_isa
    if _cl_isa is not None:
        return _cl_isa
    p = _cl_obo_path()
    text: Optional[str] = None
    if p.exists():
        try:
            text = p.read_text()
        except Exception:
            text = None
    if text is None and not _offline():
        text = _http_get(_CL_OBO_URL)                          # urllib follows the GitHub redirect
        if text:
            try:
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text(text)
            except Exception:
                pass
    g: dict[str, set] = {}
    if text:
        cur = None
        for line in text.splitlines():
            if line == "[Term]":
                cur = None
                continue
            m = re.match(r"id:\s*(CL:\d+)", line)
            if m:
                cur = m.group(1); g.setdefault(cur, set()); continue
            if cur is None:
                continue
            m = re.match(r"is_a:\s*(CL:\d+)", line)             # ignore any {is_inferred=...} qualifier
            if m:
                g[cur].add(m.group(1))
    _cl_isa = g
    return g


def cl_ancestors(cl_id: str) -> set:
    """All is-a ancestors of ``cl_id`` (including itself). Empty if the ontology is unavailable."""
    g = _cl_graph()
    out, stack = {cl_id}, [cl_id]
    while stack:
        x = stack.pop()
        for p in g.get(x, ()):
            if p not in out:
                out.add(p); stack.append(p)
    return out if g else set()


def cl_compartments(cl_id: str) -> set:
    """The major-lineage compartment(s) ``cl_id`` descends from (see :data:`_COMPARTMENT_ANCHORS`)."""
    if not cl_id:
        return set()
    anc = cl_ancestors(cl_id)
    return {name for a, name in _COMPARTMENT_ANCHORS.items() if a in anc}


# Keyword -> compartment fallback for labels that do NOT resolve to a CL id: study-specific composite
# names carry a lineage token (MFAP5_IGFBP6_fibroblast, LYVE1_macrophage, Capillary_EC, PIP_mammary_
# luminal_cell). A token/substring match on that lineage word (NOT fuzzy nearest-string) recovers the
# major compartment so the merge guard still works on RCTD / custom vocabularies. Longest-token wins.
_KEYWORD_COMPARTMENT: list[tuple[str, tuple[str, ...]]] = [
    ("immune", ("macrophage", "monocyte", "dendritic", "myeloid", "microglia", "neutrophil",
                "granulocyte", "mast", "kupffer", "lymphocyte", "leukocyte", "treg", "lymphoid",
                "t cell", "b cell", "nk cell", "plasma", "natural killer", "megakaryocyte", "basophil",
                "eosinophil", "langerhans")),
    ("endothelial", ("endothel", "ec", "capillary", "arterial", "venous", "lymphatic")),
    ("stromal", ("fibroblast", "myofibroblast", "caf", "stromal", "mesenchym", "pericyte", "adipocyte",
                 "mural", "stellate")),
    ("epithelial", ("epithel", "luminal", "basal", "keratinocyt", "carcinoma", "melanocyt", "hepatocyt",
                    "enterocyte", "goblet", "secretory", "tumor", "tumour", "club", "ciliated",
                    "alveolar", "acinar", "ductal", "podocyte", "tubule", "tubular", "nephron")),
    ("neural", ("neuron", "astrocyt", "oligodendrocyt", "opc", "ependymal", "glia", "schwann",
                "neural", "neuronal")),
    ("muscle", ("myocyte", "cardiomyocyte", "smooth muscle", "myotube")),
]


def _label_compartments(label: str) -> set:
    """Major compartment(s) for a label: the Cell Ontology first (:func:`cl_compartments`), else a
    lineage-keyword fallback (:data:`_KEYWORD_COMPARTMENT`) so composite study labels still classify."""
    cl = resolve_cl(label)
    comp = cl_compartments(cl) if cl else set()
    if comp:
        return comp
    toks = set(_norm(label).split())
    norm = _norm(label)
    for compartment, kws in _KEYWORD_COMPARTMENT:
        for k in kws:
            if (k in toks) or (" " in k and k in norm):     # short codes match a token; phrases substring
                return {compartment}
    return set()


def labels_mergeable(label_a: str, label_b: str) -> dict:
    """Should two cell-type labels be merged (coarsened) - grounded in the Cell Ontology?

    Resolves each label to a ``CL`` id (:func:`resolve_cl`) and returns whether they share a major
    lineage anchor. ``mergeable`` is ``True`` (same compartment - coarsening is defensible), ``False``
    (different major lineages - merging one into the other is biologically unsound), or ``None`` when
    the ontology cannot decide (a label did not resolve to CL, or to no known anchor, or the graph is
    offline) - the caller then falls back to its own heuristic. Never raises.
    """
    ca, cb = resolve_cl(label_a), resolve_cl(label_b)
    compa = _label_compartments(label_a)
    compb = _label_compartments(label_b)
    if not compa or not compb:
        return {"mergeable": None, "cl_a": ca, "cl_b": cb,
                "compartments_a": sorted(compa), "compartments_b": sorted(compb),
                "shared_compartments": []}
    shared = compa & compb
    return {"mergeable": bool(shared), "cl_a": ca, "cl_b": cb,
            "compartments_a": sorted(compa), "compartments_b": sorted(compb),
            "shared_compartments": sorted(shared)}


def _reset_for_tests() -> None:
    """Drop the process-level snapshot/metadata/ontology memo (unit tests seed a fresh cache per case)."""
    global _snapshot, _meta_index, _cl_isa
    _snapshot = _meta_index = _cl_isa = None
