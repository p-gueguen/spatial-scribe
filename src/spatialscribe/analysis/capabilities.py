"""The capability registry - one declarative entry per analysis operation.

What it does
------------
Collapses the three-times-repeated expression of every operation (a pure ``analysis``
function, a wizard step, a copilot tool) into ONE :class:`Capability` that both frontends
consume. Each entry binds a thin adapter over the existing pure function to:

* a human ``label`` + ``description`` (used by the wizard help AND the copilot tool schema),
* a JSON-schema ``params`` block (drives the copilot tool schema and, later, wizard widgets),
* ``requires`` / ``produces`` :class:`~spatialscribe.analysis.keys.Key`s - the explicit
  inter-stage data contract, consumed by the prerequisite check and :mod:`status`.

``run()`` is the single dispatch both the wizard and the copilot call: it checks
prerequisites (returning a structured :class:`~spatialscribe.analysis.errors.PrerequisiteError`
rather than throwing), invokes the adapter, and returns a :class:`RunResult` carrying the
JSON-able value plus a provenance ``record`` (so copilot actions land in the re-runnable
export too). ``ensure()`` adds compute-if-absent-or-forced idempotency.

The adapters are deliberately thin: the pure ``analysis`` functions keep their signatures
untouched, so this module is the one place that knows how to translate
``(adata, ctx, **params)`` into each real call. Heavy imports stay inside the adapters so
importing the registry is cheap.

How to use it
-------------
>>> from .capabilities import run, ensure, RunContext, copilot_tools
>>> ctx = RunContext(tissue="melanoma")
>>> res = run(adata, "cluster", {"resolution": 0.5}, ctx)
>>> res.ok, res.value, res.record
>>> copilot_tools()          # -> list of Anthropic tool schemas for the copilot-exposed subset
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

from .errors import PrerequisiteError, to_error_dict
from .keys import Key, Obs, Obsm, Uns, Varm


# --------------------------------------------------------------------------- #
# Run context + result
# --------------------------------------------------------------------------- #
@dataclass
class RunContext:
    """Ambient inputs an adapter may need beyond ``adata`` and its explicit params.

    Keeps the frontends from threading tissue/markers/LLM-availability through every call.
    """

    tissue: str = "melanoma"
    # Data-step checkbox: "this sample contains malignant cells". None = unanswered -> fall back to the
    # tissue-keyword heuristic (which mis-fires BOTH ways; see cnv.is_tumour_context). True/False is the
    # scientist's explicit answer and always wins, because nothing can know a priori whether a section
    # contains tumour: that is precisely what the malignant callers exist to decide.
    is_tumour: bool | None = None
    marker_sets: dict | None = None
    use_llm: bool = False
    artifacts: list = field(default_factory=list)   # figure/view artifacts a capability emits for the app
    reference: object = None          # in-memory reference AnnData chosen at the Panel-check step (optional)
    ref_label_key: str | None = None  # its detected cell-type column
    progress: Callable[[float, str], None] | None = None   # coarse progress sink (installed by the app)
    # A load_section capability appends the section it loaded here (mutable out-param, like `artifacts`):
    # the copilot loop rebinds to it for the rest of the turn, and the app swaps it onto the session.
    loaded: list = field(default_factory=list)

    def tick(self, frac: float, label: str) -> None:
        """Report coarse progress (0.0-1.0) for the running capability. No-op when no callback is
        installed, so the engine stays UI-agnostic. Clamps ``frac`` to [0,1] and never raises - the
        wrappers pass this to long engine functions, so a broken sink must not break the analysis."""
        if self.progress is None:
            return
        try:
            self.progress(max(0.0, min(1.0, float(frac))), str(label))
        except Exception:
            pass

    def markers(self, adata=None) -> dict:
        """Resolved lineage marker sets for this run.

        Explicit ``marker_sets`` override wins; else the tissue's curated set; else - for an
        unrecognised free-text tissue, when ``adata`` is given and an API key is present - a
        Claude-generated panel grounded to the section's genes (cached, so every step agrees).
        Pass ``adata`` wherever the panel is available so arbitrary tissues resolve.
        """
        if self.marker_sets is not None:
            return self.marker_sets
        from . import markers as _m
        panel = _panel_genes(adata) if adata is not None else None
        resolved = _m.resolve_markers(self.tissue, panel_genes=panel, use_llm=self.use_llm)
        # If the section carries cell_type labels the tissue-resolved sets don't cover - reference-
        # transferred / renamed / non-tissue labels the curated dictionaries never anticipated, e.g.
        # Allen mouse-brain types (Lamp5, L2_3 IT, Astro) - ground markers on the ACTUAL categories via
        # CellGuide + LLM, else verify / marker dot-plot / panel-check score against wrong-keyed markers
        # and come back empty. Curated tissues (breast/melanoma) already cover their labels, so keep them.
        if adata is not None and panel is not None and "cell_type" in getattr(adata, "obs", {}):
            from . import annotate as _an
            cats = [str(c) for c in adata.obs["cell_type"].astype("category").cat.categories
                    if not _an.is_abstention(str(c))]
            if cats and sum(1 for c in cats if c in resolved) < max(2, len(cats) // 2):
                by_type = _m.markers_for_types(cats, panel_genes=panel, tissue=self.tissue,
                                                use_llm=self.use_llm)
                # If the labels are a reference-transfer vocabulary the curated + CellGuide/LLM sets
                # STILL don't cover (fine MOCA/Allen atlas labels), ground the remaining markers on the
                # REFERENCE's own per-label DEGs. Without this a cell whose assigned type has no marker
                # set gets a NaN PMP, so the AQI's C term (C = min(median PMP, RETENTION)) collapses to
                # a coverage artifact (measured 0.08 on the E14 embryo, 14/91 covered). Gated on an
                # already-loaded reference (resolve_reference does no I/O then), so it never triggers a
                # multi-GB atlas load just to score markers. reference-DEG floor < CellGuide/LLM <
                # curated (best evidence wins per category).
                covered = sum(1 for c in cats if (c in resolved) or (by_type and by_type.get(c)))
                if covered < max(2, len(cats) // 2) and getattr(self, "reference", None) is not None:
                    from . import reference as _ref
                    ref_obj, ref_key, _ = self.resolve_reference()
                    ref_sets = _ref.reference_marker_sets(ref_obj, ref_key or "cell_type", panel) if ref_obj is not None else {}
                    if ref_sets:
                        merged = dict(ref_sets)
                        if by_type:
                            merged.update({k: v for k, v in by_type.items() if v})
                        merged.update({c: resolved[c] for c in cats if c in resolved})
                        return merged
                if by_type:
                    return by_type
        # Caveat: an UNRECOGNISED tissue fell back to the generic breast/carcinoma EPITHELIAL_LINEAGES
        # set (no curated set, no LLM/CellGuide grounding). Stamp a warning so the panel / typability /
        # composition surfaces can flag it - else a normal kidney is silently typed 'Epithelial/Tumor'
        # with no tell. Fires only on the true generic fallback (a curated or grounded set is fine).
        if adata is not None and resolved is _m.EPITHELIAL_LINEAGES and not _m.tissue_has_curated_set(self.tissue):
            try:
                adata.uns["marker_set_warning"] = {
                    "tissue": self.tissue,
                    "marker_set": "EPITHELIAL_LINEAGES (generic breast/carcinoma fallback)",
                    "message": (f"No curated marker set matches the tissue '{self.tissue}'. Falling back to "
                                "the generic breast/carcinoma lineages, so labels such as 'Epithelial/Tumor' "
                                "may be wrong-tissue. Attach a tissue-matched reference (the reference "
                                "ranking usually picks the right one) or enable the LLM/CellGuide marker path."),
                }
            except Exception:
                pass
        return resolved

    def label_context(self) -> str:
        """Free-text context string for the LLM annotator, species-aware. Adds the default
        'human' assumption for the tumour demos, but never double-prefixes when the tissue text
        already names a species (so 'mouse brain' stays 'mouse brain', not 'human mouse brain')."""
        t = (self.tissue or "").strip()
        named = ("mouse", "murine", "mus musculus", "human", "patient", "rat", "macaque", "zebrafish")
        return t if any(s in t.lower() for s in named) else f"human {t}"

    def resolve_reference(self, reference_path=None, label_key=None):
        """Resolve the reference to ``(AnnData|None, label_key|None, source_str)`` - the ONE place
        that turns (explicit path | in-memory ``ctx.reference`` | ``ctx.ref_label_key``) into a
        concrete reference. Prefers an explicit path, else the reference chosen at the Panel-check step.
        Never raises (a failed load returns ``(None, None, 'could not load ...')``) so every
        reference consumer degrades gracefully."""
        from . import reference as _ref
        if reference_path:
            try:
                ref, key = _ref.load_reference(reference_path, label_key=label_key)
                return ref, key, f"file: {reference_path}"
            except Exception as exc:
                return None, None, f"could not load reference: {exc}"
        if self.reference is not None:
            # Validate the requested key against the in-memory reference's columns (mirrors the
            # load_reference path); a key that is not a real column would KeyError downstream in
            # panel_resolvability. Fall back to the detected label column when it is missing/wrong.
            cols = set(self.reference.obs.columns)
            key = next((k for k in (label_key, self.ref_label_key) if k and k in cols), None) \
                or _ref.detect_label_key(self.reference)
            return self.reference, key, "reference chosen at the Panel-check step"
        return None, None, "no reference: pass reference_path or choose one at the Panel-check step"


@dataclass
class RunResult:
    """Outcome of a capability call: a JSON-able ``value`` + provenance ``record`` or ``error``."""

    name: str
    value: object = None
    record: dict | None = None
    error: dict | None = None

    @property
    def ok(self) -> bool:
        return self.error is None


@dataclass(frozen=True)
class Capability:
    """One registry entry: metadata + the adapter + the data contract."""

    name: str
    label: str
    description: str
    fn: Callable
    params: dict = field(default_factory=dict)
    required_params: tuple[str, ...] = ()
    requires: tuple[Key, ...] = ()
    produces: tuple[Key, ...] = ()
    copilot_exposed: bool = False
    valid: Callable | None = None   # optional predicate(adata)->bool; ensure() recomputes if it fails

    def to_tool_schema(self) -> dict:
        """Render the Anthropic tool schema for the copilot (name/description/input_schema)."""
        return {
            "name": self.name,
            "description": self.description,
            "input_schema": {
                "type": "object",
                "properties": dict(self.params),
                "required": list(self.required_params),
            },
        }


# --------------------------------------------------------------------------- #
# Adapters - thin translation from (adata, ctx, **params) to the pure function.
# --------------------------------------------------------------------------- #
def _panel_genes(adata) -> list:
    """Biological panel genes for ``adata``: controls stripped when the control mask is present,
    else all var_names. The one place adapters derive the panel from the loaded section."""
    if "control" in adata.var.columns:
        return adata.var_names[~adata.var["control"].to_numpy(dtype=bool)].tolist()
    return adata.var_names.tolist()


def _cap_describe(adata, ctx, **_):
    from . import qc as _qc
    if not Obs.TOTAL_COUNTS.present(adata):
        _qc.compute_qc(adata)
    comp = {}
    if Obs.CELL_TYPE.present(adata):
        from . import export as _ex
        counts, _abstained = _ex.composition_table(adata)     # typed cells only; not pseudo-labels
        comp = {str(k): int(v) for k, v in counts.head(12).items()}
    return {"n_cells": int(adata.n_obs), "composition": comp, "qc": _qc.qc_summary(adata)}


def _cap_panel_check(adata, ctx, **_):
    from . import panel_check as _pc
    pc = _pc.check_panel(_panel_genes(adata), marker_sets=ctx.markers(adata))
    adata.uns["panel_check"] = pc
    cov = pc["coverage"]
    merge_groups = pc.get("merge_groups", [])
    # A concrete "should I merge these?" nudge: groups the panel cannot separate (no private on-panel
    # marker). This is the panel-mismatch half of "why a cell could not be confidently labeled".
    nudge = ("Consider merging these groups - the panel has no private marker to tell them apart: "
             + "; ".join(" + ".join(g) for g in merge_groups) + ".") if merge_groups else \
            "No panel-driven merges suggested - every resolvable type has a private on-panel marker."
    return {
        "resolvable": [ct for ct, d in cov.items() if d["status"] == "green"],
        "cannot_resolve": [ct for ct, d in cov.items() if d["status"] == "red"],
        "confusable": [p["pair"] for p in pc["confusable_pairs"]],
        "merge_suggestions": merge_groups,
        "merge_nudge": nudge,
        "coverage": cov,
    }


def _valid_panel_check(adata) -> bool:
    """Invalid if an h5ad round-trip split the '/'-keyed cell types into nested groups (see is_valid)."""
    from . import panel_check as _pc
    return _pc.is_valid(adata.uns.get("panel_check"))


def _cap_reference_match(adata, ctx, reference_path=None, label_key=None, **_):
    # How well does a single-cell reference fit this panel? Resolves a reference via
    # ctx.resolve_reference (explicit reference_path, else the one chosen at the Panel-check step);
    # degrades to a 'no_reference' skip otherwise. Never raises.
    import numpy as np

    from . import reference as _ref

    panel_genes = _panel_genes(adata)
    ref, key, src = ctx.resolve_reference(reference_path, label_key)
    if ref is None or key is None:
        result = {"status": "no_reference", "global": {}, "per_type": {},
                  "clustering_nudge": None, "reference_source": src}
    else:
        depth = (float(np.median(np.asarray(adata.obs["total_counts"], dtype=float)))
                 if "total_counts" in adata.obs else None)
        result = _ref.reference_panel_match(ref, panel_genes, key, target_depth=depth, tissue=ctx.tissue)
        result["reference_source"] = src
    # When the fit is not GOOD, suggest better-matched references from the registry (the "propose a
    # new reference dataset" half of the low-match recommendation). Keyword + panel-gene-overlap
    # ranked; excludes the reference currently in use. Degrades to [] on any registry/read error.
    rec = result.get("recommendation") or {}
    alternatives: list[dict] = []
    if rec.get("verdict") in ("fair", "poor"):
        try:
            cur = (reference_path or "")
            ranked = _ref.choose_reference(ctx.tissue, panel_genes=panel_genes, top_n=4)
            alternatives = [r for r in ranked if r.get("path") != cur][:3]
        except Exception:
            alternatives = []
    result["alternatives"] = alternatives
    adata.uns["reference_match"] = result          # satisfies produces=(Uns.REFERENCE_MATCH,)
    g = result.get("global", {})
    return {"status": result.get("status"), "reference_source": result.get("reference_source", src),
            "match_score": result.get("match_score"), "recommendation": rec,
            "global": g, "clustering_nudge": result.get("clustering_nudge"),
            "n_types": g.get("n_types"), "n_resolvable": g.get("n_resolvable"),
            "alternatives": alternatives,
            "not_resolvable": [ct for ct, d in result.get("per_type", {}).items()
                               if not d.get("resolvable")]}


def _cap_auto_select_reference(adata, ctx, tissue=None, allow_fetch=False, **_):
    # Free-text tissue -> auto-CHOOSE + LOAD the best-matched precomputed reference (registry, or a
    # live CELLxGENE gget fetch when allow_fetch), then score its panel match. Sets the chosen
    # reference on ctx so a same-turn reference_match / reference_transfer uses it. Never raises;
    # degrades to a 'no_reference' status (the caller should cluster instead of forcing a bad transfer).
    from . import reference as _ref

    t = (tissue or ctx.tissue or "").strip()
    sel = _ref.auto_select_reference(t, panel_genes=_panel_genes(adata), allow_fetch=bool(allow_fetch))
    out = {"status": sel.get("status"), "chosen": sel.get("chosen"), "source": sel.get("source"),
           "message": sel.get("message"), "ranked": sel.get("ranked", [])}
    if sel.get("ref") is not None:
        ctx.reference = sel["ref"]
        ctx.ref_label_key = sel["label_key"]
        out["n_ref_cells"] = int(sel["ref"].n_obs)
        m = _cap_reference_match(adata, ctx)          # computes + stores uns['reference_match']
        out["match"] = {"match_score": m.get("match_score"), "recommendation": m.get("recommendation"),
                        "global": m.get("global"), "n_resolvable": m.get("n_resolvable"),
                        "n_types": m.get("n_types"), "clustering_nudge": m.get("clustering_nudge")}
    return out


# Optional storage-mount roots to try when a supplied path is missing its mount prefix
# (e.g. "/projects/foo/..." under a fixed mount). Empty by default; set SPATIALSCRIBE_SERVER_ROOTS
# (os.pathsep-separated) if your data lives under fixed mount roots.
import os as _os  # module-level; _resolve_server_path also imports os locally
_SERVER_ROOTS = tuple(
    r for r in _os.environ.get("SPATIALSCRIBE_SERVER_ROOTS", "").split(_os.pathsep) if r
)


def _resolve_server_path(p: str) -> str | None:
    """Return an existing filesystem path for a user/LLM-supplied ``p``, or None.

    Tries ``p`` as given, then rooted at each configured storage root (SPATIALSCRIBE_SERVER_ROOTS),
    so a path that dropped its mount prefix still loads.
    Deterministic (checks the filesystem), so it corrects the copilot's paths AND the Load tab's alike."""
    import os

    if not p:
        return None
    if os.path.exists(p):
        return p
    rel = p.lstrip("/")
    for root in _SERVER_ROOTS:
        cand = os.path.join(root, rel)
        if os.path.exists(cand):
            return cand
    return None


def _cap_load_section(adata, ctx, path=None, tissue=None, auto_reference: bool = True,
                      allow_fetch: bool = False, section_path=None, dataset_path=None,
                      data_path=None, **_):
    # Load a spatial section from a SERVER-SIDE path (a Xenium/CosMx/MERSCOPE output folder or a
    # .h5ad), make it the active section, INFER its tissue from the panel metadata, and - unless
    # auto_reference is off - auto-SELECT the best-matched reference and recommend supervised transfer
    # vs de-novo clustering. A capability cannot replace `adata` in place, so the loaded section is
    # stashed on ctx.loaded; the copilot loop rebinds to it and the app swaps it onto the session.
    import os

    from . import io as _io

    # Accept the obvious path-arg aliases: a model that emits the call as markup (the vLLM failure
    # mode) does not always follow the schema and may name the path 'section_path'/'dataset_path'.
    p = (path or section_path or dataset_path or data_path or "").strip()
    if not p:
        raise ValueError("load_section needs a 'path' (a server-side Xenium/CosMx/MERSCOPE output "
                         "folder or a .h5ad).")
    resolved = _resolve_server_path(p)      # correct a path missing its mount prefix (/projects -> /data)
    if resolved is None:
        raise ValueError(f"path not found on the server: {p!r}")
    p = resolved
    sample = _io.load(p)                    # raises a clear ValueError on a transcripts-only/unsupported dir
    new = sample.adata
    if "spatial" not in new.obsm:
        raise ValueError("loaded section has no spatial coordinates (obsm['spatial']) - not a "
                         "map-able section.")
    new.uns.setdefault("platform", getattr(sample, "platform", "spatial"))
    new.uns.setdefault("n_panel_genes", int(len(getattr(sample, "panel_genes", []) or new.var_names)))

    # The reference chooser keys off free-text tissue; a Xenium panel names its own organism/tissue,
    # so an unspecified tissue is inferred from the section rather than kept as the melanoma default.
    inferred = (tissue or "").strip() or _io.infer_tissue(new) or ctx.tissue
    entry = {"adata": new, "tissue": inferred, "path": p, "reference": None,
             "ref_label_key": None, "auto": None, "route": None}
    out = {"status": "loaded", "path": p, "platform": new.uns.get("platform"),
           "n_cells": int(new.n_obs), "n_panel_genes": int(new.uns.get("n_panel_genes", new.n_vars)),
           "panel": _io.panel_label(new.uns), "tissue": inferred}

    if auto_reference:
        # Score/select the reference on the NEW section via a sub-context (so the outer ctx.reference
        # is not clobbered until the app commits the swap). auto_select -> ranks the registry, loads
        # the best AVAILABLE tissue-matched atlas (or honestly reports 'no_reference'); the strategy
        # ladder then recommends supervised transfer vs de-novo clustering.
        sub = RunContext(tissue=inferred, use_llm=ctx.use_llm)
        try:
            auto = _cap_auto_select_reference(new, sub, tissue=inferred, allow_fetch=bool(allow_fetch))
        except Exception as exc:            # never let reference selection break the load
            auto = {"status": "error", "message": f"{type(exc).__name__}: {exc}", "ranked": []}
        entry["reference"], entry["ref_label_key"], entry["auto"] = sub.reference, sub.ref_label_key, auto
        out["reference"] = {"chosen": auto.get("chosen"), "status": auto.get("status"),
                            "source": auto.get("source"), "message": auto.get("message"),
                            "match": auto.get("match"), "ranked": (auto.get("ranked") or [])[:3]}
        try:
            route = _cap_annotation_strategy(new, sub)
        except Exception as exc:
            route = {"recommended_mode": "annotate", "reason": f"strategy failed: {exc}"}
        entry["route"] = route
        out["recommended_mode"] = route.get("recommended_mode")
        out["strategy_reason"] = route.get("reason")

    ctx.loaded.append(entry)
    return out


def _cap_compute_qc(adata, ctx, **_):
    from . import qc as _qc
    _qc.compute_qc(adata, progress=ctx.tick)
    return _qc.qc_summary(adata)


def _cap_qc_funnel(adata, ctx, **_):
    from . import markers as _m
    from . import qc as _qc
    key = "cell_type" if Obs.CELL_TYPE.present(adata) else (
        "leiden" if Obs.LEIDEN.present(adata) else "cell_type")
    # CRISP impurity (Layer 5) is a mutual-exclusivity check: pass COARSE disjoint lineages, not the fine
    # vocabulary, or a fine reference set (30+ overlapping subtypes) flags most cells as mixed-lineage
    # spillover and abstains the section. See annotate.apply_confidence.
    # Thread the reference for AQI's panel/depth CEILING ONLY when it is already loaded on the session
    # (resolve_reference returns it with no I/O then); never trigger a multi-GB atlas load in the
    # interactive funnel just for a ceiling term that rarely bites - the index is C+M driven.
    ref = ref_key = None
    if getattr(ctx, "reference", None) is not None:
        ref, ref_key, _ = ctx.resolve_reference()
    return _qc.run_funnel(adata, cluster_key=key, marker_sets=ctx.markers(adata),
                          panel_check_result=adata.uns.get("panel_check"), progress=ctx.tick,
                          lineage_markers=_m.for_tissue(ctx.tissue),
                          reference=ref, ref_key=(ref_key or "cell_type"))


def _cap_cluster(adata, ctx, resolution: float = 0.5, **_):
    from . import cluster as _cl
    _cl.cluster(adata, resolution=float(resolution), progress=ctx.tick)
    # Park the resolution on the object, not just in the return value: the UI reads it back from
    # /summary to keep its slider on the value that actually produced the clusters on screen.
    adata.uns["leiden_resolution"] = float(resolution)
    n = int(adata.obs["leiden"].nunique()) if "leiden" in adata.obs else 0
    return {"n_clusters": n, "resolution": float(resolution)}


def _cap_cluster_markers(adata, ctx, n_genes: int = 2, **_):
    # Top-N differentially-expressed genes per Leiden cluster, drawn as a dot-plot.
    #
    # This is normally a READ, not a compute: `cluster.cluster` already ran rank_genes_groups
    # (t-test, routed through backend.get_backend() so it is the GPU kernel when RAPIDS is live) and
    # parked the result on uns['rank_genes_groups']. So the button is ~instant on a 100k-cell section.
    # It only recomputes when that table is missing or STALE - i.e. its group names no longer match the
    # current leiden categories, which happens after a re-cluster at a different resolution.
    from . import cluster as _cl
    from . import views
    from .backend import get_backend, gpu_available

    cats = [str(c) for c in adata.obs["leiden"].cat.categories]
    rg = adata.uns.get("rank_genes_groups")
    have = set(rg["names"].dtype.names) if rg is not None and "names" in rg else set()
    stale = have != set(cats)
    if stale:
        ctx.tick(0.2, "ranking cluster markers")
        get_backend().rank_genes_groups(adata, groupby="leiden", method="t-test")
        rg = adata.uns["rank_genes_groups"]

    n = max(1, int(n_genes))
    per_cluster = {c: _cl.top_markers(adata, c, n=n) for c in cats}
    # Order the columns cluster-by-cluster (so a cluster's own genes sit together), de-duplicated:
    # a gene that tops two clusters is drawn once, at its first appearance.
    genes: list[str] = []
    for c in cats:
        for g in per_cluster[c]:
            if g not in genes:
                genes.append(g)
    if not genes:
        return {"rendered_view": False, "note": "no ranked genes for the current clusters"}

    ctx.tick(0.8, "drawing the dot-plot")
    # Size the figure to its content; the PNG renderer honours these. Clamped so a 60-gene x
    # 30-cluster plot stays legible without producing a multi-megabyte image.
    width = int(min(2200, max(760, 26 * len(genes) + 240)))
    height = int(min(1400, max(460, 22 * len(cats) + 170)))
    ctx.artifacts.append({
        "kind": "figure", "engine": "plotly", "title": f"top {n} DEGs per cluster",
        "fig": views.dotplot_view(adata, genes, groupby="leiden", group_order=cats,
                                  group_noun="clusters", height=height, width=width,
                                  title_suffix=f" · top {n} per cluster")})
    ctx.tick(1.0, "done")
    return {"rendered_view": True, "n_genes": n, "n_clusters": len(cats),
            "genes": genes, "per_cluster": per_cluster,
            "source": "recomputed (t-test)" if stale else "reused the ranking computed during clustering",
            "device": "GPU" if gpu_available() else "CPU"}


def _cap_cluster_confidence(adata, ctx, **_):
    from . import cluster_confidence as _cc
    return _cc.cluster_confidence(
        adata, cluster_key="leiden",
        annotation_key="cell_type" if Obs.CELL_TYPE.present(adata) else "")


def _cap_annotate(adata, ctx, reliability_weighted: bool = False, truth_key: str | None = None, **_):
    from . import annotate as _an
    from . import consensus as _cons
    from . import markers as _m

    _method_cols = ["rctd_first_type", "singler_label", "scanvi_label", "ph_fine"]
    present = [c for c in _method_cols if c in adata.obs]

    # Per-cell cross-method reconciliation is OPT-IN. Naive majority voting loses to the single best
    # method when one method dominates (annotator_bench: majority 0.858 vs RCTD alone 0.923), so when the
    # caller asks for it we weight each method's vote by its reliability - measured against `truth_key`
    # when a labeled column exists, else the documented DEFAULT_RELIABILITY prior. Off by default, which
    # leaves the cluster-only label exactly as it was.
    weights, weight_source = None, None
    if reliability_weighted and present:
        measured = _cons.reliability_from_labels(adata, present, truth_key) if truth_key else {}
        if measured:
            # The cluster label has no measured accuracy (it is created below), so it keeps its prior.
            weights = {"cell_type": _cons.DEFAULT_RELIABILITY["cell_type"], **measured}
            weight_source = "measured"
        else:
            weights = _cons.DEFAULT_RELIABILITY
            weight_source = "prior"

    df = _an.consensus_annotate(adata, cluster_key="leiden", use_llm=ctx.use_llm,
                                context=ctx.label_context(), marker_sets=ctx.markers(adata),
                                method_label_cols=present if weights else None,
                                reliability_weights=weights,
                                progress=ctx.tick)
    # popV: score cross-method consensus over any per-cell reference-method labels present BEFORE the
    # confidence call, so apply_confidence can trust it. Fires only with a diverse ensemble (a no-op on
    # the reference-free 2-method demo, which is correct - popV's calibration needs >= 3 diverse methods).
    _cons.consensus_metrics(adata, present)
    head = _an.apply_confidence(adata, cluster_key="cell_type", marker_sets=ctx.markers(adata),
                                panel_check_result=adata.uns.get("panel_check"),
                                lineage_markers=_m.for_tissue(ctx.tissue))
    return {"annotatability": head, "clusters": df.to_dict("records"),
            "consensus_methods": present, "consensus_weighting": weight_source}


def _cap_immune_exclusion(adata, ctx, type_a: str, type_b: str, **_):
    from . import annotate as _an
    from . import spatial as _sp
    # Same column as _cap_nhood: the pair's z-score here must equal the same pair's cell in the
    # enrichment heatmap shown beside it.
    return _sp.immune_exclusion(adata, type_a, type_b, cluster_key=_an.annotation_key(adata))


# Every cell-type consumer groups on the SAME column (annotate.annotation_key), so a delivered figure
# and the statistics printed beside it can never again disagree about the denominator.
def _cap_nhood(adata, ctx, **_):
    from . import annotate as _an
    from . import spatial as _sp
    return _sp.nhood_enrichment(adata, cluster_key=_an.annotation_key(adata))


def _cap_niches(adata, ctx, **_):
    from . import annotate as _an
    from . import niches as _n
    rows = _n.call_niches(adata, cluster_key=_an.annotation_key(adata),
                          progress=ctx.tick).to_dict("records")
    if "niche" in adata.obs:      # recolour the main canvas by niche (chat parity with the Spatial tab)
        ctx.artifacts.append({"kind": "map_view", "color_by": "niche"})
    return rows


def _cap_co_occurrence(adata, ctx, **_):
    from . import annotate as _an
    from . import spatial as _sp
    return _sp.co_occurrence(adata, cluster_key=_an.annotation_key(adata))


def _cap_spatial_genes(adata, ctx, n_top: int = 20, **_):
    from . import spatial as _sp
    return _sp.spatially_variable_genes(adata, n_top=n_top)


def _cap_states(adata, ctx, **_):
    from . import annotate as _an
    from . import states as _st
    mat = _st.state_by_celltype(adata, cluster_key=_an.annotation_key(adata))
    # score_fields = {state_name: obs_column} (== uns['state_columns']); lets the frontend recolour the
    # canvas by a program. Only programs with >=2 on-panel genes are scored, so it matches the matrix cols.
    fields = dict(adata.uns.get("state_columns", {}))
    if mat.empty:
        return {"cell_types": [], "states": [], "matrix": [], "score_fields": fields}
    return {"cell_types": list(map(str, mat.index)), "states": list(map(str, mat.columns)),
            "matrix": mat.values.tolist(), "score_fields": fields}


def _cap_assign_states(adata, ctx, min_z: float = 1.0, **_):
    # Type each cell with its DOMINANT cell-state program (CyteType-style): writes obs['cell_state']
    # (cycling / IFN-ISG / hypoxia / stress-HSP / EMT / T-exhaustion / T-cytotoxicity / ECM-remodeling /
    # antigen-presentation / None) and returns the state distribution + the dominant state per cell
    # type. When an LLM is available, also names each cell
    # type's state in plain language (grounded in the program z-scores + top markers).
    from . import annotate as _an
    from . import states as _st
    key = _an.annotation_key(adata)     # one column, same as the state heatmap (_cap_states)
    res = _st.assign_cell_states(adata, cluster_key=key, min_z=float(min_z))
    if ctx.use_llm and res.get("states"):
        labels = _st.name_states_llm(adata, cluster_key=key, context=f"human {ctx.tissue}")
        if labels:
            res["llm_labels"] = labels
    if "cell_state" in adata.obs:   # recolour by cell_state (chat parity with the Annotate wizard)
        ctx.artifacts.append({"kind": "map_view", "color_by": "cell_state"})
    return res


def _cap_malignant(adata, ctx, **_):
    import numpy as np
    from . import cnv as _cnv
    info = _cnv.malignant_score(adata, tissue=ctx.tissue)
    s = np.asarray(adata.obs.get("malignant_score", []))
    if "malignant_score" in adata.obs:   # "where is the tumour?" -> colour the map by the score
        ctx.artifacts.append({"kind": "map_view", "color_by": "malignant_score"})
    return {"n_markers": info.get("n_markers"),
            "pct_high_malignant": float((s > 0.6).mean()) if s.size else 0.0}


def _cap_programs(adata, ctx, n_programs: int = 8, **_):
    from . import programs as _p
    return _p.discover_programs(adata, n_programs=int(n_programs)).to_dict("records")


def _cap_name_programs(adata, ctx, **_):
    # LLM-name each de-novo NMF program from its top loading genes (grounded; degrades to "Program k"
    # with no key). Relabels obs['program'] from the STABLE dominant-program index so the map legend
    # reads by biological program, WITHOUT renaming any obsm/varm column or the program_score_<k>
    # colour field. Self-ensures discover_programs when the decomposition is absent.
    import numpy as np
    import pandas as pd
    from . import llm as _llm
    from . import programs as _p

    if not Varm.PROGRAM_LOADINGS.present(adata):
        _p.discover_programs(adata)
    load = np.asarray(adata.varm["program_loadings"])          # genes x programs
    W = np.asarray(adata.obsm["programs"])                     # cells x programs
    names = np.asarray(adata.var_names)
    n_prog = load.shape[1]
    dom = W.argmax(1)
    top = {k: [str(g) for g in names[np.argsort(-load[:, k])[:12]]] for k in range(n_prog)}
    scores = {r["program"]: r for r in _p.score_programs(adata)}   # index-based; safe pre-relabel

    named: dict = {}
    if ctx.use_llm:
        named = _llm.name_programs({f"Program {k}": top[k] for k in range(n_prog)},
                                   context=ctx.label_context())

    labelmap: dict[int, str] = {}
    for k in range(n_prog):
        lab = ((named.get(f"Program {k}") or {}).get("label"))
        if isinstance(lab, str) and lab.strip():
            labelmap[k] = lab.strip()
    disp = {k: labelmap.get(k, f"Program {k}") for k in range(n_prog)}
    adata.obs["program"] = pd.Categorical([disp[i] for i in dom])

    rows = []
    for k in range(n_prog):
        sc = scores.get(f"Program {k}", {})
        info = named.get(f"Program {k}") or {}
        rows.append({"program": disp[k], "program_id": f"Program {k}", "program_index": k,
                     "score_field": f"program_score_{k}",
                     "n_cells": int((dom == k).sum()), "top_genes": top[k],
                     "plaid_in": sc.get("plaid_in"), "plaid_out": sc.get("plaid_out"),
                     "specificity": sc.get("specificity"),
                     "label": labelmap.get(k), "confidence": info.get("confidence"),
                     "rationale": info.get("rationale")})
    if "program" in adata.obs:   # legend now reads by biological program (chat parity with wizard)
        ctx.artifacts.append({"kind": "map_view", "color_by": "program"})
    return rows


def _cap_subcluster(adata, ctx, cell_type: str, resolution: float = 0.1, **_):
    from . import subcluster as _sc
    _, rows = _sc.subcluster(adata, cell_type, resolution=float(resolution),
                             use_llm=ctx.use_llm, context=ctx.label_context())
    # Recolour the main canvas by the new subtypes (like the view tools), so "subcluster the T cells"
    # SHOWS the result on the map, not just in the chat. subcluster wrote obs['subtype']; it is in
    # _obs_fields, so the SPA's colour dropdown already lists it and switches to it on this directive.
    if "subtype" in adata.obs:
        ctx.artifacts.append({"kind": "map_view", "color_by": "subtype"})
    return {"rows": rows}


def _cap_annotation_methods(adata, ctx, **_):
    # Pure report: inspect whichever reference-method label columns are present (written by the
    # subprocesses/annotation runners) and summarize coverage + popV-style cross-method consensus.
    # Never assigns labels.
    import numpy as np
    from . import consensus as _cons
    label_cols = {"rctd": "rctd_first_type", "singler": "singler_label",
                  "scanvi": "scanvi_label", "panhumanpy": "ph_fine"}
    methods_run, coverage = [], {}
    for m, col in label_cols.items():
        if col in adata.obs:
            methods_run.append(m)
            coverage[m] = float(adata.obs[col].notna().mean())
    # popV consensus (Nat Genet 2024): cross-method AGREEMENT tracks accuracy far better than any single
    # method's self-score. Score it over the per-cell method labels present (vs the final cell_type).
    voters = [c for c in label_cols.values() if c in adata.obs]
    cm = _cons.consensus_metrics(adata, voters)
    if cm.get("status") == "ok":
        consensus = {
            "mean_agreement": round(cm["mean_agreement"], 3),
            "reliability": cm["reliability_counts"],
            "pct_trusted_ensemble": round(cm["pct_trusted"], 3),
            "note": (f"popV: cross-method agreement is the trusted uncertainty signal; only weighted "
                     f"where >= {_cons.MIN_TRUST_METHODS} methods voted (a diverse ensemble)."),
        }
    elif "consensus_agreement" in adata.obs:
        ag = np.asarray(adata.obs["consensus_agreement"], dtype=float)
        valid = ~np.isnan(ag)
        consensus = {"mean_agreement": float(ag[valid].mean()) if valid.any() else None}
    else:
        consensus = None
    return {"methods_run": methods_run, "coverage": coverage, "consensus": consensus}


def _cap_highlight(adata, ctx, criterion: str, **_):
    # Resolve the criterion to a cell mask, render a spatial view, and stash the figure on
    # ctx.artifacts so the app displays it. Returns the count + an honest note (never invents).
    import numpy as np
    from . import views
    mask, resolved, note = views.select_cells(adata, criterion, ctx)
    mask = np.asarray(mask, dtype=bool)
    n = int(mask.sum())
    done = {"criterion": criterion, "resolved": resolved, "n_matched": n,
            "n_total": int(adata.n_obs), "note": note}
    if not n:                                    # nothing matched: an honest note, no view to draw
        return {**done, "rendered_view": False}
    # If the criterion resolved to a cell-type CATEGORY ("where are the T cells"), drive the MAIN
    # specimen canvas: colour by cell_type and emphasise that category (the SPA fades the rest), rather
    # than dumping a thumbnail into the chat. The user asked the map a question; light up the map.
    cats = set(map(str, adata.obs["cell_type"].astype(str).unique())) if "cell_type" in adata.obs else set()
    if resolved in cats:
        ctx.artifacts.append({"kind": "map_view", "color_by": "cell_type", "highlight": resolved, "note": note})
        return {**done, "rendered_view": True, "drove_map": True}
    # A cross-cutting SUBSET ("low confidence", "low quality", "abstained", "malignant") - NOT a
    # cell_type, so there is no single category colour to emphasise. Collapse the mask into a transient
    # 2-class obs column and drive the MAIN canvas by colouring on it + emphasising the matched subset
    # (fade the rest) - same effect as a cell-type highlight, on the real plot instead of a chat
    # thumbnail. The column is not in _obs_fields, so it never clutters the colour dropdown, but
    # points() still colours by any obs column. Overwritten by the next highlight.
    import pandas as pd
    adata.obs["_copilot_highlight"] = pd.Categorical(
        np.where(mask, resolved, "other"), categories=[resolved, "other"])
    ctx.artifacts.append({"kind": "map_view", "color_by": "_copilot_highlight",
                          "highlight": resolved, "note": note})
    return {**done, "rendered_view": True, "drove_map": True}


def _cap_show_spatial(adata, ctx, color_by: str, **_):
    # Colour the MAIN specimen canvas by an obs field OR an on-panel gene - a gene's expression is
    # written to a persistent obs column (color_field_for) so the big canvas can colour by it. Emits a
    # 'map_view' directive (the app recolours the main spatial plot) rather than an inline figure, for
    # BOTH obs fields and genes. Grounded: an unknown key returns a note and drives no view.
    from . import views
    field, note = views.color_field_for(adata, color_by, ctx)
    if field is not None:
        ctx.artifacts.append({"kind": "map_view", "color_by": field, "note": note})
    return {"color_by": color_by, "resolved": field or color_by,
            "rendered_view": field is not None, "note": note}


def _cap_marker_dotplot(adata, ctx, cell_type: str | None = None, genes: list | None = None,
                        corrected: bool = False, **_):
    # Dot-plot of markers x cell types. genes: explicit (on-panel only) else the cell_type's canonical
    # markers else the grouped marker OVERVIEW (rows = cell types grouped so related types are adjacent,
    # cols = top-2 markers per category, ranked by expression - never LLM-chosen). Never plots off-panel
    # genes. With corrected=True and a SPLIT-purified layer present, emits a raw + purified pair.
    from . import views
    group_order = None
    groups_map: dict | None = None
    if not genes and not cell_type:                 # the no-argument overview: grouped rows + top-2/cat
        gl, group_order, note = views.category_overview(adata, ctx)
        groups_map = dict(adata.uns.get("marker_overview_groups", {}))
    else:
        gl, note = views._resolve_genes(adata, genes, cell_type, ctx)
    if not gl:
        return {"rendered_view": False, "genes": [], "note": note or "no on-panel genes to plot"}
    has_split = "split_corrected" in getattr(adata, "layers", {})
    # Once SPLIT has run, ANY dot-plot (including the auto-drawn overview, which calls this with
    # corrected unset) shows the raw + purified PAIR - driven by the layer's PRESENCE, not the
    # `corrected` flag. Otherwise "after SPLIT" silently kept rendering raw X and the correction was
    # invisible unless the user knew to re-request it.
    if has_split:
        ctx.artifacts.append({"kind": "figure", "engine": "plotly", "title": "marker dot-plot (raw)",
                              "fig": views.dotplot_view(adata, gl, groupby="cell_type",
                                                        group_order=group_order, title_suffix=" (raw)")})
        ctx.artifacts.append({"kind": "figure", "engine": "plotly", "title": "marker dot-plot (SPLIT-purified)",
                              "fig": views.dotplot_view(adata, gl, groupby="cell_type", group_order=group_order,
                                                        layer="split_corrected", title_suffix=" (SPLIT-purified)")})
        out = {"rendered_view": True, "genes": gl, "note": note, "compared": True}
    else:
        ctx.artifacts.append({"kind": "figure", "engine": "plotly", "title": "marker dot-plot",
                              "fig": views.dotplot_view(adata, gl, groupby="cell_type", group_order=group_order)})
        out = {"rendered_view": True, "genes": gl, "note": note, "split_available": has_split}
    if groups_map is not None:
        out["groups"] = groups_map
    return out


def _cap_annotate_rctd(adata, ctx, max_cells: int = 0, set_as_primary: bool = False, **_):
    # Run RCTD doublet-mode (rctd-py subprocess, GPU-optional) against the session reference and write
    # obs['rctd_first_type'/'second_type'/'spot_class'/'weight'/'singlet_score']. RCTD was the most
    # accurate annotator in the non-circular Atera recheck (0.85 lineage vs curated), which is why
    # consensus weights it highest; this makes its labels actually reach the consensus vote instead of
    # only being scored if some external process happened to write them. Needs SPATIALSCRIBE_RCTD_PYTHON
    # and an uploaded reference; degrades to 'skipped' otherwise. set_as_primary adopts rctd_first_type
    # as obs['cell_type'].
    import os

    from . import methods as _m
    from .io import SpatialSample

    ref, key, src = ctx.resolve_reference()
    env_python = os.environ.get("SPATIALSCRIBE_RCTD_PYTHON", "")
    if ref is None or key is None:
        return {"status": "skipped: no reference; upload one at the Panel-check step", "reference_source": src}
    if not env_python or not os.path.exists(env_python):
        return {"status": "skipped: SPATIALSCRIBE_RCTD_PYTHON not configured", "reference_source": src}

    import tempfile
    from pathlib import Path

    import numpy as np
    tmp = Path(tempfile.mkdtemp(prefix="ssrctd_"))
    ref_h5 = tmp / "reference.h5ad"
    ref.write_h5ad(ref_h5)
    # run_rctd only reads sample.adata; the other fields just satisfy the dataclass contract.
    sample = SpatialSample(platform=adata.uns.get("platform", "spatial"), adata=adata,
                           control_mask=np.zeros(adata.n_vars, dtype=bool),
                           panel_genes=list(adata.var_names), has_z=False)
    pq = _m.run_rctd(sample, reference_path=str(ref_h5), ref_label_key=key, env_python=env_python)
    if pq is None:
        return {"status": "skipped: RCTD subprocess produced no output", "reference_source": src}
    info = _m.join_rctd(adata, pq)
    coverage = float(adata.obs["rctd_first_type"].notna().mean()) if "rctd_first_type" in adata.obs else 0.0
    if set_as_primary and "rctd_first_type" in adata.obs:
        adata.obs["cell_type"] = adata.obs["rctd_first_type"].astype("category")
    if "rctd_first_type" in adata.obs:   # colour by RCTD's labels (or cell_type when adopted)
        ctx.artifacts.append({"kind": "map_view",
                              "color_by": "cell_type" if set_as_primary else "rctd_first_type"})
    return {"status": "ok", "method": "rctd", "coverage": round(coverage, 3),
            "set_as_primary": bool(set_as_primary), "reference_source": src, **(info or {})}


def _cap_split_purify(adata, ctx, max_cells: int = 0, weights_engine: str = "tacco", **_):
    # SPLIT reference-path spillover purification (residual-contamination removal): decontaminate
    # transcript spillover using the uploaded reference, writing layers['split_corrected']. The
    # deconvolution weights come from TACCO in-env by default (fast, no rctd-py env); pass
    # weights_engine='rctd' to use the rctd-py subprocess instead. SPLIT::rctd_free_purify is
    # annotation-agnostic, so the R side is identical either way. Degrades to 'skipped' when the
    # SPLIT-R env / reference are absent. Compare markers before/after with marker_dotplot(corrected=true).
    from . import split as _sp
    ref, key, _src = ctx.resolve_reference()
    return _sp.split_purify(adata, reference=ref, ref_label_key=key,
                            marker_sets=ctx.markers(adata), max_cells=int(max_cells),
                            weights_engine=str(weights_engine))


def _cap_expression_violin(adata, ctx, gene: str, **_):
    # Violin of a gene's expression across cell types. Grounded: off-panel gene -> note, no plot.
    from . import views
    fig, note = views.violin_view(adata, gene, groupby="cell_type")
    if fig is not None:
        ctx.artifacts.append({"kind": "figure", "engine": "plotly", "title": f"{gene} violin", "fig": fig})
    return {"gene": gene, "rendered_view": fig is not None, "note": note}


def _cap_composition_chart(adata, ctx, **_):
    from . import annotate as _an
    from . import export as _ex
    from . import views
    key = _an.annotation_key(adata)
    # The rendered artifact and the returned numbers MUST share a denominator - they did not, so the
    # chart showed "Unassigned" as a bar while the JSON silently excluded it.
    ctx.artifacts.append({"kind": "figure", "engine": "plotly", "title": "composition",
                          "fig": views.composition_view(adata, groupby=key)})
    vc, abstained = _ex.composition_table(adata, key)   # typed cells only: an abstention is not a type
    return {"rendered_view": True, "n_types": int(vc.size),
            "pct_abstained": round(100 * abstained, 1),
            "composition": {str(k): int(v) for k, v in vc.head(12).items()}}


def _cap_show_segmentation(adata, ctx, color_by: str = "cell_type", **_):
    # Draw the actual cell-segmentation polygons (Xenium cell_boundaries), colored by a field.
    # Needs a run directory with cell_boundaries.parquet (io.load records its path in uns).
    from . import views
    bpath = adata.uns.get("boundaries_path")
    if not bpath:
        return {"rendered_view": False, "note": "No cell-segmentation polygons for this section. Load a "
                "Xenium run directory (with cell_boundaries.parquet) - the processed demo ships centroids only."}
    bnd = views.load_boundaries(bpath)
    if bnd is None:
        return {"rendered_view": False, "note": f"Could not read segmentation polygons from {bpath}."}
    cb = color_by if color_by in adata.obs else ("cell_type" if "cell_type" in adata.obs else color_by)
    fig, n, note = views.segmentation_figure(adata, bnd, color_by=cb)
    if fig is not None:
        ctx.artifacts.append({"kind": "figure", "engine": "plotly", "title": "cell segmentation", "fig": fig})
    return {"rendered_view": fig is not None, "n_shown": n, "color_by": cb, "note": note}


def _cap_calibrate_confidence(adata, ctx, truth_key: str | None = None, pred_key: str = "cell_type",
                              cal_frac: float = 0.5, n_bins: int = 10, **_):
    # Post-hoc ISOTONIC calibration of the per-cell annotation confidence, fit against REAL labels. The
    # heuristic score is mis-calibrated (ECE 0.33 on a benchmarked 12k CosMx section) and non-monotonic in
    # accuracy, so the number does not mean P(correct). Fitting a monotone map from a labeled set repairs
    # the meaning. Without labels there is nothing to fit and we skip - never claim calibration we did not
    # earn. ponytail: no reference-transfer path (fit on a labeled reference, apply to the section) until
    # someone loads a labeled reference and asks; it would inherit covariate shift and need its own caveat.
    import numpy as np

    from . import calibration as _cal

    if not truth_key or truth_key not in adata.obs:
        return {"status": "skipped: no ground-truth column supplied",
                "source": None, "obs_written": None,
                "note": "Confidence stays an ordinal heuristic, not a probability. Pass truth_key=<obs "
                        "column of real labels on the SAME axis as pred_key> to fit a calibration; "
                        "calibration is only claimed when it is fit against real labels."}
    if pred_key not in adata.obs:
        return {"status": f"skipped: pred_key '{pred_key}' not in obs", "source": None,
                "obs_written": None, "note": "Run annotate first."}

    conf = np.nan_to_num(np.asarray(adata.obs[_cal.RAW_KEY], dtype=float), nan=0.0)
    correct = (adata.obs[pred_key].astype(str).to_numpy()
               == adata.obs[truth_key].astype(str).to_numpy()).astype(int)
    if correct.sum() in (0, correct.size):
        return {"status": "skipped: predictions are all-correct or all-wrong, nothing to calibrate against",
                "source": None, "obs_written": None,
                "note": f"accuracy={float(correct.mean()):.3f} against '{truth_key}' - check that the "
                        f"truth column and '{pred_key}' share a label axis."}

    # Fit on one split, report ECE on the complement: an in-sample ece_after is optimistic by construction.
    n = int(conf.size)
    idx = np.random.default_rng(0).permutation(n)
    cut = int(np.clip(round(n * float(cal_frac)), 1, n - 1)) if n > 1 else n
    fit_i, eval_i = idx[:cut], idx[cut:]
    held_out = eval_i.size > 0
    if not held_out:
        eval_i = fit_i

    cal = _cal.fit_isotonic(conf[fit_i], correct[fit_i])
    if cal is None:
        return {"status": "skipped: confidence has fewer than 2 distinct values, nothing to fit",
                "source": None, "obs_written": None, "note": ""}

    # Null model = predict the FIT half's accuracy for every cell. It is calibrated by construction, so
    # ece_after only beats it by being no worse; brier_skill is what shows real per-cell information.
    # Group by the PREDICTED label so auc_within reports discrimination free of the across-cell-type
    # Simpson confound that makes pooled AUC invert (see calibration.auc_within).
    rep = _cal.report(conf[eval_i], correct[eval_i], cal, n_bins=int(n_bins),
                      baserate=float(correct[fit_i].mean()),
                      groups=adata.obs[pred_key].astype(str).to_numpy()[eval_i])
    _cal.apply(adata, cal)                       # writes the NEW column; raw stays untouched
    adata.uns["calibration_report"] = rep

    # The honest headline: a constant "every cell is right with p = the section's accuracy" is already
    # calibrated, so a near-zero ece_after is cheap. brier_skill says whether we learned anything per-cell.
    sk = rep["brier_skill"]
    if sk < 0.01:
        skill_note = (f"But Brier skill over the base-rate null model is only {sk:+.3f}: calibration here "
                      f"reduces to reporting the section's accuracy ({rep['baserate']:.3f}) for every cell. "
                      f"The score carries no usable per-cell information.")
    else:
        skill_note = (f"Brier skill over the base-rate null model is {sk:+.3f}, so the calibrated score does "
                      f"carry per-cell information beyond the section's accuracy ({rep['baserate']:.3f}).")

    # Pooled AUC is confounded across cell types (Simpson): report the within-type term, which is what a
    # per-cell gate can actually exploit. On real sections pooled inverts to 0.42 while within is 0.48.
    a = rep["auc_after"]
    w = rep.get("auc_within")
    if a < 0.45:
        rank_note = (f"Pooled AUC {a:.2f} < 0.5: the gate is INVERTED - confidence ranks correct cells BELOW "
                     f"incorrect ones, so abstaining on low confidence discards the cells most likely right.")
    elif a < 0.55:
        rank_note = (f"Pooled AUC {a:.2f} ~ 0.5: confidence barely orders cells by correctness, so no "
                     f"abstention threshold buys much accuracy at any calibration.")
    else:
        rank_note = f"Pooled AUC {a:.2f} > 0.5: confidence ranks correct cells higher, so abstention can help."
    if w is not None:
        rank_note += (f" Within-cell-type AUC is {w:.2f}"
                      + (" - the pooled figure is a Simpson artifact of the pipeline being confident on the "
                         "cell types it labels worst, not a per-cell signal."
                         if abs(w - 0.5) < 0.05 < abs(a - 0.5) else
                         " (the term a per-cell gate can actually exploit)."))

    # Does a gate actually earn its place? Pick the threshold on the FIT half and measure it on the
    # held-out half - picking and scoring on the same cells would report a cherry-picked gain.
    gate = None
    fit_scores = np.clip(cal.predict(conf[fit_i]), 0.0, 1.0)
    candidates = [d for d in _cal.abstention_curve(fit_scores, correct[fit_i]) if d["useful"]]
    if candidates:
        t = max(candidates, key=lambda d: d["gain"])["threshold"]
        ev_scores = np.clip(cal.predict(conf[eval_i]), 0.0, 1.0)
        keep = ev_scores >= t
        base = float(correct[eval_i].mean())
        if keep.any():
            kept_acc = float(correct[eval_i][keep].mean())
            gate = {"threshold": t, "keep_frac": float(keep.mean()), "accuracy": kept_acc,
                    "gain": kept_acc - base, "selected_on": "fit half", "measured_on": "held-out half"}
    if gate and gate["gain"] > 0:
        gate_note = (f"Best gate (threshold chosen on the fit half, scored on the held-out half): keep "
                     f"calibrated >= {gate['threshold']} -> retains {gate['keep_frac']:.1%} of cells at "
                     f"accuracy {gate['accuracy']:.3f} ({gate['gain']:+.3f} vs keeping all).")
    else:
        gate_note = ("No abstention threshold both keeps >=10% of cells AND beats keeping every cell "
                     "out-of-sample: on this section the gate should be dropped, not tuned.")

    out = {"status": "ok", "source": "truth_column", "n": n,
           "n_fit": int(fit_i.size), "n_eval": int(eval_i.size),
           "obs_written": _cal.CALIBRATED_KEY, "gate": gate,
           "note": (f"Isotonic fit on {fit_i.size} cells, scored on "
                    f"{'a held-out ' + str(eval_i.size) if held_out else 'the same (in-sample, too few cells to split)'}"
                    f" cells. Calibration fixes what the score MEANS (ECE {rep['ece_before']:.3f} -> "
                    f"{rep['ece_after']:.3f}); it never reorders cells, so it cannot conjure discrimination "
                    f"that the raw score lacks. {skill_note} {rank_note} {gate_note} "
                    f"The raw '{_cal.RAW_KEY}' is unchanged.")}
    out.update({k: v for k, v in rep.items() if k != "n"})
    return out


def _cap_rejection_reasons(adata, ctx, **_):
    from . import rejection
    pc = adata.uns.get("panel_check")
    rejection.assign_rejection_reasons(adata, panel_check_result=pc)
    return {
        "breakdown": rejection.rejection_breakdown(adata, panel_check_result=pc).to_dict("records"),
        "panel_warnings": rejection.panel_resolvability_warnings(adata, pc),
    }


def _cap_self_verify(adata, ctx, neighborhood: bool = False, **_):
    from . import verify as _v
    res = _v.verify_annotation(adata, marker_sets=ctx.markers(adata), cluster_key="cell_type",
                               neighborhood=bool(neighborhood))
    res["suggestions"] = _v.suggest_reruns(res)   # advisory dry-run
    return res


def _cap_trust_ledger(adata, ctx, **_):
    # Three independent per-type verdicts (resolvable x coherent x agreed) + their disagreements - the
    # coherent-but-disputed mislabel the AQI index alone cannot see. Reuses typability + verify + the
    # per-cell consensus_agreement; advisory, never changes labels.
    from . import trust as _t
    return _t.trust_ledger(adata, marker_sets=ctx.markers(adata), cluster_key="cell_type",
                           reference_match=adata.uns.get("reference_match"))


def _cap_self_heal(adata, ctx, max_rounds: int = 2, merge_confusable: bool = False, **_):
    # Self-verify then RE-RUN what fails: the executor half of the "verify, then re-run" vision.
    # verify_annotation flags types whose cells disagree with their markers; this loops, running the
    # safe auto-fix (subcluster a heterogeneous type) and re-verifying, until nothing improves or
    # max_rounds. With merge_confusable, it ALSO applies the panel-gap coarsening first (collapse the
    # types the panel cannot separate) and reports the quality delta; relabel stays advisory. Never raises.
    from . import annotate as _an
    from . import verify as _v
    # Heal the column the report/app treat as canonical (cell_type_final once apply_confidence ran),
    # so an abstained mislabel is actually visible downstream.
    return _v.autorerun(adata, ctx, max_rounds=int(max_rounds), cluster_key=_an.annotation_key(adata),
                        merge_confusable=bool(merge_confusable))


def _cap_merge_types(adata, ctx, groups: list | None = None, dry_run: bool = False, **_):
    # APPLY the panel-driven merges (unlike self_verify, which only advises): collapse cell types the
    # panel cannot separate (no private on-panel marker) into ONE coarser label in the canonical
    # annotation column, and report the annotation-quality delta (marker-agreement + n failing types,
    # before vs after). `groups` defaults to panel_check's merge_groups, so a bare "merge the
    # confusable types" call Just Works once panel_check has run. With dry_run=True it PREVIEWS the
    # ontology-partitioned merges (`would_merge`) without mutating labels - the Curate suggest button.
    # Never raises.
    from . import annotate as _an
    from . import panel_check as _pc
    from . import verify as _v

    key = _an.annotation_key(adata)
    ms = ctx.markers(adata)
    if not groups:
        pc = adata.uns.get("panel_check")
        if not _pc.is_valid(pc):
            pc = _pc.check_panel(_panel_genes(adata), marker_sets=ms)
            adata.uns["panel_check"] = pc
        groups = pc.get("merge_groups", [])
    if not groups:
        return {"status": "no_merges", "cluster_key": key, "merged": [], "would_merge": [],
                "n_groups_merged": 0, "n_groups": 0,
                "note": "The panel resolves every type - no confusable groups to merge."}
    return _v.merge_confusable_types(adata, groups, cluster_key=key, marker_sets=ms, dry_run=bool(dry_run))


def _cap_annotation_strategy(adata, ctx, reference_path=None, label_key=None, max_rounds: int = 3, **_):
    # Self-verify the reference<->panel match and auto-rerun a remediation ladder (coarsen the labels
    # the panel cannot separate -> reselect a better-overlap tissue-matched reference), then RECOMMEND
    # supervised label transfer vs unsupervised de-novo clustering. Uses the reference chosen at Panel-check
    # (or reference_path). With NO reference, marker-based supervised annotation is the default -
    # clustering is only indicated when a reference exists but the panel cannot resolve it. The route
    # NEVER falls to clustering on a single low metric (frac_resolvable is granularity-dependent).
    # Never raises; writes uns['annotation_route'] with the full ladder for the report/UI.
    import numpy as np

    from . import reference as _ref

    panel = _panel_genes(adata)
    ref, key, src = ctx.resolve_reference(reference_path, label_key)
    if ref is None or key is None:
        route = {"status": "no_reference", "recommended_mode": "annotate",
                 "reason": ("No single-cell reference provided - use marker-based supervised annotation "
                            "(panel_check governs which types the panel resolves). Provide a tissue-"
                            "matched reference to enable label transfer, or prefer clustering if the "
                            "panel cannot resolve the expected types."),
                 "reference_used": None, "coarsen_map": None, "initial_verdict": None,
                 "final_verdict": None, "final_global": {}, "final_per_type": {},
                 "ladder": [], "clustering_nudge": None}
    else:
        depth = (float(np.median(np.asarray(adata.obs["total_counts"], dtype=float)))
                 if "total_counts" in adata.obs else None)
        try:
            route = _ref.plan_annotation_strategy(ref, panel, key, target_depth=depth, tissue=ctx.tissue,
                                                  registry=_ref._load_registry(), max_rounds=int(max_rounds))
        except Exception as exc:   # contract: ALWAYS write produces=uns['annotation_route'], never raise
            route = {"status": "error", "recommended_mode": "annotate",
                     "reason": (f"Could not evaluate the reference-panel match ({exc}); falling back to "
                                "marker-based annotation."),
                     "reference_used": None, "coarsen_map": None, "initial_verdict": None,
                     "final_verdict": None, "final_global": {}, "final_per_type": {},
                     "ladder": [], "clustering_nudge": None}
    route["reference_source"] = src
    adata.uns["annotation_route"] = route
    return {"recommended_mode": route["recommended_mode"], "reason": route["reason"],
            "initial_verdict": route.get("initial_verdict"), "final_verdict": route.get("final_verdict"),
            "reference_used": route.get("reference_used"), "coarsen_map": route.get("coarsen_map"),
            "n_ladder_steps": len(route.get("ladder", []))}


def _cap_malignant_concordance(adata, ctx, max_cells: int = 25000, cf_threshold: float = 0.5, **_):
    # Run the malignant callers (marker score always; Cancer-Finder on any tumour panel; InSituCNV
    # on >2000-probe panels) and report per-caller %-malignant + the two learned callers' pairwise
    # concordance. Gated by ctx.is_tumour (the Data-step checkbox) when answered, else by the tissue
    # keyword heuristic. The learned callers degrade to 'skipped' when their isolated envs are absent
    # (SPATIALSCRIBE_CNV_* / CANCERFINDER_*). Writes obs['malignant_score'].
    from . import cnv as _cnv
    res = _cnv.call_malignant_concordance(adata, tissue=ctx.tissue, marker_sets=ctx.markers(adata),
                                          max_cells=int(max_cells), cf_threshold=float(cf_threshold),
                                          is_tumour=ctx.is_tumour)
    if "malignant_score" in adata.obs:   # colour the map by the malignant score
        ctx.artifacts.append({"kind": "map_view", "color_by": "malignant_score"})
    return res


def _cap_reference_transfer(adata, ctx, reference_path=None, label_key=None,
                            set_as_primary: bool = False, **_):
    # Reference-anchored annotation: transfer cell-type labels from the user's single-cell reference
    # onto the section. Primary arm = a CellTypist model trained on the reference (in-env, reliable);
    # also runs TACCO OT when the optional `tacco` package is installed. Writes obs['celltypist_label']
    # (+ 'ref_label' from TACCO) and reports coverage + agreement with the current cell_type. With
    # set_as_primary, the CellTypist labels also become obs['cell_type']. Degrades to 'no_reference'.
    import numpy as np
    from . import annotate as _an

    ref, key, src = ctx.resolve_reference(reference_path, label_key)
    if ref is None or key is None:
        return {"status": "no_reference", "reference_source": src, "arms": []}
    arms = []
    ct = _an.celltypist_transfer(adata, ref, key)
    arms.append({"method": "celltypist", **ct})
    # optional TACCO OT arm (skipped cleanly when `tacco` is not installed). Feed it the section's OWN
    # estimated composition (from the CellTypist arm that just ran) as the prior - TACCO's default
    # reference-derived prior is catastrophic on a depth-capped reference (0.35 vs 0.69 accuracy on the
    # Atera breast benchmark). Truth-free and ~zero extra cost. See reference_transfer.composition_prior.
    try:
        from . import reference_transfer as _rt
        prior = _rt.composition_prior(adata, ref, key)
        tinfo = _rt.transfer_labels(adata, ref, key, annotation_prior=prior)
        arms.append({"method": "tacco", "status": "ok", "prior": "celltypist" if prior is not None
                     else "reference-default", **tinfo})
    except Exception as exc:                       # ImportError (no tacco) or a transfer error
        arms.append({"method": "tacco", "status": f"skipped: {exc}"})

    if set_as_primary and "celltypist_label" in adata.obs:
        adata.obs["cell_type"] = adata.obs["celltypist_label"].astype("category")
        # label_conflict was computed for the marker-vs-LLM naming of the LEIDEN clusters. Those
        # labels are gone; a stale flag would make apply_confidence cap unrelated cells at WARN.
        adata.obs.drop(columns=["label_conflict"], errors="ignore", inplace=True)

    # agreement of the (reliable) celltypist labels with the current cell_type, if any
    agreement = None
    if "celltypist_label" in adata.obs and "cell_type" in adata.obs:
        a = adata.obs["celltypist_label"].astype(str).to_numpy()
        b = adata.obs["cell_type"].astype(str).to_numpy()
        agreement = float(np.mean(a == b)) if len(a) else None
    if "celltypist_label" in adata.obs:   # colour by the transferred labels (or cell_type when adopted)
        ctx.artifacts.append({"kind": "map_view",
                              "color_by": "cell_type" if set_as_primary else "celltypist_label"})
    return {"status": "ok", "reference_source": src, "label_key": key, "arms": arms,
            "agreement_with_cell_type": agreement, "set_as_primary": bool(set_as_primary)}


# --------------------------------------------------------------------------- #
# The registry
# --------------------------------------------------------------------------- #
_CAPABILITIES: list[Capability] = [
    Capability(
        name="load_section", label="Load a section from a path",
        description="Load a spatial section from a SERVER-SIDE path and make it the active section. "
                    "Use this when the user asks to open/load/read a dataset at a path (a Xenium / "
                    "CosMx / MERSCOPE output folder, or a .h5ad). Auto-detects the platform, INFERS "
                    "the tissue from the panel metadata (e.g. a Xenium Mouse Brain panel -> 'mouse "
                    "brain'), and - unless auto_reference is false - automatically SELECTS the best-"
                    "matched single-cell reference for that tissue and recommends supervised label "
                    "transfer vs de-novo clustering. Returns what it loaded and which reference it "
                    "chose. After this, the section is swapped in - subsequent steps run on it.",
        fn=_cap_load_section, required_params=("path",),
        params={"path": {"type": "string", "description": "server-side path: a Xenium/CosMx/MERSCOPE output folder, or a .h5ad"},
                "tissue": {"type": "string", "description": "tissue context override (optional; else inferred from the panel)"},
                "auto_reference": {"type": "boolean", "description": "auto-select the reference + annotation strategy after loading (default true)"},
                "allow_fetch": {"type": "boolean", "description": "allow a live CELLxGENE gget fetch when no local reference fits (default false)"}},
        copilot_exposed=True,
    ),
    Capability(
        name="describe_sample", label="Describe sample",
        description="Summary of the loaded section: n cells, cell-type composition, QC medians.",
        fn=_cap_describe, produces=(), copilot_exposed=True,   # pure report (QC is a side effect)
    ),
    Capability(
        name="panel_check", label="Panel check",
        description="Which cell types the panel can/cannot resolve (panel-adequacy check), "
                    "plus pairs it cannot separate. Marker presence is necessary, not sufficient.",
        fn=_cap_panel_check, produces=(Uns.PANEL_CHECK,), copilot_exposed=True,
        valid=_valid_panel_check,
    ),
    Capability(
        name="reference_match", label="Reference match",
        description="How well does a single-cell reference match this panel? Global panel-gene "
                    "overlap plus depth-matched per-type resolvability (which cell types the panel "
                    "can confidently transfer from this reference and which it cannot), and a "
                    "clustering nudge when the reference is a poor fit. Grounded in "
                    "eval_metrics.panel_resolvability (supersedes identifiability AUC). Pass "
                    "reference_path for a .h5ad, else it uses the reference chosen at the Panel-check "
                    "step. Needs no prior clustering or annotation.",
        fn=_cap_reference_match,
        params={"reference_path": {"type": "string", "description": "path to a reference .h5ad (optional; else the one chosen at Panel-check)"},
                "label_key": {"type": "string", "description": "reference cell-type column (optional; auto-detected)"}},
        produces=(Uns.REFERENCE_MATCH,), copilot_exposed=True,
    ),
    Capability(
        name="auto_select_reference", label="Auto-select reference (free text)",
        description="Given a FREE-TEXT tissue/tumour context (e.g. 'uveal melanoma', 'mouse brain', "
                    "'lung adenocarcinoma'), automatically CHOOSE and LOAD the best-matched "
                    "pre-computed single-cell reference for this panel - not just rank them. Picks the "
                    "top registry reference whose keywords + panel-gene overlap fit the tissue and "
                    "loads it (skips a wrong-tissue atlas even if it is the only one available), or "
                    "when allow_fetch is set and nothing local fits, fetches a small CELLxGENE "
                    "reference live via gget. Then scores the reference<->panel match. Use this when "
                    "the user names a tissue instead of uploading a reference. Degrades to a "
                    "'no_reference' note (cluster instead) when nothing suitable is found.",
        fn=_cap_auto_select_reference,
        params={"tissue": {"type": "string", "description": "free-text tissue/tumour context (optional; else the Load-step tissue)"},
                "allow_fetch": {"type": "boolean", "description": "allow a live CELLxGENE gget fetch when no local reference fits (default false)"}},
        copilot_exposed=True,
    ),
    Capability(
        name="annotation_strategy", label="Annotation strategy (supervised vs clustering)",
        description="A REFERENCE tool: answers 'is my single-cell REFERENCE good enough to annotate "
                    "from, or should I cluster instead?'. Self-verifies the reference<->panel match and "
                    "reruns a ladder (coarsen the labels the panel cannot separate, then try the best-"
                    "overlap tissue-matched reference), then RECOMMENDS supervised transfer vs de-novo "
                    "clustering. NEEDS a reference (from Load or reference_path); with none it returns "
                    "status 'no_reference' (an empty ladder) and recommends marker-based annotation - "
                    "'no_reference' means NO ATLAS IS ATTACHED, it does NOT mean no section is loaded, "
                    "so never tell the user to load data because of it. To IMPROVE annotation quality by "
                    "MERGING confusable cell types with no reference (e.g. on the demo), use "
                    "`merge_types` or `self_heal` (merge_confusable=true) instead. Writes "
                    "uns['annotation_route'].",
        fn=_cap_annotation_strategy,
        params={"reference_path": {"type": "string", "description": "path to a reference .h5ad (optional; else the one chosen at Panel-check)"},
                "label_key": {"type": "string", "description": "reference cell-type column (optional; auto-detected)"},
                "max_rounds": {"type": "integer", "description": "max coarsen/reselect rounds (default 3)"}},
        produces=(Uns.ANNOTATION_ROUTE,), copilot_exposed=True,
    ),
    Capability(
        name="compute_qc", label="Compute QC",
        description="Per-cell QC metrics (counts, genes, % control) and the section summary.",
        fn=_cap_compute_qc,
        produces=(Obs.TOTAL_COUNTS, Obs.N_GENES, Obs.PCT_CONTROL),
    ),
    Capability(
        name="qc_funnel", label="QC funnel",
        description="Full six-layer QC funnel headline: segmentation, panel-indexed count floor, "
                    "purity, spatial coherence, and the annotatability breakdown "
                    "(pct pass/warn/abstain + top abstention reasons).",
        fn=_cap_qc_funnel, produces=(), copilot_exposed=True,   # pure report (QC is a side effect)
    ),
    Capability(
        name="cluster", label="Cluster",
        description="Normalize, embed and Leiden-cluster the section at a given resolution.",
        fn=_cap_cluster,
        params={"resolution": {"type": "number", "description": "Leiden resolution (0.2-2.0; default 0.5 - higher splits into more clusters)."}},
        produces=(Obs.LEIDEN, Obsm.X_UMAP, Uns.RANK_GENES),
    ),
    Capability(
        name="cluster_markers", label="Top DEGs per cluster",
        description="Top-N differentially-expressed genes for every Leiden cluster, drawn as a "
                    "dot-plot (colour = mean expression, size = % of cells expressing). Use for "
                    "'what are the marker genes of each cluster' or 'top 2 DEGs per cluster'. Reuses "
                    "the ranking computed during clustering, so it is near-instant; recomputes only "
                    "when the clusters changed.",
        fn=_cap_cluster_markers,
        params={"n_genes": {"type": "integer", "description": "genes per cluster (default 2)"}},
        requires=(Obs.LEIDEN,), produces=(), copilot_exposed=True,
    ),
    Capability(
        name="cluster_confidence", label="Cluster confidence",
        description="Data-driven over/under-clustering check: which cluster PAIRS are statistically "
                    "indistinguishable (merge) and which single clusters hide substructure (split). "
                    "Grounded in a RandomForest-vs-permutation test (sc-SHC/CHOIR-style p-value + "
                    "accuracy) and a bimodality coefficient - advisory nudges, never auto-applied. "
                    "The cluster rung of the uncertainty ladder (cluster then cell then panel).",
        fn=_cap_cluster_confidence, requires=(Obs.LEIDEN,),
        produces=(), copilot_exposed=True,
    ),
    Capability(
        name="annotate", label="Annotate",
        description="Marker + consensus cell-type annotation with per-cell confidence and abstention.",
        fn=_cap_annotate, requires=(Obs.LEIDEN,),
        params={
            "reliability_weighted": {
                "type": "boolean",
                "description": "Reconcile the per-cell labels of any reference methods that ran (RCTD / "
                               "SingleR / scANVI / panhumanpy) by a RELIABILITY-weighted vote instead of "
                               "leaving the cluster label alone. Weighted rather than naive-majority "
                               "because majority loses to the best single method when one dominates. "
                               "Default false."},
            "truth_key": {
                "type": "string",
                "description": "obs column of ground-truth labels used to MEASURE each method's reliability. "
                               "Without it the documented DEFAULT_RELIABILITY prior is used instead."},
        },
        produces=(Obs.CELL_TYPE, Obs.CELL_TYPE_FINAL, Obs.ANNOTATION_CONFIDENCE,
                  Obs.ANNOTATION_VERDICT, Obs.ANNOTATION_REASON),
    ),
    Capability(
        name="immune_exclusion", label="Immune exclusion",
        description="Neighborhood-enrichment z-score between two cell types; tells whether one is "
                    "spatially excluded from (negative) or infiltrating (positive) the other.",
        fn=_cap_immune_exclusion,
        params={"type_a": {"type": "string"}, "type_b": {"type": "string"}},
        required_params=("type_a", "type_b"), requires=(Obs.CELL_TYPE,),
        produces=(), copilot_exposed=True,   # report; the spatial graph is a shared intermediate
    ),
    Capability(
        name="neighborhood_enrichment", label="Neighborhood enrichment",
        description="Full cell-type x cell-type neighborhood-enrichment z-score matrix.",
        fn=_cap_nhood, requires=(Obs.CELL_TYPE,),
        produces=(), copilot_exposed=True,   # report; the spatial graph is a shared intermediate
    ),
    Capability(
        name="niches", label="Niches",
        description="Call TME spatial niches (neighborhood composition) and list them.",
        fn=_cap_niches, requires=(Obs.CELL_TYPE,), produces=(Obs.NICHE,), copilot_exposed=True,
    ),
    Capability(
        name="co_occurrence", label="Co-occurrence",
        description="Cell-type co-occurrence probability vs. spatial distance - which cell types are "
                    "found together, and at what radius (squidpy co_occurrence).",
        fn=_cap_co_occurrence, requires=(Obs.CELL_TYPE,),
        produces=(), copilot_exposed=True,   # report; the spatial graph is a shared intermediate
    ),
    Capability(
        name="spatial_genes", label="Spatially variable genes",
        description="Top spatially variable genes by Moran's I - genes whose expression is spatially "
                    "structured across the section rather than randomly distributed.",
        fn=_cap_spatial_genes,
        params={"n_top": {"type": "integer",
                          "description": "How many top genes to return (default 20)."}},
        requires=(), produces=(),   # label-independent; not copilot-exposed (keeps the tool count modest)
    ),
    Capability(
        name="state_by_celltype", label="Cell states",
        description="Cell-type x cell-state (cycling / IFN / hypoxia / exhaustion ...) mean-score matrix.",
        fn=_cap_states, requires=(Obs.CELL_TYPE,), produces=(Uns.STATE_COLUMNS,),
    ),
    Capability(
        name="assign_cell_states", label="Type cell states",
        description="Type each cell with its DOMINANT cell-state program (cycling / interferon / "
                    "hypoxia / stress / EMT / T-exhaustion / T-cytotoxicity / ECM-remodeling / "
                    "antigen-presentation), CyteType-style - a colourable "
                    "obs['cell_state'] label, plus the state distribution and each cell type's "
                    "dominant state. Names the states in plain language via the LLM when a key is set. "
                    "States are orthogonal to lineage identity (who a cell IS vs what it is DOING).",
        fn=_cap_assign_states,
        params={"min_z": {"type": "number", "description": "z-score a program must clear to label a cell (default 1.0; lower = more cells typed)."}},
        requires=(Obs.CELL_TYPE,), produces=(Obs.CELL_STATE,), copilot_exposed=True,
    ),
    Capability(
        name="malignant_score", label="Malignant score",
        description="Where is the tumor? Marker-based malignant score per cell; returns the "
                    "fraction of high-malignant cells.",
        fn=_cap_malignant, produces=(Obs.MALIGNANT_SCORE,), copilot_exposed=True,
    ),
    Capability(
        name="discover_programs", label="De-novo programs",
        description="Discover data-driven gene programs (NMF) beyond the fixed marker/state lists; "
                    "returns the top genes per program plus a PLAID enrichment score (how well each "
                    "program's cells express its own signature: mean in vs out + AUROC specificity).",
        fn=_cap_programs,
        params={"n_programs": {"type": "integer", "description": "Number of NMF programs."}},
        produces=(Obsm.PROGRAMS, Obsm.PROGRAM_SCORES, Obs.PROGRAM, Varm.PROGRAM_LOADINGS),
        copilot_exposed=True,
    ),
    Capability(
        name="name_programs", label="Label programs",
        description="Name each de-novo NMF program in plain language from its top loading genes "
                    "(grounded; degrades to 'Program k' with no API key). Relabels obs['program'] so "
                    "the map legend reads by biological program, and returns the program table with a "
                    "stable program_id (the program_score_<k> colour field is preserved).",
        fn=_cap_name_programs, produces=(Obs.PROGRAM,), copilot_exposed=True,
    ),
    Capability(
        name="subcluster", label="Subcluster",
        description="Subcluster ONE named cell type into finer subtypes and name them - on demand, for "
                    "the specific type the user asks about. Use this whenever the user says to subcluster "
                    "/ break down / re-cluster / find subpopulations (or subtypes / substructure) WITHIN a "
                    "named cell type. (To instead auto-subcluster only the heterogeneous or marker-failing "
                    "types across the WHOLE annotation, use self_heal.)",
        fn=_cap_subcluster,
        params={"cell_type": {"type": "string", "description": "the exact cell type label to subcluster, "
                              "taken from the section's existing cell types (e.g. 'T cell', 'Myeloid')."},
                "resolution": {"type": "number", "description": "Leiden resolution for the subset; lower = fewer subtypes (default 0.1 - a single type is fairly homogeneous, so keep it conservative)."}},
        required_params=("cell_type",), requires=(Obs.CELL_TYPE,), produces=(Obs.SUBTYPE,),
        copilot_exposed=True,
    ),
    Capability(
        name="annotation_methods", label="Annotation methods",
        description="Which reference-based annotation methods ran (RCTD / SingleR / scANVI / "
                    "panhumanpy), their per-cell coverage, and the multi-method consensus-agreement "
                    "distribution. Reports only computed coverage/agreement - it does not assign labels.",
        fn=_cap_annotation_methods, copilot_exposed=True,   # pure report over method label columns
    ),
    Capability(
        name="rejection_reasons", label="Rejection reasons",
        description="Why weren't some cells confidently typed? Granular, plain-language reasons "
                    "(too few transcripts, spatial doublet, panel lacks markers, mixed lineages, "
                    "no clear winner, ...) with the count and % of untyped cells for each.",
        fn=_cap_rejection_reasons, requires=(Obs.ANNOTATION_VERDICT,),
        produces=(Obs.REJECTION_REASON, Obs.REJECTION_DETAIL), copilot_exposed=True,
    ),
    Capability(
        name="self_verify", label="Self-verify annotation",
        description="Check whether each cell type's cells actually score highest on that type's "
                    "canonical markers (marker-argmax agreement + one-vs-rest AUC), flag the types "
                    "that fail, and suggest concrete grounded fixes (subcluster / abstain / merge a "
                    "confusable pair). Advisory - it never changes labels.",
        fn=_cap_self_verify,
        params={"neighborhood": {"type": "boolean", "description": "Also corroborate with a spatial "
                "self-consistency check (each type should neighbor itself)."}},
        requires=(Obs.CELL_TYPE,), produces=(Uns.ANNOTATION_VERIFICATION,), copilot_exposed=True,
    ),
    Capability(
        name="trust_ledger", label="Trust ledger",
        description="Three INDEPENDENT per-cell-type verdicts and their disagreements: resolvable? "
                    "(panel adequacy) x coherent? (cells express the type's markers) x agreed? (the "
                    "reference methods back the label, when >=3 voted). The informative rows are the "
                    "CONTRADICTIONS - especially 'coherent but DISPUTED', the coherent-whole-cluster-mislabel "
                    "the AQI index alone cannot see. Advisory - never changes labels.",
        fn=_cap_trust_ledger,
        requires=(Obs.CELL_TYPE,), produces=(), copilot_exposed=True,
    ),
    Capability(
        name="merge_types", label="Merge confusable types",
        description="IMPROVE / optimize annotation QUALITY by MERGING the cell types the panel cannot "
                    "separate (no private on-panel marker) into ONE coarser label, and report the "
                    "quality delta (marker-agreement + number of failing types, before vs after). "
                    "Reference-FREE: works on ANY already-annotated section (the loaded demo included) "
                    "- no atlas needed and it does NOT reload the data. This is the tool for 'merge / "
                    "coarsen the confusable cell types to raise annotation quality'. APPLIES the merges "
                    "(self_verify only advises). Defaults to the groups panel_check flagged; pass "
                    "`groups` to merge specific sets. A coin-flip between indistinguishable types "
                    "becomes one defensible call.",
        fn=_cap_merge_types,
        params={"groups": {"type": "array", "items": {"type": "array", "items": {"type": "string"}},
                           "description": "Groups of cell-type names to each collapse into one "
                                          "label; default = the confusable groups panel_check found."},
                "dry_run": {"type": "boolean", "description": "preview the ontology-partitioned merges "
                            "(would_merge + cell counts) WITHOUT applying them"}},
        requires=(Obs.CELL_TYPE,), produces=(), copilot_exposed=True,
    ),
    Capability(
        name="calibrate_confidence", label="Calibrate confidence",
        description="Fit a post-hoc ISOTONIC calibration of the per-cell annotation confidence against a "
                    "ground-truth label column, so the score finally means P(correct). Writes "
                    "annotation_confidence_calibrated plus a reliability report (ECE and Brier before/after, "
                    "reliability curves, and an AUC discrimination check); the raw heuristic is never "
                    "overwritten. Benchmarked on a 12k CosMx section the raw score had ECE 0.33 with accuracy "
                    "NON-monotonic in confidence. Skips honestly with no labels - calibration is claimed only "
                    "when it is fit against real labels. Note the AUC: calibration fixes what the number MEANS "
                    "but cannot change how it RANKS cells, and abstention depends on the ranking.",
        fn=_cap_calibrate_confidence,
        params={
            "truth_key": {"type": "string",
                          "description": "obs column of ground-truth labels, on the SAME label axis as "
                                         "pred_key. Without it the capability skips."},
            "pred_key": {"type": "string",
                         "description": "prediction column scored for correctness. Default 'cell_type'."},
            "cal_frac": {"type": "number",
                         "description": "fraction held out to FIT on; ECE is reported on the complement. "
                                        "Default 0.5."},
            "n_bins": {"type": "number", "description": "reliability/ECE bins. Default 10."},
        },
        requires=(Obs.ANNOTATION_CONFIDENCE,), produces=(), copilot_exposed=True,
    ),
    Capability(
        name="highlight_cells", label="Highlight cells",
        description="Render a NEW spatial view on the canvas with cells matching a criterion lit up. "
                    "Call this whenever the user asks to SEE / show / highlight / where-are a population "
                    "(e.g. 'highlight low-quality cells', 'show the T cells', 'where are the malignant cells'). "
                    "`criterion` accepts: a cell-type name, 'low quality', 'low confidence' / 'uncertain' / "
                    "'abstained', or 'malignant'. Returns the matched count + an honest note (it states any "
                    "fallback, e.g. when no cells were flagged) and draws the plot - then describe what it shows.",
        fn=_cap_highlight,
        params={"criterion": {"type": "string",
                              "description": "a cell type, 'low quality', 'low confidence', 'abstained', or 'malignant'"}},
        required_params=("criterion",), copilot_exposed=True,
    ),
    Capability(
        name="show_spatial", label="Show spatial",
        description="Draw the spatial map colored by a GENE's expression or a field. Use for "
                    "'color/show the tissue by <gene>', 'show <gene> expression', 'color by total counts / "
                    "malignant / niche / cell type'. `color_by` is a gene symbol or a field name. Only "
                    "on-panel genes and real fields are accepted (it says so if the key is unknown).",
        fn=_cap_show_spatial,
        params={"color_by": {"type": "string", "description": "a gene symbol, or 'cell type' / 'total counts' / 'malignant' / 'niche'"}},
        required_params=("color_by",), copilot_exposed=True,
    ),
    Capability(
        name="marker_dotplot", label="Marker dot-plot",
        description="Draw a dot-plot of marker genes across cell types (dot color = expression scaled "
                    "0-1 PER GENE i.e. relative across types, size = % of cells expressing). Use for "
                    "'dotplot of <cell type> markers' or 'dotplot of GENE1, GENE2 ...'. Off-panel genes "
                    "are dropped (reported). Once spillover purification has run (a split_corrected layer "
                    "exists) it automatically draws the raw + SPLIT-purified pair.",
        fn=_cap_marker_dotplot,
        params={"cell_type": {"type": "string", "description": "use this cell type's canonical markers (optional)"},
                "genes": {"type": "array", "items": {"type": "string"}, "description": "explicit gene list (optional)"},
                "corrected": {"type": "boolean", "description": "deprecated: the raw+purified pair is now drawn automatically whenever a split_corrected layer exists"}},
        requires=(Obs.CELL_TYPE,), copilot_exposed=True,
    ),
    Capability(
        name="annotate_rctd", label="Annotate with RCTD",
        description="Run RCTD doublet-mode deconvolution (rctd-py subprocess, GPU-optional) against the "
                    "uploaded single-cell reference and write obs['rctd_first_type'] (+ second_type, "
                    "spot_class, weight, singlet_score). RCTD was the most accurate annotator on the "
                    "non-circular benchmark, so its per-cell labels feed the consensus vote; "
                    "set_as_primary adopts them as the primary cell_type. Needs the rctd-py env "
                    "(SPATIALSCRIBE_RCTD_PYTHON) and a reference; skips honestly otherwise. Slower than "
                    "the default TACCO annotation - use when you want RCTD's doublet / contamination model.",
        fn=_cap_annotate_rctd,
        params={"max_cells": {"type": "integer", "description": "subsample cap for tractability (0 = all)."},
                "set_as_primary": {"type": "boolean",
                                   "description": "adopt rctd_first_type as obs['cell_type'] (default false)."}},
        requires=(), produces=(Obs.CELL_TYPE,), copilot_exposed=True,
    ),
    Capability(
        name="split_purify", label="SPLIT spillover purify",
        description="Decontaminate transcript spillover between neighbouring cells (SPLIT reference-path "
                    "purification with residual-contamination removal): deconvolve the section against the "
                    "reference, then SPLIT::rctd_free_purify, writing the purified counts to "
                    "layers['split_corrected']. Deconvolution weights come from TACCO by default (fast, "
                    "in-env); weights_engine='rctd' uses the rctd-py subprocess instead. Reports the median "
                    "library-size reduction (spillover removed). When the SPLIT-R env or a reference are "
                    "absent it falls back to an in-app neighbour+marker decontamination "
                    "(method='marker_neighbour'). Then compare markers with marker_dotplot(corrected=true).",
        fn=_cap_split_purify,
        params={"max_cells": {"type": "integer", "description": "subsample cap for tractability (0 = all)."},
                "weights_engine": {"type": "string", "enum": ["tacco", "rctd"],
                                   "description": "deconvolution engine for SPLIT weights (default tacco)."}},
        requires=(Obs.CELL_TYPE,), produces=(), copilot_exposed=True,
    ),
    Capability(
        name="expression_violin", label="Expression violin",
        description="Draw a violin of one gene's expression across cell types. Use for 'how is <gene> "
                    "expressed across cell types' / 'violin of <gene>'. Off-panel gene -> reported, no plot.",
        fn=_cap_expression_violin,
        params={"gene": {"type": "string", "description": "gene symbol (must be on the panel)"}},
        required_params=("gene",), requires=(Obs.CELL_TYPE,), copilot_exposed=True,
    ),
    Capability(
        name="composition_chart", label="Composition chart",
        description="Draw a bar chart of the cell-type composition (counts per type). Use for "
                    "'show the composition' / 'what are the cell-type proportions'.",
        fn=_cap_composition_chart, requires=(Obs.CELL_TYPE,), copilot_exposed=True,
    ),
    Capability(
        name="show_segmentation", label="Show segmented cells",
        description="Draw the actual cell-segmentation POLYGONS (not centroids) colored by a field. Use "
                    "for 'show the segmented cells' / 'show cell boundaries'. Needs a Xenium run directory "
                    "with cell_boundaries.parquet - says so honestly if only centroids are available.",
        fn=_cap_show_segmentation,
        params={"color_by": {"type": "string", "description": "field to color the cells by (default cell_type)"}},
        copilot_exposed=True,
    ),
    Capability(
        name="reference_transfer", label="Reference transfer",
        description="Transfer cell-type labels from a user single-cell reference onto the section "
                    "(reference-anchored annotation). Trains a CellTypist model on the reference "
                    "(the reliable in-env arm) and, when the optional `tacco` package is installed, "
                    "also runs TACCO optimal-transport transfer. Writes per-cell reference labels and "
                    "reports coverage + agreement with the current cell_type. Uses the reference "
                    "chosen at the Panel-check step (or pass reference_path); degrades to 'no_reference'. "
                    "Set set_as_primary=true to adopt the CellTypist labels as cell_type.",
        fn=_cap_reference_transfer,
        params={"reference_path": {"type": "string", "description": "path to a reference .h5ad (optional; else the one chosen at Panel-check)"},
                "label_key": {"type": "string", "description": "reference cell-type column (optional; auto-detected)"},
                "set_as_primary": {"type": "boolean", "description": "adopt the transferred labels as cell_type (default false)"}},
        produces=(), copilot_exposed=True,
    ),
    Capability(
        name="malignant_concordance", label="Malignant concordance",
        description="Where is the tumour, cross-checked? Runs the malignant callers - the always-on "
                    "marker score plus the LEARNED Cancer-Finder caller (isolated env) - and reports "
                    "each caller's %-malignant. Needs a tumour tissue context; Cancer-Finder degrades "
                    "to a 'skipped' status when its env is not configured. (InSituCNV is not run "
                    "interactively - its infercnvpy pass is minutes-long; use it offline/in batch.)",
        fn=_cap_malignant_concordance,
        params={"max_cells": {"type": "integer", "description": "subsample cap for the learned callers (default 25000)."},
                "cf_threshold": {"type": "number", "description": "Cancer-Finder probability threshold (default 0.5)."}},
        produces=(Obs.MALIGNANT_SCORE,), copilot_exposed=True,
    ),
    Capability(
        name="self_heal", label="Self-verify + re-run",
        description="Self-verify the annotation, then RE-RUN what fails: flags cell types whose cells "
                    "disagree with their canonical markers and automatically subclusters the "
                    "heterogeneous ones, re-verifying each round (up to max_rounds). Set "
                    "merge_confusable to ALSO apply the panel-gap coarsening first (collapse the "
                    "types the panel cannot separate into one label) and report the quality delta; "
                    "relabel stays advisory. Returns the per-round action log, the merge delta, and "
                    "the before/after failing types.",
        fn=_cap_self_heal,
        params={"max_rounds": {"type": "integer", "description": "max verify -> fix -> re-verify rounds (default 2)."},
                "merge_confusable": {"type": "boolean", "description": "Also APPLY the panel-gap "
                    "merges (collapse types the panel cannot separate) and report the quality delta "
                    "(default false - merging is otherwise advisory)."}},
        requires=(Obs.CELL_TYPE,), produces=(), copilot_exposed=True,
    ),
]

REGISTRY: dict[str, Capability] = {c.name: c for c in _CAPABILITIES}


# --------------------------------------------------------------------------- #
# Introspection helpers
# --------------------------------------------------------------------------- #
def get(name: str) -> Capability | None:
    return REGISTRY.get(name)


def copilot_tools() -> list[dict]:
    """Anthropic tool schemas for the copilot-exposed capabilities (drives ``agent.tools``)."""
    return [c.to_tool_schema() for c in REGISTRY.values() if c.copilot_exposed]


def copilot_names() -> set[str]:
    return {c.name for c in REGISTRY.values() if c.copilot_exposed}


def producer_of(key: Key) -> str | None:
    """Name of the first capability that produces ``key`` (used for prerequisite hints)."""
    for cap in REGISTRY.values():
        if key in cap.produces:
            return cap.name
    return None


def missing_prereqs(adata, cap: Capability) -> list[Key]:
    """The subset of ``cap.requires`` not present on ``adata``."""
    return [k for k in cap.requires if not k.present(adata)]


def _prereq_hint(missing: list[Key]) -> str:
    producers: list[str] = []
    for k in missing:
        p = producer_of(k)
        if p and p not in producers:
            producers.append(p)
    if not producers:
        return ""
    return "run " + " then ".join(f"'{p}'" for p in producers) + " first"


# --------------------------------------------------------------------------- #
# Dispatch - the single path both frontends call
# --------------------------------------------------------------------------- #
def run(adata, name: str, params: dict | None = None, ctx: RunContext | None = None) -> RunResult:
    """Run capability ``name`` on ``adata``.

    Checks prerequisites first (returning a structured ``PrerequisiteError`` dict rather than
    throwing), then invokes the adapter. On success returns the JSON-able value plus a
    provenance record ``{name, params}``; on any failure returns a structured error dict.
    """
    params = dict(params or {})
    ctx = ctx or RunContext()
    cap = REGISTRY.get(name)
    if cap is None:
        return RunResult(name, error={"error_type": "unknown_capability", "capability": name,
                                      "message": f"unknown capability {name!r}"})
    missing = missing_prereqs(adata, cap)
    if missing:
        err = PrerequisiteError(name, missing, _prereq_hint(missing))
        return RunResult(name, error=err.to_dict())
    try:
        value = cap.fn(adata, ctx, **params)
    except PrerequisiteError as exc:
        return RunResult(name, error=exc.to_dict())
    except Exception as exc:  # noqa: BLE001 - deliberately structured for the copilot
        return RunResult(name, error=to_error_dict(exc, name))
    return RunResult(name, value=value, record={"name": name, "params": params})


def ensure(adata, name: str, ctx: RunContext | None = None, force: bool = False,
           params: dict | None = None) -> RunResult:
    """Compute-if-absent-or-forced. Skips (empty ``RunResult``) when all ``produces`` are present.

    Standardizes the ``if key not in obs or button`` idiom: a step is considered done when
    every key it declares as ``produces`` is already on ``adata``.
    """
    cap = REGISTRY.get(name)
    if cap is None:
        return run(adata, name, params, ctx)  # yields the unknown-capability error
    if (not force and cap.produces and all(k.present(adata) for k in cap.produces)
            and _is_valid(cap, adata)):
        return RunResult(name, value=None, record=None)
    return run(adata, name, params, ctx)


def _is_valid(cap: Capability, adata) -> bool:
    """True if ``cap`` has no validity predicate, or it passes (a raising predicate -> invalid).

    Lets ``ensure()`` recompute a product that is present but corrupt - e.g. a ``panel_check``
    whose '/'-keyed cell types were split into nested groups by an h5ad round-trip (the source of
    ``KeyError: 'n_present'`` in the Panel-check step).
    """
    if cap.valid is None:
        return True
    try:
        return bool(cap.valid(adata))
    except Exception:
        return False
