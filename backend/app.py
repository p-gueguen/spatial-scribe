"""FastAPI backend for the SpatialScribe React SPA (branch: feat/react-spa).

Wraps the UI-agnostic ``spatialscribe.analysis`` capability registry behind HTTP. A session holds the
AnnData server-side (it can't go to the browser); endpoints run capabilities via ``cap.run`` and the
copilot runs the same Claude tool-use loop, returning a map-recolour directive. The React app never
touches AnnData - it asks for viewport points + step results + copilot answers.

Run (from the worktree root, gpu-env python):
    ANTHROPIC_API_KEY=... uvicorn backend.app:app --host 0.0.0.0 --port 8000

This is the Phase-2 target from docs/FRONTEND_STACK.md: same Python engine, React frontend.
"""
from __future__ import annotations

import os

# Cap BLAS/OMP threads BEFORE numpy/scipy/anndata import. On a compute node (64 cores) OpenBLAS otherwise
# spawns per-core threads on every call, and concurrent heavy requests (the QC funnel's spatial
# neighbours + KMeans) blow OpenBLAS's 128 memory-region limit -> the whole backend dies mid-request
# ("Program is Terminated. Because you tried to allocate too many memory regions"), which surfaced as
# the QC funnel stuck on "not run". setdefault, so the deploy env can still raise it if a host wants.
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "8")

import hashlib
import json
import threading
import uuid

import anndata as ad
import numpy as np
from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import (
    FileResponse, HTMLResponse, JSONResponse, PlainTextResponse, Response, StreamingResponse,
)

try:  # optional: orjson serializes numpy directly (float32 coords stay compact) and ~50x faster than
    import orjson as _orjson  # stdlib json on the big /points payload. Absent -> graceful stdlib fallback
except ImportError:           # (still bypasses FastAPI's jsonable_encoder), never a hard import failure.
    _orjson = None
from pydantic import BaseModel

from backend import stt
from backend.report import report_html
from spatialscribe.analysis import capabilities as cap, llm, plots, status

# Processed demo cache. `breast` = a public Xenium Prime 5K breast section (CC BY 4.0, 10x
# Genomics), shipped downsampled at data/demo_public.h5ad; override with SPATIALSCRIBE_DEMO_CACHE.
# The offline synthetic demo (/api/load_synthetic) needs no data file. See docs/DATASETS.md.
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEMO_CACHES = {
    "breast": os.environ.get("SPATIALSCRIBE_DEMO_CACHE")
    or os.path.join(_REPO, "data", "demo_public.h5ad"),
}
DEMO_CACHE = DEMO_CACHES["breast"]      # back-compat: the default section the app opens on
STEP_CAP = {"panel": "panel_check", "qc": "compute_qc", "cluster": "cluster",
            "annotate": "annotate", "spatial": "niches"}
_CAP_TO_STEP = {v: k for k, v in STEP_CAP.items()}   # capability name -> its rail step key


def _mark_ran(s: dict, cap_name: str) -> None:
    """Record that a rail STEP capability ran on this session, so the left rail greens off work the
    user actually did (not the demo's shipped pre-processing). Shared by run_step / run_cap / the
    panel auto-run so ANY of them advances the rail. Non-step capabilities (qc_funnel, subcluster,
    rank_degs, ...) are ignored. (The copilot runs caps out-of-band and is not tracked here - its map
    recolours are transient, not a rail step.)"""
    step = _CAP_TO_STEP.get(str(cap_name))
    if step:
        s.setdefault("ran", [])
        if step not in s["ran"]:
            s["ran"].append(step)

app = FastAPI(title="SpatialScribe API")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])


@app.on_event("startup")
def _warm_backend() -> None:
    """JIT-compile pynndescent + UMAP off the request path (background daemon) so the first clustering
    request doesn't eat the one-time ~30s numba compile. Best-effort; never blocks startup."""
    from spatialscribe.analysis import preload
    preload.warm_backend()
    stt.warm()   # ~5s of weights load, or the FIRST dictated question pays it

_SESSIONS: dict[str, dict] = {}   # sid -> {adata, tissue, log}


def _sess(sid: str) -> dict:
    s = _SESSIONS.get(sid)
    if s is None:
        raise HTTPException(404, "unknown session (reload the section)")
    return s


def _require_idle(s: dict) -> None:
    """409 while a background full-analysis run holds the session (single writer on adata)."""
    if s.get("pipeline_running"):
        raise HTTPException(409, "a full analysis is running on this session; wait for it to finish")


def _progress(s: dict) -> dict:
    """Per-session coarse progress the running step reports into, so GET /api/{sid}/progress can
    answer while cap.run is busy in the threadpool. Created lazily here so every session-creation site
    is covered without touching each one."""
    return s.setdefault("progress", {"running": False, "step": None, "frac": 0.0, "label": ""})


def _verdict_cache(s: dict, kind: str, inputs, compute):
    """Memoize a per-session LLM verdict so Claude is billed at most once per (section, inputs).

    ``kind`` namespaces the entry ('panel' | 'qc'); ``inputs`` is any JSON-able object whose content
    fully determines the verdict (the same values fed into the prompt); ``compute`` is a zero-arg
    callable that performs the actual (billed) LLM call. On a signature hit the cached value is
    returned WITHOUT touching Claude - so re-visiting a tab (which remounts the panel and re-fetches)
    costs nothing. Errors are not cached (they raise before the store), so a transient LLM failure
    retries next visit; the cache lives on ``s['_verdicts']`` and dies when the session is replaced.
    A fresh QC/reference genuinely changes ``inputs`` -> new signature -> one fresh call."""
    sig = hashlib.sha1(json.dumps(inputs, sort_keys=True, default=str).encode()).hexdigest()
    cache = s.setdefault("_verdicts", {})
    hit = cache.get(kind)
    if hit is not None and hit[0] == sig:
        return hit[1]
    val = compute()
    cache[kind] = (sig, val)
    return val


def _ctx(s: dict):
    # Thread the session's uploaded reference into the RunContext so run_cap('reference_match'),
    # annotate and reference_transfer see it. Both are None until a reference is uploaded/registered
    # on the session.
    # is_tumour is the Data-step checkbox; None until the scientist answers, in which case the engine
    # falls back to its tissue-keyword heuristic (which fires on 'normal breast' and misses
    # 'glioblastoma'). See cnv.is_tumour_context.
    return cap.RunContext(tissue=s["tissue"], use_llm=llm.available(),
                          is_tumour=s.get("is_tumour"),
                          reference=s.get("reference"), ref_label_key=s.get("ref_label_key"))


def _obs_fields(a) -> list[str]:
    # transcript_counts is the raw Xenium per-cell count (present before QC), so a freshly-loaded
    # section still has a meaningful default coloring; total_counts appears once QC runs. The
    # malignant / reference-transfer columns appear once those capabilities run.
    # Group all annotation outputs together: the cell-type labels, then their per-cell quality fields
    # (confidence / verdict / rejection) right next to them so the label and how much to
    # trust it read as one block; then subtype/clustering, niche, malignancy, reference transfer, raw.
    fields = [c for c in ("cell_type", "cell_type_final", "annotation_confidence", "annotation_verdict",
                          "rejection_reason", "cell_state", "subtype", "leiden",
                          "niche", "malignant_score", "cnv_score", "is_malignant", "cancerfinder_prob",
                          "rctd_first_type", "celltypist_label", "ref_label", "program")
              if c in a.obs]
    # states.score_states writes one obs['state_<name>'] z-score per program; the canvas colours by
    # them as a per-cell featureplot (points()/_colours already serve any numeric obs column). Expose
    # them in the dropdown, sorted, right after 'program'.
    fields += sorted(c for c in a.obs.columns if c.startswith("state_"))
    # Genes the copilot coloured the map by (show_spatial writes obs[<gene>] + records it here), so the
    # dropdown lists them and the <select> reflects a gene colouring instead of showing blank.
    fields += [c for c in a.uns.get("gene_color_fields", []) if c in a.obs and c not in fields]
    fields += [c for c in ("total_counts", "transcript_counts") if c in a.obs]
    return fields


def _device() -> str:
    """"GPU" if the RAPIDS path is live, else "CPU" - the honest label for the rail badge."""
    try:
        from spatialscribe.analysis import backend as _be
        return "GPU" if _be.gpu_available() else "CPU"
    except Exception:
        return "CPU"


def _state(s: dict) -> dict:
    a = s["adata"]
    done = {"load": True}
    for k, cname in STEP_CAP.items():
        done[k] = bool(status.is_done(a, cname))
    # Report is the capstone export step (no capability of its own): mark it done once every
    # compute step is complete, so the rail tracker can actually reach 7/7.
    done["report"] = all(done.values())
    from spatialscribe.analysis import io as _io
    demo = a.uns.get("spatialscribe_demo") or None
    # Reference is a SESSION resource (set in Panel-check, its first consumer); expose a small summary
    # so the Annotate step can show its status + offer supervised transfer without re-uploading it.
    ref = s.get("reference")
    reference = None
    if ref is not None:
        rk = s.get("ref_label_key")
        reference = {"name": s.get("reference_name"), "n_cells": int(ref.n_obs), "label_key": rk,
                     "n_types": int(ref.obs[rk].astype(str).nunique()) if rk and rk in ref.obs else None}
    return {"n_obs": int(a.n_obs), "n_vars": int(a.n_vars), "tissue": s["tissue"],
            "is_tumour": s.get("is_tumour"),   # None = unanswered, engine falls back to the keyword gate
            # `done` is DATA-derived (an output exists) and gates functional UI (e.g. the UMAP toggle).
            # `ran` is the steps the USER actually ran this session; the rail greens off `ran`, so a
            # bundled demo (shipped pre-processed) does not read as already-done on first load.
            "obs_fields": _obs_fields(a), "done": done, "ran": list(s.get("ran") or []),
            "device": _device(),
            "has_key": llm.available(),
            # Which model backs the copilot, so the UI can show it. provider "openai" = a local /
            # self-hosted OpenAI-compatible endpoint (e.g. our vLLM Qwen); "anthropic" = the API.
            "llm": {"available": llm.available(), "provider": llm.provider(),
                    "model": (llm.default_model() or None),
                    "providers": llm.providers_available()},   # >1 -> the UI model tag is a toggle
            # faster-whisper installed -> the ask box offers a mic button (optional dep).
            "has_stt": stt.available(),
            # `demo` describes the bundled section (source / role) so the Load panel can say which
            # one is open instead of hardcoding the breast blurb. Absent for user-loaded sections.
            "demo": ({"source": str(demo.get("source", "")), "role": str(demo.get("role", ""))}
                     if isinstance(demo, dict) else None),
            "panel_name": _io.panel_label(a.uns),   # e.g. "hMulti_100g (Human, Multi, 480 targets)" for Xenium; None otherwise
            "reference": reference}   # the session's single-cell reference summary, or None


def _hex(h: str) -> list[int]:
    h = h.lstrip("#")
    return [int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)]


def _cat_colour(i: int) -> list[int]:
    """Distinct RGB for category index ``i``: the fixed 28-entry palette, then golden-ratio HSV so
    the mapping stays INJECTIVE past the palette length. Two categories sharing an RGB would break
    the map's legend-hover and heatmap-pair emphasis, which highlight cells by exact colour match
    (a >28-cluster leiden or a many-niche map would otherwise light up the wrong cells)."""
    pal = plots._PALETTE
    if i < len(pal):
        return _hex(pal[i])
    import colorsys

    r, g, b = colorsys.hsv_to_rgb((i * 0.618033988749895) % 1.0, 0.62, 0.95)
    return [int(r * 255), int(g * 255), int(b * 255)]


# Semantic status colours (match the AnnotatePanel meter + styles.css --pass/--warn/--violet).
# annotation_verdict (apply_confidence) PASS/WARN/FAIL -> green/amber/red: the SAME semantic scale as
# the AnnotatePanel confidence meter (styles.css --pass/--warn/--fail), not the arbitrary palette.
_VERDICT_COLOURS = {"PASS": "#46E39B", "WARN": "#F2B24C", "FAIL": "#F7746E"}

# Two continuous ramps, because a continuous field is one of two different things.
#   QC / depth / confidence  -> the brand violet ramp: it reads as "more of the same measurement".
#   A feature plot (a signature score, an NMF program, a malignancy probability) -> a magma ramp:
#     a hot, wide-gamut scale that separates a sparse positive tail from a dark background, which is
#     exactly what you look for in a featureplot and what the violet ramp (dark violet -> pale violet)
#     cannot show. Both stay legible on the near-black canvas.
_RAMP_POS = np.array([0.0, 0.5, 0.8, 1.0])
_RAMP_QC = np.array([[35, 30, 56], [108, 92, 224], [168, 150, 242], [235, 230, 255]], float)
_RAMP_FEATURE_POS = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
_RAMP_FEATURE = np.array([[24, 15, 62], [113, 31, 128], [203, 70, 109], [250, 148, 88],
                          [252, 253, 191]], float)
# Fields that ARE feature plots. Prefixes cover the generated ones: `state_<program>` (states.py) and
# `program_score_<k>` (programs.py); the rest are named per-cell scores.
_FEATURE_PREFIXES = ("state_", "program_score_")
_FEATURE_FIELDS = {"malignant_score", "cnv_score", "cancerfinder_prob"}


def _is_feature_field(name: str) -> bool:
    """True when a continuous obs column is a signature/program/probability score rather than a QC
    metric, so ``_colours`` can give it the magma ramp instead of the brand violet one."""
    return str(name).startswith(_FEATURE_PREFIXES) or str(name) in _FEATURE_FIELDS


def _ramp_hex(color_by: str) -> list[str]:
    """Hex stops (low->high) of the continuous ramp ``_colours`` paints for ``color_by`` - magma for a
    feature field, brand violet otherwise. The frontend's FeatureScale key renders these SAME stops so
    the legend can never disagree with the map (feature plots are magma, not the default violet)."""
    stops = _RAMP_FEATURE if _is_feature_field(color_by) else _RAMP_QC
    return ["#%02X%02X%02X" % (int(r), int(g), int(b)) for r, g, b in stops]


def _colours(a, color_by: str, idx):
    """Per-cell RGB (flat uint8) + categorical legend, matching the Plotly/deck palette."""
    col = a.obs[color_by]
    if color_by == "annotation_verdict":
        # Fixed green/amber/red for PASS/WARN/FAIL so the spatial overlay reads the same as the
        # confidence meter, instead of the arbitrary categorical palette (cyan/magenta/green).
        vals = col.astype(str).to_numpy()[idx]
        rgb = np.array([_hex(_VERDICT_COLOURS.get(v, "#6A7080")) for v in vals],
                       dtype=np.uint8).reshape(-1).tolist()
        present = [c for c in ("PASS", "WARN", "FAIL") if (vals == c).any()]
        return rgb, [{"label": c, "color": _hex(_VERDICT_COLOURS[c])} for c in present]
    if color_by in ("cell_type", "cell_type_final"):
        # Same treatment as the delivered PNG (export._fig_spatial): the abstention pseudo-labels are
        # not five lineages, so they collapse into ONE grey "Not assigned" class. The per-cell reason
        # is still inspectable via the annotation_reason / annotation_verdict overlays.
        col = _annotate.collapse_abstention(col)
    if str(col.dtype) in ("category", "object"):
        vals = col.astype(str).to_numpy()[idx]
        cats = plots._ordered_group_labels(col, vals)
        # Drop null-like tokens (NaN leiden = low-signal cells, unassigned cell_type): they are not a
        # real category, so they must not claim a colour/legend slot - they fall through to grey below.
        cats = [c for c in cats if str(c).lower() not in ("nan", "none", "na", "<na>")]
        typed = [c for c in cats if not _annotate.is_abstention(str(c))]
        pal = {c: _cat_colour(i) for i, c in enumerate(typed)}
        for c in cats:                                  # "Not assigned" reads grey, never a lineage hue
            if _annotate.is_abstention(str(c)):
                pal[c] = [136, 146, 166]
        rgb = np.array([pal.get(c, [136, 146, 166]) for c in vals], dtype=np.uint8).reshape(-1).tolist()
        return rgb, [{"label": str(c), "color": pal[c]} for c in cats[:35]]
    v = np.asarray(col)[idx].astype(float)
    finite = v[np.isfinite(v)]
    lo, hi = (np.percentile(finite, [2, 98]) if finite.size else (0.0, 1.0))
    t = np.clip((v - lo) / (hi - lo + 1e-9), 0, 1)
    feature = _is_feature_field(color_by)
    pos = _RAMP_FEATURE_POS if feature else _RAMP_POS
    stops = _RAMP_FEATURE if feature else _RAMP_QC
    rgb = np.stack([np.interp(t, pos, stops[:, k]) for k in range(3)], 1).astype(np.uint8).reshape(-1).tolist()
    return rgb, []


def _fig_to_png_b64(fig) -> str | None:
    """Render a plotly figure (from a display capability) to an inline base64 PNG for the SPA.

    A figure that sized itself (``layout.width`` / ``layout.height``) is rendered at that size,
    clamped: a 60-gene x 30-cluster dot-plot squeezed into the 760x460 default is unreadable. Figures
    that set neither keep the default box, so every existing view renders exactly as before.
    """
    try:
        w = int(getattr(fig.layout, "width", None) or 760)
        h = int(getattr(fig.layout, "height", None) or 460)
        w, h = max(320, min(w, 2200)), max(240, min(h, 1400))
        png = fig.to_image(format="png", width=w, height=h, scale=2)  # kaleido
        import base64 as _b64
        return "data:image/png;base64," + _b64.b64encode(png).decode()
    except Exception:
        return None


def _serialize_artifacts(arts: list) -> list:
    """Turn a capability's ctx.artifacts into SPA-renderable items (map_view directives + PNGs)."""
    out: list = []
    for a in arts or []:
        k = a.get("kind")
        if k == "map_view":
            out.append({"kind": "map_view", "color_by": a.get("color_by"),
                        "highlight": a.get("highlight"), "note": a.get("note")})
        elif k == "figure":
            png = _fig_to_png_b64(a.get("fig"))
            if png:
                out.append({"kind": "figure", "title": a.get("title"), "png": png})
    return out


def _clean(obj):
    """Make an engine return value strict-JSON safe: NaN/Inf -> null, numpy scalars/arrays -> python.
    Starlette's JSON encoder uses allow_nan=False, so any NaN (e.g. a QC median over an empty subset,
    an AUC of a degenerate type, a funnel layer with no signal) would otherwise 500 the endpoint."""
    import math
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, dict):
        return {k: _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    if isinstance(obj, np.floating):
        f = float(obj); return f if math.isfinite(f) else None
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.ndarray):
        return _clean(obj.tolist())
    return obj


def _fast_json(payload: dict) -> Response:
    """Serialize a /points payload bypassing FastAPI's jsonable_encoder (which walks ~1.4M scalars
    and holds the GIL ~600ms on a 200k-point map). Returning a Response instance skips it entirely.
    orjson (OPT_SERIALIZE_NUMPY) serializes the numpy float32 coord arrays in C and emits `null` for
    NaN/Inf (the browser's JSON.parse rejects NaN); without orjson, `_clean` makes the payload strict-
    JSON-safe (numpy->list, NaN->null) and Starlette's JSONResponse serializes it - still no
    jsonable_encoder. Either way the client sees the same JSON shape (arrays of numbers)."""
    if _orjson is not None:
        return Response(_orjson.dumps(payload, option=_orjson.OPT_SERIALIZE_NUMPY),
                        media_type="application/json")
    return JSONResponse(_clean(payload))


def _confidence_pass(a, ctx) -> None:
    """Match the Streamlit annotate meter to the QC funnel: add pmp + spatial-coherence penalties and
    recompute apply_confidence with the CURRENT formula (not a cached verdict)."""
    from spatialscribe.analysis import annotate as _an, markers as _mk, purity as _pur, spatial as _sp
    ms = ctx.markers()
    coarse = _mk.for_tissue(ctx.tissue)   # coarse DISJOINT lineages for the CRISP mutual-exclusivity check
    try:
        _pur.pmp(a, assigned_label_key="cell_type", lineage_markers=ms)
        _sp.spatial_coherence(a, label_key="cell_type")
    except Exception:
        pass
    _an.apply_confidence(a, cluster_key="cell_type", marker_sets=ms,
                         panel_check_result=a.uns.get("panel_check"), lineage_markers=coarse)


@app.get("/api/demos")
def list_demos():
    """The processed demo sections this server can actually open (cache present on disk)."""
    return {"demos": [{"name": n, "available": os.path.exists(p)} for n, p in DEMO_CACHES.items()]}


@app.post("/api/load_demo")
def load_demo(name: str = "breast"):
    path = DEMO_CACHES.get(name)
    if path is None:
        raise HTTPException(404, f"unknown demo '{name}' (have: {', '.join(DEMO_CACHES)})")
    if not os.path.exists(path):
        raise HTTPException(500, f"demo cache missing: {path}")
    a = ad.read_h5ad(path)
    a.obs.drop(columns=["annotation_confidence", "annotation_verdict", "annotation_reason",
                        "cell_type_final", "pmp", "spatial_coherence", "crisp_impure"],
               errors="ignore", inplace=True)
    a.uns.pop("panel_check", None)      # the h5ad round-trip mangles its '/'-keyed coverage dict
    # Derived with the annotation above: a cached value would surface a stale usability pill in the
    # funnel's section flags before this session re-annotates.
    a.uns.pop("annotation_usability", None)
    a.uns.setdefault("platform", "xenium")
    a.uns.setdefault("n_panel_genes",
                     int((~a.var["control"]).sum()) if "control" in a.var.columns else int(a.n_vars))
    # The bundled caches were clustered before the size relabel existed, so their leiden order is the
    # backend's community-detection order, not cell count. Relabel on load (it re-keys the cached
    # rank_genes_groups in step, no recompute) so the demo opens with cluster 0 = the largest.
    if "leiden" in a.obs:
        from spatialscribe.analysis import cluster as _cl
        _cl.relabel_clusters_by_size(a, "leiden")
    tissue = a.uns.get("spatialscribe_demo", {}).get("tissue", "breast")
    sid = uuid.uuid4().hex[:12]
    _SESSIONS[sid] = {"adata": a, "tissue": tissue, "log": []}
    return {"session_id": sid, **_state(_SESSIONS[sid])}


class LoadSectionReq(BaseModel):
    path: str
    tissue: str = ""   # blank = infer from the panel metadata (see load_section); NOT a melanoma default
    is_tumour: bool | None = None   # Data-step checkbox; None = decide from the tissue keyword


@app.post("/api/load_section")
def load_section(req: LoadSectionReq):
    """Ingest a real section directly from disk: a Xenium/CosMx/MERSCOPE output folder or a .h5ad
    (server-side path). Unlike the demo cache this is raw - the user then runs the pipeline steps.
    """
    p = cap._resolve_server_path((req.path or "").strip())   # correct a missing mount prefix (/projects -> /data)
    if p is None:
        raise HTTPException(400, f"path not found on the server: {(req.path or '').strip()!r}")
    try:
        from spatialscribe.analysis import io as _io
        sample = _io.load(p)
    except Exception as exc:
        raise HTTPException(400, f"could not load section: {type(exc).__name__}: {exc}")
    a = sample.adata
    if "spatial" not in a.obsm:
        raise HTTPException(400, "loaded section has no spatial coordinates (obsm['spatial'])")
    a.uns.setdefault("platform", getattr(sample, "platform", "spatial"))
    a.uns.setdefault("n_panel_genes", int(len(getattr(sample, "panel_genes", []) or a.var_names)))
    # Tissue: use what the user typed, else INFER from the panel metadata (as the copilot's load_section
    # does), else fall back to the app default. Hardcoding 'melanoma' for any unlabelled section auto-
    # selected a skin/melanoma reference on a non-skin sample (e.g. a duodenum), which transferred foreign
    # types (Schwann cell / melanocyte / OPC) - the bogus merges + false-green typability. Inference gives
    # a panel-appropriate tissue; the user can still override in the Load tab (e.g. 'human duodenum').
    tissue = (req.tissue or "").strip() or _io.infer_tissue(a) or "melanoma"
    sid = uuid.uuid4().hex[:12]
    _SESSIONS[sid] = {"adata": a, "tissue": tissue, "log": [],
                      "is_tumour": req.is_tumour, "source_path": p}
    return {"session_id": sid, **_state(_SESSIONS[sid])}


class IsTumourReq(BaseModel):
    is_tumour: bool | None = None   # None clears the answer -> back to the tissue-keyword fallback


@app.post("/api/{sid}/is_tumour")
def set_is_tumour(sid: str, req: IsTumourReq):
    """Data-step checkbox: does this sample contain malignant cells?

    Nothing can answer this before the malignant callers run, so it must come from the scientist.
    Without it the engine substring-matches the free-text tissue, which fires on 'normal breast' and
    misses 'glioblastoma'. Setting it True makes Cancer-Finder (and InSituCNV on a >2000-probe panel)
    run at the Annotate step; False skips them; null restores the heuristic.
    """
    s = _sess(sid)
    s["is_tumour"] = req.is_tumour
    return {"ok": True, **_state(s)}


def _register_reference(s: dict, path: str, label_key: str | None = None,
                        gene_name_col: str | None = None) -> dict:
    """Load a single-cell reference .h5ad into the session (auto-detecting the label column) and run
    the panel reference-match. The same flow as the Streamlit reference chooser, lifted to HTTP. On a
    load failure returns a 400; otherwise returns the match summary + fresh state."""
    from spatialscribe.analysis import reference as _ref
    try:
        ref, key = _ref.load_reference(path, label_key=label_key, gene_name_col=gene_name_col)
    except Exception as exc:
        raise HTTPException(400, f"could not load reference: {type(exc).__name__}: {exc}")
    s["reference"] = ref
    s["ref_label_key"] = key
    s["reference_name"] = os.path.basename(str(path))
    a = s["adata"]
    a.uns.pop("reference_match", None)                 # invalidate a match cached for the old reference
    match = None
    try:                                               # best-effort: show the panel fit immediately
        out = cap.run(a, "reference_match", {}, _ctx(s))
        match = _clean(out.value) if out.ok else None
    except Exception:
        match = None
    return {"ok": True, "n_ref_cells": int(ref.n_obs), "label_key": key,
            "n_labels": int(ref.obs[key].astype(str).nunique()) if key in ref.obs else None,
            "match": match, **_state(s)}


class RefReq(BaseModel):
    path: str
    label_key: str | None = None
    gene_name_col: str | None = None


@app.post("/api/{sid}/reference")
def set_reference(sid: str, req: RefReq):
    """Register a custom single-cell reference by SERVER-SIDE path (a .h5ad)."""
    s = _sess(sid)
    p = cap._resolve_server_path((req.path or "").strip())   # correct a missing mount prefix (/projects -> /data)
    if p is None:
        raise HTTPException(400, f"reference path not found on the server: {(req.path or '').strip()!r}")
    return _register_reference(s, p, req.label_key, req.gene_name_col)


@app.post("/api/{sid}/reference/upload")
async def upload_reference(sid: str, file: UploadFile = File(...),
                          label_key: str | None = None):
    """Register a custom single-cell reference by UPLOADING a .h5ad (multipart). Streams the file to
    a temp path server-side, then loads it (the AnnData is what the session keeps)."""
    s = _sess(sid)
    import tempfile
    name = os.path.basename(file.filename or "reference.h5ad")
    if not name.endswith(".h5ad"):
        raise HTTPException(400, "reference must be a .h5ad file")
    tmp = os.path.join(tempfile.mkdtemp(prefix="ssref_"), name)
    try:
        with open(tmp, "wb") as fh:
            while True:
                chunk = await file.read(1 << 20)
                if not chunk:
                    break
                fh.write(chunk)
    finally:
        await file.close()
    return _register_reference(s, tmp, label_key)


class AutoRefReq(BaseModel):
    tissue: str | None = None
    allow_fetch: bool = False


@app.post("/api/{sid}/reference/auto")
def auto_reference(sid: str, req: AutoRefReq | None = None):
    """Free-text tissue -> auto-CHOOSE + LOAD the best-matched pre-computed reference (registry, or a
    live CELLxGENE gget fetch when ``allow_fetch``), register it on the session, and return the panel
    match plus the supervised-vs-clustering route. When nothing suitable is found, returns
    ``ok: false`` with ``recommended_mode: 'cluster'`` (do NOT force a wrong-tissue transfer)."""
    from spatialscribe.analysis import reference as _ref
    from spatialscribe.analysis.capabilities import _panel_genes
    s = _sess(sid); a = s["adata"]
    tissue = ((req.tissue if req else None) or s.get("tissue") or "").strip()
    if tissue:
        s["tissue"] = tissue
    sel = _ref.auto_select_reference(tissue, panel_genes=_panel_genes(a),
                                     allow_fetch=bool(req and req.allow_fetch))
    auto = {k: sel.get(k) for k in ("status", "chosen", "source", "message", "ranked")}
    if sel.get("ref") is None:                              # nothing loaded -> cluster instead
        return {"ok": False, "auto": auto, "recommended_mode": "cluster", **_state(s)}
    s["reference"] = sel["ref"]
    s["ref_label_key"] = sel["label_key"]
    s["reference_name"] = sel.get("chosen")
    a.uns.pop("reference_match", None)
    ctx = _ctx(s)
    match = route = None
    try:
        out = cap.run(a, "reference_match", {}, ctx)
        match = _clean(out.value) if out.ok else None
    except Exception:
        match = None
    try:
        rout = cap.run(a, "annotation_strategy", {}, ctx)   # self-verify ladder -> supervised vs cluster
        route = _clean(rout.value) if rout.ok else None
    except Exception:
        route = None
    return {"ok": True, "n_ref_cells": int(sel["ref"].n_obs), "label_key": sel["label_key"],
            "auto": auto, "match": match,
            "recommended_mode": (route or {}).get("recommended_mode"), "route": route, **_state(s)}


@app.get("/api/{sid}/state")
def get_state(sid: str):
    return _state(_sess(sid))


def _basis_coords(a, basis: str, s: dict):
    """The 2-D coordinates the map is drawn in: the spatial layout (default hero) or the UMAP embedding.

    UMAP is off the analysis critical path, so it is computed lazily HERE only when the view is actually
    requested - via the ``umap`` capability when present (deferred embedding), else it uses whatever
    ``obsm['X_umap']`` clustering already produced. If no embedding exists yet (e.g. before clustering)
    it falls back to spatial so the map never blanks. Returns ``(xy, basis_used)``.
    """
    if basis == "umap":
        if "X_umap" not in a.obsm and "umap" in cap.REGISTRY:
            cap.run(a, "umap", {}, _ctx(s))     # lazy embed - the deferred, off-critical-path step
        if "X_umap" in a.obsm:
            return np.asarray(a.obsm["X_umap"]), "umap"
    return np.asarray(a.obsm["spatial"]), "spatial"


@app.get("/api/{sid}/points")
def points(sid: str, color_by: str = "cell_type", max_points: int = 200_000, basis: str = "spatial",
           colors_only: bool = False, only_type: str | None = None, pos_min: float | None = None):
    s = _sess(sid); a = s["adata"]
    if color_by not in a.obs:
        fields = _obs_fields(a)
        color_by = "cell_type" if "cell_type" in a.obs else (fields[0] if fields else None)
    idx = plots._downsample(a.n_obs, max_points)
    xy, basis_used = _basis_coords(a, basis, s)
    if basis_used == "umap":
        # Low-signal cells are held out of the embedding (NaN UMAP coords) - drop them from the UMAP
        # view so the map fits to the real manifold instead of scattering them into a phantom disc.
        fin = np.isfinite(xy[:, 0]) & np.isfinite(xy[:, 1])
        if not fin[idx].all():
            idx = idx[fin[idx]]
    if color_by is None or color_by not in a.obs:
        # A freshly-loaded raw section with nothing to colour by yet: render a neutral monochrome
        # specimen (still shows the tissue) instead of 500-ing. Coloring appears once a step runs.
        rgb = np.tile(np.array([124, 108, 196], dtype=np.uint8), len(idx)).tolist()
        legend: list = []
        color_by = None
    else:
        rgb, legend = _colours(a, color_by, idx)
    # Restrict the colouring to ONE cell type (hover a states-heatmap cell -> the state's per-cell
    # score, but only on that cell type's cells; the rest fade to the same 14% the canvas uses for a
    # de-emphasised cell). Applied to whatever `rgb` holds so it works for a score ramp or a category.
    if only_type and "cell_type" in a.obs:
        m = a.obs["cell_type"].astype(str).to_numpy()[idx] != str(only_type)
        arr = np.asarray(rgb, dtype=np.uint8).reshape(-1, 3)
        arr[m] = (arr[m] * 0.14 + np.array([8, 9, 13], dtype=float)).astype(np.uint8)
        rgb = arr.reshape(-1).tolist()
    # pos_min: fade cells at/below a threshold on the (numeric) colour field, so a hover shows only the
    # cells POSITIVE for that criterion - e.g. the x% a malignant caller actually called (score >
    # threshold) - and the rest dim out, instead of a full-section gradient.
    if pos_min is not None and color_by in a.obs and np.issubdtype(a.obs[color_by].dtype, np.number):
        below = np.asarray(a.obs[color_by], dtype=float)[idx] <= float(pos_min)
        arr = np.asarray(rgb, dtype=np.uint8).reshape(-1, 3)
        arr[below] = (arr[below] * 0.14 + np.array([8, 9, 13], dtype=float)).astype(np.uint8)
        rgb = arr.reshape(-1).tolist()
    # A continuous field (legend empty, numeric column) is a featureplot: hand the frontend the exact
    # ramp stops it was painted with, so the FeatureScale key matches (magma vs violet), never guesses.
    ramp = (_ramp_hex(color_by) if (not legend and color_by in a.obs
                                    and np.issubdtype(a.obs[color_by].dtype, np.number)) else None)
    out = {"rgb": rgb, "n": int(len(idx)), "color_by": color_by, "legend": legend,
           "basis": basis_used, "ramp": ramp}
    # A pure RECOLOUR (hover preview / dropdown): the point positions + ids are unchanged - idx is
    # deterministic for a given (n, basis, max_points) - so the client reuses its x/y/ids and we skip
    # serialising them. x+y+ids are ~2/3 of the payload, and json.dumps of the 141k-point lists is the
    # bulk of the ~300ms response, so this roughly halves recolour latency. The client verifies `n`
    # matches its current points before merging (else it falls back to a full fetch).
    if colors_only:
        return _fast_json(out)
    # float32 coords (microns span ~0-25000 with sub-micron detail, well within float32's ~7 sig figs)
    # halve the x/y payload vs float64 string reprs; keeping them as numpy arrays lets orjson serialize
    # them in C. The client sees the same JSON arrays of numbers.
    return _fast_json({"x": xy[idx, 0].astype(np.float32), "y": xy[idx, 1].astype(np.float32),
                       "ids": np.asarray(idx, dtype=np.uint32), **out})


class RunReq(BaseModel):
    params: dict | None = None


@app.post("/api/{sid}/run/{step}")
def run_step(sid: str, step: str, req: RunReq | None = None):
    s = _sess(sid); a = s["adata"]
    _require_idle(s)
    cname = STEP_CAP.get(step)
    if not cname:
        raise HTTPException(400, f"unknown step {step}")
    ctx = _ctx(s)
    prog = _progress(s)
    ctx.progress = lambda f, l: prog.update(frac=f, label=l)   # engine's ctx.tick calls this
    prog.update(running=True, step=step, frac=0.0, label="")
    try:
        out = cap.run(a, cname, (req.params if req else None) or {}, ctx)
        if getattr(out, "record", None):
            s["log"].append(out.record)
        if step == "cluster" and out.ok:
            try:                                   # eager per-section UMAP precompute, off the critical path
                from spatialscribe.analysis import preload
                preload.start(a, "umap", ctx)
            except Exception:
                pass
        if step == "annotate" and out.ok:
            try:
                _confidence_pass(a, ctx)   # pmp + spatial-coherence, so React confidence == Streamlit/funnel
            except Exception:
                pass
        if out.ok:                             # record that the USER ran this step this session
            _mark_ran(s, cname)
        return {"ok": bool(out.ok), "error": (out.error or None), "state": _state(s)}
    finally:
        prog.update(running=False, frac=1.0)   # clear `running` even if the step raised


class CapReq(BaseModel):
    params: dict | None = None


@app.post("/api/{sid}/run_cap/{name}")
def run_cap(sid: str, name: str, req: CapReq | None = None):
    """Run ANY registered capability by name (params passed through), returning its value plus any
    figure/map artifacts it emitted. This is how the wizard panels reach the full scientific surface
    (qc_funnel, immune_exclusion, neighborhood_enrichment, state_by_celltype, malignant_score,
    discover_programs, subcluster, rejection_reasons, annotation_methods, and the display views)."""
    s = _sess(sid); a = s["adata"]
    _require_idle(s)
    if name not in cap.REGISTRY:
        raise HTTPException(400, f"unknown capability {name}")
    ctx = _ctx(s)
    prog = _progress(s)
    ctx.progress = lambda f, l: prog.update(frac=f, label=l)   # engine's ctx.tick calls this
    prog.update(running=True, step=name, frac=0.0, label="")
    try:
        out = cap.run(a, name, (req.params if req else None) or {}, ctx)
        if getattr(out, "record", None):
            s["log"].append(out.record)
        if name == "annotate" and out.ok:
            try:
                _confidence_pass(a, ctx)
            except Exception:
                pass
        if out.ok:                             # a panel that runs a rail-step cap advances the rail too
            _mark_ran(s, name)
        return {"ok": bool(out.ok), "error": (out.error or None), "value": _clean(out.value),
                "artifacts": _serialize_artifacts(ctx.artifacts), "state": _state(s)}
    finally:
        prog.update(running=False, frac=1.0)   # clear `running` even if the capability raised


@app.post("/api/{sid}/ran/{step}")
def mark_step_ran(sid: str, step: str):
    """Green a rail STEP once its panel finished the first-visit auto-compute. The QC/Panel-check
    panels compute via sub-caps (qc_funnel) that ``_mark_ran`` skips, or via a GET (panel_report)
    whose response the SPA never folds into its state - so the rail needs this explicit nudge.
    Idempotent; only the five compute steps are markable. Returns fresh state so the SPA greens it."""
    s = _sess(sid)
    if step in STEP_CAP:
        s.setdefault("ran", [])
        if step not in s["ran"]:
            s["ran"].append(step)
    return _state(s)


@app.get("/api/{sid}/progress")
def get_progress(sid: str):
    """Coarse progress of the step currently running in this session (see run_step / run_cap). Poll it
    from the UI while a step runs; both run endpoints are sync ``def`` so this GET is served from the
    threadpool and answerable mid-run. ``frac`` is 0.0..1.0; ``label`` is the engine's checkpoint text
    (empty until the first ctx.tick, so the UI can keep its own placeholder)."""
    return _progress(_sess(sid))


class PipelineReq(BaseModel):
    tumour: bool | None = None   # the malignant gate; None keeps the session's current answer
    rctd: bool = False
    split: bool = False
    # Match the interactive `cluster` cap default (0.5): a coarse first pass yields a few clean clusters
    # that annotate cleanly. Higher resolutions fragment (e.g. a targeted panel's neurons into many
    # unresolvable sub-clusters), so the user can still raise it on the Cluster tab. (Subcluster stays 0.1.)
    resolution: float = 0.5


def _mark_pipeline_ran(s: dict, a) -> None:
    """Green the rail for the steps a full "Run full analysis" produced (it bypasses run_step/run_cap,
    which are what normally mark `ran`). Data-derived: a step's output being present == it ran this run.

    The `annotate` rail step means "the section is typed", true for ANY annotate route (consensus
    `annotate`, reference_transfer, or RCTD). is_done("annotate") only fires for the `annotate` cap's
    full key set (confidence/verdict/reason), so a reference-transfer route (e.g. a mouse-brain section
    typed from an Allen atlas) would leave the step - and, downstream, the Report step (which needs every
    compute step) - grey even though cells WERE typed. Green it off the outcome (`cell_type` present)."""
    for _cap in STEP_CAP.values():
        if status.is_done(a, _cap):
            _mark_ran(s, _cap)
    if "cell_type" in a.obs and "annotate" not in s.get("ran", []):
        s.setdefault("ran", []).append("annotate")


def _run_pipeline_thread(sid: str, opts_kwargs: dict) -> None:
    """Run the full spine in the background and record its outcome on the session. Never raises."""
    s = _SESSIONS.get(sid)
    if s is None:
        return
    from spatialscribe.analysis import pipeline as _pl
    a = s["adata"]
    ctx = _ctx(s)
    prog = _progress(s)
    ctx.progress = lambda f, l: prog.update(frac=f, label=l)   # engine's ctx.tick -> pollable progress
    prog.update(running=True, step="pipeline", frac=0.0, label="starting")
    try:
        res = _pl.run_pipeline(a, ctx, _pl.PipelineOptions(export=False, **opts_kwargs))
        s["log"].extend(res.get("action_log", []))
        s["pipeline"] = {"done": True, "error": None, "route": res["route"],
                         "stages": res["stages"], "summary": res["summary"]}
        # The SPA reloads state via onMutate when the poll reports done, so the rail greens then.
        _mark_pipeline_ran(s, a)
    except Exception as exc:   # noqa: BLE001 - a background job must record its failure, never crash
        s["pipeline"] = {"done": True, "error": str(exc), "route": None, "stages": [], "summary": {}}
    finally:
        s["pipeline_running"] = False
        prog.update(running=False, frac=1.0)


@app.post("/api/{sid}/run_pipeline", status_code=202)
def run_pipeline_ep(sid: str, req: PipelineReq | None = None):
    """Run the FULL analysis spine (analysis.pipeline.run_pipeline) as a BACKGROUND job on this
    session, so the app becomes a viewer of precomputed results. Returns 202 immediately; poll
    GET .../pipeline_status for progress + the per-stage outcome, then reload state when done."""
    s = _sess(sid)
    if s.get("pipeline_running"):
        raise HTTPException(409, "a full analysis is already running on this session")
    req = req or PipelineReq()
    if req.tumour is not None:
        s["is_tumour"] = req.tumour           # apply the Data-step checkbox before the run
    s["pipeline_running"] = True
    s["pipeline"] = {"done": False, "error": None, "route": None, "stages": [], "summary": {}}
    opts = {"rctd": req.rctd, "split": req.split, "resolution": req.resolution}
    threading.Thread(target=_run_pipeline_thread, args=(sid, opts), daemon=True).start()
    return {"started": True}


@app.get("/api/{sid}/pipeline_status")
def pipeline_status_ep(sid: str):
    """Progress + outcome of the background full-analysis run (see POST run_pipeline)."""
    s = _sess(sid)
    prog = _progress(s)
    pl = s.get("pipeline") or {"done": False, "error": None, "route": None, "stages": [], "summary": {}}
    return {"running": bool(s.get("pipeline_running")), "frac": prog.get("frac", 0.0),
            "label": prog.get("label", ""), "stage": prog.get("step"),
            "done": bool(pl["done"]), "error": pl["error"], "route": pl["route"],
            "stages": pl["stages"], "summary": pl.get("summary", {}),
            "state": _state(s) if pl["done"] and not s.get("pipeline_running") else None}


class IdxReq(BaseModel):
    indices: list[int]


@app.post("/api/{sid}/region_qc")
def region_qc_ep(sid: str, req: IdxReq):
    """QC summary for a lasso/box-selected region (H1) vs the whole section."""
    s = _sess(sid); a = s["adata"]
    from spatialscribe.analysis import qc as _qc
    idx = [int(i) for i in req.indices if 0 <= int(i) < a.n_obs]
    if not idx:
        raise HTTPException(400, "empty region selection")
    if "total_counts" not in a.obs:   # QC not computed yet (e.g. lasso before the QC step) -> ensure it
        try:
            cap.run(a, "compute_qc", {}, _ctx(s))
        except Exception:
            pass
    return _clean({"n": len(idx), "region": _qc.region_qc(a, idx), "section": _qc.qc_summary(a)})


class RegionFilterReq(BaseModel):
    indices: list[int]
    mode: str = "exclude"   # "exclude" the selected cells, or "keep" only them


@app.post("/api/{sid}/region_filter")
def region_filter_ep(sid: str, req: RegionFilterReq):
    """Region filtering: drop (exclude) or crop-to (keep) the selected cells, in place on the session."""
    s = _sess(sid); a = s["adata"]
    mask = np.zeros(a.n_obs, dtype=bool)
    for i in req.indices:
        if 0 <= int(i) < a.n_obs:
            mask[int(i)] = True
    keep = mask if req.mode == "keep" else ~mask
    n_keep = int(keep.sum())
    if n_keep == 0 or n_keep == a.n_obs:
        raise HTTPException(400, f"filter would keep {n_keep}/{a.n_obs} cells (no-op or empties the section)")
    s["adata"] = a[keep].copy()
    s["log"].append({"name": "region_filter", "params": {"mode": req.mode, "n_removed": int((~keep).sum())}})
    return {"session_id": sid, **_state(s)}


class RenameReq(BaseModel):
    old: str
    new: str


@app.post("/api/{sid}/rename_celltype")
def rename_ep(sid: str, req: RenameReq):
    """Rename a cell type (reuse an existing name to MERGE), then recompute the confidence pass."""
    s = _sess(sid); a = s["adata"]
    if "cell_type" not in a.obs:
        raise HTTPException(400, "run annotate first")
    ct = a.obs["cell_type"].astype(str)
    ct[ct == req.old] = req.new
    a.obs["cell_type"] = ct.astype("category")
    a.obs.drop(columns=["annotation_confidence", "annotation_verdict", "annotation_reason",
                        "cell_type_final", "pmp", "spatial_coherence"], errors="ignore", inplace=True)
    try:
        _confidence_pass(a, _ctx(s))
    except Exception:
        pass
    s["log"].append({"name": "rename_celltype", "params": {"old": req.old, "new": req.new}})
    return {"session_id": sid, **_state(s)}


def _merge_nudge(pc: dict) -> str:
    groups = [g for g in (pc.get("merge_groups") or []) if len(g) > 1]
    if not groups:
        return "No panel-driven merges suggested - every resolvable type has a private on-panel marker."
    return ("Consider merging these groups - the panel has no private marker to tell them apart: "
            + "; ".join(" + ".join(g) for g in groups) + ".")


# Abstention labels (annotate.apply_confidence) are not cell types - never ask the LLM for markers
# for them, and never let their presence trigger the LLM path. Re-exported from the ONE definition in
# analysis.annotate; keeping a second private copy here is how the two vocabularies drifted apart.
from spatialscribe.analysis import annotate as _annotate  # noqa: E402

_ABSTENTION_PREFIXES = _annotate.ABSTENTION_PREFIXES


def _apply_llm_markers_for_categories(a, ctx) -> None:
    """Panel-check only: AUGMENT the per-type marker sets for the section's CURRENT annotated
    categories with CellGuide's Cell-Ontology-grounded canonical markers, unioned ONTO the curated
    baseline, and pin the result on ``ctx.marker_sets`` so coverage / identifiability / typability
    reflect the labels actually on the section (incl. renamed / subclustered / reference-transferred
    types). No-op when unannotated. The ctx is per-request (``_ctx``) - never leaks to other steps.

    **Why augment, not just gap-fill.** The curated dictionaries list only a handful of markers per
    lineage (Mast 4, NK 4, Endothelial 4), so a rich panel reads 'not typable' even when it carries
    good markers the curated set never named. CellGuide knows those literature markers (Mast +HPGDS/
    KIT/GATA2, NK +GZMA/KLRF1/EOMES, Endothelial +FLT1/CD34/SOX17); unioning them lets more land on
    the panel -> higher identifiability AUC and more typable lineages. CellGuide is a STATIC, versioned,
    key-free DB, so this stays deterministic across runs (the old worry was an LLM's per-run marker
    choices, not a fixed DB). Ordering keeps curated flagship genes FIRST, then CellGuide extras; the
    LLM is a last resort that only fires for a label CellGuide cannot resolve (truly novel names), so
    known-lineage scoring never depends on a model's choices.
    """
    if ctx.marker_sets is not None or "cell_type" not in a.obs:
        return
    from spatialscribe.analysis import markers as _markers
    cats = [c for c in map(str, a.obs["cell_type"].astype("category").cat.categories)
            if not _annotate.is_abstention(c)]
    if not cats:
        return
    curated = _markers.for_tissue(ctx.tissue)
    ground = _markers.markers_for_types(cats, a.var_names, ctx.tissue, use_llm=llm.available())
    merged: dict[str, list[str]] = {}
    for c in cats:
        base = list(curated.get(c, []))
        seen = {g.upper() for g in base}
        combined = base + [g for g in ground.get(c, []) if g.upper() not in seen]
        if combined:
            merged[c] = combined
    ctx.marker_sets = merged or None


def _typability(a, ctx, s) -> list[dict]:
    """Computed per-type 'confidently typable?' table - depth-matched F1 when a reference is loaded,
    else one-vs-rest identifiability AUC, else marker coverage (flagged weak). Never raises."""
    from spatialscribe.analysis import panel_check as _pc
    ref_match = None
    if s.get("reference") is not None:
        try:
            cap.run(a, "reference_match", {}, ctx)
            ref_match = a.uns.get("reference_match")
        except Exception:
            ref_match = None
    try:
        return _pc.typability_table(a, cluster_key="cell_type", marker_sets=ctx.markers(),
                                    reference_match=ref_match, panel_check_result=a.uns.get("panel_check"))
    except Exception:
        return []


@app.get("/api/{sid}/panel_report")
def panel_report_ep(sid: str):
    """Comprehensive panel adequacy + reference-marker concordance (global + per cell type):
    marker coverage traffic lights, one-vs-rest identifiability AUC, confusable pairs, merge nudges."""
    s = _sess(sid); a = s["adata"]; ctx = _ctx(s)
    from spatialscribe.analysis import panel_check as _pc
    _apply_llm_markers_for_categories(a, ctx)
    cap.run(a, "panel_check", {}, ctx)
    _mark_ran(s, "panel_check")             # opening the Panel step runs it -> advance the rail
    pc = a.uns.get("panel_check") or {}
    cov = pc.get("coverage", {}) or {}
    auc: dict = {}
    if "cell_type" in a.obs:
        try:
            auc = _pc.identifiability_auc(a, cluster_key="cell_type", marker_sets=ctx.markers())
        except Exception:
            auc = {}
    per_type = []
    for ct, d in cov.items():
        if not isinstance(d, dict):
            continue
        per_type.append({
            "cell_type": ct, "status": d.get("status"),
            "n_present": d.get("n_present"), "n_markers": d.get("n_markers"),
            "coverage_frac": round(float(d.get("coverage_frac", 0) or 0), 3),
            "missing": (d.get("missing") or [])[:6],
            "auc": (auc.get(ct, {}) or {}).get("auc"),
        })
    # Computed per-type "confidently typable?" decision (F1 / AUC, not gene presence). Merge its
    # verdict onto each coverage row so the table can show a computed can/can't-type column.
    typ = {t["cell_type"]: t for t in _typability(a, ctx, s)}
    from spatialscribe.analysis.panel_check import MIN_TYPABLE_CELLS as _MIN_TYP
    for r in per_type:
        t = typ.get(r["cell_type"], {})
        r["confidently_typable"] = t.get("confidently_typable")
        r["typability_basis"] = t.get("basis")
        r["typability_score"] = t.get("score")
        r["n_cells"] = t.get("n_cells")
        # WHY a type is not typable, so the table explains 'no' instead of contradicting the markers /
        # AUC columns - the most common surprise is a high AUC on very few cells (deliberately discounted
        # as noise; see panel_check.typability_table), which otherwise looks like "good AUC yet not typable".
        reason = None
        if t.get("confidently_typable") is False:
            nc = t.get("n_cells")
            if t.get("zero_on_panel"):
                reason = "none of its markers are on this panel"
            elif nc is not None and nc < _MIN_TYP:
                reason = f"too few cells (n={nc} < {_MIN_TYP})"
            elif t.get("basis") == "identifiability_auc" and t.get("score") is not None \
                    and t.get("threshold") is not None and t["score"] < t["threshold"]:
                reason = f"identifiability below {t['threshold']}"
            elif t.get("basis") == "depth_matched_f1" and t.get("confused_with"):
                reason = f"confused with {t['confused_with']}"
            elif t.get("basis") == "coverage_only":
                reason = "too few markers on this panel"
        r["not_typable_reason"] = reason
    per_type.sort(key=lambda r: (r["auc"] is None, r["auc"] if r["auc"] is not None else 0))
    n_present_tot = sum(int(d.get("n_present") or 0) for d in cov.values() if isinstance(d, dict))
    n_markers_tot = sum(int(d.get("n_markers") or 0) for d in cov.values() if isinstance(d, dict))
    aucs = [v.get("auc") for v in auc.values() if isinstance(v, dict) and v.get("auc") is not None]
    typ_vals = list(typ.values())
    bases = {t.get("basis") for t in typ_vals}
    glob = {
        "n_types": len(cov),
        "overall_marker_coverage": round(n_present_tot / n_markers_tot, 3) if n_markers_tot else 0.0,
        "n_resolvable": sum(1 for d in cov.values() if isinstance(d, dict) and d.get("status") == "green"),
        "n_weak": sum(1 for d in cov.values() if isinstance(d, dict) and d.get("status") == "amber"),
        "n_cannot": sum(1 for d in cov.values() if isinstance(d, dict) and d.get("status") == "red"),
        # the cell types behind each marker-coverage count, so the UI can list them on hover instead of
        # spending vertical space. Same `cov` source as the counts above, so they always agree.
        "resolvable_types": sorted(ct for ct, d in cov.items() if isinstance(d, dict) and d.get("status") == "green"),
        "weak_types": sorted(ct for ct, d in cov.items() if isinstance(d, dict) and d.get("status") == "amber"),
        "cannot_types": sorted(ct for ct, d in cov.items() if isinstance(d, dict) and d.get("status") == "red"),
        "mean_identifiability_auc": round(sum(aucs) / len(aucs), 3) if aucs else None,
        # the headline computed decision: how many types are confidently typable, and on what basis.
        # Count only STRONG evidence (weak_evidence=False) in the headline, so a gene-presence-only or
        # tiny-n row is not tallied identically to a computed AUC/F1 separability (WTA panels read
        # "10/10" while half the rows were absent lineages scored on marker presence).
        "n_confidently_typable": sum(1 for t in typ_vals
                                     if t.get("confidently_typable") and not t.get("weak_evidence")),
        "n_confidently_typable_weak": sum(1 for t in typ_vals
                                          if t.get("confidently_typable") and t.get("weak_evidence")),
        "n_typability_assessed": len(typ_vals),
        "typability_basis": ("depth_matched_f1" if "depth_matched_f1" in bases
                             else "identifiability_auc" if "identifiability_auc" in bases
                             else "coverage_only"),
        "n_confusable_pairs": len(pc.get("confusable_pairs", []) or []),
        "marker_db": pc.get("marker_db"),
    }
    return _clean({"global": glob, "per_type": per_type,
                   "confusable": [p.get("pair") for p in (pc.get("confusable_pairs") or []) if isinstance(p, dict)],
                   "merge_nudge": _merge_nudge(pc)})


@app.get("/api/{sid}/panel_verdict")
def panel_verdict_ep(sid: str):
    """Auto LLM plain-language panel verdict (Haiku - cheap, prompt-cached system prompt)."""
    s = _sess(sid); a = s["adata"]
    if not llm.available():
        return {"verdict": None, "note": "No AI model connected (set ANTHROPIC_API_KEY or SPATIALSCRIBE_LLM_BASE_URL)"}
    from spatialscribe.analysis import llm as _llm
    ctx = _ctx(s)
    _apply_llm_markers_for_categories(a, ctx)   # same LLM-per-category markers as the panel report
    pc = a.uns.get("panel_check")
    if not pc:
        cap.run(a, "panel_check", {}, ctx); pc = a.uns.get("panel_check")
    typ = _typability(a, ctx, s)                # computed F1/AUC decision - grounds the verdict
    try:
        verdict = _verdict_cache(s, "panel", {"pc": pc, "typ": typ},
                                 lambda: _llm.panel_verdict(pc, typability=typ))
        return {"verdict": verdict,
                "typability_basis": (typ[0]["basis"] if typ else None)}
    except Exception as exc:
        return {"verdict": None, "note": f"verdict failed: {type(exc).__name__}: {exc}"}


@app.get("/api/{sid}/qc_verdict")
def qc_verdict_ep(sid: str):
    """Auto LLM plain-language QC verdict: can this section's signal be trusted for annotation?
    Grounded in the six-layer funnel + top rejection reasons + (if a reference is loaded) the
    reference-to-dataset match. Haiku by default, prompt-cached (cheap)."""
    s = _sess(sid); a = s["adata"]
    ctx = _ctx(s)
    cap.ensure(a, "compute_qc", ctx)   # ground the verdict on real QC metrics even if the funnel fails
    funnel = cap.run(a, "qc_funnel", {}, ctx)
    headline = funnel.value if funnel.ok else {}
    rej = None
    if "annotation_verdict" in a.obs:                          # rejection needs annotation to have run
        rr = cap.run(a, "rejection_reasons", {}, ctx)
        rej = (rr.value or {}).get("breakdown") if rr.ok else None
    ref_match = None
    if s.get("reference") is not None:                         # fold the reference-to-dataset match in
        rm = cap.run(a, "reference_match", {}, ctx)
        if rm.ok:
            ref_match = {k: rm.value.get(k) for k in ("match_score", "recommendation", "global")}
    if not llm.available():
        return _clean({"verdict": None, "note": "No AI model connected (set ANTHROPIC_API_KEY or SPATIALSCRIBE_LLM_BASE_URL)",
                       "funnel": headline, "reference_match": ref_match})
    from spatialscribe.analysis import llm as _llm
    payload = dict(headline)
    if ref_match:
        payload["reference_to_dataset_match"] = ref_match
    section_metrics = _summary(s).get("qc")
    try:
        verdict = _verdict_cache(s, "qc", {"payload": payload, "rej": rej, "metrics": section_metrics},
                                 lambda: _llm.qc_verdict(payload, rejection_breakdown=rej,
                                                         section_metrics=section_metrics))
        return _clean({"verdict": verdict, "reference_match": ref_match})
    except Exception as exc:
        return _clean({"verdict": None, "note": f"verdict failed: {type(exc).__name__}: {exc}",
                       "reference_match": ref_match})


@app.get("/api/{sid}/verify_report")
def verify_report_ep(sid: str, neighborhood: bool = False):
    """Curated self-verify payload for the dedicated React panel: per-type marker-argmax agreement +
    one-vs-rest AUC + pass/fail, the failing types with their cause, and the advisory remediation."""
    s = _sess(sid); a = s["adata"]
    if "cell_type" not in a.obs:
        return {"status": "no_annotation", "note": "Run Annotate first - self-verify audits the labels."}
    out = cap.run(a, "self_verify", {"neighborhood": bool(neighborhood)}, _ctx(s))
    if not out.ok:
        return _clean({"status": "error", "error": out.error})
    v = out.value or {}
    per_type = []
    for ct, d in (v.get("per_type", {}) or {}).items():
        if not isinstance(d, dict):
            continue
        per_type.append({"cell_type": ct, "status": d.get("status"),
                         "argmax_agreement": d.get("argmax_agreement"), "auc": d.get("auc"),
                         "confused_with": d.get("confused_with"), "cause": d.get("cause"),
                         "n_cells": d.get("n_cells")})
    per_type.sort(key=lambda r: (r["status"] != "fail",
                                 r["argmax_agreement"] if r["argmax_agreement"] is not None else 1.0))
    return _clean({"status": "ok", "section_agreement": v.get("section_agreement"),
                   "mean_auc": v.get("mean_auc"), "n_types_scored": v.get("n_types_scored"),
                   "failed": v.get("failed", []), "per_type": per_type,
                   "suggestions": v.get("suggestions", [])})


@app.post("/api/load_synthetic")
def load_synthetic():
    """Load the synthetic melanoma demo (a melanoma tissue) - matches Streamlit's option."""
    from spatialscribe.analysis import demo as _demo
    sample = _demo.load_demo()
    a = sample.adata
    a.uns.setdefault("platform", getattr(sample, "platform", "xenium"))
    sid = uuid.uuid4().hex[:12]
    _SESSIONS[sid] = {"adata": a, "tissue": "melanoma", "log": []}
    return {"session_id": sid, **_state(_SESSIONS[sid])}


class CopilotReq(BaseModel):
    prompt: str


def _commit_load(s: dict, loaded: list) -> None:
    """Swap a section the copilot's ``load_section`` loaded onto the session (adata + inferred tissue
    + the auto-selected reference), so the next request runs on it. No-op when nothing was loaded. The
    load's provenance record is already appended to ``s['log']`` by the copilot loop, so this does not
    re-log it; it just clears the previous section's derived caches (LLM verdicts)."""
    if not loaded:
        return
    latest = loaded[-1]
    s["adata"] = latest["adata"]
    s["tissue"] = latest.get("tissue") or s.get("tissue")
    if latest.get("reference") is not None:
        s["reference"] = latest["reference"]
        s["ref_label_key"] = latest.get("ref_label_key")
    s["source_path"] = latest.get("path")
    s.pop("_verdicts", None)
    # A new section = a fresh run: clear the steps the user "ran" on the old one so the rail greens
    # reset to only Load (a manual load already gets this free via a brand-new session; the copilot
    # reuses the session, so reset it here). `done` is data-derived and recomputes from the new adata.
    s["ran"] = []


# Annotation columns whose CONTENT change means the copilot mutated the section in place (a merge /
# relabel / subcluster / cell-state / niche run), so the SPA must repaint the map + drop cached panels.
# Materialising a gene into obs for colouring is NOT here, so a recolour never reads as a mutation.
_ANNOTATION_STATE_COLS = ("cell_type", "cell_type_final", "subtype", "cell_state", "niche",
                          "annotation_verdict")


def _annotation_sig(a) -> tuple:
    """A cheap signature of the section's annotation state - changes iff a copilot capability mutated
    the labels/subtypes/states/niches in place (used to tell the SPA to refresh)."""
    parts = []
    obs = getattr(a, "obs", None)
    if obs is None:
        return ()
    for c in _ANNOTATION_STATE_COLS:
        if c in obs.columns:
            vc = obs[c].astype(str).value_counts()
            parts.append((c, tuple(sorted((str(k), int(v)) for k, v in vc.items()))))
    return tuple(parts)


@app.post("/api/{sid}/copilot")
def copilot(sid: str, req: CopilotReq):
    s = _sess(sid); a = s["adata"]
    if not llm.available():
        raise HTTPException(400, "no LLM configured (set ANTHROPIC_API_KEY or SPATIALSCRIBE_LLM_BASE_URL)")
    from spatialscribe.agent.tools import run_copilot
    arts: list = []
    loaded: list = []
    sig_before = _annotation_sig(a)
    reply = run_copilot(a, req.prompt, s["tissue"], action_log=s["log"], use_llm=True, artifacts=arts,
                        loaded=loaded)
    _commit_load(s, loaded)
    mutated = (not loaded) and (_annotation_sig(s["adata"]) != sig_before)
    mv = next((x for x in reversed(arts) if x.get("kind") == "map_view"), None)
    figs = _serialize_artifacts([x for x in arts if x.get("kind") == "figure"])
    return {"reply": reply,
            "map_view": ({"color_by": mv["color_by"], "highlight": mv.get("highlight")} if mv else None),
            "figures": figs, "mutated": bool(mutated), "state": _state(s)}


@app.get("/api/{sid}/summary")
def summary(sid: str):
    return _clean(_summary(_sess(sid)))


def _summary(s: dict) -> dict:
    """Per-step panel data computed from the current adata state (only what's present)."""
    a = s["adata"]
    out: dict = {}
    pc = a.uns.get("panel_check")
    if isinstance(pc, dict) and "coverage" in pc:
        cov = pc["coverage"]
        cp = pc.get("confusable_pairs", [])
        out["panel"] = {
            "resolvable": [ct for ct, d in cov.items() if d.get("status") == "green"],
            "weak": [ct for ct, d in cov.items() if d.get("status") in ("amber", "yellow")],
            "cannot": [ct for ct, d in cov.items() if d.get("status") == "red"],
            "confusable": [p.get("pair") for p in cp] if isinstance(cp, list) else [],
        }
    try:
        from spatialscribe.analysis import qc as _qc
        q = _qc.qc_summary(a)
        out["qc"] = {"n_cells": int(a.n_obs),
                     "median_genes": round(float(q.get("median_genes_per_cell", 0)), 1),
                     "median_counts": round(float(q.get("median_transcripts_per_cell", 0)), 1),
                     "empty_frac": round(float(q.get("fraction_empty_cells", 0)), 4),
                     # surfaced so the neg-control flag chip has a visible number + threshold behind it
                     "pct_counts_control": round(float(q.get("pct_counts_control", 0)), 2)}
    except Exception:
        pass
    if "leiden" in a.obs:
        # `resolution` is what the cluster step stamped on uns; absent for a cache clustered before
        # this was recorded, and the panel then just keeps its own slider value.
        out["cluster"] = {"n_clusters": int(a.obs["leiden"].nunique()),
                          "resolution": (float(a.uns["leiden_resolution"])
                                         if "leiden_resolution" in a.uns else None)}
    ann: dict = {}
    if "annotation_verdict" in a.obs:
        v = a.obs["annotation_verdict"].value_counts(); n = max(1, int(a.n_obs))
        ann = {"confident_pct": round(100 * int(v.get("PASS", 0)) / n),
               "tentative_pct": round(100 * int(v.get("WARN", 0)) / n),
               "abstained_pct": round(100 * int(v.get("FAIL", 0)) / n)}
        if "rejection_reason" in a.obs:
            rr = a.obs["rejection_reason"].astype(str).value_counts().head(8)
            ann["rejections"] = [{"reason": str(k), "count": int(x)} for k, x in rr.items() if k and k != "nan"]
    if "cell_type" in a.obs:
        ct = a.obs["cell_type"].astype(str).value_counts().head(14)
        ann["cell_types"] = [{"name": str(k), "count": int(x)} for k, x in ct.items()]
    if ann:
        out["annotate"] = ann
    if "niche" in a.obs:
        nc = a.obs["niche"].astype(str).value_counts()
        out["spatial"] = {"n_niches": int(nc.size),
                          "niches": [{"id": str(k), "count": int(x)} for k, x in nc.head(35).items()]}
    out["report"] = {"log": [r for r in s["log"][-30:]]}
    return out


@app.get("/api/{sid}/export/script")
def export_script_ep(sid: str):
    s = _sess(sid)
    from spatialscribe.analysis import export as _ex
    body = _ex.build_runnable_script(s["log"], adata=s["adata"], tissue=s["tissue"],
                                     source_path=s.get("source_path"))
    return PlainTextResponse(body,
                             headers={"Content-Disposition": 'attachment; filename="analysis.py"'})


@app.get("/api/{sid}/export/h5ad")
def export_h5ad_ep(sid: str):
    s = _sess(sid)
    import os as _os
    import tempfile
    from spatialscribe.analysis import export as _ex
    p = _os.path.join(tempfile.mkdtemp(prefix="ss_export_"), "annotated.h5ad")
    _ex.export_h5ad(s["adata"], p)
    return FileResponse(p, filename="annotated.h5ad", media_type="application/octet-stream")


@app.get("/api/{sid}/export/report")
def export_report_ep(sid: str):
    # The self-contained HTML report lives in backend/report.py, kept out of this thin HTTP layer.
    # It needs app-side helpers (summary / colours / obs-fields / device); pass them in so report.py
    # does not import app.py back (circular).
    return HTMLResponse(
        report_html(_sess(sid), summary=_summary, colours=_colours,
                    obs_fields=_obs_fields, device=_device),
        headers={"Content-Disposition": 'attachment; filename="spatialscribe_report.html"'})


@app.post("/api/{sid}/copilot/stream")
def copilot_stream(sid: str, req: CopilotReq):
    """Server-sent events: a status tick, the map-recolour directive (as soon as it's known), then
    the answer word-streamed for a typing effect, then a done event with fresh state."""
    s = _sess(sid); a = s["adata"]
    if not llm.available():
        raise HTTPException(400, "no LLM configured (set ANTHROPIC_API_KEY or SPATIALSCRIBE_LLM_BASE_URL)")
    import json as _json
    import time as _time

    def sse(obj):
        return "data: " + _json.dumps(obj) + "\n\n"

    def gen():
        yield sse({"type": "status", "text": "thinking …"})
        from spatialscribe.agent.tools import run_copilot
        import queue as _queue
        import threading as _threading
        arts: list = []
        loaded: list = []
        # Run the (blocking) tool-use loop in a thread and stream its per-tool status here, so a slow
        # step (loading a section + selecting a reference is ~10s) shows "loading the section …" /
        # "selecting a reference …" instead of a static "thinking …" for the whole time.
        q: "_queue.Queue" = _queue.Queue()
        result: dict = {}
        sig_before = _annotation_sig(a)

        def _run():
            try:
                result["reply"] = run_copilot(a, req.prompt, s["tissue"], action_log=s["log"],
                                              use_llm=True, artifacts=arts, loaded=loaded,
                                              on_status=lambda m: q.put(m))
            except Exception as exc:                       # never leave the stream hanging
                result["reply"] = f"copilot error: {exc}"
            finally:
                q.put(None)                                # sentinel: the loop finished

        th = _threading.Thread(target=_run, daemon=True)
        th.start()
        while True:
            msg = q.get()
            if msg is None:
                break
            yield sse({"type": "status", "text": msg})
        th.join()
        reply = result.get("reply", "")
        _commit_load(s, loaded)
        mv = next((x for x in reversed(arts) if x.get("kind") == "map_view"), None)
        if mv:
            yield sse({"type": "map_view", "color_by": mv["color_by"], "highlight": mv.get("highlight")})
        figs = _serialize_artifacts([x for x in arts if x.get("kind") == "figure"])
        if figs:
            yield sse({"type": "figures", "figures": figs})
        acc = ""
        for w in (reply or "").split(" "):
            acc = (acc + " " + w).strip()
            yield sse({"type": "token", "text": acc})
            _time.sleep(0.012)   # visible typing (sync generator runs in uvicorn's threadpool)
        # `loaded` tells the SPA the copilot swapped in a NEW section this turn, so it drops the old
        # section's cached panels and repaints the map from the new one (else the demo points linger).
        # `mutated` tells it a capability changed the labels/subtypes/states IN PLACE (merge / relabel /
        # subcluster / self_heal) - the SPA must then invalidate cached panels + repaint the map, or the
        # change is applied server-side but invisible (the live "optimize merging did nothing" report).
        mutated = (not loaded) and (_annotation_sig(s["adata"]) != sig_before)
        yield sse({"type": "done", "state": _state(s), "loaded": bool(loaded), "mutated": bool(mutated)})

    return StreamingResponse(gen(), media_type="text/event-stream")


@app.post("/api/stt")
async def stt_ep(audio: UploadFile = File(...)):
    """Transcribe a dictated question (faster-whisper, on this CPU - the audio never leaves the box).

    Session-free: the recording is all the transcriber needs, so the SPA posts it straight from
    MediaRecorder. Returns ``{"text": ...}``. An empty string means the user said nothing, which the
    UI treats as a no-op rather than an error."""
    if not stt.available():
        raise HTTPException(503, "speech-to-text is not installed on the server (see the internal docs)")
    data = await audio.read()
    if not data:
        raise HTTPException(400, "empty recording")
    try:
        return {"text": stt.transcribe(data)}
    except stt.RecordingTooLarge as e:
        raise HTTPException(413, str(e)) from e
    # Order matters: PyAV's decode errors are ValueError subclasses, so this must come SECOND or a
    # 9-byte junk upload gets reported as "recording too large".
    except Exception as e:                       # undecodable blob: wrong codec, truncated upload
        raise HTTPException(422, f"could not decode the recording: {e}") from e


@app.get("/api/health")
def health():
    return {"ok": True, "sessions": len(_SESSIONS)}


class ProviderReq(BaseModel):
    provider: str


@app.post("/api/llm/provider")
def set_llm_provider(req: ProviderReq):
    """Switch the copilot LLM backend at RUNTIME (process-wide) - the UI's clickable model tag. Only
    switches to a configured provider; returns the resulting llm descriptor."""
    ok = llm.set_provider(req.provider)
    return {"ok": ok, "llm": {"available": llm.available(), "provider": llm.provider(),
                              "model": (llm.default_model() or None), "providers": llm.providers_available()}}


# Serve the built SPA at "/" (single origin) if present - so the whole app is ONE service.
# Mounted last, so the /api/* routes above take precedence. In dev, run vite instead (proxy /api).
from pathlib import Path as _Path
from fastapi.staticfiles import StaticFiles

_DIST = _Path(__file__).resolve().parents[1] / "webapp" / "dist"
if (_DIST / "index.html").exists():
    app.mount("/", StaticFiles(directory=str(_DIST), html=True), name="spa")
