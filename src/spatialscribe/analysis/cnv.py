"""Tumor calling (gap #1) - "which cells are the malignant tumor?"

Two paths:

- `malignant_score` (dependency-free): a marker-based malignant probability from the
  tissue's malignant/tumor signature. Cheap, always available; a first-pass "where is the
  tumor" layer that does NOT rely on CNV.
- `call_malignant_cnv` (`infercnvpy` via an isolated subprocess): the real, marker-independent
  call from copy-number aberration - catches tumor cells that down-regulate lineage markers and
  confirms malignancy by aneuploidy. Runs ``subprocesses/cnv/run_infercnv.py`` in the
  ``cnv_env`` env (infercnvpy will not import in the main interpreter) with Jensen
  neighbor-smoothing + a non-tumor reference, and joins the per-cell CNV burden back. Degrades
  gracefully (returns a 'skipped' status, never raises) so the core stays lean.

The CNV path mirrors the `insitucnv-analysis` skill conventions.
"""

from __future__ import annotations

import os

from . import markers as _m


def malignant_score(adata, tissue: str = "melanoma", key_added: str = "malignant_score") -> dict:
    """Marker-based malignant probability in ``adata.obs[key_added]`` (0-1, min-max).

    Uses the tissue's malignant/tumor lineage signature (melanoma: MLANA/MITF/SOX10…;
    epithelial tumors: EPCAM/KRT8/18…). A cheap, CNV-free "where is the tumor" layer.
    """
    import numpy as np

    sets = _m.for_tissue(tissue)
    mal_key = next((k for k in sets if any(t in k for t in ("Malignant", "Tumor", "Epithelial"))), None)
    genes = _m.on_panel(adata.var_names, sets.get(mal_key, [])) if mal_key else []
    if not genes:
        adata.obs[key_added] = 0.0
        return {"status": "no malignant markers on panel", "n_markers": 0}
    from .backend import get_backend
    get_backend().score_genes(adata, gene_list=genes, score_name="_mal_raw", ctrl_size=min(50, adata.n_vars))
    raw = np.asarray(adata.obs["_mal_raw"], dtype=float)
    lo, hi = np.nanpercentile(raw, [1, 99])
    adata.obs[key_added] = np.clip((raw - lo) / (hi - lo + 1e-9), 0, 1)
    del adata.obs["_mal_raw"]
    return {"status": "ok", "n_markers": len(genes), "markers": genes}


# infercnvpy 0.6.0 + its numba/llvmlite stack will not import in the main SpatialScribe
# interpreter, so CNV runs in an isolated subprocess env. All three are site-specific with NO
# public default: configure them via env vars for a real CNV run. Unset -> call_malignant_cnv
# skips gracefully and the app falls back to the marker-based malignant_score.
_GI_INSITUCNV_PY = os.environ.get("SPATIALSCRIBE_CNV_PYTHON", "")   # a python with infercnvpy
_GI_INSITUCNV_LIB = os.environ.get("SPATIALSCRIBE_CNV_LIB", "")     # its lib dir (LD_LIBRARY_PATH)
_GTF_DEFAULT = os.environ.get("SPATIALSCRIBE_CNV_GTF", "")          # a GTF with gene coordinates


def _join_cnv(adata, parquet_path) -> int:
    """Join a ``{cell_id, cnv_score, is_malignant}`` parquet onto ``obs`` (reindexed). Returns n covered."""
    import pandas as pd

    df = pd.read_parquet(parquet_path)
    df["cell_id"] = df["cell_id"].astype(str)
    aligned = df.set_index("cell_id").reindex(adata.obs_names.astype(str))
    adata.obs["cnv_score"] = aligned["cnv_score"].to_numpy(dtype=float)
    adata.obs["is_malignant"] = aligned["is_malignant"].fillna(False).to_numpy(dtype=bool)
    return int(aligned["cnv_score"].notna().sum())


def _build_reference(adata, normal_key, normal_categories, marker_sets,
                     reference_purity: float, min_reference: int):
    """Build the diploid-reference labels for infercnv. Returns ``(labels, categories)`` where
    ``labels`` is an ndarray of ``'reference'`` / ``'query'``, or ``(None, error_msg)``.

    Reference = cells whose type is a non-tumor category. Optionally restricted to **marker-PURE**
    reference cells (``pmp >= reference_purity``) when ``marker_sets`` is given and enough pure
    cells remain: on a noisy-label section the immune/stromal "reference" is contaminated with
    mislabeled tumour, which flattens the CNV baseline; keeping only high-purity reference cells
    sharpens the tumour-vs-normal separation.
    """
    import numpy as np

    types = adata.obs[normal_key].astype(str)
    if not normal_categories:
        normal_categories = sorted({t for t in types.unique()
                                    if not any(k in t for k in ("Malignant", "Tumor", "Epithelial"))})
    normal_categories = [c for c in normal_categories if (types == c).sum() >= 20]
    if not normal_categories:
        return None, "no usable non-tumor reference (>=20 cells)"

    ref_mask = types.isin(normal_categories).to_numpy()
    if marker_sets:
        from . import purity as _purity
        _purity.pmp(adata, assigned_label_key=normal_key, lineage_markers=marker_sets)
        pmp = np.nan_to_num(np.asarray(adata.obs["pmp"], dtype=float), nan=0.0)
        pure = ref_mask & (pmp >= reference_purity)
        if int(pure.sum()) >= min_reference:               # only if the pure subset is big enough
            ref_mask = pure
    if int(ref_mask.sum()) < 20:
        return None, "reference too small after purity filtering (<20 cells)"
    return np.where(ref_mask, "reference", "query"), normal_categories


_TUMOUR_KEYWORDS = ("melanoma", "tumor", "tumour", "carcinoma", "cancer", "breast",
                    "skin", "uveal", "cervical", "lung", "ovarian", "glioma")
# Non-malignant lineages that make the diploid CNV reference (immune / stromal / endothelial).
_NORMAL_SUBSTRINGS = ("T cell", "T Lymph", "Myeloid", "B/", "B Cell", "Plasma", "Stromal",
                      "CAF", "Fibro", "Endothelial", "Pericyte", "NK", "Mast", "Dendritic", "Macrophage")


def is_tumour_context(tissue: str, is_tumour: bool | None = None) -> tuple[bool, str]:
    """Should the malignant callers run? Returns ``(run, how_it_was_decided)``.

    An EXPLICIT ``is_tumour`` (the Data-step checkbox) always wins. Only when it is ``None`` do we
    fall back to substring-matching ``tissue`` against :data:`_TUMOUR_KEYWORDS`, which is wrong in
    BOTH directions and must never be the only gate. Measured: 'normal breast', 'healthy skin
    biopsy' and 'lung (healthy donor)' all MATCH, so malignant calling runs on normal tissue; while
    'glioblastoma', 'sarcoma', 'lymphoma', 'leukemia', 'neuroblastoma' and 'myeloma' all MISS, so it
    is silently skipped on real tumours. Callers must surface the returned reason, so a report can
    never read as "no malignant cells" when in truth the gate simply never tripped.
    """
    if is_tumour is not None:
        return bool(is_tumour), "explicit"
    return any(t in str(tissue).lower() for t in _TUMOUR_KEYWORDS), "tissue_keyword"


def call_malignant_concordance(adata, tissue: str = "melanoma", marker_sets: dict | None = None,
                               max_cells: int = 25000, cf_threshold: float = 0.5,
                               is_tumour: bool | None = None) -> dict:
    """Run the malignant callers and report per-caller %-malignant + their pairwise concordance.

    Lifts the Streamlit "tumor calling" panel into the shared engine so BOTH frontends (and the
    copilot) reach it. Always runs the cheap marker ``malignant_score``; runs Cancer-Finder on any
    tumour panel and InSituCNV only on a >2000-probe panel (both in isolated envs, both degrade to a
    'skipped: ...' status when unconfigured). Reports the fraction of jointly-scored cells where the
    two LEARNED callers agree. Never raises - returns a JSON-able dict; off a tumour context it
    returns ``{'status': 'skipped: non-tumour context'}``.

    ``is_tumour`` is the Data-step checkbox ("this sample contains malignant cells"). When set it
    overrides the tissue-keyword heuristic entirely (see :func:`is_tumour_context`). The returned
    dict carries ``gate`` so a caller can tell which path decided, and never mistake a gate that
    did not fire for a section with no tumour in it.
    """
    import numpy as np

    from . import cancerfinder as _cf

    tissue = str(tissue)
    n_genes = int((~adata.var["control"]).sum()) if "control" in adata.var.columns else adata.n_vars
    run, gate = is_tumour_context(tissue, is_tumour)
    if not run:
        why = ("the 'contains malignant cells' box is unticked" if gate == "explicit" else
               f"'{tissue}' is not a tumour context; tick 'contains malignant cells' at the Data "
               f"step (or set a tumour tissue) to enable malignant calling")
        return {"status": "skipped: non-tumour context", "tumour_context": False, "gate": gate,
                "n_genes": n_genes, "callers": [], "concordance": None, "notes": [why]}

    callers, notes = [], []
    # 1) cheap marker score (always available)
    malignant_score(adata, tissue=tissue)
    if "malignant_score" in adata.obs:
        s = np.asarray(adata.obs["malignant_score"], dtype=float)
        callers.append({"name": "marker score", "status": "ok",
                        "pct_malignant": float((s > 0.6).mean()), "threshold": 0.6})
    # 2) Cancer-Finder (any tumour panel)
    cf = _cf.call_cancerfinder(adata, threshold=cf_threshold, max_cells=max_cells)
    callers.append({"name": "Cancer-Finder", "status": cf.get("status", "ok"),
                    "pct_malignant": float(cf.get("pct_malignant", 0.0)),
                    "threshold": cf.get("threshold")})
    if cf.get("status") != "ok":
        notes.append(f"Cancer-Finder: {cf.get('status')}")
    # InSituCNV is intentionally NOT run: the isolated-env infercnvpy pass (Jensen smoothing + infercnv
    # over a full GTF) takes many minutes even on a >2000-probe panel, too slow for the interactive
    # app. `call_malignant_cnv` stays available for offline/batch use. With it off there is only one
    # learned caller (Cancer-Finder), so the concordance below stays None (no `is_malignant`).

    # concordance between the two LEARNED callers over cells both scored (only if InSituCNV was run)
    concordance, n_joint = None, 0
    if "cancerfinder_malignant" in adata.obs and "is_malignant" in adata.obs:
        cfp = np.asarray(adata.obs.get("cancerfinder_prob", np.full(adata.n_obs, np.nan)), dtype=float)
        cnvs = np.asarray(adata.obs.get("cnv_score", np.full(adata.n_obs, np.nan)), dtype=float)
        m = ~np.isnan(cfp) & ~np.isnan(cnvs)
        n_joint = int(m.sum())
        if m.any():
            concordance = float((adata.obs["cancerfinder_malignant"].to_numpy(dtype=bool)[m]
                                 == adata.obs["is_malignant"].to_numpy(dtype=bool)[m]).mean())
    return {"status": "ok", "tumour_context": True, "gate": gate, "n_genes": n_genes, "callers": callers,
            "concordance": concordance, "n_jointly_scored": n_joint, "notes": notes}


def call_malignant_cnv(adata, normal_key: str = "cell_type", normal_categories=None,
                       gtf: str | None = None, env_python: str | None = None,
                       max_cells: int = 0, marker_sets: dict | None = None,
                       reference_purity: float = 0.5, min_reference: int = 200) -> dict:
    """CNV-based malignant calling via infercnvpy, run in the isolated ``cnv_env`` env.

    infercnvpy cannot be imported in the main interpreter, so this exports the section and runs
    ``subprocesses/cnv/run_infercnv.py`` in ``cnv_env`` (Jensen neighbor-smoothing + infercnv +
    per-cell RMS CNV burden), then joins ``obs['cnv_score']`` + ``obs['is_malignant']`` back.

    Parameters
    ----------
    normal_key: the ``obs`` column holding cell types (the CNV reference axis).
    normal_categories: the non-tumor reference cell types (immune / stromal / endothelial). If
        None, inferred as every type NOT matching Malignant / Tumor / Epithelial.
    marker_sets: if given, restrict the diploid reference to **marker-pure** cells
        (``pmp >= reference_purity``) so a noisy-label reference does not flatten the CNV baseline.
    reference_purity / min_reference: purity cutoff and the minimum pure-reference count below
        which the filter is skipped (keep all reference-type cells).
    gtf / env_python: overrides for the GTF and the infercnvpy python (site-configured via
        SPATIALSCRIBE_CNV_GTF / SPATIALSCRIBE_CNV_PYTHON; no public default).
    max_cells: uniform subsample cap for tractability (0 = all; CNV separation stays strong at
        ~15-25k on a 5K panel). infercnv on 10^5 cells needs a big-memory node.

    Returns a summary with ``pct_malignant`` + ``cnv_threshold``; on ANY failure (env missing,
    subprocess error, too few positioned genes) returns ``{'status': 'skipped: ...'}`` and never
    raises, so the app degrades to the marker-based ``malignant_score``. CNV quality scales with
    panel size - reliable on >=5K Prime / WTA, unreliable on a 480-gene panel.
    """
    import os
    import subprocess
    import tempfile
    from pathlib import Path

    import numpy as np

    def _skip(msg: str) -> dict:
        return {"status": f"skipped: {msg}", "pct_malignant": 0.0, "cnv_threshold": float("nan")}

    env_python = env_python or _GI_INSITUCNV_PY
    gtf = gtf or _GTF_DEFAULT
    if not env_python:
        return _skip("no CNV env configured; set SPATIALSCRIBE_CNV_PYTHON (a python with infercnvpy)")
    if not Path(env_python).exists():
        return _skip(f"CNV env not found ({env_python})")
    if not gtf:
        return _skip("no GTF configured; set SPATIALSCRIBE_CNV_GTF (a GTF with gene coordinates)")
    if not Path(gtf).exists():
        return _skip(f"GTF not found ({gtf})")
    if normal_key not in adata.obs:
        return _skip(f"no '{normal_key}' column for the CNV reference")

    ref_labels, normal_categories = _build_reference(
        adata, normal_key, normal_categories, marker_sets, reference_purity, min_reference)
    if ref_labels is None:
        return _skip(normal_categories)   # holds the error message on failure

    import anndata as ad
    import pandas as pd

    tmp = Path(tempfile.mkdtemp(prefix="sscnv_"))
    h5, out = tmp / "section.h5ad", tmp / "cnv.parquet"
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    # Export a single purpose-built reference column ('reference'/'query') so infercnv's baseline
    # is exactly the (optionally purity-filtered) diploid reference, independent of the label noise.
    a_exp = ad.AnnData(X=counts.copy(),
                       obs=pd.DataFrame({"_cnv_ref": ref_labels}, index=adata.obs_names.astype(str)))
    a_exp.var_names = adata.var_names.astype(str)
    a_exp.write_h5ad(h5)

    script = str(Path(__file__).resolve().parents[3] / "subprocesses" / "cnv" / "run_infercnv.py")
    ld = os.environ.get("LD_LIBRARY_PATH", "")
    if _GI_INSITUCNV_LIB:
        ld = f"{_GI_INSITUCNV_LIB}:{ld}"
    env = dict(os.environ, PYTHONNOUSERSITE="1", LD_LIBRARY_PATH=ld)
    cmd = [env_python, script, "--h5ad", str(h5), "--out", str(out), "--gtf", gtf,
           "--cell-type-key", "_cnv_ref", "--reference", "reference",
           "--max-cells", str(int(max_cells))]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=7200, env=env)
    except Exception as exc:
        return _skip(f"subprocess error ({exc})")
    if proc.returncode != 0 or not out.exists():
        tail = (proc.stderr or "").strip().splitlines()
        return _skip(f"infercnv subprocess failed ({tail[-1] if tail else 'no output'})")

    covered = _join_cnv(adata, out)
    score = np.asarray(adata.obs["cnv_score"], dtype=float)
    ismal = np.asarray(adata.obs["is_malignant"], dtype=bool)
    valid = ~np.isnan(score)
    ref_valid = (ref_labels == "reference") & valid
    thr = float(np.nanpercentile(score[ref_valid], 95)) if ref_valid.any() else float("nan")
    pct = float(ismal[valid].mean()) if covered else 0.0
    return {"status": "ok", "cnv_threshold": thr, "pct_malignant": pct,
            "coverage": covered / max(1, adata.n_obs),
            "n_reference": int((ref_labels == "reference").sum()),
            "reference_types": normal_categories,
            "reference_purity_filtered": bool(marker_sets)}
