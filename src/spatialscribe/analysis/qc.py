"""Quality control for imaging-based spatial data - global and region-scoped (H1).

What it does
------------
`compute_qc` adds per-cell QC columns (counts, genes, % control probes). `qc_summary`
reduces them to section-level medians and flags them against the Xenium thresholds.
`region_qc` computes the same summary for a box-select-selected subset of cells (H1), and
`suggest_filter` returns sensible default cell filters.

Depends on
----------
scanpy, numpy. Reads the ``control`` boolean var column set by ``io.load``.

Thresholds are section-level heuristics transcribed from
``the 10x section-level QC thresholds`` and are context-dependent
(fresh/FFPE, panel size, discovery-vs-validation) - advisory, not hard cutoffs.
"""

from __future__ import annotations

import logging

import numpy as np

# Section-level warn/error thresholds (metric -> (direction, warn, error)).
# direction "low" = smaller is worse; "high" = larger is worse.
XENIUM_THRESHOLDS: dict[str, tuple[str, float, float]] = {
    "median_genes_per_cell": ("low", 20, 10),
    "median_transcripts_per_cell": ("low", 50, 25),
    "fraction_empty_cells": ("high", 0.15, 0.25),
    "pct_counts_control": ("high", 2.0, 5.0),  # per-cell % control probe counts
}

# Per-cell filter defaults. minCounts=10 keeps low-RNA T cells/neutrophils that
# spillover hits hardest (repo convention; lowered from 20).
DEFAULT_FILTER = {"min_counts": 10, "min_genes": 5, "max_pct_control": 5.0}


def compute_qc(adata, progress=None) -> None:
    """Add per-cell QC metrics in place.

    Produces ``total_counts``, ``n_genes_by_counts`` and ``pct_counts_control``. ``progress(frac,
    label)`` (optional) reports coarse checkpoints for the app's progress bar.

    ``pct_counts_control`` is computed over the STRICT negative-control features
    (``var['neg_control']``, see :func:`io.build_neg_control_mask`), NOT the broad ``var['control']``
    panel mask. The broad mask includes deprecated codewords, which carry real transcripts: on Xenium
    Prime 5K they hold ~15% of every cell's counts while true negative controls hold 0.000%, so the
    old metric saturated the contamination penalty and abstained 94.8% of the section.
    """
    from .backend import get_backend
    from .io import build_neg_control_mask

    if "control" not in adata.var.columns:
        adata.var["control"] = False
    if "neg_control" not in adata.var.columns:
        adata.var["neg_control"] = build_neg_control_mask(adata)
    if progress:
        progress(0.3, "per-cell QC metrics")
    # QC counts/genes/%-control must come from RAW counts. If X has been normalized (log1p/scale) but a
    # raw 'counts' layer is present (e.g. the pre-processed demo caches, or a section clustered before QC
    # re-runs), compute the metrics on the counts layer - else median transcripts/cell can fall BELOW
    # median genes/cell (biologically impossible) and a shallow section can cross the depth warn/error
    # thresholds on log-space values. One guard here; every caller (the compute_qc capability,
    # run_pipeline stage 1, qc_summary's lazy back-fill) inherits it.
    use_counts = ("counts" in getattr(adata, "layers", {})) and _x_is_normalized(adata)
    _saved_X = None
    if use_counts:
        _saved_X = adata.X
        adata.X = adata.layers["counts"]
    try:
        get_backend().calculate_qc_metrics(adata, qc_vars=["neg_control"])
    finally:
        if use_counts:
            adata.X = _saved_X
    # Keep the long-standing column name every caller reads (annotate, rejection, qc_summary, the
    # frontend), now with corrected semantics: background counts, not "non-gene" counts.
    adata.obs["pct_counts_control"] = adata.obs["pct_counts_neg_control"]
    # MERSCOPE keeps its Blank negative-control probes in obsm['blank_genes'] (var is panel-only, so the
    # var-based neg_control mask is all-False). Without this, pct_counts_control is a STRUCTURAL 0 on
    # MERSCOPE and the contamination / background-rate QC (+ the apply_confidence contamination penalty)
    # are blind to a technology whose negative controls are right there - compute the control percentage
    # from the blanks, matching the Xenium neg-control semantics (% of a cell's counts on a control).
    blank = adata.obsm["blank_genes"] if ("blank_genes" in getattr(adata, "obsm", {})) else None
    if blank is not None:
        b = blank.values if hasattr(blank, "values") else blank   # DataFrame -> ndarray; sparse/ndarray as-is
        blank_pc = np.asarray(b.sum(axis=1)).ravel().astype(float)
        gene_pc = np.asarray(adata.obs.get("total_counts", np.zeros(adata.n_obs)), dtype=float)
        denom = gene_pc + blank_pc
        adata.obs["pct_counts_control"] = np.where(denom > 0, 100.0 * blank_pc / denom, 0.0)
    if progress:
        progress(1.0, "done")


def _x_is_normalized(adata) -> bool:
    """True if ``adata.X`` looks normalized (log1p / scaled), so QC must read ``layers['counts']``.

    Signals: scanpy stamps ``uns['log1p']`` after ``sc.pp.log1p``; and raw imaging counts are
    non-negative integers, so any negative or fractional value means ``X`` is not raw counts. Cheap:
    samples the sparse ``.data`` (or a small dense slice), never densifies the whole matrix.
    """
    if "log1p" in getattr(adata, "uns", {}):
        return True
    X = getattr(adata, "X", None)
    if X is None:
        return False
    data = getattr(X, "data", None)
    sample = np.asarray(data[:2000]) if data is not None else np.asarray(X[:50]).ravel()
    finite = sample[np.isfinite(sample)] if sample.size else sample
    if finite.size == 0:
        return False
    return bool((finite < 0).any() or not np.allclose(finite, np.round(finite)))


def qc_summary(adata, subset: np.ndarray | None = None, platform: str | None = None) -> dict:
    """Section (or subset) QC summary with platform-aware warn/error flags.

    Parameters
    ----------
    subset: optional boolean/int index selecting cells (used by ``region_qc``).
    platform: overrides ``uns['platform']`` (stamped by ``io.load``); drives the per-platform
        threshold profile (looser neg-control on CosMx, no retention hard-fail on MERSCOPE).
    """
    # Never silently zero-default the QC columns: a missing total_counts would read as "0 genes,
    # 0 counts, 100% empty" and drive a false 'section unusable' verdict. Compute them (backend-aware,
    # from the raw counts) when a step has not populated them yet - e.g. the LLM QC verdict runs
    # before the QC step, or a GPU path did not persist the columns.
    if "total_counts" not in adata.obs or "n_genes_by_counts" not in adata.obs:
        compute_qc(adata)
    a = adata[subset] if subset is not None else adata
    platform = platform or adata.uns.get("platform", "xenium")
    n = int(a.n_obs)
    counts = np.asarray(a.obs.get("total_counts", np.zeros(n)))
    genes = np.asarray(a.obs.get("n_genes_by_counts", np.zeros(n)))
    pct_ctrl = np.asarray(a.obs.get("pct_counts_control", np.zeros(n)), dtype=float)

    # NaN-safe: an empty cell has pct_counts_control = 0/0 = NaN, and np.median of an array holding
    # ONE NaN is NaN - which makes every `metric > threshold` comparison below False, silently
    # passing the neg-control flag. Measured on 4 of 5 real benchmark sections (2026-07).
    all_nan = bool(n) and bool(np.isnan(pct_ctrl).all())
    metrics = {
        "n_cells": n,
        "platform": platform,
        "median_genes_per_cell": float(np.median(genes)) if n else 0.0,
        "median_transcripts_per_cell": float(np.median(counts)) if n else 0.0,
        "fraction_empty_cells": float(np.mean(counts < 1)) if n else 0.0,
        "pct_counts_control": float(np.nanmedian(pct_ctrl)) if n and not all_nan else 0.0,
    }
    metrics["flags"] = _flag(metrics, platform)
    # Surface the annotation output-usability gate (stamped by annotate.apply_confidence) next to
    # the input-QC flags, so "ran to completion, QC ok" cannot hide a section that abstained ~95%.
    usability = adata.uns.get("annotation_usability") if hasattr(adata, "uns") else None
    if isinstance(usability, dict) and usability.get("flag") in ("warn", "fail"):
        metrics["flags"]["annotation_usability"] = str(usability["flag"])
    # Completeness (informational, not a hard flag - a green-coverage lineage may be genuinely absent):
    # surface how many panel-resolvable lineages annotated to 0 cells, so a "usability ok" headline
    # cannot hide an annotation that recovered only some of the lineages the panel can resolve.
    comp = adata.uns.get("annotation_completeness") if hasattr(adata, "uns") else None
    if isinstance(comp, dict) and comp.get("resolvable_absent"):
        metrics["resolvable_absent_lineages"] = list(comp["resolvable_absent"])
    return metrics


def _platform_flag_modifiers(platform: str) -> tuple[float, set[str]]:
    """Per-platform section-flag modifiers, read from the ``cross_platform`` YAML profile.

    Returns ``(neg_control_factor, no_hard_error_metrics)``:
      * ``neg_control_factor`` scales the ``pct_counts_control`` warn/fail thresholds - CosMx's
        higher background + spatial aggregation of negative-control signal needs a looser margin
        (``neg_control_threshold_scaling: cosmx = loosen`` -> 2x);
      * ``no_hard_error_metrics`` are metrics that must never hard-``error`` on a platform whose
        typical FFPE retention floor is < 50% (MERSCOPE legitimately loses most cells), derived
        from ``typical_qc_retention_ffpe`` so a newly-added low-retention platform is auto-covered.
    """
    from . import config

    scaling = config.get("cross_platform", "neg_control_threshold_scaling", default={}) or {}
    raw = scaling.get(platform, 1.0)
    if isinstance(raw, (int, float)):
        factor = float(raw)
    elif str(raw).strip().lower() == "loosen":
        factor = 2.0
    else:
        try:
            factor = float(raw)
        except (TypeError, ValueError):
            factor = 1.0

    retention = config.get("cross_platform", "typical_qc_retention_ffpe", default={}) or {}
    no_error: set[str] = set()
    lo = retention.get(platform)
    if isinstance(lo, (list, tuple)) and lo and float(lo[0]) < 0.5:
        no_error.add("fraction_empty_cells")
    return factor, no_error


def _flag(metrics: dict, platform: str = "xenium") -> dict[str, str]:
    """Return {metric: 'ok'|'warn'|'error'} against XENIUM_THRESHOLDS, per-platform-adjusted."""
    factor, no_error = _platform_flag_modifiers(platform)
    out = {}
    for key, (direction, warn, err) in XENIUM_THRESHOLDS.items():
        if key not in metrics:
            continue
        w, e = (warn * factor, err * factor) if key == "pct_counts_control" else (warn, err)
        val = metrics[key]
        worse = (lambda v, t: v < t) if direction == "low" else (lambda v, t: v > t)
        flag = "error" if worse(val, e) else "warn" if worse(val, w) else "ok"
        if flag == "error" and key in no_error:
            flag = "warn"  # low-retention platform: down-rank, don't hard-fail the section
        out[key] = flag
    return out


def region_qc(adata, cell_index) -> dict:
    """QC summary for a box-select-selected region (H1). ``cell_index`` = bool/int index."""
    return qc_summary(adata, subset=cell_index)


def apply_ovrlpy_vsi(adata, vsi_parquet, threshold: float = 0.7) -> dict:
    """Join per-cell VSI (from ``subprocesses/ovrlpy/run_ovrlpy.py``) onto ``adata.obs`` (H5).

    Adds ``obs['vsi']`` (0-1 confidence) and ``obs['vsi_confident']``. Flags, never drops
    (low-VSI cells are spatially non-random; dropping biases the analysis).
    """
    import pandas as pd

    df = pd.read_parquet(vsi_parquet)
    df["cell_id"] = df["cell_id"].astype(str)
    aligned = df.set_index("cell_id").reindex(adata.obs_names.astype(str))
    adata.obs["vsi"] = aligned["vsi"].to_numpy()
    adata.obs["vsi_confident"] = adata.obs["vsi"] >= threshold
    n_low = int((adata.obs["vsi"] < threshold).sum())
    return {"n_low_vsi": n_low, "pct_low_vsi": float(n_low / max(1, adata.n_obs)), "threshold": threshold}


def suggest_filter(adata) -> dict:
    """Default per-cell filter thresholds (see ``DEFAULT_FILTER``)."""
    return dict(DEFAULT_FILTER)


def segmentation_qc(adata) -> dict:
    """Layer 1 segmentation QC from Xenium ``cells.parquet`` columns (loaded by ``io``).

    Writes ``obs['nucleus_present']``, ``obs['nucleus_to_cell_ratio']`` and a dataset-relative
    ``obs['seg_area_flag']`` (small/large/ok). No-ops (status='skipped') when the platform
    lacks segmentation area columns.
    """
    import numpy as np

    if "cell_area" not in adata.obs.columns:
        return {"status": "skipped", "reason": "no cell_area column (platform lacks segmentation metrics)"}

    area = np.asarray(adata.obs["cell_area"], dtype=float)
    has_nuc = "nucleus_area" in adata.obs.columns
    if has_nuc:
        nuc = np.asarray(adata.obs["nucleus_area"], dtype=float)
        adata.obs["nucleus_present"] = nuc > 0
        with np.errstate(divide="ignore", invalid="ignore"):
            adata.obs["nucleus_to_cell_ratio"] = np.where(area > 0, nuc / area, np.nan)

    lo, hi = np.nanpercentile(area, [0.5, 99.5])
    # Guard a degenerate area distribution (constant/placeholder column, or ties at the extremes). If
    # the 0.5/99.5 span collapses (hi <= lo), a naive `area < lo` / `area > hi` either flags nothing or
    # half the section. When the span is real, use INCLUSIVE boundaries so a block tied at the minimum
    # (e.g. 0-area fragments) is flagged small, not skipped by a strict `<`. Same boundary + degeneracy
    # lesson as the count floor (a distributional cutoff that lands on a tie removes nothing).
    if not (np.isfinite(lo) and np.isfinite(hi)) or hi <= lo:
        flag = np.full(area.shape, "ok", dtype=object)
    else:
        # Strict percentile tails (so a MODAL boundary value - most cells sharing one area - is not
        # all flagged), but ALWAYS flag 0-area empties as small: a block of 0-area junk ties the low
        # percentile so a strict `< lo` alone misses it. The count-floor lesson (catch the empties)
        # without the over-flagging an inclusive `<=`/`>=` on both tails would cause.
        flag = np.where((area <= 0) | (area < lo), "small", np.where(area > hi, "large", "ok"))
    adata.obs["seg_area_flag"] = flag
    pct_outlier = float(np.mean(flag != "ok"))
    # The 0.5/99.5 tails put the nominal outlier fraction near ~1%; a much larger fraction means a
    # heavy tail of implausibly tiny/huge segments (fragments or merges, or a block of 0-area junk).
    # Warn so the segmentation dot REFLECTS pct_area_outlier instead of always reading green.
    return {
        "status": "warn" if pct_outlier > 0.05 else "ok",
        "median_cell_area": float(np.nanmedian(area)),
        "pct_no_nucleus": float(np.mean(~adata.obs["nucleus_present"].to_numpy(dtype=bool))) if has_nuc else None,
        "pct_area_outlier": pct_outlier,
    }


def resolve_panel_size(adata, n_panel_genes: int | None = None) -> int:
    """Panel gene count, in priority order: explicit arg > ``uns['n_panel_genes']``
    (stamped by ``io.load``) > non-control ``var`` count > ``n_vars``.

    The ``uns`` hop is what lets a small unit-test AnnData (few genes in ``X``) declare itself
    a 5K/WTA section, and what keeps every adata-only downstream function panel-aware.
    """
    if n_panel_genes is not None:
        return int(n_panel_genes)
    if "n_panel_genes" in adata.uns:
        return int(adata.uns["n_panel_genes"])
    if "control" in adata.var.columns:
        return int((~adata.var["control"]).sum())
    return int(adata.n_vars)


# Minimum detected genes for a cell to carry reliable identity signal on a RICH panel (>=1000 genes).
# ABSOLUTE, not a fraction of the section median: a deep section (e.g. atera5k, 552 counts/cell) has
# absolutely-fine cells well below its own median, and a fraction-of-median floor wrongly dropped 9% of
# that positive control. Below this a cell is near-empty (<0.3% of a 5K panel): Leiden splits such cells
# into single-gene artefact micro-clusters (the breast 5K demo made 16 of them at res 0.1, each a
# housekeeping/lncRNA "marker") and annotate cannot type them. Tuned on breast 5K (junk tail 6-13 genes
# vs real populations at 47+; a floor of 15 collapsed 19 clusters -> 3 while excluding 0.0% of atera5k).
# Targeted panels (<1000 genes) keep their own, more permissive floor - a flat 15 would gut low-plex cells.
RICH_PANEL_MIN_GENES = 15


def panel_indexed_floor(counts, n_panel_genes: int, percentile: float | None = None) -> tuple[float, str]:
    """THE single source of truth for the panel-size-indexed count floor (docs Layer 2).

    Returns ``(floor, mode)``. A fixed floor is not panel-transferable: ``<10`` removes ~0.2%
    on a targeted panel but ~60% on Xenium Prime 5K. So ``<1000``-gene (targeted) panels get a
    fixed floor; ``>=1000``-gene (Prime 5K / WTA) panels get a section-relative percentile floor.

    Used by BOTH ``suggest_count_floor`` (Layer 2 filtering) AND ``annotate.apply_confidence``
    (the Layer 5 low-signal gate) so the two can never diverge into the exact inconsistency the
    docs warn about - a distributional filter followed by a fixed ``counts<10`` abstention gate.
    """
    import numpy as np

    from . import config

    if int(n_panel_genes) < 1000:
        sched = config.get("layer2_counts", "total_counts_panel_schedule", "targeted_lt1000_genes", default={}) or {}
        return float(sched.get("fail", 10)), "fixed"
    pct = percentile if percentile is not None else float(
        config.get("layer2_counts", "large_panel_percentile", default=0.02) or 0.02)
    counts = np.asarray(counts, dtype=float)
    # Floor at 1 so a SHALLOW section whose pct-th percentile is 0 (>= pct% of cells have 0 counts -
    # common on a low-depth 5K panel) still drops its truly-empty cells. A raw quantile of 0 gives a
    # floor of 0, and `counts < 0` removes NOTHING - not even the 0-count segments - so the count-floor
    # layer silently does nothing ("count floor 0, 0.0% removed"). Clamping to 1 only ever removes
    # 0-count cells; it never reaches the aggressive fixed <10 floor that over-removes on 5K.
    floor = max(float(np.quantile(counts, pct)) if counts.size else 0.0, 1.0)
    return floor, "distributional"


def suggest_count_floor(adata, n_panel_genes: int | None = None, platform: str | None = None) -> dict:
    """Panel-size-indexed minimum-count floor (docs Layer 2).

    Keys the floor off the *declared* panel size (see ``resolve_panel_size``: honours
    ``uns['n_panel_genes']`` stamped by ``io.load``), not the number of genes that survived into
    ``X``. Returns the chosen floor, its mode, the schedule bucket and its impact so the UI can
    preview before applying.
    """
    import numpy as np

    n_panel_genes = resolve_panel_size(adata, n_panel_genes)
    platform = platform or adata.uns.get("platform", "xenium")
    counts = np.asarray(adata.obs.get("total_counts", np.zeros(adata.n_obs)), dtype=float)

    floor, mode = panel_indexed_floor(counts, n_panel_genes)
    bucket = ("targeted_lt1000_genes" if n_panel_genes < 1000
              else "xenium_prime_5k" if n_panel_genes < 12000 else "whole_transcriptome")
    n_removed = int((counts < floor).sum())
    return {
        "mode": mode,
        "floor": floor,
        "n_panel_genes": int(n_panel_genes),
        "platform": platform,
        "schedule_bucket": bucket,
        "n_removed": n_removed,
        "pct_removed": float(n_removed / max(1, adata.n_obs)),
    }


def apply_filter(adata, min_counts=None, min_genes=None, max_pct_control=None,
                 exclude_cells: np.ndarray | None = None):
    """Return a filtered copy. ``exclude_cells`` drops region-excluded cells (H1)."""
    import scanpy as sc

    f = {**DEFAULT_FILTER}
    if min_counts is not None:
        f["min_counts"] = min_counts
    if min_genes is not None:
        f["min_genes"] = min_genes
    if max_pct_control is not None:
        f["max_pct_control"] = max_pct_control

    a = adata.copy()
    if exclude_cells is not None:
        keep = np.ones(a.n_obs, dtype=bool)
        keep[exclude_cells] = False
        a = a[keep].copy()
    sc.pp.filter_cells(a, min_counts=f["min_counts"])
    sc.pp.filter_cells(a, min_genes=f["min_genes"])
    if "pct_counts_control" in a.obs:
        a = a[a.obs["pct_counts_control"] <= f["max_pct_control"]].copy()
    return a


def run_funnel(adata, cluster_key: str = "cell_type", marker_sets=None,
               panel_check_result: dict | None = None, progress=None,
               lineage_markers=None, reference=None, ref_key: str = "cell_type") -> dict:
    """Run the QC funnel end-to-end over whatever data is available and return the section headline.

    Always runs Layer 0 (section summary), Layer 1 (segmentation, no-ops without columns) and the
    Layer 2 panel-indexed count floor. If ``cluster_key`` is present, also runs Layer 3 PMP, Layer 6
    spatial coherence, and the Layer 5 fusion (``apply_confidence``) to produce the annotatability
    headline. Purely orchestration; each metric owns its own obs columns. ``progress(frac, label)``
    (optional) reports coarse checkpoints for the app's progress bar.
    """
    from . import annotate as _annotate
    from . import purity as _purity
    from . import spatial as _spatial
    from . import markers as _markers

    # Curated marker symbols are human-UPPERCASE, but the panel may ship mouse title-case (Gfap vs
    # GFAP). Re-key every set to the panel's ACTUAL casing via markers.on_panel (the species/case-
    # adaptive resolver panel_check already uses) BEFORE the marker-based metrics (PMP purity, per-cell
    # confidence, and the AQI's C/M) look genes up in var_names - otherwise every marker reads as
    # off-panel on a mouse section, no_dict_frac hits 1.0, C and M collapse to None, and the AQI
    # silently degrades to the internal-validity fallback. No-op on a human panel (UPPERCASE==UPPERCASE).
    if marker_sets:
        _panel = list(map(str, adata.var_names))
        marker_sets = {t: _markers.on_panel(_panel, g) for t, g in marker_sets.items()}

    if "total_counts" not in adata.obs.columns:
        compute_qc(adata)                            # no progress: its 0.3/1.0 would break the funnel's monotone sequence

    if progress:
        progress(0.10, "section metrics")
    section = qc_summary(adata)
    if progress:
        progress(0.25, "segmentation")
    seg = segmentation_qc(adata)
    if progress:
        progress(0.40, "count floor")
    headline: dict = {
        "n_cells": int(adata.n_obs),
        "section": section,
        "segmentation": seg,
        "count_floor": suggest_count_floor(adata),
        "layers_run": ["0_section", "1_segmentation", "2_count_floor"],
    }

    if cluster_key in adata.obs.columns:
        if progress:
            progress(0.55, "purity (PMP)")
        _purity.pmp(adata, assigned_label_key=cluster_key, lineage_markers=marker_sets)
        if progress:
            progress(0.70, "spatial coherence")
        headline["spatial_coherence"] = _spatial.spatial_coherence(adata, label_key=cluster_key)
        if progress:
            progress(0.85, "per-cell confidence")
        headline["annotatability"] = _annotate.apply_confidence(
            adata, cluster_key=cluster_key, marker_sets=marker_sets,
            panel_check_result=panel_check_result, lineage_markers=lineage_markers,
        )
        headline["layers_run"] += ["3_pmp", "6_spatial_coherence", "5_confidence"]

        # Reference-free annotation-quality battery (internal validity + marker-program fidelity;
        # docs/research/cell-annotation-quality-metrics.md). Guarded: never breaks the funnel.
        # Subsample large sections: annotation_quality's marker-fidelity + neighborhood-purity parts
        # touch ALL cells (only its silhouette self-subsamples), so on a 100k-cell section they run
        # for minutes and the funnel appears "stuck". The dominant term is an O(N^2) Ward linkage
        # inside eval_metrics.internal_validity (~85% of the capability's wall clock), so the working
        # set governs cost quadratically: a 4k random subsample keeps the section-level estimate
        # representative (measured AQI delta < 1e-3 vs 15k) while cutting annotation_quality ~7x
        # (8.1s -> 1.1s at 20k). (try/except only catches errors, NOT a slow-but-succeeding call, so
        # the cap - not the guard - is what prevents the hang.)
        if progress:
            progress(0.95, "annotation-quality battery")
        try:
            import numpy as np

            from . import annotate as _an
            from . import eval_metrics as _em
            aq = adata
            if adata.n_obs > 4000:
                idx = np.sort(np.random.default_rng(0).choice(adata.n_obs, 4000, replace=False))
                aq = adata[idx].copy()
            # AQI inputs: reference-free purity+fidelity always; the panel/depth CEILING when a reference
            # is threaded; the WITHIN-section abstention signal when >=3 reference methods have voted.
            depth = (float(np.median(np.asarray(aq.obs["total_counts"])))
                     if "total_counts" in aq.obs else None)
            mcols = [c for c in ("rctd_first_type", "singler_label", "scanvi_label", "ph_fine")
                     if c in aq.obs.columns]
            headline["annotation_quality"] = _em.annotation_quality(
                aq, label_key=cluster_key, marker_sets=marker_sets,
                reference=reference, ref_key=ref_key, panel_genes=list(map(str, aq.var_names)),
                median_depth=depth, method_label_cols=(mcols if len(mcols) >= 3 else None),
                abstention_labels=list(_an.ABSTENTION_LABELS),
                platform=aq.uns.get("platform"), n_panel_genes=aq.uns.get("n_panel_genes"),
                normalization="log1p-median")
            headline["layers_run"].append("quality_metrics")
        except Exception:  # degrade-graceful, but never SILENT (a deploy version-skew or a bad section)
            # The AQI battery must never break the funnel, but swallowing bare-`pass` turned a
            # one-line-diagnosable failure (e.g. an in-memory spatial_anno_metrics older than on disk,
            # rejecting the current kwargs) into a mystery empty AQI tile behind an HTTP 200. Log the
            # traceback so the next failure is loud; the funnel still returns (tile just degrades empty).
            logging.getLogger(__name__).warning(
                "annotation_quality failed - AQI tile will be empty", exc_info=True)

    if progress:
        progress(1.0, "done")
    return headline
