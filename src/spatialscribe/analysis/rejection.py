"""Granular, plain-language *reasons* for why a cell could not be confidently typed.

What it does
------------
``apply_confidence`` (see :mod:`spatialscribe.analysis.annotate`) already writes a coarse
verdict/reason per cell (``annotation_verdict`` PASS/WARN/FAIL and one of five coarse
``annotation_reason`` classes). This module ENRICHES that into a SPECIFIC, actionable
reason - "too few transcripts (7)", "your panel lacks the markers to resolve NK cells",
"markers from multiple lineages (spillover)", "no cell type scored clearly above the
others" - by inspecting whichever QC columns the funnel produced. It is robust to missing
columns: a demo cache without ``vsi``/``seg_area_flag``/``spatial_coherence`` still gets a
useful reason inferred from counts / genes / purity / confidence, and it falls back to the
coarse ``annotation_reason`` when present.

How to use
----------
    from spatialscribe.analysis import rejection
    rejection.assign_rejection_reasons(adata)           # writes obs columns, returns summary
    tab = rejection.rejection_breakdown(adata)          # per-reason table (the section headline)
    warn = rejection.panel_resolvability_warnings(adata, adata.uns.get("panel_check"))

Writes ``obs['rejection_reason']`` (granular code, ``""`` for PASS cells) and
``obs['rejection_detail']`` (a plain-language string with the specific type / numbers).

Precedence
----------
When several reasons apply, the most upstream / decisive one wins (a broken segment or a
near-empty cell is reported before contamination, which is reported before a panel gap,
which is reported before a low-margin / spatial-fit call). ``spatially_incoherent`` is the
gentlest reason: it is advisory (a down-weight, never a hard reject) because rare-but-real
infiltrating cells are legitimately surrounded by other types.

Depends on: numpy, pandas; reads columns written by :mod:`annotate`, :mod:`qc`,
:mod:`purity`, :mod:`spatial` and the ``panel_check`` result. All heavy imports are inside
functions (module import stays cheap). Grounded strictly in real computed columns - no
invented numbers.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .annotate import CONF_WARN

if TYPE_CHECKING:
    import pandas

# Thresholds (mirrored from the layers that own each metric so the reason agrees with the
# funnel that produced it).
GENE_FLOOR = 5          # n_genes_by_counts hard gate (annotate.apply_confidence)
COUNT_HARD_FLOOR = 10   # total_counts hard gate (annotate.apply_confidence)
VSI_FLOOR = 0.7         # qc.apply_ovrlpy_vsi default confident threshold
CONTROL_PCT_MAX = 5.0   # qc.DEFAULT_FILTER max_pct_control
COHERENCE_FLOOR = 0.2   # spatial.spatial_coherence pas_threshold default


# (label, generic plain-language description). The per-cell ``rejection_detail`` fills in the
# specific cell type and numbers; these are the human-readable names/blurbs for the UI table.
REASON_LIBRARY: dict[str, tuple[str, str]] = {
    "too_few_counts": (
        "Too few transcripts",
        "Too few transcripts to identify this cell.",
    ),
    "too_few_genes": (
        "Too few genes detected",
        "Too few distinct genes were detected to identify this cell.",
    ),
    "poor_segmentation": (
        "Implausible segment",
        "Implausible segment (fragment or merged).",
    ),
    "spatial_doublet": (
        "Spatial doublet / vertical overlap",
        "Overlaps another cell vertically (low signal integrity).",
    ),
    "mixed_lineages": (
        "Mixed lineages (spillover)",
        "Markers from multiple lineages (spillover or contamination).",
    ),
    "high_background": (
        "High control-probe background",
        "A large fraction of counts are control-probe background.",
    ),
    "panel_cannot_separate": (
        "Panel cannot separate two types",
        "Panel cannot separate two candidate types; reported as the coarser group.",
    ),
    "panel_cannot_resolve": (
        "Panel lacks the markers",
        "Panel lacks the markers to confidently call this cell type.",
    ),
    "novel_unknown": (
        "Novel / unknown",
        "Does not match any reference type.",
    ),
    "low_margin": (
        "No clear winner",
        "No cell type scored clearly above the others.",
    ),
    "method_disagreement": (
        "Annotation methods disagree",
        "Annotation methods disagree on this cell's type (low cross-method consensus).",
    ),
    "spatially_incoherent": (
        "Label does not fit its neighborhood",
        "Label does not fit its spatial neighborhood (advisory, not a hard reject).",
    ),
}

# Evaluated in order; the first matching reason for a non-PASS cell wins (see "Precedence").
_PRECEDENCE = [
    "poor_segmentation",
    "too_few_counts",
    "too_few_genes",
    "spatial_doublet",
    "mixed_lineages",
    "high_background",
    "panel_cannot_separate",
    "panel_cannot_resolve",
    "novel_unknown",
    "method_disagreement",
    "spatially_incoherent",
    "low_margin",
]

# Coarse annotation_reason -> a sensible granular default, used only when no computed column
# lets us be more specific (keeps every non-PASS cell covered).
_COARSE_DEFAULT = {
    "low_quality": "too_few_counts",
    "ambiguous_mixed": "mixed_lineages",
    "unresolvable_panel": "panel_cannot_resolve",
    "uncertain_lowconf": "low_margin",
    "novel": "novel_unknown",
}


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _counts_and_genes(adata):
    """Per-cell (total_counts, n_genes) - read the raw ``counts`` layer when present so the
    gate agrees with ``apply_confidence`` (obs['total_counts'] may be normalized)."""
    import numpy as np

    n = adata.n_obs
    if "counts" in adata.layers:
        raw = adata.layers["counts"]
        counts = np.asarray(raw.sum(1)).ravel().astype(float)
        genes = np.asarray((raw > 0).sum(1)).ravel().astype(float)
        return counts, genes
    counts = np.asarray(adata.obs.get("total_counts", np.zeros(n)), dtype=float)
    genes = np.asarray(adata.obs.get("n_genes_by_counts", np.zeros(n)), dtype=float)
    return counts, genes


def _count_floor(adata) -> float:
    """Panel-indexed minimum-count floor (>= the hard floor of 10). Falls back to 10."""
    try:
        from . import qc

        floor = float(qc.suggest_count_floor(adata).get("floor", COUNT_HARD_FLOOR))
    except Exception:
        floor = float(COUNT_HARD_FLOOR)
    return max(floor, float(COUNT_HARD_FLOOR))


def _non_pass_mask(adata):
    """Boolean array: cells that are NOT confidently typed (WARN or FAIL).

    Prefers ``annotation_verdict``; falls back to ``annotation_confidence`` thresholds, else
    to any hard QC flag - so the function still works on a cache that lacks the verdict.
    """
    import numpy as np

    n = adata.n_obs
    if "annotation_verdict" in adata.obs:
        return adata.obs["annotation_verdict"].astype(str).to_numpy() != "PASS"
    if "annotation_confidence" in adata.obs:
        conf = np.asarray(adata.obs["annotation_confidence"], dtype=float)
        return conf < CONF_WARN
    counts, genes = _counts_and_genes(adata)
    impure = np.asarray(adata.obs.get("crisp_impure", np.zeros(n)), dtype=bool)
    return (counts < _count_floor(adata)) | (genes < GENE_FLOOR) | impure


def _panel_maps(panel_check_result):
    """From a check_panel result build (status_by_type, coverage_by_type, confusable_partner).

    ``confusable_partner`` maps a cell type to (other_type, [A, B]) for the first pair it is
    in. Guards the h5ad round-trip corruption (missing 'status') via ``panel_check.is_valid``.
    """
    from . import panel_check as _pc

    if not panel_check_result or not _pc.is_valid(panel_check_result):
        return {}, {}, {}
    cov = panel_check_result.get("coverage", {})
    status = {ct: d["status"] for ct, d in cov.items() if isinstance(d, dict) and "status" in d}
    partner: dict[str, tuple[str, list]] = {}
    for pair in panel_check_result.get("confusable_pairs", []):
        p = pair.get("pair", [])
        if len(p) == 2:
            a, b = p
            partner.setdefault(a, (b, [a, b]))
            partner.setdefault(b, (a, [a, b]))
    return status, cov, partner


# --------------------------------------------------------------------------- #
# main entry point
# --------------------------------------------------------------------------- #
def assign_rejection_reasons(adata, panel_check_result: dict | None = None) -> dict:
    """Enrich the coarse verdict into a granular per-cell rejection reason + detail.

    Writes ``obs['rejection_reason']`` (a code from :data:`REASON_LIBRARY`, ``""`` for
    PASS/confident cells) and ``obs['rejection_detail']`` (plain-language string with the
    specific cell type and numbers filled in). Robust to missing QC columns. Returns a small
    summary ``{n_non_pass, n_reasoned, breakdown}``.
    """
    import numpy as np

    n = adata.n_obs
    non_pass = _non_pass_mask(adata)

    # --- gather whatever signals are available (missing column -> all-False predicate) ---
    counts, genes = _counts_and_genes(adata)
    floor = _count_floor(adata)

    seg_flag = adata.obs["seg_area_flag"].astype(str).to_numpy() if "seg_area_flag" in adata.obs \
        else np.array([""] * n)
    # No-nucleus is NOT a segmentation problem: in a thin tissue section the nucleus often sits in an
    # adjacent cut, so an in-plane anucleate cell is EXPECTED, not implausible. Only size outliers
    # (small = fragment, large = merge/doublet) mark a poor segment. The section-level no-nucleus
    # fraction is still reported by qc.segmentation_qc (pct_no_nucleus) as an informational metric.
    seg_bad = np.isin(seg_flag, ("small", "large"))

    vsi = np.asarray(adata.obs["vsi"], dtype=float) if "vsi" in adata.obs else np.full(n, np.nan)
    doublet = np.nan_to_num(vsi, nan=1.0) < VSI_FLOOR

    impure = adata.obs["crisp_impure"].to_numpy(dtype=bool) if "crisp_impure" in adata.obs \
        else np.zeros(n, dtype=bool)

    pct_ctrl = np.asarray(adata.obs.get("pct_counts_control", np.zeros(n)), dtype=float)
    high_bg = pct_ctrl > CONTROL_PCT_MAX

    coh = (np.asarray(adata.obs["spatial_coherence"], dtype=float)
           if "spatial_coherence" in adata.obs else np.full(n, np.nan))
    incoherent = np.nan_to_num(coh, nan=1.0) < COHERENCE_FLOOR

    reason_col = (adata.obs["annotation_reason"].astype(str).to_numpy()
                  if "annotation_reason" in adata.obs else np.array([""] * n))
    novel = reason_col == "novel"

    ct = (adata.obs["cell_type"].astype(str).to_numpy()
          if "cell_type" in adata.obs else np.array(["this cell type"] * n))

    pc_in = panel_check_result or adata.uns.get("panel_check")
    status_by_ct, cov_by_ct, partner = _panel_maps(pc_in)
    cell_status = np.array([status_by_ct.get(c, "green") for c in ct])
    panel_weak = np.isin(cell_status, ("red", "amber"))
    confusable = np.array([c in partner for c in ct])

    # popV consensus disagreement (Nat Genet 2024): the annotation methods disagree AND a diverse
    # ensemble voted (below MIN_TRUST_METHODS the signal is noise). Threshold + min-methods are tunable
    # config; columns absent (written by :mod:`consensus`) -> all-False, so this reason is a no-op.
    from . import config as _cfg
    from . import consensus as _cons
    agr = (np.asarray(adata.obs["consensus_agreement"], dtype=float)
           if "consensus_agreement" in adata.obs else np.full(n, 1.0))
    n_meth = (np.asarray(adata.obs["consensus_n_methods"], dtype=float)
              if "consensus_n_methods" in adata.obs else np.zeros(n))
    dis_thr = float(_cfg.get("consensus_popv", "disagreement_agreement_max", default=0.5))
    method_disagree = (np.nan_to_num(agr, nan=1.0) < dis_thr) & (n_meth >= _cons.MIN_TRUST_METHODS)

    # --- predicate mask per code (only over non-PASS cells) ---
    masks = {
        "poor_segmentation": seg_bad,
        "too_few_counts": counts < floor,
        "too_few_genes": genes < GENE_FLOOR,
        "spatial_doublet": doublet,
        "mixed_lineages": impure,
        "high_background": high_bg,
        "panel_cannot_separate": confusable,
        "panel_cannot_resolve": panel_weak,
        "novel_unknown": novel,
        "method_disagreement": method_disagree,
        "spatially_incoherent": incoherent,
        "low_margin": np.ones(n, dtype=bool),   # catch-all for any remaining non-PASS cell
    }

    code = np.array([""] * n, dtype=object)
    unset = non_pass.copy()
    for name in _PRECEDENCE:
        take = unset & masks[name]
        code[take] = name
        unset &= ~take

    # Anything still unset (e.g. non-PASS but no computed column matched) -> coarse default.
    if unset.any():
        for i in np.nonzero(unset)[0]:
            code[i] = _COARSE_DEFAULT.get(reason_col[i], "low_margin")

    detail = _fill_details(code, counts, genes, seg_flag, vsi, pct_ctrl,
                           coh, ct, cov_by_ct, partner)

    import pandas as pd

    adata.obs["rejection_reason"] = pd.Categorical(code)
    adata.obs["rejection_detail"] = detail

    reasoned = code[code != ""]
    breakdown = {k: int(v) for k, v in pd.Series(reasoned).value_counts().items()}
    return {
        "n_non_pass": int(non_pass.sum()),
        "n_reasoned": int((code != "").sum()),
        "breakdown": breakdown,
    }


def _fill_details(code, counts, genes, seg_flag, vsi, pct_ctrl, coh, ct,
                  cov_by_ct, partner):
    """Build the per-cell plain-language ``rejection_detail`` string array (grounded numbers)."""
    import numpy as np

    n = len(code)
    detail = np.array([""] * n, dtype=object)

    for i in np.nonzero(code != "")[0]:
        c = code[i]
        if c == "too_few_counts":
            detail[i] = f"Too few transcripts ({int(counts[i])}) to identify this cell."
        elif c == "too_few_genes":
            detail[i] = f"Too few genes detected ({int(genes[i])})."
        elif c == "poor_segmentation":
            bits = []
            if seg_flag[i] == "small":
                bits.append("unusually small segment (possible fragment)")
            elif seg_flag[i] == "large":
                bits.append("unusually large segment (possible merge / doublet)")
            detail[i] = "Implausible segment: " + (", ".join(bits) or "fragment / merged") + "."
        elif c == "spatial_doublet":
            v = vsi[i]
            vtxt = f" (signal integrity {v:.2f} < {VSI_FLOOR:.2f})" if np.isfinite(v) else ""
            detail[i] = f"Overlaps another cell vertically{vtxt} - low signal integrity."
        elif c == "mixed_lineages":
            detail[i] = ("Markers from multiple lineages detected in this cell "
                         "(spillover or contamination from neighbors).")
        elif c == "high_background":
            detail[i] = (f"{pct_ctrl[i]:.1f}% of counts are control-probe background "
                         f"(> {CONTROL_PCT_MAX:.0f}%).")
        elif c == "panel_cannot_separate":
            other, pair = partner.get(ct[i], (None, [ct[i], "another type"]))
            detail[i] = (f"Panel cannot separate {pair[0]} vs {pair[1]}; "
                         "reported as the coarser group.")
        elif c == "panel_cannot_resolve":
            d = cov_by_ct.get(ct[i], {})
            k, m = d.get("n_present", 0), d.get("n_markers", 0)
            detail[i] = (f"Panel lacks markers to confidently call {ct[i]} "
                         f"(only {k} of {m} canonical markers on panel).")
        elif c == "novel_unknown":
            detail[i] = (f"{ct[i]} does not match any reference type - "
                         "possibly a state or type not in the reference.")
        elif c == "spatially_incoherent":
            frac = coh[i]
            ftxt = f" ({frac:.0%} same-type neighbors)" if np.isfinite(frac) else ""
            detail[i] = (f"Label ({ct[i]}) does not fit its neighborhood{ftxt} - "
                         "advisory, not a hard reject.")
        elif c == "low_margin":
            detail[i] = ("No cell type scored clearly above the others "
                         "(low margin between the top candidates).")
        else:
            detail[i] = REASON_LIBRARY.get(c, ("", ""))[1]
    return detail


# --------------------------------------------------------------------------- #
# reporting
# --------------------------------------------------------------------------- #
def rejection_breakdown(adata, panel_check_result: dict | None = None) -> "pandas.DataFrame":  # noqa: F821
    """Per-reason table (the section headline): count, % of non-PASS cells, label, description.

    Sorted by count descending. Runs :func:`assign_rejection_reasons` first when the column
    is missing so it can be called standalone (e.g. from the copilot tool).
    """
    import pandas as pd

    if "rejection_reason" not in adata.obs:
        assign_rejection_reasons(adata, panel_check_result=panel_check_result)

    codes = adata.obs["rejection_reason"].astype(str)
    nz = codes[codes != ""]
    total = int(len(nz))
    counts = nz.value_counts()
    rows = []
    for code, n in counts.items():
        label, desc = REASON_LIBRARY.get(code, (code, ""))
        rows.append({
            "reason": code,
            "label": label,
            "n_cells": int(n),
            "pct_of_untyped": round(100.0 * n / total, 1) if total else 0.0,
            "description": desc,
        })
    df = pd.DataFrame(rows, columns=["reason", "label", "n_cells", "pct_of_untyped", "description"])
    return df.sort_values("n_cells", ascending=False).reset_index(drop=True)


def panel_resolvability_warnings(adata, panel_check_result: dict | None = None) -> list[dict]:
    """Per cell type the panel cannot confidently resolve, with the number of cells affected.

    Returns a list of dicts (weak/absent coverage and confusable pairs), each with a plain-
    language ``message`` and ``n_cells`` = how many cells carry that (candidate) type.
    """
    from . import panel_check as _pc

    pc = panel_check_result or adata.uns.get("panel_check")
    if not pc or not _pc.is_valid(pc):
        return []

    ct = adata.obs["cell_type"].astype(str) if "cell_type" in adata.obs else None
    counts = ct.value_counts().to_dict() if ct is not None else {}

    # Marker COUNT is a weak proxy: a type with few on-panel markers can still be separable
    # (2 highly-discriminative markers beat 8 redundant ones). Defer to the same typability
    # verdict the Panel-check "resolved" column reads, so the two surfaces cannot disagree -
    # suppress the coverage warning for any confidently-typable type. Reconstruct marker sets
    # from pc's on-panel `present` genes (identical identifiability AUC to the full set; the
    # depth-matched-F1 path uses panel genes, not marker sets, so a loaded reference is honoured).
    typable: set[str] = set()
    if ct is not None:
        try:
            ms = {t: d.get("present", []) for t, d in pc.get("coverage", {}).items()}
            typable = {r["cell_type"] for r in _pc.typability_table(
                adata, cluster_key="cell_type", marker_sets=ms,
                reference_match=adata.uns.get("reference_match"),
                panel_check_result=pc) if r.get("confidently_typable")}
        except Exception:
            typable = set()   # fall back to coverage-only (conservative: keep the warning)

    out: list[dict] = []
    for typ, d in pc.get("coverage", {}).items():
        status = d.get("status")
        if status not in ("red", "amber"):
            continue
        if typ in typable:
            continue   # separable despite thin coverage - Panel-check calls it "resolved"
        k, m = d.get("n_present", 0), d.get("n_markers", 0)
        adj = "cannot be resolved" if status == "red" else "is weakly resolved"
        out.append({
            "kind": "coverage",
            "cell_type": typ,
            "status": status,
            "n_markers_present": int(k),
            "n_markers_total": int(m),
            "n_cells": int(counts.get(typ, 0)),
            "message": (f"{typ} {adj} on this panel (only {k} of {m} canonical markers present)."),
        })

    seen_pairs = set()
    for pair in pc.get("confusable_pairs", []):
        p = pair.get("pair", [])
        if len(p) != 2:
            continue
        key = tuple(sorted(p))
        if key in seen_pairs:
            continue
        seen_pairs.add(key)
        a, b = p
        n_affected = int(counts.get(a, 0) + counts.get(b, 0))
        out.append({
            "kind": "confusable",
            "pair": [a, b],
            "n_cells": n_affected,
            "message": (f"Panel cannot separate {a} vs {b} - they share their on-panel markers; "
                        "reported as the coarser group."),
        })

    out.sort(key=lambda r: r["n_cells"], reverse=True)
    return out
