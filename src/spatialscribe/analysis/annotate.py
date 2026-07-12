"""Multi-method cell-type annotation, reconciled to a consensus.

Signals (independent):
  * marker scoring - ``sc.tl.score_genes`` per lineage, argmax per cluster;
  * Claude-as-annotator - top markers per cluster -> reasoned label (``llm``);
  * (optional) reference-transfer arms, when the user supplies a single-cell reference: a
    CellTypist model TRAINED on that reference (:func:`celltypist_transfer`, the reliable in-env
    arm) and/or TACCO OT (:mod:`reference_transfer`). Their per-cell labels reconcile in through
    ``consensus_annotate(method_label_cols=[...])``.

`consensus_annotate` merges them into one label per cluster, flags disagreements, and
writes ``adata.obs['cell_type']``. All scoring on Xenium is panel-restricted first and
uses a small control set (the maxRank=150 lesson: default gene ranks dilute to noise on
<500-gene cells).

Depends on: scanpy, numpy, pandas; :mod:`markers`, :mod:`cluster`, :mod:`llm`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from . import cluster as _cluster
from . import markers as _markers

if TYPE_CHECKING:
    import pandas


def _safe_col(name: str) -> str:
    """h5-safe score column name (no '/' - disallowed as an h5 key)."""
    return "score_" + name.replace("/", "_").replace(" ", "_")


def score_marker_sets(adata, marker_sets: dict[str, list[str]] | None = None) -> dict[str, str]:
    """Add one score column per marker set (panel-restricted). Returns {label: column}."""
    import scanpy as sc

    from .markers import present

    import hashlib

    from .backend import get_backend

    marker_sets = marker_sets or _markers.LINEAGE_MARKERS
    panel = set(adata.var_names)
    present_sets = {name: genes for name, genes in present(panel, marker_sets).items() if genes}
    colmap: dict[str, str] = {name: _safe_col(name) for name in present_sets}
    # score_genes is the single hottest call on the annotate path (several times per step: consensus +
    # confidence, each recomputing). Cache it: skip when the SAME marker sets were already scored on
    # this adata and their score columns are still present. The backend routes the compute to rapids on
    # GPU; ctrl_size stays modest for small panels (the maxRank=150 spirit).
    sig = hashlib.sha1(
        repr(sorted((n, tuple(sorted(g))) for n, g in present_sets.items())).encode()
    ).hexdigest()
    if adata.uns.get("_score_marker_sig") == sig and all(c in adata.obs for c in colmap.values()):
        return colmap
    be = get_backend()
    for name, genes in present_sets.items():
        be.score_genes(adata, gene_list=genes, score_name=colmap[name], ctrl_size=min(50, adata.n_vars))
    adata.uns["_score_marker_sig"] = sig
    return colmap


def marker_labels(adata, cluster_key: str = "leiden",
                  marker_sets: dict[str, list[str]] | None = None) -> dict[str, str]:
    """Argmax marker-score per cluster -> {cluster: lineage label} (original names kept)."""
    import numpy as np

    colmap = score_marker_sets(adata, marker_sets)
    if not colmap:
        return {}
    labels = {}
    for cl, sub in adata.obs.groupby(cluster_key, observed=True):
        means = {name: float(np.mean(sub[col])) for name, col in colmap.items()}
        labels[str(cl)] = max(means, key=lambda k: means[k])
    return labels


def claude_labels(adata, cluster_key: str = "leiden", n_markers: int = 15,
                  context: str = "human skin / melanoma",
                  allowed_labels: list[str] | None = None) -> dict[str, dict]:
    """Claude-as-annotator over each cluster's top marker genes.

    ``allowed_labels`` closes the vocabulary in the prompt (see :func:`llm.annotate_clusters`).
    """
    from . import llm

    groups = list(adata.obs[cluster_key].cat.categories) if hasattr(
        adata.obs[cluster_key], "cat") else sorted(adata.obs[cluster_key].unique())
    per_cluster = {str(g): _cluster.top_markers(adata, str(g), n=n_markers) for g in groups}
    return llm.annotate_clusters(per_cluster, context=context, allowed_labels=allowed_labels)


def _canonicalise(label, allowed: dict[str, str]) -> str | None:
    """Map a model's label onto the allowed vocabulary, or None when it is off-vocabulary.

    Case- and whitespace-tolerant EXACT matching only. Deliberately no fuzzy/nearest-string snapping:
    an invented label ("pDC", "RBC") is near nothing in a ~10-lineage space and would snap to a wrong
    lineage, converting an honest "I don't know" into a confident error.
    """
    if label is None:
        return None
    return allowed.get(str(label).strip().lower())


def consensus_annotate(adata, cluster_key: str = "leiden", use_llm: bool = True,
                       context: str = "human skin / melanoma",
                       marker_sets: dict[str, list[str]] | None = None,
                       method_label_cols: list[str] | None = None,
                       reliability_weights: dict | None = None,
                       progress=None) -> "pandas.DataFrame":  # noqa: F821
    """Reconcile marker + Claude labels into ``adata.obs['cell_type']``.

    Returns a per-cluster table (marker_label, claude_label, confidence, agreement, conflict,
    off_vocabulary, final) for the UI.

    The LLM's job is to NAME a cluster within a closed vocabulary (the curated lineages on this
    section plus ``llm.NOVEL_LABEL``), not to adjudicate biology:

      * agreeing label      -> used;
      * ``Novel / unknown`` -> used (an honest "no lineage fits this panel");
      * off-vocabulary      -> discarded, marker argmax used, ``off_vocabulary=True``
        (never snapped to the nearest curated string - see :func:`_canonicalise`);
      * a DIFFERENT valid lineage -> ``conflict=True``, marker argmax used, and every cell in that
        cluster is flagged in ``obs['label_conflict']`` so :func:`apply_confidence` caps its verdict
        below PASS.

    Rationale: on the 2026-07 five-section benchmark (two LLMs naming the SAME leiden clusters) 44 of
    59 disagreements were naming-only, which the closed vocabulary removes outright; the remaining 15
    were real biology and nothing establishes the LLM as more accurate than the marker argmax there.

    SCOPE OF THE CLOSED VOCABULARY: it constrains the LLM naming step only. When ``method_label_cols``
    is given, the per-cell consensus vote below OVERWRITES ``cell_type`` with the winning reference
    label (RCTD/SingleR/scANVI/panhumanpy), which carries the reference's own vocabulary - deliberately,
    since a reference's label space is itself controlled and usually finer than these ~10 lineages.
    So ``cell_type`` is closed over the curated keys only on the marker+LLM path.

    When ``method_label_cols`` (a list of ``obs`` column names) is given, after the cluster
    labelling above ALSO reconciles per-cell method labels (e.g. RCTD/SingleR/scANVI/
    panhumanpy/TACCO): for each cell, takes the majority label across
    ``[cluster cell_type] + method_label_cols`` (ignoring NaN/None), writes it to
    ``obs['cell_type']``, and records the winning fraction in ``obs['consensus_agreement']``.
    Ties keep the cluster label if it is among the tied winners, else the alphabetically
    first tied label (deterministic). Columns absent from ``obs`` are skipped, not an error.
    Default (``method_label_cols=None``) leaves this entire block a no-op - the pre-existing
    cluster-only behavior is unchanged.

    ``reliability_weights`` (opt-in) switches that per-cell vote from a naive majority to a
    RELIABILITY-weighted one - see :func:`consensus.weighted_vote`. Pass measured weights from
    :func:`consensus.reliability_from_labels`, or the documented :data:`consensus.DEFAULT_RELIABILITY`
    prior. Needed because a naive majority loses to the single best method when one method dominates:
    on the annotator_bench, majority scored 0.858 against RCTD's 0.923. ``None`` (the default) keeps
    the majority path exactly as it was.
    """
    import pandas as pd

    from . import llm as _llm

    if progress:
        progress(0.15, "scoring lineage markers")
    m_lab = marker_labels(adata, cluster_key, marker_sets=marker_sets)

    # Closed vocabulary: every lineage that was SCORED on this section, plus the escape hatch.
    # Not the argmax winners - a lineage that no cluster happened to win must still be nameable, or
    # the LLM is forbidden from correcting an argmax that missed it.
    vocab = list(marker_sets or {}) or list(score_marker_sets(adata, marker_sets))
    allowed = list(dict.fromkeys([*vocab, _llm.NOVEL_LABEL]))
    canon = {a.strip().lower(): a for a in allowed}

    if progress:
        progress(0.45, "reconciling methods")
        progress(0.75, "asking Claude to name the clusters" if use_llm
                 else "assigning marker-argmax labels")
    c_lab = claude_labels(adata, cluster_key, context=context,
                          allowed_labels=allowed) if use_llm else {}

    rows, mapping, conflicted = [], {}, set()
    groups = [str(g) for g in (adata.obs[cluster_key].cat.categories
              if hasattr(adata.obs[cluster_key], "cat") else sorted(adata.obs[cluster_key].unique()))]
    for cl in groups:
        marker = m_lab.get(cl, "Unknown")
        cl_info = c_lab.get(cl, {})
        raw = cl_info.get("label") if use_llm else None
        claude = _canonicalise(raw, canon)               # None => off-vocabulary (or no LLM)
        off_vocab = bool(use_llm and raw is not None and claude is None)

        # The LLM NAMES clusters; it does not adjudicate biology. It can (a) confirm the marker
        # argmax, (b) say "Novel / unknown", or (c) be ignored. A genuine lineage disagreement is
        # FLAGGED and resolved toward the computed marker evidence, because nothing shows the LLM is
        # more accurate than the argmax: on the 2026-07 benchmark 44/59 cross-model disagreements
        # were naming-only and the remaining 15 were never adjudicated against truth.
        real_marker = marker not in ("Unknown",)
        conflict = bool(claude and real_marker and claude != marker and claude != _llm.NOVEL_LABEL)
        if conflict or off_vocab or not claude:
            final = marker                                # computed evidence wins / is the fallback
        else:
            final = claude
        if conflict:
            conflicted.add(cl)

        mapping[cl] = final
        rows.append({
            "cluster": cl, "marker_label": marker, "claude_label": raw if raw is not None else marker,
            "confidence": cl_info.get("confidence", "n/a"),
            "agreement": bool(claude is not None and claude == marker),
            "conflict": conflict, "off_vocabulary": off_vocab, "final": final,
            "rationale": cl_info.get("rationale", ""),
        })

    # Cells with no cluster label (leiden NaN = low-signal cells held out of the embedding) map to
    # NaN; label them "Unassigned" rather than leaving a NaN category, which would surface as a "nan"
    # cell type in the legend/barplot and crash squidpy neighborhood-enrichment. The low-signal gate in
    # apply_confidence still overrides them to "Unassigned: low quality" in cell_type_final.
    adata.obs["cell_type"] = (
        adata.obs[cluster_key].astype(str).map(mapping).fillna("Unassigned").astype("category")
    )
    # Per-cell flag for clusters where the LLM named a DIFFERENT lineage than the marker argmax.
    # apply_confidence reads it to cap the verdict below PASS: a conflict is a reason to look, not a
    # reason to trust either side.
    adata.obs["label_conflict"] = (
        adata.obs[cluster_key].astype(str).isin(conflicted).to_numpy()
    )

    if method_label_cols:
        import numpy as np

        from . import consensus as _cons
        cols = [c for c in method_label_cols if c in adata.obs.columns]
        voters = ["cell_type"] + cols
        votes = adata.obs[voters].astype("object")
        winners, agree = [], []
        for _, row in votes.iterrows():
            vals = [v for v in row.to_numpy() if v is not None and not (isinstance(v, float) and np.isnan(v))]
            if not vals:
                winners.append(None); agree.append(np.nan); continue
            cluster_lab = str(row["cell_type"]) if row["cell_type"] is not None else None
            if reliability_weights:
                win, frac = _cons.weighted_vote(row.to_dict(), reliability_weights, prefer=cluster_lab)
            else:
                counts: dict[str, int] = {}
                for v in vals:
                    counts[str(v)] = counts.get(str(v), 0) + 1
                top = max(counts.values())
                tied = sorted(k for k, c in counts.items() if c == top)
                win = cluster_lab if cluster_lab in tied else tied[0]
                frac = counts[win] / len(vals)
            winners.append(win)
            agree.append(frac)
        adata.obs["cell_type"] = pd.Categorical([w if w is not None else "Unassigned" for w in winners])
        adata.obs["consensus_agreement"] = np.asarray(agree, dtype=float)

    if progress:
        progress(1.0, "done")
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Per-cell confidence + abstention (QC Layer 5; see docs/research/cell-annotation-qc.md).
# What this score is WORTH is measured in docs/research/confidence-calibration.md: it is an ordinal
# heuristic whose within-cell-type AUC is 0.48-0.57 across three ground-truthed sections, so it barely
# ranks cells by correctness. Calibrate it against real labels (analysis/calibration.py) before reading
# it as a probability, and never quote an ECE improvement without the base-rate null.
# --------------------------------------------------------------------------- #
# Fallback verdict cutoffs. The AUTHORITATIVE values live in the tunable thresholds YAML,
# read via config.get("layer5_confidence", "composite_confidence", ...); these constants are
# used only when the YAML / PyYAML is unavailable (keeps the demo offline-deterministic).
CONF_WARN = 0.5     # below -> tentative (greyed)
CONF_FAIL = 0.25    # below -> abstain
# Output-usability gate: abstention rates above these fractions make an annotation warn / unusable.
USABILITY_WARN = 0.25
USABILITY_FAIL = 0.50
_ABSTAIN = {
    "low_quality": "Unassigned: low quality",
    "ambiguous_mixed": "Ambiguous: mixed",
    "unresolvable_panel": "Unresolvable: panel",
    "uncertain_lowconf": "Uncertain: low confidence",
    "label_conflict": "Ambiguous: label conflict",
    "novel": "Novel / unknown",
}

# --------------------------------------------------------------------------- #
# Abstention vocabulary - THE single source of truth.
#
# These strings live in the same pd.Categorical as real cell types (`cell_type_final`), so every
# consumer that groups by cell type must exclude or collapse them. They are NOT lineages: grouping the
# DELIVERED composition figure on `cell_type_final` printed "Ambiguous: mixed" in the legend as if it
# were biology, while `niches`/`nhood_enrichment` grouped on `cell_type` - two different denominators
# in one report. Match EXACTLY, never by prefix: "Novel epithelial subtype X" is a real label.
# --------------------------------------------------------------------------- #
NOT_ASSIGNED = "Not assigned"           # the one category abstained cells collapse into
ABSTENTION_LABELS: frozenset = frozenset(_ABSTAIN.values()) | {
    "Unassigned",       # leiden-NaN cells (consensus_annotate)
    "Unknown",          # marker_labels fallback for a cluster with no marker signal
    NOT_ASSIGNED,       # the collapsed category itself
}
# Kept for the backend's category filter; DERIVED, never re-typed by hand.
ABSTENTION_PREFIXES: tuple = tuple(sorted({s.split(":")[0] for s in ABSTENTION_LABELS}))


def is_abstention(label) -> bool:
    """True when ``label`` is an abstention pseudo-label rather than a cell type (exact match)."""
    return str(label) in ABSTENTION_LABELS


def abstention_mask(labels):
    """Boolean numpy mask over an array/Series of label strings."""
    import numpy as np

    return np.asarray([is_abstention(v) for v in np.asarray(labels, dtype=object)], dtype=bool)


def annotation_key(adata) -> str:
    """The ONE column every cell-type consumer must group on.

    ``cell_type_final`` (labels + the abstention overlay) once ``apply_confidence`` has run, else the
    raw ``cell_type``. Consumers used to disagree, which is the bug this exists to prevent.
    """
    return "cell_type_final" if "cell_type_final" in adata.obs else "cell_type"


def collapse_abstention(labels, other: str = NOT_ASSIGNED):
    """Merge every abstention pseudo-label into ONE ``other`` category.

    For consumers where abstained cells must stay in the data (they occupy space, so the spatial
    neighbour graph needs them) but must not read as several distinct lineages.
    """
    import pandas as pd

    s = pd.Series(labels).astype(str)
    vals = [other if is_abstention(v) else v for v in s]
    return pd.Series(pd.Categorical(vals), index=s.index)


def typed_mask(adata, key: str | None = None):
    """Boolean mask selecting cells carrying a REAL cell-type label (abstained cells excluded)."""
    key = key or annotation_key(adata)
    return ~abstention_mask(adata.obs[key].astype(str).to_numpy())


def _combine_soft(factors):
    """Combine soft down-weight factors into ONE multiplier by AVERAGING (neutral 1.0 if none).

    Averaging, not multiplying, is deliberate. These factors (marker purity, spatial coherence,
    label stability, reference posteriors) are independent, individually mild corroborators, each
    in ~[0.65, 1]. Multiplying them compounds - a stack of three 0.7-0.9 factors craters an
    otherwise-good call (0.7^3 ~= 0.34) - which, together with a structurally-broken PMP, made
    68% of a real 5K breast section abstain. Averaging preserves each signal's pull without letting
    a pile of mild penalties collapse a confident call.
    """
    import numpy as np

    if not factors:
        return 1.0
    return np.mean(np.vstack([np.asarray(f, dtype=float) for f in factors]), axis=0)


def usability_flag(pct_abstain: float) -> str:
    """Output-usability gate: ``"ok" | "warn" | "fail"`` from the abstention rate.

    Input QC answers "is this section loadable?", never "did we actually label it?". A Prime 5K run
    completed with ``qc_flags = ok`` while abstaining 94.8% of its cells - completed is not usable.
    """
    if pct_abstain > USABILITY_FAIL:
        return "fail"
    return "warn" if pct_abstain > USABILITY_WARN else "ok"


def annotation_completeness(adata, cluster_key: str = "cell_type",
                            panel_check_result: dict | None = None) -> dict | None:
    """Reconcile which lineages the PANEL can resolve against which the annotation actually PRODUCED.

    A lineage with green marker coverage (the panel carries its markers) that annotates to ZERO cells
    is worth flagging: on a deep section it can mean an annotation collapse (e.g. all T cells folding
    into Myeloid), which nothing else in the pipeline surfaces - "10/10 typable, usability ok" can sit
    on top of an annotation that recovered only 4 of 10 lineages. Informational only (a green-coverage
    lineage may also be genuinely absent), so it lists the gap; it never mutates or hard-fails.

    Returns ``{"resolvable_absent": [...], "n": int, "note": str}`` or ``None`` (no valid panel_check).
    """
    from . import panel_check as _pc

    pc = panel_check_result if _pc.is_valid(panel_check_result) else adata.uns.get("panel_check")
    if not _pc.is_valid(pc) or cluster_key not in getattr(adata, "obs", {}):
        return None
    cov = pc.get("coverage", {})
    present = {str(x) for x in adata.obs[cluster_key].astype(str).unique()}
    resolvable_absent = sorted(
        str(ct) for ct, d in cov.items()
        if isinstance(d, dict) and d.get("status") == "green"
        and str(ct) not in present and not is_abstention(str(ct)))
    return {"resolvable_absent": resolvable_absent, "n": len(resolvable_absent),
            "note": ("These lineages have green panel coverage (their markers are on-panel) yet "
                     "annotated to 0 cells. On a deep section this can indicate an annotation collapse "
                     "(a lineage folded into another); on a shallow one they may be genuinely absent. "
                     "Verify the resolvable-but-absent lineages against the composition.")}


def apply_confidence(adata, cluster_key: str = "cell_type",
                     marker_sets: dict[str, list[str]] | None = None,
                     panel_check_result: dict | None = None,
                     n_panel_genes: int | None = None,
                     use_calibrated: bool = False,
                     lineage_markers: dict[str, list[str]] | None = None) -> dict:
    """Attach a calibrated per-cell confidence + verdict + abstention reason.

    Penalized-posterior recipe (docs section 10.2): a **margin-based** posterior over the
    lineage marker scores, multiplied by a contamination penalty (control % + CRISP
    impurity), a panel-coverage ceiling (green 1.0 / amber 0.6 / red -> abstain), and a
    **panel-indexed** low-signal gate (``counts < panel_indexed_floor`` or ``genes <
    qc.RICH_PANEL_MIN_GENES`` on a rich panel, ``genes < 5`` on a targeted one).
    The count floor comes from ``qc.panel_indexed_floor`` - the SAME source Layer 2 filtering
    uses - so a targeted panel gates at a fixed 10 while a Prime 5K / WTA section uses its
    section-relative distributional floor (never the fixed 10 that would abstain ~60% of 5K
    cells). Writes ``annotation_confidence``, ``annotation_verdict`` (PASS/WARN/FAIL),
    ``annotation_reason`` and ``cell_type_final`` (abstention labels override on FAIL). Returns
    the section annotatability headline.
    """
    import numpy as np
    import pandas as pd

    from . import purity

    marker_sets = marker_sets or _markers.LINEAGE_MARKERS
    colmap = score_marker_sets(adata, marker_sets)
    names = list(colmap)
    if not names:
        # No on-panel lineage markers scored: abstain on every cell (the honest call), and still
        # write the four annotation columns so annotate's declared products hold on a barren panel.
        n = adata.n_obs
        adata.obs["annotation_confidence"] = np.zeros(n)
        adata.obs["annotation_verdict"] = np.full(n, "FAIL", dtype=object)
        adata.obs["annotation_reason"] = np.full(n, "unresolvable_panel", dtype=object)
        adata.obs["cell_type_final"] = pd.Categorical([_ABSTAIN["unresolvable_panel"]] * n)
        return {"pct_pass": 0.0, "pct_warn": 0.0, "pct_abstain": 1.0,
                "top_abstention_reasons": {"unresolvable_panel": n}, "mean_confidence": 0.0}
    S = adata.obs[[colmap[n] for n in names]].to_numpy(dtype=float)
    n = S.shape[0]

    order = np.argsort(-S, axis=1)
    rows = np.arange(n)
    top1 = S[rows, order[:, 0]]
    top2 = S[rows, order[:, 1]] if S.shape[1] > 1 else np.zeros(n)
    margin = top1 - top2
    # Margin -> posterior via a logistic scaled by the margin spread (arm-relative).
    spread = float(np.std(margin)) + 1e-9
    posterior = 1.0 / (1.0 + np.exp(-margin / spread))

    # Contamination penalty: negative-control % and CRISP impurity. `pct_counts_control` MUST come
    # from the strict negative-control mask (qc.compute_qc): computed over the broad panel mask it
    # counts deprecated-codeword SIGNAL, hits this 5% cut on every Prime 5K cell, and abstains the
    # section. See io.build_neg_control_mask.
    #
    # CRISP impurity is a MUTUAL-EXCLUSIVITY spillover check, so it needs COARSE, DISJOINT lineages -
    # NOT the fine annotation vocabulary. Feeding it 30 fine subtypes whose markers overlap (e.g. a gut
    # atlas' Stromal_1/2/3, Arterial/Venous endothelial, Stem/TA) makes almost every cell "co-express
    # multiple lineages" and flags ~80% of the section as spillover (the observed artifact). Use the
    # curated coarse lineage set (`lineage_markers`, e.g. markers.for_tissue) when provided.
    # CRISP impurity ("mixed lineages / spillover"). The vendored crisp_purity flags a cell on >= 1
    # DETECTED marker from >= 2 lineages - far too loose for single-molecule imaging: one stray ambient
    # transcript trips it, so ~2/3 of a dense Xenium section is called spillover regardless of how coarse
    # the lineages are (the observed 73-89% duodenum artifact). Recompute it SCORE-based on the COARSE
    # disjoint lineages: a cell is mixed only when it is STRONGLY positive (z > 1 across cells) for >= 2
    # lineages at once - a genuine co-expression, not a lone transcript. Falls back to the loose flag if
    # the coarse lineages cannot be scored (no on-panel markers). Robust to score sign (z-standardized).
    _lm = lineage_markers or marker_sets
    purity.crisp_purity(adata, _lm)
    try:
        _lcols = score_marker_sets(adata, _lm)
        if len(_lcols) >= 2:
            _M = adata.obs[list(_lcols.values())].to_numpy(dtype=float)
            _Z = (_M - _M.mean(0)) / (_M.std(0) + 1e-9)
            adata.obs["crisp_impure"] = (_Z > 1.0).sum(1) >= 2
    except Exception:
        pass
    pct_ctrl = np.nan_to_num(np.asarray(adata.obs.get("pct_counts_control", np.zeros(n)), dtype=float), nan=0.0)
    contam = np.clip(pct_ctrl / 5.0, 0, 1)
    contam = np.maximum(contam, adata.obs["crisp_impure"].to_numpy(dtype=float) * 0.5)

    # Panel-coverage ceiling for each cell's argmax lineage. Guard against a panel_check that was
    # corrupted by an h5ad round-trip (the '/' in cell-type keys splits the nested dict, so the
    # coverage values lose their 'status' field) - in that case skip the ceiling (default 1.0).
    ceiling = np.ones(n)
    if panel_check_result:
        cov = panel_check_result.get("coverage", {})
        status = {ct: d["status"] for ct, d in cov.items() if isinstance(d, dict) and "status" in d}
        if status:
            ceil_map = {"green": 1.0, "amber": 0.6, "red": 0.0}
            top_names = [names[i] for i in order[:, 0]]
            ceiling = np.array([ceil_map.get(status.get(tn, "green"), 1.0) for tn in top_names])

    # Base confidence: the HARD evidence (marker margin x contamination x panel-coverage ceiling)
    # stays multiplicative - a genuinely amber panel or a real contaminant SHOULD cap the call.
    confidence = np.clip(posterior * (1 - contam) * ceiling, 0, 1)

    # Soft corroborating penalties (Layers 3/5/6 + optional reference posteriors), each in ~[0.65,1]
    # and applied only when its column exists. They are AVERAGED into ONE factor via _combine_soft,
    # NOT multiplied: multiplying independent mild penalties compounds and made PASS nearly
    # unreachable on real 5K panels (see _combine_soft). A NaN (e.g. PMP with no marker signal) is
    # treated as neutral.
    soft = []
    if "pmp" in adata.obs:
        pmp = np.nan_to_num(np.asarray(adata.obs["pmp"], dtype=float), nan=1.0)
        soft.append(0.7 + 0.3 * np.clip(pmp, 0, 1))                          # marker purity
    if "spatial_coherence" in adata.obs:
        coh = np.nan_to_num(np.asarray(adata.obs["spatial_coherence"], dtype=float), nan=1.0)
        # Spatial incoherence is a CORROBORATOR, not a primary penalty. A lone cell is only suspect when
        # an INDEPENDENT mixed-signal cue says its transcripts bled in from neighbours (CRISP impurity, or
        # an ovrlpy VSI vertical-overlap doublet - the actual spillover mechanism). A clean, isolated cell
        # is real biology (an infiltrating / rare cell) and keeps full confidence, so honest QC never greys
        # out the immune infiltration the app is built to surface (cf. spatial.immune_exclusion). Location
        # is a prior, not evidence about identity. Kept in the averaged soft-factor list (neutral 1.0 for
        # clean cells) rather than a standalone multiply, per the _combine_soft design.
        mixed = (adata.obs["crisp_impure"].to_numpy(dtype=bool)
                 if "crisp_impure" in adata.obs else np.zeros(n, dtype=bool))
        if "vsi" in adata.obs:
            mixed = mixed | (np.nan_to_num(np.asarray(adata.obs["vsi"], dtype=float), nan=1.0) < 0.7)
        soft.append(np.where(mixed, 0.7 + 0.3 * np.clip(coh, 0, 1), 1.0))    # coherence penalty ONLY when corroborated
    if "annotation_stability" in adata.obs:
        stab = np.nan_to_num(np.asarray(adata.obs["annotation_stability"], dtype=float), nan=0.0)
        soft.append(1.0 - 0.35 * np.clip(stab, 0, 1))                       # label-flip stability
    if "scanvi_confidence" in adata.obs:
        scf = np.nan_to_num(np.asarray(adata.obs["scanvi_confidence"], dtype=float), nan=1.0)
        soft.append(0.7 + 0.3 * np.clip(scf, 0, 1))                          # ref-posterior max-prob
    if "rctd_weight" in adata.obs:
        rw = np.nan_to_num(np.asarray(adata.obs["rctd_weight"], dtype=float), nan=1.0)
        soft.append(0.7 + 0.3 * np.clip(rw, 0, 1))                           # RCTD doublet weight
    confidence = np.nan_to_num(np.clip(confidence * _combine_soft(soft), 0, 1), nan=0.0)

    # popV consensus (Kimmel/Ergen/Yosef, Nat Genet 2024): cross-method AGREEMENT is the trusted,
    # well-calibrated uncertainty signal - far better than any single method's self-score (which popV
    # showed is "futile" to weight by). Applied as a PRIMARY multiplicative factor (steeper than the
    # soft corroborators above), but ONLY when a diverse ensemble voted (>= consensus.MIN_TRUST_METHODS);
    # popV's calibration depends on ensemble diversity, so a 2-method demo consensus stays display-only.
    ensemble_frac = 0.0     # fraction of cells whose abstention is backed by a trusted (>=3) ensemble
    if "consensus_agreement" in adata.obs:
        from . import consensus as _cons
        agr = np.nan_to_num(np.asarray(adata.obs["consensus_agreement"], dtype=float), nan=1.0)
        n_meth = np.asarray(adata.obs.get("consensus_n_methods", np.zeros(n)), dtype=float)
        trusted = n_meth >= _cons.MIN_TRUST_METHODS
        ensemble_frac = float(trusted.mean())
        cfac = np.where(trusted, 0.5 + 0.5 * np.clip(agr, 0, 1), 1.0)      # agreement 1->1.0, 0->0.5
        confidence = np.nan_to_num(np.clip(confidence * cfac, 0, 1), nan=0.0)

    # Low-signal hard gate on RAW counts. On a processed section ``adata.X`` (and thus
    # ``obs['total_counts']``) is normalized, so read the raw counts layer when present -
    # otherwise a normalized total mis-fires the ``counts < 10`` gate and abstains everything.
    if "counts" in adata.layers:
        raw = adata.layers["counts"]
        counts = np.asarray(raw.sum(1)).ravel().astype(float)
        genes = np.asarray((raw > 0).sum(1)).ravel().astype(float)
    else:
        counts = np.asarray(adata.obs.get("total_counts", np.zeros(n)), dtype=float)
        genes = np.asarray(adata.obs.get("n_genes_by_counts", np.zeros(n)), dtype=float)
    # Panel-indexed low-signal floor (NOT a fixed counts<10): shared with qc.suggest_count_floor
    # so Layer 2 and Layer 5 can't diverge. Fixed 10 on targeted panels; distributional on 5K/WTA.
    from . import qc as _qc
    n_panel = _qc.resolve_panel_size(adata, n_panel_genes)
    floor, _floor_mode = _qc.panel_indexed_floor(counts, n_panel)
    # Gene floor matches cluster._low_signal_mask on rich panels, so a cell held OUT of the embedding
    # (NaN leiden) is consistently abstained here as low_quality rather than getting a spurious label.
    gene_floor = _qc.RICH_PANEL_MIN_GENES if n_panel >= 1000 else 5
    low = (counts < floor) | (genes < gene_floor)

    from . import calibration, config
    warn_thr = float(config.get("layer5_confidence", "composite_confidence", "warn", default=CONF_WARN))
    fail_thr = float(config.get("layer5_confidence", "composite_confidence", "fail", default=CONF_FAIL))

    # The score the PASS/WARN/FAIL gate is cut on. By default this is the raw heuristic, whose thresholds
    # are tuned for it. `use_calibrated=True` cuts on the isotonic-calibrated score instead, so the cutoffs
    # read as real probabilities - but only once `calibration.calibrate_confidence` has fit one against real
    # labels. The raw score is still what gets written to obs; it is never overwritten by the calibrated one.
    gate = confidence
    if use_calibrated and calibration.CALIBRATED_KEY in adata.obs:
        gate = np.nan_to_num(np.asarray(adata.obs[calibration.CALIBRATED_KEY], dtype=float), nan=0.0)

    verdict = np.where(
        low | (gate < fail_thr), "FAIL",
        np.where(gate < warn_thr, "WARN", "PASS"),
    )

    # A cluster where the LLM named a different lineage than the marker argmax is UNRESOLVED, not
    # confidently typed. Cap it at WARN (tentative/greyed) rather than abstain the whole cluster: the
    # marker evidence still stands, it just has not been corroborated. Never silently PASS.
    conflict = (np.asarray(adata.obs["label_conflict"], dtype=bool)
                if "label_conflict" in adata.obs else np.zeros(n, dtype=bool))
    verdict = np.where(conflict & (verdict == "PASS"), "WARN", verdict)

    reason = np.full(n, "", dtype=object)
    fail = verdict == "FAIL"
    reason[low] = "low_quality"
    reason[fail & ~low & (ceiling == 0)] = "unresolvable_panel"
    reason[fail & ~low & (ceiling > 0) & (contam >= 0.5)] = "ambiguous_mixed"
    reason[fail & ~low & (ceiling > 0) & (contam < 0.5)] = "uncertain_lowconf"
    # A conflicted cell that ALSO fails on its own merits reports the conflict, which is the more
    # actionable reason (re-cluster / check the markers) than a generic low-confidence note.
    reason[fail & ~low & conflict] = "label_conflict"

    adata.obs["annotation_confidence"] = confidence
    adata.obs["annotation_verdict"] = verdict
    adata.obs["annotation_reason"] = reason

    base = adata.obs[cluster_key].astype(str).to_numpy()
    final = base.copy().astype(object)
    final[fail] = [_ABSTAIN.get(r, "Unassigned") for r in reason[fail]]
    adata.obs["cell_type_final"] = pd.Categorical(final)

    reasons = pd.Series(reason[fail]).value_counts().head(3).to_dict()
    pct_abstain = float(fail.mean())
    # Output-usability gate. Stamped on uns so qc.qc_summary can surface it next to the input-QC
    # flags - a run that "completed with acceptable QC" but abstained most of its cells is a failure.
    usability = usability_flag(pct_abstain)
    adata.uns["annotation_usability"] = {"flag": usability, "pct_abstain": pct_abstain,
                                         "top_abstention_reasons": reasons}
    # Completeness: reconcile panel-resolvable lineages vs the ones actually annotated (a resolvable
    # lineage that produced 0 cells can be an annotation collapse). Informational; never mutates.
    try:
        comp = annotation_completeness(adata, cluster_key, panel_check_result)
        if comp is not None:
            adata.uns["annotation_completeness"] = comp
    except Exception:
        pass
    # What the abstention actually RESTS on. Measured (2026-07-10 GT benchmark, the internal docs):
    # the per-cell confidence heuristic barely ranks correct cells (within-lineage AUC ~0.54), so
    # abstention gated on it alone is near-random; cross-method AGREEMENT does rank them (~0.77). So a
    # section with no reference ensemble has HEURISTIC (untrustworthy) abstention - say so, don't imply
    # the greying is reliable. Running >=3 reference annotators (the `annotation_methods` capability)
    # upgrades it. `ensemble_frac` > 0 means a trusted ensemble weighted at least some cells.
    ensemble_backed = ensemble_frac > 0.0
    adata.uns["abstention_basis"] = {
        "basis": "ensemble_agreement" if ensemble_backed else "confidence_heuristic",
        "frac_ensemble_backed": ensemble_frac,
        "trustworthy": ensemble_backed,
        "note": ("Abstention is gated on cross-method agreement (>=3 reference annotators), which ranks "
                 "correct cells within lineage (AUC ~0.77, 2026-07-10 GT benchmark)."
                 if ensemble_backed else
                 "Abstention rests on the confidence heuristic alone - no reference ensemble voted. That "
                 "heuristic barely ranks correct cells (within-lineage AUC ~0.54), so treat the greying as "
                 "advisory, not reliable. Run >=3 reference annotators (annotation_methods) for "
                 "agreement-based abstention."),
    }
    return {
        "pct_pass": float((verdict == "PASS").mean()),
        "pct_warn": float((verdict == "WARN").mean()),
        "pct_abstain": pct_abstain,
        "top_abstention_reasons": reasons,
        "mean_confidence": float(confidence.mean()),
        "usability": usability,
        # Whether abstention is ensemble-backed (trustworthy) or heuristic-only (near-random) - see uns.
        "abstention_basis": adata.uns["abstention_basis"]["basis"],
        # Fraction of cells in clusters where the LLM and the marker argmax named different lineages.
        "pct_label_conflict": float(conflict.mean()),
    }


# --------------------------------------------------------------------------- #
# Reference-free confidence proxy: label-flip rate under subsampling (Layer 5)
# --------------------------------------------------------------------------- #
def _subsample_counts(counts, drop_frac, rng):
    """Binomial-thin a counts matrix (drop ~drop_frac of transcripts). Returns a dense float array."""
    import numpy as np

    arr = counts.toarray() if hasattr(counts, "toarray") else np.asarray(counts, dtype=float)
    keep_p = max(0.0, 1.0 - float(drop_frac))
    return rng.binomial(arr.astype(int), keep_p).astype("float32")


def _argmax_marker_label(adata, sets, names):
    """Per-cell argmax lineage over score_genes columns (log-normalized adata). Returns str array.

    ``n_bins`` is capped to the gene pool size (scanpy's default of 25 leaves ~1 gene per
    bin on small panels, so no non-list gene ever shares a bin with a marker gene and
    ``score_genes`` raises "No control genes found in any cut").
    """
    import numpy as np

    from .backend import get_backend

    be = get_backend()
    n_bins = min(25, max(2, adata.n_vars // 2))
    S = np.zeros((adata.n_obs, len(names)), dtype=float)
    for j, name in enumerate(names):
        be.score_genes(adata, gene_list=sets[name], score_name="_stab",
                       ctrl_size=min(50, adata.n_vars), n_bins=n_bins)
        S[:, j] = adata.obs["_stab"].to_numpy()
    if "_stab" in adata.obs:
        del adata.obs["_stab"]
    return np.array(names)[np.argmax(S, axis=1)]


def annotation_stability(adata, marker_sets: dict | None = None, drop_frac: float = 0.2,
                         reps: int = 5, seed: int = 0) -> None:
    """Layer 5 reference-free confidence proxy: label-flip rate under transcript subsampling.

    Drop ~``drop_frac`` of each cell's transcripts, re-score the marker argmax over ``reps``
    reps, and record how often the label changes from the full-data call. Writes
    ``obs['annotation_stability']`` in [0,1] (higher = less stable). For large sections, run on
    a subsample of cells (this copies + re-normalizes per rep).
    """
    import numpy as np
    import scanpy as sc

    from . import markers as _m
    from .markers import present

    marker_sets = marker_sets or _m.LINEAGE_MARKERS
    panel = set(adata.var_names)
    sets = {k: v for k, v in present(panel, marker_sets).items() if v}
    names = list(sets)
    if len(names) < 2:
        adata.obs["annotation_stability"] = 0.0
        return

    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X

    # Baseline label on a log-normalized copy of the full counts.
    base_ad = adata.copy()
    base_ad.X = counts.toarray().astype("float32") if hasattr(counts, "toarray") else np.asarray(counts, dtype="float32")
    sc.pp.normalize_total(base_ad)
    sc.pp.log1p(base_ad)
    base = _argmax_marker_label(base_ad, sets, names)

    rng = np.random.default_rng(seed)
    flips = np.zeros(adata.n_obs, dtype=float)
    for _ in range(reps):
        tmp = adata.copy()
        tmp.X = _subsample_counts(counts, drop_frac, rng)
        sc.pp.normalize_total(tmp)
        sc.pp.log1p(tmp)
        lab = _argmax_marker_label(tmp, sets, names)
        flips += (lab != base).astype(float)
    adata.obs["annotation_stability"] = flips / reps


# --------------------------------------------------------------------------- #
# Reference-transfer arm: a CellTypist model trained on the user's own reference.
# This is the in-env (no subprocess, no optional dep) reference-transfer method - the
# reliable counterpart to the TACCO / RCTD / SingleR / scANVI paths.
# --------------------------------------------------------------------------- #
def _lognorm_10k(a):
    """A log1p-CPM(1e4)-normalized copy of ``a`` (CellTypist's expected input), from the raw
    counts layer when present, else from ``X``. Never mutates the input."""
    import scanpy as sc

    b = a.copy()
    if "counts" in b.layers:
        b.X = b.layers["counts"].copy()
    sc.pp.normalize_total(b, target_sum=1e4)
    sc.pp.log1p(b)
    return b


def celltypist_transfer(adata, reference, ref_label_key: str, *, max_ref_cells: int = 20000,
                        key_added: str = "celltypist_label") -> dict:
    """Train CellTypist on a custom reference (panel-restricted) and predict per-cell labels.

    The in-env reference-transfer arm: restricts the reference to the section's shared genes,
    subsamples it for tractable training, fits a CellTypist logistic model, and writes the
    predicted per-cell label to ``adata.obs[key_added]`` (+ ``obs[key_added + '_prob']`` when a
    probability matrix is available). Degrades gracefully - returns ``{'status': 'skipped: ...'}``
    and never raises (missing celltypist, too few shared genes, a training error), so it is safe to
    call unconditionally. Feed ``key_added`` into ``consensus_annotate(method_label_cols=[...])`` to
    reconcile it with the marker/Claude labels.
    """
    def _skip(msg: str) -> dict:
        return {"status": f"skipped: {msg}"}

    try:
        import celltypist
        import numpy as np
        import pandas as pd
    except Exception as exc:                                   # pragma: no cover - env-dependent
        return _skip(f"celltypist unavailable ({exc})")

    try:
        if reference is None:
            return _skip("no reference supplied")
        if ref_label_key not in reference.obs:
            return _skip(f"reference has no '{ref_label_key}' label column")
        shared = [g for g in adata.var_names if g in set(map(str, reference.var_names))]
        if len(shared) < 5:
            return _skip(f"only {len(shared)} genes shared with the reference")
        ref = reference
        if ref.n_obs > max_ref_cells:                          # subsample for tractable training
            idx = np.sort(np.random.default_rng(0).choice(ref.n_obs, size=max_ref_cells, replace=False))
            ref = ref[idx]
        ref = ref[:, shared].copy()
        ref.obs[ref_label_key] = ref.obs[ref_label_key].astype(str)
        # drop labels with <2 cells (a classifier cannot learn a singleton class)
        vc = ref.obs[ref_label_key].value_counts()
        keep = ref.obs[ref_label_key].isin(vc[vc >= 2].index).to_numpy()
        ref = ref[keep].copy()
        if ref.obs[ref_label_key].nunique() < 2:
            return _skip("reference has fewer than 2 usable label classes on the shared panel")

        model = celltypist.train(_lognorm_10k(ref), labels=ref_label_key, n_jobs=1,
                                 use_SGD=True, feature_selection=False, check_expression=False)
        sec = _lognorm_10k(adata[:, shared].copy())
        pred = celltypist.annotate(sec, model=model, majority_voting=False)
        col = pred.predicted_labels
        col = col["predicted_labels"] if hasattr(col, "columns") else col
        labels = np.asarray(col).astype(str).ravel()
        adata.obs[key_added] = pd.Categorical(labels)
        prob = getattr(pred, "probability_matrix", None)
        if prob is not None:
            adata.obs[key_added + "_prob"] = np.asarray(prob).max(1)
        return {"status": "ok", "n_labels": int(len(set(labels))),
                "n_shared_genes": len(shared), "n_ref_cells": int(ref.n_obs),
                "coverage": 1.0}
    except Exception as exc:                                   # pragma: no cover - defensive
        return _skip(str(exc))
