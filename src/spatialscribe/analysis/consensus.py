"""popV-style consensus uncertainty - trust cross-method AGREEMENT, not any single method's self-score.

What it does
------------
Ports the load-bearing idea of **popV** (Kimmel/Ergen/Yosef, *Nature Genetics* 2024): across an ensemble
of annotation methods, the number that AGREE on a cell's label (the "consensus score") tracks annotation
accuracy far better than any individual method's self-reported certainty - which popV showed are
"calibrated differently ... futile" to weight by. A LOW consensus is the interpretable flag for the three
cases that need review: an ambiguous continuum state, an out-of-reference / novel type, or a wrong
reference label.

SpatialScribe runs a natural popV ensemble (the marker+Claude cluster label plus per-cell RCTD / SingleR /
scANVI / panhumanpy / TACCO). :func:`consensus_metrics` scores agreement across whichever per-cell
method-label columns are present and writes:

* ``obs['consensus_score']``       - integer count of methods agreeing with the winning label,
* ``obs['consensus_n_methods']``   - number of methods that voted (non-null) for the cell,
* ``obs['consensus_agreement']``   - ``score / n_methods`` in [0, 1],
* ``obs['consensus_reliability']`` - a plain-language bin (very_high / high / moderate / low).

``annotate.apply_confidence`` reads ``consensus_agreement`` as the TRUSTED confidence factor, but only when
a *diverse* ensemble voted (``>= MIN_TRUST_METHODS``) - popV's calibration depends on ensemble diversity,
so a 2-method demo consensus is treated as display-only, not weighted.

Caveats (deliberate)
--------------------
* popV's agreement is **ontology-aware**; we use **exact-label matching**. This is the same simplification
  LatchBio adopted when productionizing popV (simple majority voting performed on par with the ontology
  consensus in their PBMC benchmark), and it avoids needing a cell-ontology graph.
* Running popV *itself* needs a reference + GPU (scANVI/scArches); this module reuses the labels the
  existing subprocess annotators already produce, so it stays reference-free and CPU.

Depends on: numpy, pandas; :mod:`spatialscribe.analysis.config` for tunable bin thresholds. Heavy imports
stay inside functions.
"""
from __future__ import annotations

# popV's calibration needs a diverse ensemble; below this many voting methods we surface the consensus
# for transparency but do NOT weight the confidence by it (see apply_confidence).
MIN_TRUST_METHODS = 3

# Per-method RELIABILITY prior (accuracy in (0,1)), keyed by the obs label column each method writes.
# MEASURED as the MEDIAN lineage accuracy across the three independent-GT sections (CosMx breast 12k,
# Atera WTA breast, Atera WTA cervical) in the 2026-07-10 quality benchmark - superseding the earlier
# n=1 CosMx-only numbers and the two guesses the old prior carried. Per-section acc was:
#   rctd 0.927/0.352/0.834 | singler 0.671/0.568/0.708 | scanvi 0.686/0.307/0.422
#   panhuman 0.274/0.446/0.445 | cell_type(app) 0.694/0.901/0.737
# MEDIAN (robust to the one bad-reference section - Atera WTA breast, whose capped reference tanks every
# reference method) is used as the prior; the mean would down-weight RCTD to 0.70 on that single outlier.
# Two numbers this corrects: scanvi_label was a 0.35 PLACEHOLDER (scANVI had never finished on CPU) - it
# now ran on GPU and measures 0.42; cell_type was a 0.60 knob - the marker+LLM app label is in fact the
# MOST ROBUST method (never catastrophic; median 0.74). Benchmark: quality_jobs_2026-07-10/annotators.md.
# RCTD is bimodal (0.83 median, but 0.35 on a bad reference) - a reference-quality dependency, not a method
# flaw; reference QC should catch a bad reference, not this prior. Pass measured weights from
# `reliability_from_labels` whenever a labeled calibration set exists; the prior is only the no-labels fallback.
DEFAULT_RELIABILITY = {
    "rctd_first_type": 0.83,    # median of 0.927/0.352/0.834 (bimodal: reference-quality dependent)
    "singler_label": 0.67,      # median of 0.671/0.568/0.708 (unchanged - the n=1 value held up)
    "scanvi_label": 0.42,       # MEASURED (was a 0.35 placeholder): median of 0.686/0.307/0.422
    "ph_fine": 0.45,            # median of 0.274/0.446/0.445 (reference-free FM; weak on targeted panels)
    "cell_type": 0.74,          # MEASURED (was a 0.60 knob): median of 0.694/0.901/0.737 (most robust)
}

# Reliabilities are squashed away from 0/1 before the log-odds transform so one voter cannot get infinite
# weight from a rounded 1.0 accuracy.
_RELIABILITY_CLAMP = (0.02, 0.98)


def reliability_bin(agreement: float, n_methods: int) -> str:
    """Map a consensus agreement fraction to a plain-language reliability bin (popV/Latch style).

    Thresholds are tunable via ``config['consensus_popv']['reliability_bins']`` (defaults very_high>=0.9,
    high>=0.7, moderate>=0.5, else low; tightened from the initial 0.875/0.6/0.4 guesses after a real
    reference/query benchmark showed 0.6 agreement was only ~44% accurate - see docs/POPV_UNCERTAINTY.md).
    ``n_methods < 2`` (no ensemble) returns ``'single'``.
    """
    from . import config

    if n_methods is not None and int(n_methods) < 2:
        return "single"
    vh = float(config.get("consensus_popv", "reliability_bins", "very_high", default=0.9))
    hi = float(config.get("consensus_popv", "reliability_bins", "high", default=0.7))
    mo = float(config.get("consensus_popv", "reliability_bins", "moderate", default=0.5))
    if agreement >= vh:
        return "very_high"
    if agreement >= hi:
        return "high"
    if agreement >= mo:
        return "moderate"
    return "low"


def _is_missing(v) -> bool:
    """True for a non-vote (None / NaN / empty / 'nan')."""
    if v is None:
        return True
    try:
        import math

        if isinstance(v, float) and math.isnan(v):
            return True
    except Exception:
        pass
    s = str(v).strip().lower()
    return s in ("", "nan", "none")


def _logit_weight(p: float) -> float:
    """Reliability -> vote weight: ``max(0, log(p / (1-p)))``, the log-odds of the method being right.

    The log-odds transform is load-bearing, not decoration. It is the optimal weight for combining
    independent noisy voters (Nitzan-Paroush): a voter's evidence scales with its log-odds, not its
    accuracy. Weighting by raw accuracy still gets the benchmark case backwards - RCTD 0.923 loses to
    SingleR 0.667 + panhumanpy 0.272 = 0.939 - while log-odds gets it right: 2.44 > 0.71 + 0.
    Clamping at 0 means a method no better than a coin flip is IGNORED rather than allowed to vote
    against the field.

    ponytail: binary log-odds on a multiclass vote is an approximation (the exact multiclass form adds a
    log(K-1) term). Upgrade if reliabilities ever cluster near 0.5 where the term stops being a constant.
    """
    import math

    lo, hi = _RELIABILITY_CLAMP
    p = min(max(float(p), lo), hi)
    return max(0.0, math.log(p / (1.0 - p)))


def weighted_vote(labels: dict, weights: dict, prefer: str | None = None) -> tuple:
    """One cell's RELIABILITY-weighted vote. Returns ``(winner, agreement)``.

    ``labels`` maps voter (an obs column name) -> that voter's label for this cell; missing votes are
    dropped. ``weights`` maps voter -> reliability in (0,1); an absent voter defaults to 0.5, i.e. it is
    ignored. Ties keep ``prefer`` (the cluster label) when it is among the tied winners, else the
    alphabetically first - the same deterministic rule the unweighted path uses.

    ``agreement`` is the UNWEIGHTED fraction of voters backing the winner. Weighting moves the *winner*;
    it never silently redefines ``consensus_agreement``, which keeps its popV meaning.

    Why weight at all, when popV says not to: popV showed weighting by each method's own SELF-reported
    certainty is futile, because those scores are calibrated differently and are not comparable. A
    method's *reliability* is a different quantity - its accuracy measured against labels, on one shared
    scale. The annotator_bench is the evidence that this matters: naive majority scored 0.858, WORSE than
    RCTD alone (0.923), because SingleR and panhumanpy diluted the one method that was right.
    """
    tally: dict[str, float] = {}
    voted: list[str] = []
    for voter, lab in labels.items():
        if _is_missing(lab):
            continue
        lab = str(lab)
        voted.append(lab)
        tally[lab] = tally.get(lab, 0.0) + _logit_weight(weights.get(voter, 0.5))
    if not voted:
        return None, float("nan")

    top = max(tally.values())
    tied = sorted(k for k, v in tally.items() if v >= top - 1e-12)
    win = str(prefer) if (prefer is not None and str(prefer) in tied) else tied[0]
    return win, voted.count(win) / len(voted)


def reliability_from_labels(adata, label_cols, gt_col: str) -> dict:
    """Measure each method's accuracy against a ground-truth obs column -> a weights dict.

    This is the ONLY path allowed to call the weighting calibrated: it is fit against real labels, on the
    same label axis as ``gt_col``. Returns ``{}`` when ``gt_col`` is absent, so the caller falls back to
    :data:`DEFAULT_RELIABILITY` (a documented prior) and says so.
    """
    if gt_col not in adata.obs:
        return {}
    import numpy as np

    gt = adata.obs[gt_col].astype(str).to_numpy()
    gt_ok = np.array([not _is_missing(v) for v in gt])
    out = {}
    for c in label_cols:
        if c not in adata.obs:
            continue
        pred = adata.obs[c].astype(str).to_numpy()
        m = gt_ok & np.array([not _is_missing(v) for v in pred])
        if m.any():
            out[c] = float((pred[m] == gt[m]).mean())
    return out


def consensus_metrics(adata, label_cols, winner_col: str | None = "cell_type") -> dict:
    """popV consensus scoring over the per-cell method-label columns present in ``label_cols``.

    For each cell: the winning label is ``obs[winner_col]`` when present, else the per-cell majority of the
    voters (ties broken alphabetically, deterministic). ``consensus_score`` = number of voters equal to the
    winner; ``consensus_agreement`` = score / (number of non-null voters). Writes the four ``consensus_*``
    columns and returns a summary (mean agreement, reliability-bin counts). No-ops with
    ``{"status": "skipped"}`` when none of ``label_cols`` are present.
    """
    import numpy as np
    import pandas as pd

    present = [c for c in label_cols if c in adata.obs]
    if not present:
        return {"status": "skipped", "reason": "no annotation-method label columns present",
                "n_methods_available": 0}

    n = adata.n_obs
    votes = adata.obs[present].astype("object").to_numpy()          # n x n_present, may contain None/NaN
    have_winner = bool(winner_col) and winner_col in adata.obs
    winner_in = adata.obs[winner_col].astype(str).to_numpy() if have_winner else None

    score = np.zeros(n, dtype=int)
    n_meth = np.zeros(n, dtype=int)
    agree = np.zeros(n, dtype=float)
    winners = np.empty(n, dtype=object)

    for i in range(n):
        vals = [str(v) for v in votes[i] if not _is_missing(v)]
        n_meth[i] = len(vals)
        if not vals:
            winners[i] = winner_in[i] if have_winner else None
            continue
        if have_winner:
            win = winner_in[i]
        else:
            counts: dict[str, int] = {}
            for v in vals:
                counts[v] = counts.get(v, 0) + 1
            top = max(counts.values())
            win = sorted(k for k, c in counts.items() if c == top)[0]   # deterministic tie-break
        winners[i] = win
        score[i] = sum(1 for v in vals if v == win)
        agree[i] = score[i] / len(vals)

    reliab = np.array([reliability_bin(float(agree[i]), int(n_meth[i])) if n_meth[i] else "single"
                       for i in range(n)], dtype=object)

    adata.obs["consensus_score"] = score
    adata.obs["consensus_n_methods"] = n_meth
    adata.obs["consensus_agreement"] = agree
    adata.obs["consensus_reliability"] = pd.Categorical(reliab)

    scored = n_meth > 0
    bins = pd.Series(reliab[scored]).value_counts().to_dict() if scored.any() else {}
    return {
        "status": "ok",
        "n_methods_available": len(present),
        "methods": present,
        "mean_agreement": float(agree[scored].mean()) if scored.any() else 0.0,
        "reliability_counts": {str(k): int(v) for k, v in bins.items()},
        "pct_trusted": float((n_meth >= MIN_TRUST_METHODS).mean()),
    }
