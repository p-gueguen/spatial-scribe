"""Post-hoc calibration of the per-cell annotation confidence - fit against REAL labels, or not at all.

Why this module exists
----------------------
``annotate.apply_confidence`` produces a heuristic score: a marker-margin posterior times a contamination
penalty times a panel-coverage ceiling, times an average of soft factors. A benchmark against independent
ground truth (49 fine types coarsened to a shared lineage axis, 12k-cell CosMx breast section) measured two
distinct defects:

* **Mis-calibration** - ECE 0.33 (0.45 on the 100k demo). The number does not mean ``P(correct)``.
* **Non-monotonicity** - lineage accuracy peaks (~0.84) near confidence 0.44 and *falls* to 0.26-0.46
  through the 0.5-0.75 band. Keeping only ``confidence >= 0.5`` scored 0.46 versus 0.69 for keeping every
  cell, so the abstention gate was effectively inverted.

Only the first defect is repairable after the fact, and the distinction is the whole point:

* **Calibration** asks *does the number mean P(correct)?* A monotone isotonic fit on ``(confidence,
  correct)`` from a labeled set repairs this - it relabels each score with the empirical accuracy of the
  cells that scored it.
* **Discrimination** asks *does the number RANK correct cells above incorrect ones?* No monotone map can
  change a ranking, so calibration cannot fix this, and abstention depends entirely on it. :func:`auc`
  measures it so the report can state the limitation rather than imply it away.

This is deliberately fit-only. With no labels there is nothing to fit and nothing is written: the standing
project stance is that the confidence is an ordinal heuristic and that calibration is unsolved. We claim
calibration exactly when real labels were supplied, and never otherwise.

Isotonic (not Platt/temperature) because the input is non-monotonic: a sigmoid fit assumes the score is
already correctly ordered, which is the assumption the benchmark falsified. Where the score carries no
usable ordering, an increasing isotonic fit correctly collapses toward the base rate - the honest answer.

Depends on: numpy, scikit-learn (isotonic), scipy (rank statistics) - all already transitive deps.
Heavy imports stay inside functions.
"""
from __future__ import annotations

RAW_KEY = "annotation_confidence"
CALIBRATED_KEY = "annotation_confidence_calibrated"


def _arrays(conf, correct):
    import numpy as np

    c = np.asarray(conf, dtype=float).ravel()
    y = np.asarray(correct).ravel().astype(float)
    return c, y


# --------------------------------------------------------------------------- #
# metrics
# --------------------------------------------------------------------------- #
def ece(conf, correct, n_bins: int = 10) -> float:
    """Expected Calibration Error: count-weighted mean ``|mean(conf) - mean(correct)|`` over equal-width
    bins. 0.0 = the score means exactly what it says. Empty bins contribute nothing."""
    import numpy as np

    c, y = _arrays(conf, correct)
    if c.size == 0:
        return 0.0
    n_bins = int(n_bins)
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    idx = np.clip(np.digitize(c, edges[1:-1], right=False), 0, n_bins - 1)
    total = 0.0
    for b in range(n_bins):
        m = idx == b
        if m.any():
            total += m.sum() * abs(c[m].mean() - y[m].mean())
    return float(total / c.size)


def brier(conf, correct) -> float:
    """Brier score ``mean((conf - correct)^2)`` - calibration and sharpness together. Lower is better."""
    import numpy as np

    c, y = _arrays(conf, correct)
    return float(np.mean((c - y) ** 2)) if c.size else 0.0


def auc(conf, correct) -> float:
    """Probability a random CORRECT cell scores above a random INCORRECT one (rank AUC, ties = 0.5).

    This is the discrimination check, and it is the number that decides whether abstention can ever work:

    * ``auc > 0.5`` - confidence ranks correct cells higher; dropping low-confidence cells raises accuracy.
    * ``auc ~ 0.5`` - the score is noise with respect to correctness; NO threshold helps, at any calibration.
    * ``auc < 0.5`` - the gate is inverted; abstaining on low confidence discards the cells most likely right.

    Isotonic calibration can never REORDER two cells, but it is non-*strictly* increasing: its flat regions
    merge distinct raw scores into ties. Where the raw score was locally anti-correlated with correctness,
    erasing that ordering moves those pairs from a losing contribution to a neutral one, so ``auc`` after
    calibration can be strictly HIGHER than before. Report both rather than assuming invariance.

    Returns 0.5 for a degenerate single-class input.
    """
    import numpy as np

    c, y = _arrays(conf, correct)
    pos = y > 0.5
    n1, n0 = int(pos.sum()), int((~pos).sum())
    if n1 == 0 or n0 == 0:
        return 0.5
    from scipy.stats import rankdata

    r = rankdata(c)                       # average ranks -> ties contribute 0.5 as they should
    return float((r[pos].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def auc_within(conf, correct, groups, min_n: int = 20) -> dict:
    """Discrimination measured INSIDE each group, free of across-group confounding.

    Pooled :func:`auc` mixes two very different things: whether confidence ranks correct cells above
    incorrect ones *within* a cell type, and whether the cell types the pipeline is confident about happen
    to be the ones it labels well. The second term produces Simpson reversals. Measured on real sections:

    ==================  ==========  =============  =============
    section             pooled AUC  within-lineage across-lineage
    ==================  ==========  =============  =============
    CosMx breast 1k         0.4985         0.4910         0.5746
    Atera breast 18k        0.4191         0.4790         0.1340
    Atera cervical 18k      0.5780         0.5710         0.4270
    ==================  ==========  =============  =============

    Atera breast inverts *only* because the pipeline is most confident on Endothelial (52.7% accurate) and
    least confident on Epithelial (98.6% accurate, 58% of cells). The within-lineage signal is flat.
    A per-cell gate can only exploit the within-group term, so that is the number to report.

    Confidence is rank-normalised inside each group (removing level offsets), then scored pooled. Groups
    with fewer than ``min_n`` cells, or with only one outcome class, are dropped from ``per_group`` but
    still contribute their rank-normalised cells to ``auc_within``.
    """
    import numpy as np

    c, y = _arrays(conf, correct)
    g = np.asarray(groups).astype(str)
    if c.size == 0:
        return {"auc_within": 0.5, "n_groups": 0, "per_group": {}}

    ranked = np.full(c.size, 0.5)
    per: dict[str, dict] = {}
    for lab in np.unique(g):
        m = g == lab
        n = int(m.sum())
        if n > 1:
            # rank -> [0, 1] inside the group; ties get distinct ranks, which is fine for a pooled AUC
            ranked[m] = c[m].argsort().argsort() / (n - 1)
        if n >= min_n and 0 < int(y[m].sum()) < n:
            per[lab] = {"n": n, "accuracy": float(y[m].mean()), "mean_conf": float(c[m].mean()),
                        "auc": auc(c[m], y[m])}
    return {"auc_within": auc(ranked, y), "n_groups": len(per), "per_group": per}


def reliability_curve(conf, correct, n_bins: int = 10) -> dict:
    """Per-bin mean confidence vs empirical accuracy (the reliability diagram), JSON-serializable.

    Empty bins are dropped, so all four lists stay aligned and plottable as-is.
    """
    import numpy as np

    c, y = _arrays(conf, correct)
    n_bins = int(n_bins)
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    out: dict[str, list] = {"bin_mid": [], "bin_conf": [], "bin_acc": [], "bin_count": []}
    if c.size == 0:
        return out
    idx = np.clip(np.digitize(c, edges[1:-1], right=False), 0, n_bins - 1)
    for b in range(n_bins):
        m = idx == b
        if not m.any():
            continue
        out["bin_mid"].append(float((edges[b] + edges[b + 1]) / 2))
        out["bin_conf"].append(float(c[m].mean()))
        out["bin_acc"].append(float(y[m].mean()))
        out["bin_count"].append(int(m.sum()))
    return out


# --------------------------------------------------------------------------- #
# fit / apply / report
# --------------------------------------------------------------------------- #
def fit_isotonic(conf, correct):
    """Fit a monotone-increasing isotonic map ``confidence -> P(correct)`` on a LABELED set.

    Returns ``None`` (never raises) when there is nothing to fit - fewer than two distinct confidence
    values, or scikit-learn unavailable - so the caller skips honestly instead of faking a calibrator.
    """
    import numpy as np

    c, y = _arrays(conf, correct)
    if c.size < 2 or np.unique(c).size < 2:
        return None
    try:
        from sklearn.isotonic import IsotonicRegression
    except Exception:
        return None
    return IsotonicRegression(increasing=True, out_of_bounds="clip", y_min=0.0, y_max=1.0).fit(c, y)


def apply(adata, calibrator, src: str = RAW_KEY, dst: str = CALIBRATED_KEY) -> int:
    """Write ``obs[dst] = calibrator(obs[src])``, clipped to [0, 1]. Returns the number of cells written.

    Refuses ``src == dst``: the raw heuristic is the audit trail and is never overwritten. A ``None``
    calibrator is a no-op returning 0 (nothing was fit, so nothing is claimed).
    """
    import numpy as np

    if src == dst:
        raise ValueError(f"refusing to overwrite the raw heuristic {src!r} with a calibrated score; "
                         f"write to a separate column")
    if calibrator is None or src not in adata.obs:
        return 0
    vals = np.asarray(adata.obs[src], dtype=float)
    adata.obs[dst] = np.clip(calibrator.predict(vals), 0.0, 1.0)
    return int(vals.size)


def abstention_curve(conf, correct, thresholds=None, min_keep: float = 0.10) -> list:
    """Accuracy and kept-fraction for the gate "keep cells with ``conf >= t``", per threshold.

    This is the number that decides whether an abstention gate earns its place: ``gain`` is the accuracy
    on the kept cells minus the accuracy on all cells. A gate that only wins by discarding almost every
    cell has not bought anything, so entries keeping less than ``min_keep`` of the section are marked
    ``useful=False`` rather than silently reported as a win.
    """
    import numpy as np

    c, y = _arrays(conf, correct)
    if c.size == 0:
        return []
    base = float(y.mean())
    ts = np.arange(0.05, 1.0, 0.05) if thresholds is None else np.asarray(thresholds, dtype=float)
    out = []
    for t in ts:
        keep = c >= t
        if not keep.any():
            continue
        kf = float(keep.mean())
        acc = float(y[keep].mean())
        out.append({"threshold": round(float(t), 3), "keep_frac": kf, "accuracy": acc,
                    "gain": acc - base, "useful": bool(kf >= min_keep and acc > base)})
    return out


def report(conf, correct, calibrator, n_bins: int = 10, baserate: float | None = None,
           groups=None) -> dict:
    """Before/after calibration metrics, the AUC discrimination check, and the abstention curve.

    Evaluate this on a HELD-OUT split (the calibrator must not be fit on the cells it is scored against),
    otherwise ``ece_after`` is optimistic by construction. All values are JSON-serializable.

    ``auc_before`` vs ``auc_after`` differ only through the ties isotonic introduces - see :func:`auc`.

    **A low ``ece_after`` proves almost nothing on its own.** The constant predictor "every cell is correct
    with probability = the section's accuracy" is already near-perfectly calibrated (``ece_baserate`` ~ 0).
    Any calibrator that collapses toward that constant inherits its near-zero ECE without having learned
    anything per-cell. ``brier_skill`` = ``1 - brier_after / brier_baserate`` is the honest headline: it is
    the fraction of squared error the calibrated score removes *beyond* just knowing the base rate, and it
    is ~0 exactly when the confidence carries no usable information about which cells are right.

    ``baserate`` should be the accuracy of the FIT split (the null model may not peek at the eval labels).
    Defaults to the eval accuracy, which flatters the null model slightly and is disclosed here.
    """
    import numpy as np

    c, y = _arrays(conf, correct)
    after = np.clip(calibrator.predict(c), 0.0, 1.0) if calibrator is not None else c
    br = float(y.mean()) if baserate is None else float(baserate)
    const = np.full(c.size, br)
    brier_base = brier(const, y)
    within = auc_within(c, y, groups) if groups is not None else None
    return {
        "n": int(c.size),
        # pooled AUC is Simpson-confounded across cell types; this is the term a per-cell gate can use
        "auc_within": (within or {}).get("auc_within"),
        "per_group": (within or {}).get("per_group"),
        "ece_before": ece(c, y, n_bins), "ece_after": ece(after, y, n_bins),
        "brier_before": brier(c, y), "brier_after": brier(after, y),
        "auc_before": auc(c, y), "auc_after": auc(after, y),      # the abstention diagnostic
        "accuracy": float(y.mean()) if y.size else 0.0,
        # the null model: predict the base rate for every cell. Calibrated by construction.
        "baserate": br,
        "ece_baserate": ece(const, y, n_bins), "brier_baserate": brier_base,
        "brier_skill": float(1.0 - brier(after, y) / brier_base) if brier_base > 0 else 0.0,
        "reliability_before": reliability_curve(c, y, n_bins),
        "reliability_after": reliability_curve(after, y, n_bins),
        "abstention_after": abstention_curve(after, y),
    }
