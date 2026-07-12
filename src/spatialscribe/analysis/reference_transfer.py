"""Reference-anchored annotation via optimal-transport label transfer (TACCO).

Anchors annotation to a user reference atlas (e.g. a user reference atlas) without a
GPU. TACCO's ``tc.tl.annotate`` returns a per-cell **compositional** assignment over
reference types; the dominant fraction is the label, its value a confidence, and the
entropy an ambiguity flag - which feeds the same abstention machinery as the marker path.

Dependency: ``tacco``, declared as the ``transfer`` extra in ``pyproject.toml`` (pixi pulls it
via ``extras = ["transfer"]``). Use when a reference ``.h5ad`` is supplied; the marker + Claude
path (:mod:`annotate`) is the no-reference default.

TACCO needs RAW INTEGER COUNTS. ``cluster.preprocess`` normalizes and log1p's ``adata.X`` in
place (stashing the originals in ``adata.layers['counts']``), and reference transfer runs AFTER
clustering, so ``adata.X`` is always log-transformed by the time we get here. Passing it to
``tc.tl.annotate`` raises "Some of the counts dont look like integers", ``capabilities.py``
catches that and records ``status="skipped: ..."``, and the celltypist arm silently wins. That
is why the TACCO arm never actually ran. :func:`_as_counts` restores the counts view for the
call. Never "fix" this with ``assume_valid_counts=True``: that runs optimal transport on
log-transformed data and returns confident nonsense.
"""

from __future__ import annotations

from contextlib import contextmanager


def _looks_like_counts(X) -> bool:
    """True when X holds non-negative integers (a cheap sample, not a full scan)."""
    import numpy as np

    data = X.data if hasattr(X, "data") else np.asarray(X).ravel()
    if data.size == 0:
        return True
    s = data[:10_000]
    return bool((s >= 0).all() and np.allclose(s, np.rint(s)))


@contextmanager
def _as_counts(a, what: str):
    """Temporarily expose raw integer counts as ``a.X`` for the duration of the block.

    Prefers ``a.layers['counts']`` (written by :func:`cluster.preprocess`); falls back to ``a.X``
    when that already looks like counts. Raises with an actionable message otherwise rather than
    letting TACCO fail deep inside its own preprocessing. Always restores the original ``.X``.
    """
    original = a.X
    if "counts" in getattr(a, "layers", {}):
        a.X = a.layers["counts"]
    elif not _looks_like_counts(original):
        raise RuntimeError(
            f"TACCO needs raw integer counts, but {what}.X is not integral and has no "
            f"'counts' layer. Keep the raw matrix in {what}.layers['counts']."
        )
    try:
        yield a
    finally:
        a.X = original


def composition_prior(adata, reference, ref_label_key, label_col="celltypist_label"):
    """A CELL-SPACE composition prior for TACCO, estimated from an already-computed label column
    (default the in-env CellTypist prediction). Returns a pandas Series over the reference's type
    names summing to 1, or ``None`` when the column is absent.

    Why this matters: ``tc.tl.annotate(annotation_prior=None)`` derives its prior by count-weighting
    the REFERENCE's own composition. When the reference is depth-capped / subsampled (a flat prior),
    that prior is wildly wrong for a real section - on the Atera breast benchmark it drove TACCO to
    0.35 accuracy vs CellTypist's 0.69, purely because it under-called the 54%-tumour section as
    ~11% tumour. Feeding the section's OWN estimated composition (truth-free, from CellTypist which
    the app runs anyway) lifts TACCO back to 0.69/0.82 (WTA/5K) - measured, validation_2026-07-10.
    Empirically the prior is CELL-space (a cell-fraction beat a read-weighted one on both panels)."""
    import numpy as np
    import pandas as pd

    if label_col not in adata.obs:
        return None
    ref_types = pd.Index(pd.unique(reference.obs[ref_label_key].astype(str)))
    vc = adata.obs[label_col].astype(str).value_counts()
    prior = vc.reindex(ref_types).fillna(0.0)
    total = float(prior.sum())
    if total <= 0:
        return None
    # a small floor so no reference type is hard-zeroed out of the transport, then renormalize
    prior = prior + 1.0
    return prior / float(prior.sum()) if np.isfinite(prior.sum()) else None


def transfer_labels(adata, reference, ref_label_key: str,
                    result_key: str = "ref_transfer", annotation_prior=None) -> dict:
    """OT label transfer from ``reference`` (AnnData) onto ``adata``.

    Writes ``adata.obs['ref_label']`` (argmax), ``adata.obs['ref_confidence']`` (dominant
    fraction) and ``adata.obs['ref_ambiguous']`` (bool, entropy-based), plus the full
    composition in ``adata.obsm[result_key]``. Returns a small summary.

    ``annotation_prior`` (a pandas Series over reference types, cell-space) overrides TACCO's default
    reference-derived prior. Pass :func:`composition_prior` for the deployable, truth-free estimate;
    ``None`` reproduces the old (often bad on a capped reference) behaviour. See that function's note.
    """
    try:
        import tacco as tc
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError(
            "reference transfer needs the `tacco` package (install the 'transfer' extra)."
        ) from exc
    import numpy as np

    with _as_counts(adata, "adata"), _as_counts(reference, "reference"):
        tc.tl.annotate(adata, reference, annotation_key=ref_label_key, result_key=result_key,
                       annotation_prior=annotation_prior)
    comp = np.asarray(adata.obsm[result_key])
    types = list(adata.obsm[result_key].columns) if hasattr(adata.obsm[result_key], "columns") \
        else [f"t{i}" for i in range(comp.shape[1])]
    top = comp.argmax(1)
    conf = comp.max(1)
    # Normalized entropy over the composition -> ambiguity.
    p = np.clip(comp, 1e-9, 1)
    ent = -(p * np.log(p)).sum(1) / np.log(comp.shape[1])
    adata.obs["ref_label"] = [types[i] for i in top]
    adata.obs["ref_confidence"] = conf
    adata.obs["ref_ambiguous"] = ent > 0.7
    return {
        "n_types": len(types),
        "mean_confidence": float(conf.mean()),
        "pct_ambiguous": float((ent > 0.7).mean()),
    }
