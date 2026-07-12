"""Regression: ensure() must not skip a present-but-INVALID product.

The demo cache bakes ``uns['panel_check']`` into an ``.h5ad``. AnnData's writer splits the ``/``
in cell-type keys (e.g. ``Malignant/Melanocyte``) into nested h5 groups, so on reload every
coverage entry loses its ``n_present`` / ``status`` fields. ``ensure()`` saw the key present and
skipped recompute, so the corrupt result reached the wizard's Panel-check step and raised
``KeyError: 'n_present'``. ``ensure()`` now honours an optional capability validity predicate and
recomputes an invalid product.
"""

from __future__ import annotations


def _corrupt_panel_check() -> dict:
    """What an h5ad round-trip makes of ``{'Malignant/Melanocyte': {...}}``: a nested group whose
    values are missing ``n_present`` / ``status``."""
    return {"coverage": {"Malignant": {"Melanocyte": {"n_markers": 5}}},
            "confusable_pairs": [], "merge_groups": []}


def test_ensure_recomputes_invalid_panel_check(raw_adata, ctx):
    from spatialscribe.analysis import capabilities as cap
    from spatialscribe.analysis import panel_check

    # A fresh, valid panel_check.
    assert cap.run(raw_adata, "panel_check", {}, ctx).ok
    assert panel_check.is_valid(raw_adata.uns["panel_check"])

    # Simulate the corrupt result an h5ad round-trip bakes into a cache.
    raw_adata.uns["panel_check"] = _corrupt_panel_check()
    assert not panel_check.is_valid(raw_adata.uns["panel_check"])

    # ensure() must NOT treat the present-but-corrupt result as done - it recomputes.
    cap.ensure(raw_adata, "panel_check", ctx)
    assert panel_check.is_valid(raw_adata.uns["panel_check"]), \
        "ensure() skipped an invalid panel_check - the KeyError:'n_present' bug"
    # Every coverage entry has the field the UI reads (the one that KeyError'd).
    assert all("n_present" in d for d in raw_adata.uns["panel_check"]["coverage"].values())


def test_ensure_still_skips_a_valid_product(raw_adata, ctx):
    """The validity check must not defeat normal caching: a valid product is still skipped."""
    from spatialscribe.analysis import capabilities as cap

    assert cap.run(raw_adata, "panel_check", {}, ctx).ok
    res = cap.ensure(raw_adata, "panel_check", ctx)   # valid + present -> cached skip
    assert res.record is None


def test_array_confusable_pairs_is_invalid_and_recomputed(raw_adata, ctx):
    """An h5ad round-trip turns confusable_pairs (a list) into a numpy array, which crashed the UI
    at ``if pc['confusable_pairs']:``. is_valid must catch it and ensure must recompute a list."""
    import numpy as np

    from spatialscribe.analysis import capabilities as cap
    from spatialscribe.analysis import panel_check

    assert cap.run(raw_adata, "panel_check", {}, ctx).ok
    pc = raw_adata.uns["panel_check"]
    pc["confusable_pairs"] = np.array([{"pair": ["A", "B"]}, {"pair": ["C", "D"]}], dtype=object)
    raw_adata.uns["panel_check"] = pc
    assert not panel_check.is_valid(pc)                     # array -> caught as invalid

    cap.ensure(raw_adata, "panel_check", ctx)              # recomputes
    assert isinstance(raw_adata.uns["panel_check"]["confusable_pairs"], list)
    # a plain truthiness check on the recomputed list no longer raises
    assert len(raw_adata.uns["panel_check"]["confusable_pairs"]) >= 0
