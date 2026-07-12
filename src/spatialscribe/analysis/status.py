"""Derived pipeline state - what has been computed, and what can run next.

What it does
------------
Answers "which stage is this AnnData at?" by inspecting which capability ``produces`` keys
are actually present, instead of trusting a hand-maintained ``done`` set that can lie (e.g.
marking every step done on cache-load when nothing was computed). Both frontends read the
same derived truth: the wizard for its checkmarks and greying-out, the copilot for "run X
first" guidance.

How to use it
-------------
>>> from .status import pipeline_status, is_done, available, missing_prereqs
>>> pipeline_status(adata)        # {'cluster': True, 'annotate': False, ...}
>>> is_done(adata, "annotate")    # False
>>> available(adata)              # capabilities whose prerequisites are all met
>>> missing_prereqs(adata, "annotate")   # [Key('obs','leiden')] if clustering hasn't run

Depends on
----------
:mod:`spatialscribe.analysis.capabilities` (the registry) and :mod:`keys`.
"""

from __future__ import annotations

from . import capabilities as _cap
from .keys import Key


def is_done(adata, name: str) -> bool:
    """True iff every key ``name`` declares as ``produces`` is present on ``adata``.

    A capability with no declared products (a pure report, e.g. ``describe_sample``) is
    never "done" in this sense - it is re-runnable on demand - so this returns False for it.
    """
    cap = _cap.REGISTRY.get(name)
    if cap is None or not cap.produces:
        return False
    return all(k.present(adata) for k in cap.produces)


def pipeline_status(adata) -> dict[str, bool]:
    """``{capability_name: is_done}`` for every registered capability."""
    return {name: is_done(adata, name) for name in _cap.REGISTRY}


def missing_prereqs(adata, name: str) -> list[Key]:
    """The prerequisite keys ``name`` needs that are not yet present (empty if runnable)."""
    cap = _cap.REGISTRY.get(name)
    if cap is None:
        return []
    return _cap.missing_prereqs(adata, cap)


def can_run(adata, name: str) -> bool:
    """True iff all of ``name``'s prerequisites are satisfied on ``adata``."""
    return not missing_prereqs(adata, name)


def available(adata) -> list[str]:
    """Capability names whose prerequisites are all met on ``adata`` (runnable now)."""
    return [name for name in _cap.REGISTRY if can_run(adata, name)]
