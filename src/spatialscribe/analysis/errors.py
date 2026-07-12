"""Structured errors for the analysis engine - recoverable vs genuine failure.

What it does
------------
Gives the two frontends a *typed* error contract instead of raw exception strings. A
:class:`PrerequisiteError` means "run an earlier step first" (recoverable - the wizard
greys the step out, the copilot is told which capability to call). Any other exception is
a genuine failure. ``to_error_dict`` renders either into a small JSON-able dict that the
copilot can reason over and the UI can render.

How to use it
-------------
>>> raise PrerequisiteError("annotate", missing=[Obs.CELL_TYPE], hint="run 'cluster' then 'annotate'")
>>> to_error_dict(exc)   # -> {'error_type': 'prerequisite_missing', 'capability': ..., ...}

Depends on
----------
:mod:`spatialscribe.analysis.keys` (for the :class:`Key` type only).
"""

from __future__ import annotations

from .keys import Key


class SpatialScribeError(Exception):
    """Base class for engine errors surfaced through the capability layer."""


class PrerequisiteError(SpatialScribeError):
    """A capability was invoked before an upstream step produced what it needs.

    Recoverable: the caller should run the producing capability first. ``missing`` lists
    the absent :class:`Key`s; ``hint`` is a plain-language remediation string.
    """

    error_type = "prerequisite_missing"

    def __init__(self, capability: str, missing: list[Key], hint: str = "") -> None:
        self.capability = capability
        self.missing = list(missing)
        self.hint = hint
        keys = ", ".join(str(k) for k in self.missing) or "(unknown)"
        super().__init__(
            f"'{capability}' needs {keys} which is not present yet."
            + (f" {hint}" if hint else "")
        )

    def to_dict(self) -> dict:
        return {
            "error_type": self.error_type,
            "capability": self.capability,
            "missing": [str(k) for k in self.missing],
            "hint": self.hint,
        }


class CapabilityError(SpatialScribeError):
    """A capability ran but failed (bad params, empty result, upstream exception)."""

    error_type = "capability_failed"

    def __init__(self, capability: str, message: str) -> None:
        self.capability = capability
        self.message = message
        super().__init__(f"'{capability}' failed: {message}")

    def to_dict(self) -> dict:
        return {"error_type": self.error_type, "capability": self.capability,
                "message": self.message}


def to_error_dict(exc: Exception, capability: str | None = None) -> dict:
    """Render any exception into the JSON-able error contract.

    Typed engine errors keep their structured form; anything else becomes a generic
    ``capability_failed`` so the copilot never receives an opaque traceback string.
    """
    if isinstance(exc, (PrerequisiteError, CapabilityError)):
        return exc.to_dict()
    return {"error_type": "capability_failed",
            "capability": capability or "", "message": str(exc)}
