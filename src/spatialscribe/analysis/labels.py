"""Tolerant cell-type-label resolution, shared by every capability that takes a label from the LLM.

The copilot echoes the user's wording ("T cells", "tumor") rather than the exact category
("T cell", "Epithelial/Tumor"). A capability that matches the label EXACTLY then errors on the
near-miss, and the copilot loop silently recovers into a different tool - the user gets the wrong
result with no error. Route every LLM/user-supplied label through :func:`match_cell_type` first.
"""
from __future__ import annotations


def norm_label(s: str) -> str:
    """Lowercase, unify hyphens/underscores to spaces, collapse whitespace, drop a trailing 's'."""
    s = " ".join(str(s).lower().replace("-", " ").replace("_", " ").split())
    return s[:-1] if s.endswith("s") else s


def match_cell_type(requested: str, categories) -> str | None:
    """Resolve an LLM/user-supplied label to an actual category, tolerant of case and plural/singular
    (exact -> case-insensitive -> normalized -> unambiguous substring), or None if nothing matches
    unambiguously (the caller should then raise a helpful "Available: [...]" error)."""
    cats = [str(c) for c in categories]
    if requested in cats:
        return requested
    low = {c.lower(): c for c in cats}
    if str(requested).lower() in low:
        return low[str(requested).lower()]
    nreq = norm_label(requested)
    hits = [c for c in cats if norm_label(c) == nreq]
    if len(hits) == 1:
        return hits[0]
    subs = [c for c in cats if nreq and (nreq in norm_label(c) or norm_label(c) in nreq)]
    return subs[0] if len(subs) == 1 else None


def demo() -> None:
    cats = ["T cell", "B cell", "Myeloid", "Epithelial/Tumor"]
    assert match_cell_type("T cells", cats) == "T cell"      # plural
    assert match_cell_type("t cell", cats) == "T cell"       # case
    assert match_cell_type("myeloid cells", cats) == "Myeloid"
    assert match_cell_type("tumor", cats) == "Epithelial/Tumor"  # substring
    assert match_cell_type("neuron", cats) is None           # absent
    assert match_cell_type("cell", cats) is None             # ambiguous
    print("labels demo OK")


if __name__ == "__main__":
    demo()
