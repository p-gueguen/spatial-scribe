"""Typed return shapes for the analysis engine (annotation-only, non-breaking).

What it does
------------
Documents, as ``TypedDict``s, the JSON-able dicts the copilot-facing capabilities return,
so both frontends index a named shape instead of guessing string keys independently. These
are pure annotations - they add no runtime behaviour and never wrap or validate the real
returns - but they make the contract greppable and give the conformance suite a shape to
assert against.

Depends on
----------
Nothing at runtime (``typing`` only).
"""

from __future__ import annotations

from typing import TypedDict


class QcSummary(TypedDict, total=False):
    n_cells: int
    median_genes_per_cell: float
    median_transcripts_per_cell: float
    fraction_empty_cells: float
    pct_counts_control: float
    flags: dict[str, str]


class CountFloor(TypedDict, total=False):
    mode: str
    floor: float
    n_panel_genes: int
    n_removed: int
    pct_removed: float


class ImmuneExclusion(TypedDict, total=False):
    tumor: str
    tcell: str
    zscore: float
    verdict: str
    error: str


class NhoodEnrichment(TypedDict, total=False):
    cluster_key: str
    categories: list[str]
    zscore: list[list[float]]
    dropped_sparse_types: list[str]


class Annotatability(TypedDict, total=False):
    pct_pass: float
    pct_warn: float
    pct_abstain: float
    top_abstention_reasons: dict[str, int]
    mean_confidence: float
    usability: str          # output-usability gate: "ok" | "warn" | "fail"


class PanelSummary(TypedDict, total=False):
    resolvable: list[str]
    cannot_resolve: list[str]
    confusable: list[list[str]]


class DescribeSample(TypedDict, total=False):
    n_cells: int
    composition: dict[str, int]
    qc: QcSummary
