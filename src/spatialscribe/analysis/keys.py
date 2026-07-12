"""The AnnData data contract - every obs/uns/obsm/varm/layers key the pipeline uses.

What it does
------------
Names, in ONE place, every ``adata`` slot the analysis engine reads or writes. Each
:class:`Key` carries its container (``obs`` / ``uns`` / ``obsm`` / ``varm`` / ``layers`` /
``var`` / ``obsp``) and its name, and knows how to test its own presence on an object.
This turns the inter-stage contract - previously a scatter of string literals across the
UI, the copilot, and ``analysis/`` - into greppable symbols with a single definition.

How to use it
-------------
>>> from .keys import Obs
>>> Obs.CELL_TYPE.present(adata)          # True iff adata.obs has 'cell_type'
>>> str(Obs.CELL_TYPE)                     # 'obs:cell_type'
>>> Obs.CELL_TYPE.name                     # 'cell_type'  (the bare adata key)

Depends on
----------
Nothing (pure). ``present`` accepts any object exposing ``.obs``/``.uns``/... mappings,
so it works on a real AnnData without importing anndata here.

Notes
-----
* Keys whose exact name is computed at runtime (per-signature score columns) are exposed
  as *prefixes* (:data:`SCORE_PREFIX`, :data:`STATE_PREFIX`), not individual keys.
* ``all_keys()`` enumerates every declared Key; the registry-conformance test asserts that
  every capability's ``requires``/``produces`` references only keys declared here.
"""

from __future__ import annotations

from dataclasses import dataclass

# Valid AnnData containers a Key can live in.
_SPACES = ("obs", "var", "uns", "obsm", "varm", "layers", "obsp")


@dataclass(frozen=True)
class Key:
    """One named slot on an AnnData, e.g. ``Key('obs', 'cell_type')``.

    ``space`` is the AnnData attribute ('obs', 'uns', 'obsm', ...); ``name`` is the key
    within it. ``present(adata)`` is the single, container-aware existence check used by
    the capability prerequisite logic.
    """

    space: str
    name: str

    def __post_init__(self) -> None:
        if self.space not in _SPACES:
            raise ValueError(f"Key space must be one of {_SPACES}, got {self.space!r}")

    def present(self, adata) -> bool:
        """True iff this slot exists on ``adata`` (columns for obs/var, keys otherwise)."""
        container = getattr(adata, self.space, None)
        if container is None:
            return False
        try:
            return self.name in container
        except TypeError:  # pragma: no cover - defensive for exotic containers
            return False

    def __str__(self) -> str:
        return f"{self.space}:{self.name}"


# --------------------------------------------------------------------------- #
# obs (per-cell) columns
# --------------------------------------------------------------------------- #
class Obs:
    # QC (qc.compute_qc via scanpy.calculate_qc_metrics)
    TOTAL_COUNTS = Key("obs", "total_counts")
    N_GENES = Key("obs", "n_genes_by_counts")
    PCT_CONTROL = Key("obs", "pct_counts_control")
    # Segmentation (qc.segmentation_qc)
    NUCLEUS_PRESENT = Key("obs", "nucleus_present")
    NUCLEUS_TO_CELL_RATIO = Key("obs", "nucleus_to_cell_ratio")
    SEG_AREA_FLAG = Key("obs", "seg_area_flag")
    # ovrlpy VSI (qc.apply_ovrlpy_vsi)
    VSI = Key("obs", "vsi")
    VSI_CONFIDENT = Key("obs", "vsi_confident")
    # Clustering (cluster.cluster)
    LEIDEN = Key("obs", "leiden")
    # Annotation (annotate.consensus_annotate / apply_confidence)
    CELL_TYPE = Key("obs", "cell_type")
    CELL_TYPE_FINAL = Key("obs", "cell_type_final")
    ANNOTATION_CONFIDENCE = Key("obs", "annotation_confidence")
    # Post-hoc isotonic calibration of the above, fit against real labels (calibration.py). Written only
    # by the calibrate_confidence capability; the raw score is never overwritten.
    ANNOTATION_CONFIDENCE_CALIBRATED = Key("obs", "annotation_confidence_calibrated")
    ANNOTATION_VERDICT = Key("obs", "annotation_verdict")
    ANNOTATION_REASON = Key("obs", "annotation_reason")
    ANNOTATION_STABILITY = Key("obs", "annotation_stability")
    # Purity (purity.pmp / crisp_purity)
    PMP = Key("obs", "pmp")
    CRISP_IMPURE = Key("obs", "crisp_impure")
    # Spatial (spatial.spatial_coherence)
    SPATIAL_COHERENCE = Key("obs", "spatial_coherence")
    # Niches (niches.call_niches)
    NICHE = Key("obs", "niche")
    # Cell-state typing (states.assign_cell_states) - the dominant program per cell
    CELL_STATE = Key("obs", "cell_state")
    # Programs (programs.discover_programs)
    PROGRAM = Key("obs", "program")
    # Tumour (cnv.malignant_score / call_malignant_cnv)
    MALIGNANT_SCORE = Key("obs", "malignant_score")
    CNV_SCORE = Key("obs", "cnv_score")
    IS_MALIGNANT = Key("obs", "is_malignant")
    # Subclustering (subcluster.subcluster)
    SUBTYPE = Key("obs", "subtype")
    # Rejection reasons (rejection.assign_rejection_reasons)
    REJECTION_REASON = Key("obs", "rejection_reason")
    REJECTION_DETAIL = Key("obs", "rejection_detail")
    # Xenium segmentation extras carried by io.load
    TRANSCRIPT_COUNTS = Key("obs", "transcript_counts")
    CELL_AREA = Key("obs", "cell_area")
    NUCLEUS_AREA = Key("obs", "nucleus_area")


# --------------------------------------------------------------------------- #
# var (per-gene) columns
# --------------------------------------------------------------------------- #
class Var:
    CONTROL = Key("var", "control")
    FEATURE_TYPES = Key("var", "feature_types")


# --------------------------------------------------------------------------- #
# uns (unstructured) entries
# --------------------------------------------------------------------------- #
class Uns:
    PANEL_CHECK = Key("uns", "panel_check")
    LOG1P = Key("uns", "log1p")
    NEIGHBORS = Key("uns", "neighbors")
    RANK_GENES = Key("uns", "rank_genes_groups")
    MORAN_I = Key("uns", "moranI")
    STATE_COLUMNS = Key("uns", "state_columns")
    DEMO_META = Key("uns", "spatialscribe_demo")
    REFERENCE_MATCH = Key("uns", "reference_match")
    ANNOTATION_VERIFICATION = Key("uns", "annotation_verification")
    ANNOTATION_ROUTE = Key("uns", "annotation_route")
    CALIBRATION_REPORT = Key("uns", "calibration_report")


# --------------------------------------------------------------------------- #
# obsm / varm / layers / obsp
# --------------------------------------------------------------------------- #
class Obsm:
    SPATIAL = Key("obsm", "spatial")
    X_PCA = Key("obsm", "X_pca")
    X_UMAP = Key("obsm", "X_umap")
    PROGRAMS = Key("obsm", "programs")
    PROGRAM_SCORES = Key("obsm", "program_scores")


class Varm:
    PROGRAM_LOADINGS = Key("varm", "program_loadings")


class Layers:
    COUNTS = Key("layers", "counts")


class Obsp:
    SPATIAL_CONNECTIVITIES = Key("obsp", "spatial_connectivities")
    CONNECTIVITIES = Key("obsp", "connectivities")


# --------------------------------------------------------------------------- #
# Runtime-computed key families (per-signature score columns): prefixes, not keys.
# --------------------------------------------------------------------------- #
SCORE_PREFIX = "score_"      # annotate.score_marker_sets -> obs['score_<lineage>']
STATE_PREFIX = "state_"      # states.score_states       -> obs['state_<program>']


def nhood_enrichment_key(cluster_key: str) -> Key:
    """The uns slot squidpy writes for neighborhood enrichment of ``cluster_key``."""
    return Key("uns", f"{cluster_key}_nhood_enrichment")


def co_occurrence_key(cluster_key: str) -> Key:
    """The uns slot squidpy writes for co-occurrence of ``cluster_key``."""
    return Key("uns", f"{cluster_key}_co_occurrence")


def all_keys() -> list[Key]:
    """Every explicitly-declared :class:`Key` (used by the conformance test).

    Excludes the runtime-computed families (score_/state_/<key>_nhood_enrichment), which
    are validated by prefix, not by identity.
    """
    out: list[Key] = []
    for cls in (Obs, Var, Uns, Obsm, Varm, Layers, Obsp):
        for attr in vars(cls).values():
            if isinstance(attr, Key):
                out.append(attr)
    return out


def is_declared(key: Key) -> bool:
    """True iff ``key`` is one of the explicitly-declared keys (identity by value)."""
    return key in set(all_keys())
