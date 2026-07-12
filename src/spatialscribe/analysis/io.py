"""Platform-agnostic ingestion for imaging-based spatial transcriptomics.

What it does
------------
`load(path)` detects the platform (Xenium / CosMx / MERSCOPE / Atera) from the
on-disk directory signature and returns a :class:`SpatialSample` - one common
in-memory contract so every downstream step is platform-agnostic by construction.

How to use it
-------------
>>> sample = load("/path/to/xenium_output")          # auto-detects
>>> sample.adata                                      # cells x genes AnnData, obsm['spatial'] set
>>> sample.panel_genes                                # biological genes only (controls stripped)
>>> sample.control_mask                               # bool over adata.var_names
>>> if sample.has_z: run_ovrlpy(sample.transcripts_path)   # VSI needs the raw molecule table

Depends on
----------
scanpy, anndata, pandas, numpy (+ pyarrow to read parquet). ``spatialdata_io`` is
optional and only used for the non-demo readers; the Xenium demo path deliberately
avoids it (see ``_read_xenium``).

Notes / verified gotchas
------------------------
* There is **no** ``squidpy.read.xenium``. Xenium is read via a direct
  ``cell_feature_matrix.h5`` + ``cells.parquet`` path here.
* ``spatialdata_io.xenium`` returns a table **without** ``obsm['spatial']`` (centroids
  live in the shapes element), so the direct path is preferred for the app.
* z-coordinates live only in the transcript table, never on the cell x gene matrix;
  ``has_z`` reflects the transcript schema, and ovrlpy (VSI) requires it.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Optional

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import anndata

# Non-biological probe name patterns, unioned across platforms. Preferred signal is
# ``var['feature_types']`` when present (Xenium); these regexes are the fallback.
_CONTROL_PATTERNS = re.compile(
    r"^(NegControlProbe|NegControlCodeword|UnassignedCodeword|DeprecatedCodeword|"
    r"Genomic|GenomicControl|IntergenicRegion|antisense|BLANK|Blank|"
    r"Negative|NegPrb|SystemControl|FalseCode)",
    flags=re.IGNORECASE,
)

# feature_types values that are real signal (everything else is treated as control).
_SIGNAL_FEATURE_TYPES = {"Gene Expression", "Protein Expression"}

# TRUE negative controls: features that should capture (near) zero real transcripts, so the share of
# a cell's counts landing on them measures background/contamination.
#
# Deliberately EXCLUDES 'Deprecated Codeword' and 'Unassigned Codeword'. Deprecated codewords are
# decommissioned GENE probes that still decode real transcripts: on a Xenium Prime 5K breast section
# their 3294 features carry a median 14.96% of every cell's counts, while the true negative controls
# carry 0.000%. Counting them as controls made ``pct_counts_control`` ~15%, saturated the
# contamination penalty (``clip(pct/5, 0, 1)`` -> 1) and abstained 94.8% of the section.
_NEG_CONTROL_FEATURE_TYPES = {
    "Negative Control Probe", "Negative Control Codeword", "Genomic Control", "Blank Codeword",
}
_NEG_CONTROL_PATTERNS = re.compile(
    r"^(NegControlProbe|NegControlCodeword|Genomic|GenomicControl|BLANK|Blank|"
    r"Negative|NegPrb|SystemControl|FalseCode)",
    flags=re.IGNORECASE,
)


@dataclass
class SpatialSample:
    """Common contract for one imaging-based spatial section.

    Attributes
    ----------
    platform:      "xenium" | "cosmx" | "merscope" | "atera" | "h5ad"
    adata:         cells x genes AnnData, ``obsm['spatial']`` = 2D centroids (microns).
    control_mask:  bool array over ``adata.var_names``; True = non-biological probe.
    panel_genes:   list of biological gene symbols (controls removed).
    has_z:         True if a per-transcript z-coordinate is available (needed for ovrlpy).
    transcripts_path: path to the raw molecule table (parquet/csv), or None. Kept as a
                   path (not an eager DataFrame) so multi-million-row tables are read
                   lazily by the consumer (e.g. the ovrlpy subprocess) with polars.
    sdata:         optional SpatialData handle when a reader produced one; else None.
    """

    platform: str
    adata: "anndata.AnnData"  # noqa: F821  (imported lazily in readers)
    control_mask: np.ndarray
    panel_genes: list[str]
    has_z: bool
    transcripts_path: Optional[Path] = None
    sdata: object = None
    meta: dict = field(default_factory=dict)

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return (
            f"SpatialSample(platform={self.platform!r}, n_cells={self.adata.n_obs:,}, "
            f"n_panel_genes={len(self.panel_genes):,}, has_z={self.has_z}, "
            f"transcripts={'yes' if self.transcripts_path else 'no'})"
        )


# --------------------------------------------------------------------------- #
# Platform detection
# --------------------------------------------------------------------------- #
def detect_platform(path: Path) -> str:
    """Infer the platform from a directory's file signature.

    Globs tolerantly because export filename suffixes drift across pipeline versions.
    Atera preview data is Xenium-v4-shaped, so it is detected as ``xenium`` and
    re-stamped ``atera`` only when an ~18k-gene panel or Atera metadata is seen.
    """
    path = Path(path)
    if path.is_file() and path.suffix in {".h5ad", ".h5"}:
        return "h5ad" if path.suffix == ".h5ad" else "xenium"

    names = {p.name for p in path.iterdir()} if path.is_dir() else set()

    def has(*needles: str) -> bool:
        return any(any(n in name for name in names) for n in needles)

    if "experiment.xenium" in names or ("cell_feature_matrix.h5" in names and has("cells.parquet")):
        return "xenium"
    if has("_exprMat_file", "exprMat") and has("_metadata_file", "metadata_file"):
        return "cosmx"
    if has("cell_by_gene") and has("cell_metadata"):
        return "merscope"
    # VisiumHD (SpaceRanger) ships binned_outputs/ (square_NNNum) and/or segmented_outputs/.
    if "binned_outputs" in names or "segmented_outputs" in names:
        return "visiumhd"
    # A molecule-table-only export cannot be loaded as cells: it needs segmentation / binning first.
    if names == {"transcripts.parquet"}:
        raise ValueError(
            "transcripts-only export (a molecule table, no cell x gene matrix) - needs cell "
            "segmentation or binning before it is a loadable section; not a direct-load format."
        )
    # Classic (spot) Visium: a filtered matrix alongside a top-level spatial/ dir, but no
    # binned_outputs/segmented_outputs (those are Visium HD, matched above). Unsupported by design -
    # SpatialScribe targets single-cell imaging platforms + Visium HD - so REJECT it clearly here
    # instead of letting it fall through to the non-spatial Flex reader below, which would silently
    # drop the spot coordinates and load it as dissociated scRNA (a confidently-wrong load).
    if has("filtered_feature_bc_matrix.h5") and "spatial" in names:
        raise ValueError(
            "classic (spot) Visium is not supported - SpatialScribe targets single-cell imaging "
            "platforms (Xenium / CosMx / MERSCOPE) and Visium HD. Re-run with Visium HD output "
            "(binned_outputs/ or segmented_outputs/), or load a segmented cell x gene .h5ad."
        )
    # 10x Flex / dissociated scRNA: a single filtered_feature_bc_matrix.h5 (or matrix.h5), no spatial.
    if has("filtered_feature_bc_matrix.h5") or "matrix.h5" in names:
        return "flex"
    raise ValueError(
        f"Could not detect a supported spatial platform in {path}. "
        f"Found: {sorted(names)[:12]}"
    )


# --------------------------------------------------------------------------- #
# Control-probe mask
# --------------------------------------------------------------------------- #
def build_control_mask(adata) -> np.ndarray:
    """Boolean mask over ``adata.var_names``; True = non-biological control probe.

    Prefers the ``feature_types`` column (Xenium onboard analysis) and falls back to
    a name-pattern match for platforms/exports that lack it.
    """
    if "feature_types" in adata.var.columns:
        ft = adata.var["feature_types"].astype(str)
        return (~ft.isin(_SIGNAL_FEATURE_TYPES)).to_numpy()
    return np.array([bool(_CONTROL_PATTERNS.match(str(g))) for g in adata.var_names])


def build_neg_control_mask(adata) -> np.ndarray:
    """Boolean mask over ``adata.var_names``; True = TRUE negative-control feature.

    A strict subset of :func:`build_control_mask`. Use this - never the broad control mask - for
    ``pct_counts_control``: the broad mask also covers deprecated / unassigned codewords, which carry
    real transcript signal (see :data:`_NEG_CONTROL_FEATURE_TYPES`) and make the metric mean
    "non-gene counts" instead of "background counts".
    """
    if "feature_types" in adata.var.columns:
        ft = adata.var["feature_types"].astype(str)
        return ft.isin(_NEG_CONTROL_FEATURE_TYPES).to_numpy()
    return np.array([bool(_NEG_CONTROL_PATTERNS.match(str(g))) for g in adata.var_names])


# --------------------------------------------------------------------------- #
# Readers
# --------------------------------------------------------------------------- #
# Columns the app actually consumes from cells.parquet (id + centroids + segmentation).
_CELL_COLS = ("cell_id", "x_centroid", "y_centroid", "x", "y", "center_x", "center_y",
              "transcript_counts", "cell_area", "nucleus_area")


def _read_cells_table(cells: Path):
    """Read cells.parquet fast via polars column projection (pandas fallback).

    Only the id + centroid + segmentation columns the pipeline uses are materialized, which
    is markedly faster than a full pandas read on 10^5-10^6-cell sections (polars is a
    multithreaded Rust engine and skips the unused columns entirely). Returns a pandas
    DataFrame indexed by the cell-id column; falls back to a full pandas read on any error.
    """
    try:
        import polars as pl

        avail = pl.scan_parquet(cells).collect_schema().names()
        id_col = "cell_id" if "cell_id" in avail else avail[0]
        cols = list(dict.fromkeys([id_col] + [c for c in _CELL_COLS if c in avail]))
        return pl.read_parquet(cells, columns=cols).to_pandas().set_index(id_col)
    except Exception:
        cdf = pd.read_parquet(cells)
        return cdf.set_index(cdf.columns[0]) if "cell_id" not in cdf.columns else cdf.set_index("cell_id")


def _xenium_panel_meta(path: Path) -> dict:
    """Panel metadata from a Xenium run's ``experiment.xenium`` (name / organism / tissue / target
    count). Returns ``{}`` when the file is absent or unreadable - never raises.

    The total target count is predesigned + custom (a 380-gene predesigned panel with a 100-gene
    add-on reports as 480). Only non-empty fields are returned.
    """
    meta = Path(path) / "experiment.xenium"
    if not meta.exists():
        return {}
    try:
        import json
        m = json.loads(meta.read_text())
    except Exception:
        return {}
    n_tot = int(m.get("panel_num_targets_predesigned") or 0) + int(m.get("panel_num_targets_custom") or 0)
    out = {"name": m.get("panel_name"), "organism": m.get("panel_organism"),
           "tissue_type": m.get("panel_tissue_type"), "n_targets": (n_tot or None),
           "predesigned_id": m.get("panel_predesigned_id"),
           "n_custom": (m.get("panel_num_targets_custom") or None)}
    return {k: v for k, v in out.items() if v not in (None, "")}


def panel_label(uns) -> Optional[str]:
    """A one-line human panel label from ``adata.uns`` (e.g. ``'hMulti_100g (Human, Multi, 480 targets)'``),
    or ``None`` if the section carries no panel name (non-Xenium / processed inputs)."""
    name = uns.get("panel_name")
    if not name:
        return None
    p = uns.get("panel", {}) if isinstance(uns.get("panel"), dict) else {}
    extras = [str(p[k]) for k in ("organism", "tissue_type") if p.get(k)]
    if p.get("n_targets"):
        extras.append(f"{int(p['n_targets'])} targets")
    return f"{name} ({', '.join(extras)})" if extras else str(name)


def infer_tissue(adata) -> Optional[str]:
    """Free-text tissue context inferred from a Xenium panel's ``uns['panel']`` (organism + tissue
    type), e.g. ``'mouse brain'`` for the Mouse Brain panel or ``'human breast'``. This is what the
    reference chooser keys off (``reference.choose_reference`` matches 'mouse'/'brain' keywords), so a
    freshly-loaded section can auto-select its own reference instead of defaulting to the melanoma
    tissue. Returns ``None`` when the section carries no panel organism/tissue (non-Xenium / processed
    inputs) - the caller then keeps whatever tissue it already had.
    """
    p = adata.uns.get("panel", {}) if isinstance(adata.uns.get("panel"), dict) else {}
    parts = [str(p[k]).strip() for k in ("organism", "tissue_type") if p.get(k)]
    tissue = " ".join(parts).strip().lower()
    return tissue or None


def _read_xenium(path: Path):
    """Direct Xenium reader (demo path): cell_feature_matrix.h5 + cells.parquet.

    Avoids ``spatialdata_io`` so ``obsm['spatial']`` is populated straight away. Also stamps the
    panel metadata from ``experiment.xenium`` (``uns['panel']`` + ``uns['panel_name']``) when present.
    """
    import scanpy as sc

    h5 = path / "cell_feature_matrix.h5"
    cells = path / "cells.parquet"
    adata = sc.read_10x_h5(h5, gex_only=False)
    adata.var_names_make_unique()

    if cells.exists():
        cdf = _read_cells_table(cells).reindex(adata.obs_names)
        xcol = next((c for c in ("x_centroid", "x", "center_x") if c in cdf.columns), None)
        ycol = next((c for c in ("y_centroid", "y", "center_y") if c in cdf.columns), None)
        if xcol and ycol:
            adata.obsm["spatial"] = cdf[[xcol, ycol]].to_numpy(dtype=float)
        for extra in ("transcript_counts", "cell_area", "nucleus_area"):
            if extra in cdf.columns:
                adata.obs[extra] = cdf[extra].to_numpy()

    pm = _xenium_panel_meta(path)
    if pm:
        adata.uns["panel"] = pm
        if pm.get("name"):
            adata.uns["panel_name"] = str(pm["name"])

    tx = path / "transcripts.parquet"
    return adata, (tx if tx.exists() else None)


def _read_merscope(path: Path):
    """Direct MERSCOPE / Vizgen reader: resolve the CSV filenames, then delegate to squidpy.

    ``squidpy.read.vizgen`` takes the counts + metadata filenames as **required
    keyword-only** args (it will not infer them from the directory), so we glob the
    Vizgen naming convention - ``*cell_by_gene*.csv`` / ``*cell_metadata*.csv`` - rather
    than assuming fixed names (the export prefix drifts across pipeline versions, the same
    reason ``detect_platform`` matches on substrings). The reader populates
    ``obsm['spatial']`` from the metadata ``center_x`` / ``center_y`` columns and moves the
    platform's ``Blank-*`` negative-control probes into ``obsm['blank_genes']`` (so ``var``
    is already panel-only and the returned ``control_mask`` is all-False by construction).
    """
    import squidpy as sq

    def _find(token: str) -> str:
        hits = sorted(p.name for p in path.iterdir() if token in p.name and p.name.endswith(".csv"))
        if not hits:
            raise FileNotFoundError(f"MERSCOPE reader: no '*{token}*.csv' found in {path}")
        return hits[0]

    return sq.read.vizgen(path, counts_file=_find("cell_by_gene"), meta_file=_find("cell_metadata"))


def _fov_file_is_indexable(p: Path) -> bool:
    """True when ``squidpy`` can index this FOV-positions file.

    ``squidpy.read.nanostring`` hardcodes ``index_col="fov"`` for that file (it is NOT a parameter),
    while real AtoMx exports head the column ``FOV`` - so passing one raised
    ``ValueError: Index fov invalid`` and NO real CosMx export could be loaded at all.

    Skipping an unindexable one is lossless: squidpy sources ``obsm['spatial']``, ``obsm['spatial_fov']``
    and ``obs['fov']`` from the META file, and the fov_file only decorates
    ``uns['spatial'][<fov>]['metadata']`` - which nothing in this package reads.
    """
    import pandas as pd

    try:
        return "fov" in pd.read_csv(p, nrows=0).columns
    except Exception:
        return False


def _read_cosmx(path: Path):
    """Direct CosMx / NanoString reader: resolve the flat-file names, then delegate to squidpy.

    ``squidpy.read.nanostring`` takes the counts + metadata (+ optional fov) filenames as **required
    keyword-only** args (it will NOT infer them from the directory), so - exactly like the MERSCOPE
    reader - we glob the NanoString naming convention rather than assuming fixed names: the export
    prefix drifts across AtoMx versions (``<slide>_exprMat_file.csv`` / ``_metadata_file`` /
    ``_fov_positions_file``, plain or ``.csv.gz``). A bare ``sq.read.nanostring(path)`` raised
    TypeError on every real CosMx export; this fixes it. squidpy populates ``obsm['spatial']`` from
    the metadata ``CenterX_global_px`` / ``CenterY_global_px`` columns.

    ponytail: squidpy reads the counts CSV densely, so a ~1M-cell 1000-plex export needs >64 GB.
    Raise the job's memory rather than expecting this reader to chunk.
    """
    import squidpy as sq

    def _find(token: str, required: bool = True) -> Optional[str]:
        hits = sorted(p.name for p in path.iterdir()
                      if token in p.name and (p.name.endswith(".csv") or p.name.endswith(".csv.gz")))
        if not hits:
            if required:
                raise FileNotFoundError(f"CosMx reader: no '*{token}*.csv[.gz]' found in {path}")
            return None
        return hits[0]

    fov_file = _find("fov_positions", required=False)
    if fov_file is not None and not _fov_file_is_indexable(path / fov_file):
        fov_file = None

    return sq.read.nanostring(path, counts_file=_find("exprMat"), meta_file=_find("metadata_file"),
                              fov_file=fov_file)


def _read_flex(path: Path):
    """10x Flex / dissociated scRNA reader: a single ``filtered_feature_bc_matrix.h5`` (or ``matrix.h5``).

    Flex is not imaging-based, so there is no ``obsm['spatial']`` - QC / clustering / annotation /
    reference-transfer all run, only the spatial canvas and neighbourhood statistics are unavailable
    (they read ``obsm['spatial']`` and simply do not fire). Returns ``adata`` (no transcript table).
    """
    import scanpy as sc

    h5 = next((p for p in sorted(path.iterdir(), key=lambda q: q.name)
               if p.name.endswith(".h5") and ("filtered_feature_bc_matrix" in p.name or p.name == "matrix.h5")),
              None)
    if h5 is None:
        raise FileNotFoundError(f"Flex reader: no '*filtered_feature_bc_matrix.h5' / 'matrix.h5' in {path}")
    adata = sc.read_10x_h5(h5, gex_only=False)
    adata.var_names_make_unique()
    return adata


def _read_visiumhd(path: Path):
    """VisiumHD (SpaceRanger) reader.

    Prefers the **8um binned** matrix (``binned_outputs/square_008um``) - bins act as pseudo-cells and
    ``spatial/tissue_positions.parquet`` gives ``obsm['spatial']`` (full-res pixel centroids). Falls back
    to the **segmented cell** matrix (``segmented_outputs/filtered_feature_cell_matrix.h5``, real cells)
    when a section shipped only the segmentation; centroids are set only if a positions parquet is present,
    otherwise the section still runs QC / cluster / annotate without a canvas. The 8um bin is the standard
    downstream unit; a finer 2um bin would be ~10x the cells for no annotation gain on a targeted question.
    """
    import scanpy as sc

    binned = path / "binned_outputs" / "square_008um"
    bh5 = binned / "filtered_feature_bc_matrix.h5"
    if bh5.exists():
        adata = sc.read_10x_h5(bh5, gex_only=False)
        adata.var_names_make_unique()
        tp = binned / "spatial" / "tissue_positions.parquet"
        if tp.exists():
            pos = pd.read_parquet(tp)
            if "barcode" in pos.columns:
                pos = pos.set_index(pos["barcode"].astype(str))
            xy = pos.reindex(adata.obs_names.astype(str))
            if {"pxl_col_in_fullres", "pxl_row_in_fullres"} <= set(xy.columns):
                adata.obsm["spatial"] = xy[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(dtype=float)
        return adata

    sh5 = path / "segmented_outputs" / "filtered_feature_cell_matrix.h5"
    if sh5.exists():
        adata = sc.read_10x_h5(sh5, gex_only=False)
        adata.var_names_make_unique()
        return adata

    raise FileNotFoundError(
        f"VisiumHD reader: neither binned_outputs/square_008um/filtered_feature_bc_matrix.h5 nor "
        f"segmented_outputs/filtered_feature_cell_matrix.h5 found under {path}"
    )


def _transcripts_have_z(transcripts_path: Optional[Path]) -> bool:
    """Cheaply peek the transcript parquet schema for a z-coordinate column."""
    if transcripts_path is None or not Path(transcripts_path).exists():
        return False
    try:
        import pyarrow.parquet as pq

        cols = {c.lower() for c in pq.ParquetFile(transcripts_path).schema.names}
        return any(z in cols for z in ("z", "z_location", "global_z"))
    except Exception:
        return False


def load(path: str | Path, platform: str = "auto") -> SpatialSample:
    """Load one spatial section into a :class:`SpatialSample`.

    Parameters
    ----------
    path: directory of platform output, or a ``.h5ad`` file.
    platform: "auto" (default) or one of xenium/cosmx/merscope/atera/h5ad.
    """
    path = Path(path)
    if platform == "auto":
        platform = detect_platform(path)

    transcripts_path: Optional[Path] = None
    sdata = None

    if platform in ("xenium", "atera"):
        adata, transcripts_path = _read_xenium(path)
        # Record the cell-segmentation polygon file if the run dir ships one (Xenium
        # cell_boundaries.parquet) so `views.segmentation_figure` can draw real cell shapes.
        if path.is_dir() and (path / "cell_boundaries.parquet").exists():
            adata.uns["boundaries_path"] = str(path / "cell_boundaries.parquet")
        # Re-stamp Atera by panel size (~18k WTA) when routed through the Xenium path.
        if platform == "xenium" and adata.n_vars > 12000:
            platform = "atera"
    elif platform == "cosmx":
        adata = _read_cosmx(path)
    elif platform == "merscope":
        adata = _read_merscope(path)
    elif platform == "visiumhd":
        adata = _read_visiumhd(path)
    elif platform == "flex":
        adata = _read_flex(path)
    elif platform == "h5ad":
        import anndata as ad

        adata = ad.read_h5ad(path)
    else:  # pragma: no cover
        raise ValueError(f"Unsupported platform: {platform}")

    control_mask = build_control_mask(adata)
    panel_genes = adata.var_names[~control_mask].tolist()
    adata.var["control"] = control_mask       # broad: everything that is not a panel gene
    adata.var["neg_control"] = build_neg_control_mask(adata)   # strict: qc_vars for pct_counts_control
    # Stamp platform + panel size on uns so every adata-only downstream step (count floor,
    # confidence gate, section-QC flagging) is panel/platform aware without re-threading `sample`.
    adata.uns["platform"] = platform
    adata.uns["n_panel_genes"] = int(len(panel_genes))

    return SpatialSample(
        platform=platform,
        adata=adata,
        control_mask=control_mask,
        panel_genes=panel_genes,
        has_z=_transcripts_have_z(transcripts_path),
        transcripts_path=transcripts_path,
        sdata=sdata,
        meta={"source_path": str(path)},
    )
