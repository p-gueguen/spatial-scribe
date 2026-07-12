"""Tests for the platform-agnostic ingestion layer (run in the main pixi env)."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from spatialscribe.analysis.io import (
    _xenium_panel_meta, build_control_mask, detect_platform, infer_tissue, panel_label,
)


def _touch(d, *names):
    for n in names:
        (d / n).write_bytes(b"")


def test_detect_xenium(tmp_path):
    _touch(tmp_path, "experiment.xenium", "cell_feature_matrix.h5", "cells.parquet")
    assert detect_platform(tmp_path) == "xenium"


def test_xenium_panel_meta_and_label(tmp_path):
    import json
    (tmp_path / "experiment.xenium").write_text(json.dumps({
        "panel_name": "hMulti_100g", "panel_organism": "Human", "panel_tissue_type": "Multi",
        "panel_num_targets_predesigned": 380, "panel_num_targets_custom": 100,
        "panel_predesigned_id": "hImmune_v1"}))
    pm = _xenium_panel_meta(tmp_path)
    assert pm["name"] == "hMulti_100g" and pm["organism"] == "Human"
    assert pm["n_targets"] == 480                           # predesigned + custom
    label = panel_label({"panel_name": pm["name"], "panel": pm})
    assert label == "hMulti_100g (Human, Multi, 480 targets)"


def test_panel_meta_absent_is_empty(tmp_path):
    assert _xenium_panel_meta(tmp_path) == {}               # no experiment.xenium -> {}
    assert panel_label({}) is None                          # no panel name -> no label


def test_infer_tissue_from_panel():
    # A Xenium Mouse Brain panel names its own organism + tissue -> the reference chooser can key off
    # "mouse brain" instead of the melanoma default (this is what auto-select-on-load rides on).
    a = ad.AnnData(np.zeros((3, 2), dtype="float32"))
    a.uns["panel"] = {"name": "Xenium Mouse Brain Gene Expression", "organism": "Mouse",
                      "tissue_type": "Brain", "n_targets": 247}
    assert infer_tissue(a) == "mouse brain"
    # Organism alone still yields a usable context; no panel -> None (caller keeps its tissue).
    a.uns["panel"] = {"organism": "Human"}
    assert infer_tissue(a) == "human"
    assert infer_tissue(ad.AnnData(np.zeros((3, 2), dtype="float32"))) is None
    assert panel_label({"panel_name": "P1"}) == "P1"        # bare name, no extras


def test_detect_cosmx(tmp_path):
    _touch(tmp_path, "run_exprMat_file.csv", "run_metadata_file.csv")
    assert detect_platform(tmp_path) == "cosmx"


def test_detect_merscope(tmp_path):
    _touch(tmp_path, "cell_by_gene.csv", "cell_metadata.csv")
    assert detect_platform(tmp_path) == "merscope"


def test_detect_classic_visium_rejected(tmp_path):
    # Classic (spot) Visium ships filtered_feature_bc_matrix.h5 + a spatial/ dir but NO binned/
    # segmented outputs. It must be rejected CLEARLY, never silently misrouted to the non-spatial
    # flex reader (which would drop the spot coordinates and load it as dissociated scRNA).
    _touch(tmp_path, "filtered_feature_bc_matrix.h5")
    (tmp_path / "spatial").mkdir()
    with pytest.raises(ValueError, match="classic .spot. Visium"):
        detect_platform(tmp_path)


def test_detect_flex_without_spatial_dir(tmp_path):
    # A bare filtered matrix with NO spatial/ dir is dissociated scRNA (Flex) - still detected as flex,
    # so the Visium guard above does not accidentally swallow real Flex input.
    _touch(tmp_path, "filtered_feature_bc_matrix.h5")
    assert detect_platform(tmp_path) == "flex"


def test_detect_unknown_raises(tmp_path):
    _touch(tmp_path, "random.txt")
    with pytest.raises(ValueError):
        detect_platform(tmp_path)


def _write_merscope(d, n_cells=40, seed=0):
    """Write a minimal Vizgen/MERSCOPE run dir (cell_by_gene.csv + cell_metadata.csv).

    Mirrors the real export contract squidpy's vizgen reader expects: the counts
    CSV is cells x genes with the cell id as the first column, and the metadata
    CSV carries `center_x`/`center_y` (-> obsm['spatial']). Two `Blank-*` probes
    stand in for the platform's negative controls.
    """
    rng = np.random.default_rng(seed)
    genes = ["CD3D", "MLANA", "PECAM1", "EPCAM", "Blank-1", "Blank-2"]
    cells = [str(i) for i in range(n_cells)]
    cbg = pd.DataFrame(rng.poisson(3, size=(n_cells, len(genes))), index=cells, columns=genes)
    cbg.index.name = "cell"
    cbg.to_csv(d / "cell_by_gene.csv")
    meta = pd.DataFrame(
        {
            "fov": 0,
            "volume": rng.uniform(80, 200, n_cells),
            "center_x": rng.uniform(0, 1000, n_cells),
            "center_y": rng.uniform(0, 1000, n_cells),
        },
        index=cells,
    )
    meta.index.name = "cell"
    meta.to_csv(d / "cell_metadata.csv")
    return genes, cells


def test_load_merscope_synthetic(tmp_path):
    """End-to-end MERSCOPE read: detect -> squidpy vizgen -> SpatialSample contract."""
    from spatialscribe.analysis.io import load

    genes, cells = _write_merscope(tmp_path)
    s = load(tmp_path)

    assert s.platform == "merscope"
    assert s.adata.n_obs == len(cells)
    # spatial centroids populated from center_x/center_y
    assert "spatial" in s.adata.obsm
    assert s.adata.obsm["spatial"].shape == (len(cells), 2)
    # squidpy's vizgen reader separates the Blank probes into obsm['blank_genes'],
    # leaving var panel-only (control_mask all-False, panel = the 4 real genes)
    assert "blank_genes" in s.adata.obsm
    assert s.adata.obsm["blank_genes"].shape[1] == 2
    assert s.control_mask.sum() == 0
    assert all("Blank" not in g for g in s.panel_genes)
    assert len(s.panel_genes) == len(genes) - 2
    # platform + panel size stamped on uns for downstream platform-aware QC
    assert s.adata.uns["platform"] == "merscope"
    assert s.adata.uns["n_panel_genes"] == len(genes) - 2


def test_control_mask_from_feature_types():
    var = pd.DataFrame(
        {"feature_types": ["Gene Expression", "Gene Expression", "Negative Control Probe"]},
        index=["CD3D", "MLANA", "NegControlProbe_00001"],
    )
    a = ad.AnnData(X=np.zeros((2, 3), dtype="float32"), var=var)
    mask = build_control_mask(a)
    assert mask.tolist() == [False, False, True]


def test_control_mask_regex_fallback():
    var = pd.DataFrame(index=["CD3D", "BLANK_0001", "NegPrb1", "MLANA"])
    a = ad.AnnData(X=np.zeros((2, 4), dtype="float32"), var=var)
    mask = build_control_mask(a)
    assert mask.tolist() == [False, True, True, False]
