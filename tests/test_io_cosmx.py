"""CosMx / NanoString flat-file reader.

Reproduces (and guards against) the bug where ``io.load`` called ``sq.read.nanostring(path)`` with
no ``counts_file``/``meta_file`` - squidpy makes those REQUIRED keyword-only args, so it raised
TypeError on any real CosMx flat-file directory (the same class of bug that was fixed for MERSCOPE;
prior "CosMx validated" runs used a prebuilt .h5ad, never this reader).

The synthetic section uses MOUSE (Title-case) gene symbols, so it also exercises the mouse code path
(no real mouse CosMx dataset exists on the system).
"""
from __future__ import annotations

import gzip

import numpy as np
import pandas as pd

from spatialscribe.analysis import io


# CosMx Mouse Neuroscience-style symbols (Title-case) + two negative probes.
_MOUSE_GENES = ["Snap25", "Gfap", "Aqp4", "Cx3cr1", "Pdgfra", "Mbp", "Slc17a7", "Gad1",
                "Cldn5", "Flt1", "Aif1", "Rbfox3", "Meg3", "Plp1", "Sox10", "Vip",
                "Sst", "Pvalb", "Olig1", "Olig2", "NegPrb1", "NegPrb2"]


def _write_cosmx_flatfiles(d, prefix="MouseBrain", n_fov=2, per_fov=40, gzip_out=False):
    """Write a minimal but valid NanoString CosMx flat-file trio into directory ``d``.

    Genes = the named mouse markers + ~70 filler genes whose per-gene rates SPAN the expression range
    (0.3-25), so ``sc.tl.score_genes`` (used by annotate/states) can draw control genes from every
    expression bin - a hand-built fixture with only uniformly-low genes raises 'No control genes
    found' (real sections have hundreds of genes and never hit this; see project CLAUDE.md)."""
    rng = np.random.default_rng(0)
    fillers = [f"Gm{1000 + i}" for i in range(70)]              # mouse-style predicted-gene symbols
    genes = _MOUSE_GENES + fillers
    rates = np.concatenate([np.full(len(_MOUSE_GENES), 4.0), rng.uniform(0.3, 25.0, len(fillers))])
    rows_expr, rows_meta = [], []
    for fov in range(1, n_fov + 1):
        for cid in range(1, per_fov + 1):
            counts = rng.poisson(rates).astype(int)
            rows_expr.append({"fov": fov, "cell_ID": cid, **dict(zip(genes, counts))})
            lx, ly = float(rng.uniform(0, 5000)), float(rng.uniform(0, 5000))
            rows_meta.append({"fov": fov, "cell_ID": cid,
                              "CenterX_local_px": lx, "CenterY_local_px": ly,
                              "CenterX_global_px": lx + fov * 6000.0, "CenterY_global_px": ly,
                              "Area": float(rng.uniform(50, 400)),
                              "nCount_RNA": int(counts.sum())})
    expr = pd.DataFrame(rows_expr)
    meta = pd.DataFrame(rows_meta)
    fov_pos = pd.DataFrame({"fov": list(range(1, n_fov + 1)),
                            "x_global_px": [i * 6000.0 for i in range(n_fov)],
                            "y_global_px": [0.0] * n_fov})
    ext = ".csv.gz" if gzip_out else ".csv"

    def _write(df, name):
        p = d / f"{prefix}_{name}{ext}"
        if gzip_out:
            with gzip.open(p, "wt") as fh:
                df.to_csv(fh, index=False)
        else:
            df.to_csv(p, index=False)
        return p

    _write(expr, "exprMat_file")
    _write(meta, "metadata_file")
    _write(fov_pos, "fov_positions_file")
    return d


def test_detect_and_load_cosmx_flatfiles(tmp_path):
    _write_cosmx_flatfiles(tmp_path)
    assert io.detect_platform(tmp_path) == "cosmx"
    sample = io.load(tmp_path)
    assert sample.platform == "cosmx"
    a = sample.adata
    assert a.n_obs == 80                                     # 2 fov x 40 cells
    assert "spatial" in a.obsm and a.obsm["spatial"].shape == (80, 2)
    # mouse Title-case symbols preserved (biological genes present; a real gene, not a control)
    assert "Snap25" in set(map(str, a.var_names))
    # control probes flagged, not treated as signal
    assert bool(sample.control_mask.any())


def test_load_cosmx_flatfiles_gzipped(tmp_path):
    # Real CosMx exports ship .csv.gz - the reader must resolve those too.
    _write_cosmx_flatfiles(tmp_path, gzip_out=True)
    sample = io.load(tmp_path)
    assert sample.platform == "cosmx"
    assert sample.adata.n_obs == 80
    assert "spatial" in sample.adata.obsm


def test_load_cosmx_uppercase_fov_positions_header(tmp_path):
    """Real AtoMx exports head the FOV-positions file ``FOV``, not ``fov``.

    ``squidpy.read.nanostring`` hardcodes ``index_col="fov"`` for that file (it is not a parameter),
    so passing it raised ``ValueError: Index fov invalid`` and NO real CosMx export could be loaded.
    Both header casings must load. Guards the 2026-07 overnight-benchmark failure.
    """
    _write_cosmx_flatfiles(tmp_path)
    p = tmp_path / "MouseBrain_fov_positions_file.csv"
    pd.read_csv(p).rename(columns={"fov": "FOV"}).to_csv(p, index=False)

    sample = io.load(tmp_path)
    assert sample.platform == "cosmx"
    assert sample.adata.n_obs == 80
    assert "spatial" in sample.adata.obsm


def test_load_cosmx_without_fov_positions_file(tmp_path):
    """The FOV-positions file is optional - an export missing it must still load."""
    _write_cosmx_flatfiles(tmp_path)
    (tmp_path / "MouseBrain_fov_positions_file.csv").unlink()
    assert io.load(tmp_path).adata.n_obs == 80
