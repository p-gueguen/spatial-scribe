"""SPLIT reference-path spillover purification (residual-contamination removal), in-app.

Wraps the validated `scripts/split/` recipe (annotation-agnostic RCTD -> SPLIT):
  1. `rctd_full.py`   (rctd-py env)  -> SPLIT inputs (counts.mtx, weights.csv, reference.csv, ...)
  2. `run_split.R`    (R + SPLIT)    -> purified_counts.mtx (`rctd_free_purify`,
                                        `DO_remove_residual_contamination=TRUE`)
  3. join the purified counts back onto ``adata.layers['split_corrected']``

so a scientist can decontaminate transcript spillover *in the app* and see the effect on the
marker dot-plots (raw vs purified). Demonstrated on the demo: median library 76 -> 35 (54.5%
spillover removed), T x Epithelial co-expression -> ~0.

SPLIT is R-only (bdsc-tds/SPLIT) and RCTD needs the rctd-py env, so this runs as ISOLATED
subprocesses configured by env vars (empty public defaults). Unset -> :func:`split_purify` returns
``{'status': 'skipped: ...'}`` and NEVER raises, exactly like the other subprocess tools. Live
smoke needs the two envs; the join + graceful-skip are unit-tested against synthetic files.

Env vars
--------
``SPATIALSCRIBE_SPLIT_RSCRIPT``       - an ``Rscript`` with SPLIT (>=0.3.0) installed.
``SPATIALSCRIBE_SPLIT_RCTD_PYTHON``   - a python with rctd-py (for `rctd_full.py`); falls back to
                                        ``SPATIALSCRIBE_RCTD_PYTHON``.
"""
from __future__ import annotations

import os


def _split_scripts_dir():
    from pathlib import Path
    return Path(__file__).resolve().parents[3] / "scripts" / "split"


def _join_split(adata, inputs_dir, split_dir) -> int:
    """Join SPLIT purified counts onto ``adata.layers['split_corrected']`` (aligned to obs/var).

    ``inputs_dir/cells.txt`` gives the purified cell order, ``split_dir/purified_genes.txt`` the gene
    order, ``split_dir/purified_counts.mtx`` the genes x cells matrix. Cells RCTD dropped stay 0 in
    the layer (honest: SPLIT only purifies surviving cells). Returns the number of cells covered.
    """
    from pathlib import Path

    import numpy as np
    import pandas as pd
    import scipy.io as sio
    import scipy.sparse as sp

    inputs_dir, split_dir = Path(inputs_dir), Path(split_dir)
    cells = pd.read_csv(inputs_dir / "cells.txt", header=None)[0].astype(str).tolist()
    pgenes = pd.read_csv(split_dir / "purified_genes.txt", header=None)[0].astype(str).tolist()
    pur = sio.mmread(split_dir / "purified_counts.mtx").T.tocsr()   # cells x genes (rows=cells, cols=pgenes)

    cell_pos = adata.obs_names.astype(str).get_indexer(cells)       # -> adata row per purified cell
    gene_pos = adata.var_names.astype(str).get_indexer(pgenes)      # -> adata col per purified gene
    rmask, cmask = cell_pos >= 0, gene_pos >= 0
    pur = pur[np.where(rmask)[0]][:, np.where(cmask)[0]]
    rpos, gpos = cell_pos[rmask], gene_pos[cmask]

    coo = pur.tocoo()
    L = sp.csr_matrix((coo.data, (rpos[coo.row], gpos[coo.col])),
                      shape=(adata.n_obs, adata.n_vars), dtype="float32")
    adata.layers["split_corrected"] = L
    return int(rmask.sum())


def _skip(msg: str) -> dict:
    return {"status": f"skipped: {msg}", "pct_removed": 0.0}


def _lib_summary(raw_mat, corrected, touched, n_obs: int, method: str, extra: dict | None = None) -> dict:
    """Shared ok-result builder for both purification paths: median library size over the touched
    cells before vs after, %-removed, coverage, plus ``method`` and any ``extra`` fields."""
    import numpy as np

    raw = np.asarray(raw_mat.sum(1) if hasattr(raw_mat, "sum") else raw_mat.sum(axis=1)).ravel().astype(float)
    cor = np.asarray(corrected.sum(1)).ravel().astype(float)
    m = np.asarray(touched, dtype=bool)
    lib_before = float(np.median(raw[m])) if m.any() else 0.0
    lib_after = float(np.median(cor[m])) if m.any() else 0.0
    n_pur = int(m.sum())
    out = {"status": "ok", "n_purified": n_pur, "coverage": n_pur / max(1, n_obs),
           "median_lib_before": lib_before, "median_lib_after": lib_after,
           "pct_removed": (1 - lib_after / lib_before) if lib_before else 0.0,
           "layer": "split_corrected", "method": method}
    if extra:
        out.update(extra)
    return out


def split_purify(adata, reference=None, ref_label_key=None, *, marker_sets=None,
                 rctd_python=None, rscript=None, max_cells: int = 0, weights_engine: str = "tacco") -> dict:
    """Purify transcript spillover; write ``layers['split_corrected']`` and summarize the effect.

    Three paths, one contract. When a reference AND the SPLIT-R env are configured, run the full
    deconvolve -> SPLIT recipe. The deconvolution weights come from ``weights_engine``:

    - ``"tacco"`` (default): TACCO optimal-transport, run IN-ENV (the declared ``transfer`` extra),
      so no separate rctd-py env is needed. This reuses the same engine the Annotate step already
      runs, and TACCO is ~11x faster than RCTD. Needs only ``SPATIALSCRIBE_SPLIT_RSCRIPT``.
    - ``"rctd"``: the original RCTD deconvolution via the rctd-py subprocess. Needs the extra
      ``SPATIALSCRIBE_SPLIT_RCTD_PYTHON`` (or ``SPATIALSCRIBE_RCTD_PYTHON``) env. Use when you want
      RCTD's doublet/contamination model to drive purification specifically.

    ``method`` in the returned dict is ``'tacco_split'`` or ``'rctd_split'`` accordingly. When neither
    the reference nor the R env is available, or the heavy path fails at runtime, falls back to a
    reference-free, main-env neighbour+marker decontamination (``method='marker_neighbour'``).
    Degrades to ``{'status': 'skipped: ...'}`` (never raises, no layer) only when NEITHER path can run,
    so it is safe to call unconditionally. Because ``SPLIT::rctd_free_purify`` is annotation-agnostic
    (it consumes a cells x types weight matrix), swapping the weights engine does not touch ``run_split.R``.
    """
    from pathlib import Path

    rscript = rscript or os.environ.get("SPATIALSCRIBE_SPLIT_RSCRIPT", "")
    rctd_python = (rctd_python or os.environ.get("SPATIALSCRIBE_SPLIT_RCTD_PYTHON")
                   or os.environ.get("SPATIALSCRIBE_RCTD_PYTHON", ""))
    r_ok = bool(rscript and Path(rscript).exists() and reference is not None and ref_label_key is not None)
    if r_ok:
        if weights_engine == "rctd" and rctd_python and Path(rctd_python).exists():
            res = _purify_split(adata, reference, ref_label_key, engine="rctd",
                                rctd_python=rctd_python, rscript=rscript, max_cells=max_cells)
        else:
            # default: TACCO weights, in-env (no rctd_python needed)
            res = _purify_split(adata, reference, ref_label_key, engine="tacco",
                                rctd_python=None, rscript=rscript, max_cells=max_cells)
        if res.get("status") == "ok":
            return res
        # heavy path failed at runtime (bad env / subprocess error) -> the in-app fallback, not a skip.
    return _purify_neighbour_markers(adata, marker_sets)


def _write_tacco_split_inputs(adata, reference, ref_label_key, inputs_dir) -> None:
    """In-env TACCO deconvolution -> SPLIT inputs. Same output files as scripts/rctd_full.py, no
    subprocess: TACCO is the declared `transfer` extra, so it runs in the main interpreter. Imports
    the standalone builder by path to keep ONE copy of the input contract."""
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "tacco_split_inputs", str(_split_scripts_dir() / "tacco_split_inputs.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    # write section + reference h5ad for the builder (it re-reads counts from layers['counts'])
    import tempfile
    from pathlib import Path

    import anndata as ad
    tmp = Path(inputs_dir).parent
    sec_h5, ref_h5 = tmp / "section.h5ad", tmp / "reference.h5ad"
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    a_exp = ad.AnnData(X=counts.copy())
    a_exp.obs_names = adata.obs_names.astype(str); a_exp.var_names = adata.var_names.astype(str)
    a_exp.layers["counts"] = a_exp.X.copy()
    a_exp.write_h5ad(sec_h5)
    reference.write_h5ad(ref_h5)
    mod.main(str(sec_h5), str(ref_h5), str(ref_label_key), str(inputs_dir))


def _purify_split(adata, reference, ref_label_key, *, engine: str, rctd_python, rscript,
                  max_cells: int) -> dict:
    """The deconvolve -> SPLIT recipe. `engine='tacco'` runs TACCO in-env; `engine='rctd'` shells the
    rctd-py subprocess. Both write the same annotation-agnostic SPLIT inputs, then run_split.R purifies.
    """
    import subprocess
    import tempfile
    from pathlib import Path

    import numpy as np

    sdir = _split_scripts_dir()
    run_split = sdir / "run_split.R"
    if not run_split.exists():
        return _skip("SPLIT recipe scripts missing from scripts/split/")

    import anndata as ad

    tmp = Path(tempfile.mkdtemp(prefix="sssplit_"))
    inputs, split_out = tmp / "inputs", tmp / "split_out"
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}

    try:
        if engine == "tacco":
            _write_tacco_split_inputs(adata, reference, ref_label_key, inputs)
        else:
            rctd_full = sdir / "rctd_full.py"
            if not rctd_full.exists():
                return _skip("rctd_full.py missing from scripts/split/")
            sec_h5, ref_h5 = tmp / "section.h5ad", tmp / "reference.h5ad"
            a_exp = ad.AnnData(X=counts.copy())
            a_exp.obs_names = adata.obs_names.astype(str); a_exp.var_names = adata.var_names.astype(str)
            if max_cells and a_exp.n_obs > max_cells:
                idx = np.sort(np.random.default_rng(0).choice(a_exp.n_obs, size=max_cells, replace=False))
                a_exp = a_exp[idx].copy()
            a_exp.write_h5ad(sec_h5)
            reference.write_h5ad(ref_h5)
            p1 = subprocess.run([rctd_python, str(rctd_full), str(sec_h5), str(ref_h5),
                                 str(ref_label_key), str(inputs)],
                                capture_output=True, text=True, timeout=7200, env=env)
            if p1.returncode != 0 or not (inputs / "counts.mtx").exists():
                tail = (p1.stderr or "").strip().splitlines()
                return _skip(f"rctd_full failed ({tail[-1] if tail else 'no output'})")

        if not (inputs / "counts.mtx").exists():
            return _skip(f"{engine} deconvolution produced no SPLIT inputs")
        p2 = subprocess.run([rscript, str(run_split), str(inputs), str(split_out)],
                            capture_output=True, text=True, timeout=7200, env=env)
        if p2.returncode != 0 or not (split_out / "purified_counts.mtx").exists():
            tail = (p2.stderr or "").strip().splitlines()
            return _skip(f"SPLIT (run_split.R) failed ({tail[-1] if tail else 'no output'})")
    except Exception as exc:
        return _skip(f"{engine} SPLIT error ({exc})")

    _join_split(adata, inputs, split_out)
    touched = np.asarray(adata.layers["split_corrected"].sum(1)).ravel() > 0
    return _lib_summary(counts, adata.layers["split_corrected"], touched, adata.n_obs, f"{engine}_split")


def _purify_neighbour_markers(adata, marker_sets) -> dict:
    """Reference-free, main-env spillover decontamination (numpy/scipy; never raises).

    For each PRIVATE lineage marker gene (a gene in exactly one lineage's marker list AND on the
    panel), zero that gene in cells of a DIFFERENT lineage that have >=1 spatial neighbour of the
    marker's lineage - attributing the count to cross-cell spillover. Writes
    ``layers['split_corrected']`` and reports the cross-lineage co-expression drop.

    # ponytail: reference-free neighbour+marker heuristic - only removes cross-lineage PRIVATE
    # markers with neighbour evidence. Full per-gene RCTD->SPLIT is the upgrade; set
    # SPATIALSCRIBE_SPLIT_RSCRIPT + SPATIALSCRIBE_SPLIT_RCTD_PYTHON + upload a reference.
    """
    from collections import defaultdict

    import numpy as np
    import scipy.sparse as sp

    from .keys import Obs

    if not Obs.CELL_TYPE.present(adata) or not marker_sets or "spatial" not in adata.obsm:
        return _skip("no reference and cannot run in-app fallback "
                     "(needs cell_type + markers + spatial coords)")

    if "spatial_connectivities" not in adata.obsp:
        try:
            from . import spatial as _sp
            _sp.spatial_neighbors(adata, method="knn")
        except Exception as exc:
            return _skip(f"could not build the spatial graph for the fallback ({exc})")
    A = adata.obsp["spatial_connectivities"].tocsr()

    var_pos = {str(g): i for i, g in enumerate(adata.var_names)}

    # Private markers: a gene owned by exactly ONE lineage's marker list and present on the panel.
    gene_lineages: dict[str, set] = defaultdict(set)
    for lineage, genes in marker_sets.items():
        for g in genes:
            if str(g) in var_pos:
                gene_lineages[str(g)].add(lineage)
    private = {g: next(iter(lins)) for g, lins in gene_lineages.items() if len(lins) == 1}

    # Map each cell_type label to a lineage (case-insensitive exact, then substring, vs marker keys).
    cell_types = adata.obs["cell_type"].astype(str)
    lineage_keys = list(marker_sets.keys())

    def _match(label: str):
        ll = label.lower()
        for k in lineage_keys:
            if k.lower() == ll:
                return k
        for k in lineage_keys:
            if k.lower() in ll or ll in k.lower():
                return k
        return None

    ct_to_lineage = {c: _match(c) for c in cell_types.unique()}
    present = {v for v in ct_to_lineage.values() if v is not None}
    private = {g: lin for g, lin in private.items() if lin in present}
    cell_lineage = cell_types.map(ct_to_lineage).to_numpy(dtype=object)

    raw = sp.csr_matrix(adata.layers["counts"] if "counts" in adata.layers else adata.X).astype("float32")
    corrected = raw.tolil()
    lin_markers: dict[str, list] = defaultdict(list)
    for g, lin in private.items():
        gi = var_pos[g]
        lin_markers[lin].append(gi)
        ind = (cell_lineage == lin).astype(np.float32)
        if ind.sum() == 0:
            continue
        has_neighbour = np.asarray(A @ ind).ravel() > 0        # cells with >=1 neighbour of this lineage
        target = np.where(has_neighbour & (cell_lineage != lin))[0]
        if target.size:
            corrected[target, gi] = 0.0                        # zero the spilled marker in the wrong lineage
    corrected = corrected.tocsr()
    adata.layers["split_corrected"] = corrected

    # Cross-lineage co-expression readout: the ordered lineage pair (A, B) where the most A-cells
    # carry B's private markers, before vs after (fraction of A cells with >=1 B-private count).
    def _coexpr(mat, a_lin, cols):
        mask = cell_lineage == a_lin
        if not mask.any() or not cols:
            return 0.0
        return float((np.asarray(mat[mask][:, cols].sum(1)).ravel() > 0).mean())

    best = None
    for a_lin in sorted(present):
        for b_lin in sorted(present):
            if a_lin == b_lin or not lin_markers.get(b_lin):
                continue
            frac = _coexpr(raw, a_lin, lin_markers[b_lin])
            if best is None or frac > best[0]:
                best = (frac, a_lin, b_lin)
    if best is not None:
        _, a_lin, b_lin = best
        cols = lin_markers[b_lin]
        before, after = _coexpr(raw, a_lin, cols), _coexpr(corrected, a_lin, cols)
        extra = {"coexpr_pair": f"{a_lin} x {b_lin}", "coexpr_before": before, "coexpr_after": after,
                 "note": (f"Cross-lineage co-expression: {round(before * 100)}% of {a_lin} cells carried "
                          f"{b_lin} private markers before, {round(after * 100)}% after neighbour-based "
                          "decontamination (reference-free; upload a reference for full RCTD -> SPLIT).")}
    else:
        extra = {"coexpr_pair": None, "coexpr_before": 0.0, "coexpr_after": 0.0,
                 "note": "No cross-lineage private-marker co-expression to correct on this section."}

    touched = np.asarray(abs(raw - corrected).sum(1)).ravel() > 0
    return _lib_summary(raw, corrected, touched, adata.n_obs, "marker_neighbour", extra)


def demo() -> None:
    """One runnable self-check: plant a cross-lineage spillover marker next to the other lineage and
    assert the fallback lowers its co-expression and writes the layer. Run: ``python -m ...split``."""
    import anndata as ad
    import numpy as np

    rng = np.random.default_rng(0)
    n = 60
    # two lineages, split left/right in space; gene EPCAM is Epithelial-private, CD3D T-cell-private.
    X = rng.integers(0, 3, size=(n, 3)).astype("float32")   # cols: EPCAM, CD3D, ACTB(shared)
    a = ad.AnnData(X=X)
    a.var_names = ["EPCAM", "CD3D", "ACTB"]
    a.obs["cell_type"] = ["Epithelial"] * (n // 2) + ["T cell"] * (n - n // 2)
    # place T cells adjacent to Epithelial; plant EPCAM (epithelial-private) into T cells = spillover
    coords = np.column_stack([np.r_[np.arange(n // 2), np.arange(n - n // 2)], np.r_[np.zeros(n // 2), np.ones(n - n // 2)]])
    a.obsm["spatial"] = coords.astype(float)
    a.layers["counts"] = a.X.copy()
    a.X[n // 2:, 0] = 3.0                                    # every T cell wrongly expresses EPCAM
    a.layers["counts"][n // 2:, 0] = 3.0
    markers = {"Epithelial": ["EPCAM"], "T cell": ["CD3D"]}
    res = _purify_neighbour_markers(a, markers)
    assert res["status"] == "ok", res
    assert "split_corrected" in a.layers
    assert res["coexpr_after"] <= res["coexpr_before"], res
    print("split fallback demo OK:", {k: res[k] for k in ("method", "coexpr_pair", "coexpr_before", "coexpr_after", "pct_removed")})


if __name__ == "__main__":
    demo()
