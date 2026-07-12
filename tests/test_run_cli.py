"""The headless CLI (``scripts/run.py``) runs the SAME spine as the app and emits regeneration artifacts.

Contract: on a real (synthetic) section it exits 0, writes ``annotated.h5ad`` + ``rerun.py`` +
``run.json``, types the cells (``cell_type`` present), and reports NO failed stage - an optional arm
that reports ``skipped:`` is fine. A non-zero exit iff a stage failed is what lets an agent assert
success at the process level, per the pipeline's honest per-stage status contract. A second test
pins the clean usage-error exit (2) for an unloadable path.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")
pytest.importorskip("scanpy")

# scripts/run.py is a standalone module (run as `python -m scripts.run`), not part of the package.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))


def _synthetic_h5ad(path: str) -> str:
    """Two spatially-separated populations (T vs melanocyte) + noise + controls, written to .h5ad.

    Mirrors test_pipeline_synthetic so clustering + marker labels resolve deterministically offline.
    """
    rng = np.random.default_rng(0)
    t = ["CD3D", "CD3E", "CD2", "CD8A", "TRAC"]
    mel = ["MLANA", "PMEL", "TYR", "DCT", "MITF", "SOX10"]
    noise = [f"NOISE{i}" for i in range(30)]
    controls = ["NegControlProbe_00001", "BLANK_0001"]
    genes = t + mel + noise + controls
    n = 400
    grp = np.array([0] * (n // 2) + [1] * (n - n // 2))
    X = rng.poisson(0.3, size=(n, len(genes))).astype("float32")
    gi = {g: i for i, g in enumerate(genes)}
    for g in t:
        X[grp == 0, gi[g]] += rng.poisson(8, (grp == 0).sum())
    for g in mel:
        X[grp == 1, gi[g]] += rng.poisson(8, (grp == 1).sum())
    var = pd.DataFrame(index=genes)
    var["feature_types"] = ["Gene Expression"] * (len(genes) - 2) + ["Negative Control Probe"] * 2
    a = anndata.AnnData(X=X, var=var)
    a.obs_names = [f"cell{i}" for i in range(n)]
    a.obsm["spatial"] = np.column_stack([rng.normal(np.where(grp == 0, 0, 50), 4),
                                         rng.normal(0, 4, n)])
    a.write_h5ad(path)
    return path


def test_cli_runs_spine_and_emits_regeneration_artifacts(tmp_path):
    import run as cli  # scripts/run.py

    section = _synthetic_h5ad(str(tmp_path / "section.h5ad"))
    out = str(tmp_path / "out")
    rc = cli.main(["--section", section, "--tissue", "melanoma", "--out", out, "--resolution", "1.0"])

    assert rc == 0
    for f in ("annotated.h5ad", "rerun.py", "run.json"):
        assert os.path.exists(os.path.join(out, f)), f"missing output: {f}"

    rec = json.loads(open(os.path.join(out, "run.json")).read())
    # No stage may FAIL; skips (unconfigured optional arms, e.g. RCTD/SPLIT) are fine.
    assert rec["summary"].get("failed", []) == [], rec["summary"]
    assert all(s["status"] == "ok" or s["status"].startswith("skipped:")
               for s in rec["stages"]), rec["stages"]

    # The point of the run: the section is typed, and rerun.py actually reproduces the pipeline.
    typed = anndata.read_h5ad(os.path.join(out, "annotated.h5ad"))
    assert "cell_type" in typed.obs
    rerun = open(os.path.join(out, "rerun.py")).read()
    assert "cluster" in rerun and "annotate" in rerun


def test_cli_bad_path_returns_nonzero(tmp_path):
    import run as cli

    rc = cli.main(["--section", str(tmp_path / "does_not_exist.h5ad"), "--out", str(tmp_path / "o")])
    assert rc == 2
