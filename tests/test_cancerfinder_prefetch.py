"""Cancer-Finder subprocess lifecycle: temp dirs never leak, a prefetched run is reused only for the
same input. Uses a fake 'python' that writes a parquet, so no torch/checkpoint is needed."""
from __future__ import annotations

import glob
import os
import stat
import sys
import tempfile

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp

from spatialscribe.analysis import cancerfinder as cf


@pytest.fixture
def fake_cf(tmp_path, monkeypatch):
    """A stand-in CANCERFINDER_PYTHON: reads --h5ad, writes prob = 0.9 for every cell to --out."""
    py = tmp_path / "fake_python"
    py.write_text(f"""#!{sys.executable}
import sys, anndata as ad, pandas as pd
a = sys.argv[sys.argv.index("--h5ad") + 1]; o = sys.argv[sys.argv.index("--out") + 1]
names = ad.read_h5ad(a).obs_names
pd.DataFrame({{"cell_id": names, "cancerfinder_prob": 0.9}}).to_parquet(o)
""")
    py.chmod(py.stat().st_mode | stat.S_IEXEC)
    ckpt = tmp_path / "repo" / "checkpoints" / "sc_pretrain_article.pkl"   # derived from the repo
    ckpt.parent.mkdir(parents=True)
    ckpt.write_text("")
    monkeypatch.setattr(cf, "_CF_PY", str(py))
    monkeypatch.setattr(cf, "_CF_REPO", str(tmp_path / "repo"))
    monkeypatch.setattr(cf, "_CF_CKPT", str(ckpt))
    monkeypatch.setattr(cf, "_PENDING", {})
    launches = []
    real = cf._launch
    monkeypatch.setattr(cf, "_launch", lambda *a, **k: launches.append(1) or real(*a, **k))
    return launches


def _adata(seed=0):
    rng = np.random.default_rng(seed)
    a = ad.AnnData(X=sp.csr_matrix(rng.poisson(1.0, (50, 20)).astype(np.float32)))
    a.obs_names = [f"c{i}" for i in range(50)]
    return a


def _sscf_dirs():
    return set(glob.glob(os.path.join(tempfile.gettempdir(), "sscf_*")))


def test_temp_dir_removed_after_run(fake_cf):
    before = _sscf_dirs()
    r = cf.call_cancerfinder(_adata())
    assert r["status"] == "ok" and r["coverage"] == 1.0
    assert _sscf_dirs() == before          # the 9.6 GB of leaked sscf_* dirs on the server


def test_prefetch_is_reused_for_identical_input(fake_cf):
    a = _adata()
    h = cf.prefetch(a, threshold=0.5, max_cells=25)
    assert h is not None and len(fake_cf) == 1
    r = cf.call_cancerfinder(a, threshold=0.5, max_cells=25)
    assert r["status"] == "ok" and len(fake_cf) == 1      # joined, not relaunched
    assert not cf._PENDING


def test_prefetch_not_reused_when_input_changed(fake_cf):
    a = _adata()
    h = cf.prefetch(a, threshold=0.5, max_cells=25)
    b = a[:40].copy()                                    # different cells -> different answer
    assert cf.call_cancerfinder(b, threshold=0.5, max_cells=25)["status"] == "ok"
    assert len(fake_cf) == 2
    before = _sscf_dirs()
    cf.discard(h)                                        # the unjoined run is killed + cleaned
    assert not cf._PENDING and len(_sscf_dirs()) <= len(before)
