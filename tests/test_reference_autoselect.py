"""auto_select_reference: free-text tissue -> the best-matched reference is CHOSEN and LOADED
(the 'auto' half choose_reference was missing), and the live gget fetch degrades cleanly.

Network-free: the registry points at tiny on-disk .h5ad fixtures, and the fetch path is exercised
by monkeypatching the runner / fetch so no CELLxGENE call is made.
"""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

from spatialscribe.analysis import reference as ref


def _tiny_ref(path, types=("T_cell", "B_cell", "Macrophage"), genes=("CD3D", "CD19", "CD68", "EPCAM")):
    rng = np.random.default_rng(0)
    n = 60
    X = rng.poisson(2.0, size=(n, len(genes))).astype("float32")
    a = ad.AnnData(X)
    a.var_names = list(genes)
    a.obs["cell_type"] = pd.Categorical([types[i % len(types)] for i in range(n)])
    a.write_h5ad(path)
    return str(path)


def _registry(tmp_path, available=("lung",)):
    """A synthetic registry; only the tissues in `available` get a real on-disk path."""
    reg = {
        "lung": {"path": "", "label_key": "cell_type", "gene_name_col": None,
                 "description": "lung reference", "keywords": ["lung", "pulmonary", "nsclc"]},
        "kidney": {"path": "", "label_key": "cell_type", "gene_name_col": None,
                   "description": "kidney reference", "keywords": ["kidney", "renal", "nephron"]},
    }
    for key in available:
        reg[key]["path"] = _tiny_ref(tmp_path / f"{key}.h5ad")
    return reg


def test_auto_selects_and_loads_matching_reference(tmp_path):
    reg = _registry(tmp_path, available=("lung", "kidney"))
    sel = ref.auto_select_reference("lung adenocarcinoma", registry=reg)
    assert sel["status"] == "registry"
    assert sel["chosen"] == "lung"
    assert sel["ref"] is not None and sel["label_key"] == "cell_type"
    assert sel["ref"].n_obs == 60


def test_skips_wrong_tissue_reference(tmp_path):
    # Only a LUNG reference is on disk, but the query is KIDNEY -> must NOT auto-load the lung atlas
    # (wrong-tissue transfer annotates nonsense); it falls through to no_reference / cluster.
    reg = _registry(tmp_path, available=("lung",))
    sel = ref.auto_select_reference("kidney", registry=reg)
    assert sel["status"] == "no_reference"
    assert sel["ref"] is None


def test_coverage_tiebreak_picks_organ_specific_over_generic(tmp_path):
    """Among references that TIE on keyword score (e.g. two melanoma-organ atlases both hitting 1.0 on
    the 'melanoma' token), the one matching MORE of the query tokens (organ AND disease) wins - even
    when it is listed later. Guards the RCTD/UM registry, where a 'uveal melanoma liver metastasis'
    section must select the liver reference, not the (first-listed) eye reference."""
    # um_eye listed FIRST; both match uveal+melanoma, only um_liver also matches liver+hepatic.
    reg = {
        "um_eye": {"path": _tiny_ref(tmp_path / "eye.h5ad"), "label_key": "cell_type",
                   "gene_name_col": None, "description": "UM eye/retina reference (Melanoma added)",
                   "keywords": ["eye", "retina", "ocular", "choroid", "uveal", "melanoma", "um"]},
        "um_liver": {"path": _tiny_ref(tmp_path / "liver.h5ad"), "label_key": "cell_type",
                     "gene_name_col": None, "description": "UM liver reference (Melanoma added)",
                     "keywords": ["liver", "hepatic", "hepatocyte", "metastasis", "uveal", "melanoma", "um"]},
    }
    ranked = ref.choose_reference("uveal melanoma liver metastasis", registry=reg, top_n=2)
    assert ranked[0]["tissue_key"] == "um_liver"                 # more query tokens covered
    assert ranked[0]["match_coverage"] > ranked[1]["match_coverage"]
    assert ranked[0]["keyword_score"] == ranked[1]["keyword_score"] == 1.0   # tie on primary score
    # ...and it actually loads that one
    sel = ref.auto_select_reference("uveal melanoma liver metastasis", registry=reg)
    assert sel["status"] == "registry" and sel["chosen"] == "um_liver"

    # A section in the SAME organ but NO disease token still resolves to that organ's reference.
    assert ref.choose_reference("liver", registry=reg, top_n=1)[0]["tissue_key"] == "um_liver"


def test_fetch_path_when_no_local_reference(tmp_path, monkeypatch):
    reg = _registry(tmp_path, available=())         # nothing available locally
    fetched_h5ad = _tiny_ref(tmp_path / "cellxgene.h5ad")
    monkeypatch.setattr(ref, "fetch_cellxgene_reference",
                        lambda *a, **k: {"status": "ok", "path": fetched_h5ad, "label_key": "cell_type"})
    sel = ref.auto_select_reference("lung", registry=reg, allow_fetch=True)
    assert sel["status"] == "fetched"
    assert sel["ref"] is not None and sel["chosen"] == "cellxgene"


def test_fetch_degrades_without_gget_runner(monkeypatch):
    # No gget env and no uv -> the live fetch returns a clean 'skipped' status, never raises.
    monkeypatch.setattr(ref, "_resolve_gget_runner", lambda: None)
    out = ref.fetch_cellxgene_reference("skin of body")
    assert out["status"] == "skipped"
    assert "gget" in out["message"].lower()


def test_no_fetch_without_local_returns_no_reference(tmp_path):
    reg = _registry(tmp_path, available=())
    sel = ref.auto_select_reference("lung", registry=reg, allow_fetch=False)
    assert sel["status"] == "no_reference"
    assert sel["ref"] is None
    assert "clustering" in sel["message"].lower()
