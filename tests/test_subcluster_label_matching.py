"""subcluster must tolerate the label the copilot actually passes.

Regression: `subcluster.subcluster` matched cell_type EXACTLY, so when the copilot echoed the user's
wording ("T cells") instead of the exact category ("T cell") it raised "Too few (0)", and the
multi-turn copilot loop silently recovered into a different tool (assign_cell_states) - the user asked
to subcluster T cells and got cell states instead. This pins the tolerant resolver + helpful error.
"""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from spatialscribe.analysis import subcluster as sc

CATS = ["T cell", "B cell", "Myeloid", "Mast", "Epithelial/Tumor"]


@pytest.mark.parametrize("requested,expected", [
    ("T cell", "T cell"),        # exact
    ("T cells", "T cell"),       # plural (the actual failure)
    ("t cell", "T cell"),        # case
    ("T-cell", "T cell"),        # hyphen
    ("myeloid", "Myeloid"),      # case
    ("Epithelial/Tumor", "Epithelial/Tumor"),
])
def test_match_resolves_near_misses(requested, expected):
    assert sc._match_cell_type(requested, CATS) == expected


@pytest.mark.parametrize("requested", ["Neuron", "", "cell", "xyzzy"])
def test_match_returns_none_when_absent_or_ambiguous(requested):
    # "cell" is a substring of several categories -> ambiguous -> None (not a wrong lucky pick)
    assert sc._match_cell_type(requested, CATS) is None


def test_subcluster_raises_helpful_error_listing_types():
    a = ad.AnnData(X=np.zeros((40, 3), dtype="float32"))
    a.obs["cell_type"] = pd.Categorical(["T cell"] * 20 + ["B cell"] * 20)
    with pytest.raises(ValueError) as e:
        sc.subcluster(a, "Neuron", use_llm=False)
    msg = str(e.value)
    assert "not a cell type" in msg and "T cell" in msg and "B cell" in msg  # tells the LLM what IS valid


def test_subcluster_emits_map_view_recolor(processed_adata, ctx):
    """The copilot's subcluster should ALSO recolour the map by the new subtypes (a map_view directive,
    like the view tools), not just answer in chat - so 'subcluster the T cells' switches the canvas to
    obs['subtype'] instead of leaving it on cell_type."""
    from spatialscribe.analysis import capabilities as cap
    top = str(processed_adata.obs["cell_type"].value_counts().index[0])
    res = cap.run(processed_adata, "subcluster", {"cell_type": top}, ctx)
    assert res.ok, res.error
    assert "subtype" in processed_adata.obs
    assert any(a.get("kind") == "map_view" and a.get("color_by") == "subtype" for a in ctx.artifacts)
