"""Synthetic melanoma section - an instant, reproducible fallback demo sample.

Used by the app's "Load melanoma example" when the real 10x bundle isn't present
(set ``SPATIALSCRIBE_DEMO`` to use the real FFPE Human Skin Melanoma Xenium Prime data).

The spatial layout is deliberately realistic for the demo: a central malignant tumor
nest, T cells pushed to a ring at the margin (an **immune-excluded** phenotype - the
copilot hero moment), TAMs infiltrating the nest, and stroma/endothelium/keratinocytes
in the periphery.
"""

from __future__ import annotations

from . import markers as _m


def make_demo_adata(n_cells: int = 3000, seed: int = 0):
    """Return a Xenium-like AnnData (counts, feature_types, obsm['spatial'])."""
    import anndata as ad
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(seed)

    # Cell-type composition and their signature genes (panel-present subset).
    comp = {
        "Malignant/Melanocyte": 0.40,
        "T cell": 0.15,
        "Myeloid": 0.12,
        "Endothelial": 0.08,
        "Stromal/CAF": 0.12,
        "B/Plasma": 0.05,
        "Keratinocyte": 0.08,
    }
    sig = {k: _m.LINEAGE_MARKERS[k] for k in comp}
    tam = _m.TAM_STATES  # add TAM-state genes to the panel

    genes = sorted({g for v in sig.values() for g in v} | {g for v in tam.values() for g in v})
    noise = [f"BG{i:03d}" for i in range(120)]
    controls = [f"NegControlProbe_{i:05d}" for i in range(10)] + ["BLANK_0001", "BLANK_0002"]
    all_genes = genes + noise + controls
    gidx = {g: i for i, g in enumerate(all_genes)}

    labels = rng.choice(list(comp), size=n_cells, p=list(comp.values()))
    X = rng.poisson(0.2, size=(n_cells, len(all_genes))).astype("float32")
    for ct, gset in sig.items():
        m = labels == ct
        for g in gset:
            X[m, gidx[g]] += rng.poisson(10, m.sum())
    # TAMs (subset of Myeloid) express an IFN-TAM program - gives H6 something to find.
    mye = labels == "Myeloid"
    for g in tam["IFN-TAM"]:
        X[mye, gidx[g]] += rng.poisson(6, mye.sum())

    coords = _immune_excluded_layout(labels, rng)

    var = pd.DataFrame(index=all_genes)
    var["feature_types"] = (
        ["Gene Expression"] * (len(all_genes) - len(controls))
        + ["Negative Control Probe"] * len(controls)
    )
    a = ad.AnnData(X=X, var=var)
    a.obs_names = [f"cell_{i:05d}" for i in range(n_cells)]
    a.obs["true_type"] = pd.Categorical(labels)
    a.obsm["spatial"] = coords
    return a


def _immune_excluded_layout(labels, rng):
    """Place a central tumor nest with T cells excluded to the margin ring."""
    import numpy as np

    n = len(labels)
    coords = np.zeros((n, 2))
    cx, cy, r_nest, r_margin = 0.0, 0.0, 300.0, 380.0
    for i, ct in enumerate(labels):
        if ct == "Malignant/Melanocyte":
            r = r_nest * np.sqrt(rng.uniform(0, 1))          # filled tumor nest
        elif ct == "Myeloid":
            r = rng.uniform(0, r_margin)                      # infiltrating
        elif ct == "T cell":
            r = rng.normal(r_margin, 25)                      # excluded ring at the margin
        else:
            r = rng.uniform(r_margin + 20, 700)               # peripheral stroma/vasc/epi
        theta = rng.uniform(0, 2 * np.pi)
        coords[i] = [cx + r * np.cos(theta), cy + r * np.sin(theta)]
    return coords


def load_demo():
    """Return a :class:`~spatialscribe.analysis.io.SpatialSample` for the synthetic demo."""
    from .io import build_control_mask, SpatialSample

    a = make_demo_adata()
    mask = build_control_mask(a)
    a.var["control"] = mask
    a.uns["platform"] = "xenium"          # a targeted-panel Xenium section, for downstream QC
    a.uns["n_panel_genes"] = int((~mask).sum())
    return SpatialSample(
        platform="xenium-demo",
        adata=a,
        control_mask=mask,
        panel_genes=a.var_names[~mask].tolist(),
        has_z=False,
        transcripts_path=None,
        meta={"source": "synthetic-melanoma-demo"},
    )
