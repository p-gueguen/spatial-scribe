"""Plotting - the interactive spatial canvas (plotly) + static report figures.

`spatial_figure` returns a WebGL scatter (handles 10^5-10^6 cells) plus the index array
mapping plotted points back to ``adata`` positions - so a box-select selection's
``point_indices`` resolve to real cells (see docs: stable-index requirement). Shared by
the Streamlit app and the HTML report.

The canvas is styled as a dark-field microscopy view (the native look of imaging-based
spatial transcriptomics): a near-black field where cells glow as bright fluorescent points.
"""

from __future__ import annotations

# A qualitative palette tuned to glow on the near-black (#070a10) canvas - soft fluorescence hues,
# NOT saturated primaries. The first 28 are the original Tailwind-400/300 set (restored: the numerically
# packed 41-colour luminance-grid that briefly replaced it maximised CIELAB separation but read as
# garish - pure #f80700 red, acid #6eff24 green, hot #f9118c magenta). The last 12 are additional soft
# pastels for sections with >28 categories (fine cell types / many niches), greedily chosen to sit as
# far as possible from the first 28 in CIELAB so extending the list does NOT reduce distinctness.
# Verified properties (see tests/test_palette.py, recomputed from the hex strings):
#   min relative luminance >= 0.18 -> every entry reads on #070a10 (measured 0.33)
#   min pairwise CIE76     >= 8    -> some light-blues are close BY DESIGN (aesthetic continuity was
#     chosen over maximal separation; the readability floor is what actually matters on black). Measured 8.5.
# Deterministic literal so a cluster keeps its colour across runs. `_cat_colour` in the backend
# extends past this length with golden-ratio HSV; do not shorten below 40 or that fallback kicks in early.
_PALETTE = [
    "#22d3ee", "#e879f9", "#34d399", "#fbbf24", "#f87171", "#a78bfa", "#4ade80",
    "#f472b6", "#38bdf8", "#facc15", "#fb923c", "#2dd4bf", "#c084fc", "#60a5fa",
    "#f9a8d4", "#86efac", "#fde047", "#fca5a5", "#5eead4", "#d8b4fe", "#93c5fd",
    "#fdba74", "#bef264", "#67e8f9", "#f0abfc", "#a5f3fc", "#fef08a", "#fbcfe8",
    "#fbf8cc", "#adc178", "#94d2bd", "#caffbf", "#ffd6a5", "#c3aed6", "#d4a373",
    "#b5e48c", "#ffc6ff", "#7dd3fc", "#98f5e1", "#ffd8be",
]

# Verdict colors for the confidence overlay (glow on black).
VERDICT_COLORS = {"PASS": "#34d399", "WARN": "#fbbf24", "FAIL": "#fb7185"}

PLOT_BG = "#070a10"
_FONT = "JetBrains Mono, ui-monospace, SF Mono, monospace"


def _apply_dark(fig, hide_axes: bool = False):
    """Give a plotly figure the dark-field viewport look (near-black, light mono type)."""
    fig.update_layout(
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor=PLOT_BG,
        font=dict(color="#c7d0dc", family=_FONT, size=11),
        legend=dict(bgcolor="rgba(10,14,20,.55)", bordercolor="#1d2531", borderwidth=1,
                    font=dict(size=10, color="#c7d0dc")),
    )
    if hide_axes:
        opt = dict(showgrid=False, zeroline=False, showticklabels=False, showline=False, ticks="")
        fig.update_xaxes(**opt)
        fig.update_yaxes(**opt)
    else:
        grid = dict(showgrid=True, gridcolor="rgba(120,140,170,.09)", zeroline=False,
                    linecolor="#1d2531", tickfont=dict(color="#7d8a9c", size=9))
        fig.update_xaxes(**grid)
        fig.update_yaxes(**grid)
    return fig


def _downsample(n: int, max_points: int, seed: int = 0):
    import numpy as np

    if n <= max_points:
        return np.arange(n)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n, max_points, replace=False))


def _ordered_group_labels(col, present) -> list[str]:
    """Legend/trace order for a discrete obs column, returned as strings.

    Honors a pandas ``Categorical``'s declared category order - so integer cluster labels
    read 0, 1, 2, ..., 10 rather than the lexical 0, 1, 10, 2 that a plain
    ``sorted(set(...))`` produces (see ``cluster.order_clusters_numeric``) - intersected
    with the labels actually ``present`` in the plotted subset. Any present label outside
    the declared categories is appended lexically so nothing is dropped. Plain
    (non-categorical) columns fall back to a lexical sort, the prior behavior.
    """
    present = {str(v) for v in present}
    if hasattr(col, "cat"):                       # a categorical pandas Series (the real call sites)
        cats = col.cat.categories
    else:                                         # a raw pd.Categorical, else None -> lexical fallback
        cats = getattr(col, "categories", None)
    if cats is None:
        return sorted(present)
    ordered = [str(c) for c in cats if str(c) in present]
    return ordered + sorted(present - set(ordered))


def spatial_figure(adata, color_by: str, max_points: int = 400_000, point_size: float | None = None,
                   height: int = 780):
    """Interactive dark-field spatial scatter colored by an obs column. Returns (fig, idx).

    ``idx`` maps each plotted point to its position in ``adata`` (for selection mapping).
    Point size auto-scales to the cell density so cells read as a dense tissue without heavy
    overlap (a coarse proxy for the mean nearest-neighbour spacing of the section).
    """
    import numpy as np
    import plotly.graph_objects as go

    xy = adata.obsm["spatial"]
    idx = _downsample(adata.n_obs, max_points)
    x, y = xy[idx, 0], xy[idx, 1]
    if point_size is None:
        # Target the on-screen nearest-neighbour spacing so cells read as a dense tissue but do
        # not overlap: ~ plot_pixels / sqrt(n), scaled down for the section's partial fill.
        point_size = float(np.clip(height / (max(1, len(idx)) ** 0.5) * 0.42, 0.5, 2.2))
    fig = go.Figure()

    if color_by in adata.obs and str(adata.obs[color_by].dtype) in ("category", "object"):
        vals = adata.obs[color_by].astype(str).to_numpy()[idx]
        cats = _ordered_group_labels(adata.obs[color_by], vals)
        for i, c in enumerate(cats):
            m = vals == c
            fig.add_trace(go.Scattergl(
                x=x[m], y=y[m], mode="markers", name=str(c),
                marker=dict(size=point_size, color=_PALETTE[i % len(_PALETTE)], opacity=0.9,
                            line=dict(width=0)),
                customdata=idx[m], hovertext=[str(c)] * int(m.sum()), hoverinfo="text+name",
            ))
    else:
        vals = np.asarray(adata.obs[color_by])[idx] if color_by in adata.obs else np.zeros(len(idx))
        fig.add_trace(go.Scattergl(
            x=x, y=y, mode="markers",
            marker=dict(size=point_size, color=vals, colorscale="Turbo", showscale=True, opacity=0.9,
                        colorbar=dict(title=dict(text=color_by, font=dict(color="#c7d0dc")),
                                      tickfont=dict(color="#7d8a9c"), outlinewidth=0, thickness=12)),
            customdata=idx,
        ))
    fig.update_layout(
        height=height, dragmode="box-select", showlegend=True,
        margin=dict(l=0, r=0, t=0, b=0), yaxis=dict(scaleanchor="x", scaleratio=1),
        legend=dict(itemsizing="constant", x=0.997, y=0.997, xanchor="right", yanchor="top",
                    bgcolor="rgba(8,11,16,.66)", bordercolor="#1d2531", borderwidth=1,
                    font=dict(size=10)),
    )
    _apply_dark(fig, hide_axes=True)
    fig.update_yaxes(autorange="reversed")   # image convention
    return fig, idx


def umap_figure(adata, color_by: str):
    """UMAP scatter colored by an obs column (plotly), dark-field styled."""
    import plotly.express as px

    if "X_umap" not in adata.obsm:
        return None
    um = adata.obsm["X_umap"]
    df = {"UMAP1": um[:, 0], "UMAP2": um[:, 1], color_by: adata.obs[color_by].astype(str).to_numpy()}
    fig = px.scatter(df, x="UMAP1", y="UMAP2", color=color_by,
                     color_discrete_sequence=_PALETTE, render_mode="webgl")
    fig.update_traces(marker=dict(size=3, opacity=0.9, line=dict(width=0)))
    fig.update_layout(height=500, margin=dict(l=0, r=0, t=10, b=0))
    # UMAP coordinates are not meaningful - hide the axis ticks/values (and titles) for the
    # clean dark-field look, matching spatial_figure.
    _apply_dark(fig, hide_axes=True)
    fig.update_layout(xaxis_title=None, yaxis_title=None)
    return fig
