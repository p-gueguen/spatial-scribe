"""Copilot-driven spatial VIEWS: resolve a plain-language cell criterion to a mask, and render
a spatial plot that highlights those cells.

This is what lets the copilot *act on the canvas* rather than only return numbers: a request like
"highlight low-quality cells" resolves to a boolean mask (``select_cells``) and a dark-field spatial
figure with the matching cells lit up over a greyed background (``highlight_figure``). The capability
``highlight_cells`` wires these into the copilot loop and stashes the figure on ``RunContext.artifacts``
so the app renders it.

Criteria understood (case-insensitive, fuzzy): a cell-type name; ``low quality``; ``low confidence`` /
``uncertain`` / ``abstained``; ``malignant`` / ``tumor``. Honest by construction - e.g. "low quality"
when the QC layers flagged nothing falls back to the lowest-count decile with an explicit note.
"""
from __future__ import annotations


def select_cells(adata, criterion: str, ctx=None):
    """Resolve ``criterion`` to ``(mask, resolved_label, note)``.

    ``mask`` is a boolean array over ``adata`` cells; ``resolved_label`` is the concrete rule applied;
    ``note`` is a one-line plain-language explanation (states honest fallbacks). Never raises - an
    unresolvable criterion returns an all-False mask + an explanatory note.
    """
    import numpy as np

    c = str(criterion).strip().lower()
    n = adata.n_obs
    obs = adata.obs
    tissue = getattr(ctx, "tissue", "melanoma") if ctx is not None else "melanoma"

    cats = {}
    if "cell_type" in obs:
        cats = {str(x).lower(): str(x) for x in obs["cell_type"].astype(str).unique()}
        if c in cats:
            m = (obs["cell_type"].astype(str) == cats[c]).to_numpy()
            return m, cats[c], f"{int(m.sum()):,} cells annotated {cats[c]}."

    # malignant / tumor
    if any(k in c for k in ("malignant", "tumor", "tumour")):
        if "malignant_score" not in obs:
            from . import cnv
            cnv.malignant_score(adata, tissue=tissue)
        s = np.nan_to_num(np.asarray(obs.get("malignant_score", np.zeros(n)), dtype=float))
        m = s > 0.6
        return m, "high malignant score (>0.6)", f"{int(m.sum()):,} cells with malignant score > 0.6."

    # low confidence / uncertain / abstained
    if any(k in c for k in ("confiden", "uncertain", "tentative", "abstain", "unassigned", "reject", "ambiguous")):
        if "annotation_verdict" in obs:
            v = obs["annotation_verdict"].astype(str)
            if any(k in c for k in ("abstain", "unassigned", "reject", "ambiguous")):
                m = (v == "FAIL").to_numpy()
                return m, "abstained (verdict FAIL)", f"{int(m.sum()):,} cells the annotator abstained on (verdict FAIL)."
            m = v.isin(["WARN", "FAIL"]).to_numpy()
            return m, "low-confidence (WARN or FAIL)", f"{int(m.sum()):,} low-confidence cells (WARN/FAIL)."
        if "annotation_confidence" in obs:
            conf = np.nan_to_num(np.asarray(obs["annotation_confidence"], dtype=float))
            m = conf < 0.5
            return m, "confidence < 0.5", f"{int(m.sum()):,} cells with annotation confidence < 0.5."
        return np.zeros(n, bool), criterion, "No annotation confidence yet - run Annotate first, then ask again."

    # low quality (QC): under the panel-indexed count floor, or segmentation/annotation-quality flags
    if any(k in c for k in ("quality", "low_qual", "lowqual", "poor", "bad", "junk")):
        from . import qc as _qc
        if "total_counts" not in obs:
            _qc.compute_qc(adata)
        tc = np.nan_to_num(np.asarray(adata.obs["total_counts"], dtype=float))
        floor = float(_qc.suggest_count_floor(adata).get("floor", 10) or 10)
        flagged = tc < floor
        if "seg_area_flag" in obs:
            flagged = flagged | (obs["seg_area_flag"].astype(str) != "ok").to_numpy()
        if "annotation_reason" in obs:
            flagged = flagged | obs["annotation_reason"].astype(str).str.contains("low_quality", case=False).to_numpy()
        if flagged.sum() > 0:
            return flagged, f"low quality (< count floor {floor:.0f} or flagged)", \
                f"{int(flagged.sum()):,} cells flagged low-quality (under the count floor or a QC flag)."
        # honest fallback: nothing was flagged -> show the lowest-count decile as the nearest population
        thr = float(np.quantile(tc, 0.10))
        m = tc <= thr
        return m, f"lowest-count decile (<= {thr:.0f} counts)", (
            f"No cells were FLAGGED as low-quality by the QC layers that ran (count floor {floor:.0f} removed 0). "
            f"Highlighting the lowest-count 10% ({int(m.sum()):,} cells, <= {thr:.0f} counts) as the "
            f"nearest-to-low-quality population - not a failing set.")

    # fuzzy cell-type substring match as a last resort
    for lc, orig in cats.items():
        if c and (c in lc or lc in c):
            m = (obs["cell_type"].astype(str) == orig).to_numpy()
            return m, orig, f"{int(m.sum()):,} cells annotated {orig} (matched '{criterion}')."

    return np.zeros(n, bool), criterion, (
        f"Could not resolve '{criterion}'. Try a cell type, 'low quality', 'low confidence', or 'malignant'.")


def highlight_figure(adata, mask, title: str, height: int = 620):
    """Dark-field spatial scatter with ``mask`` cells lit (magenta) over a greyed background."""
    import numpy as np
    import plotly.graph_objects as go

    from .plots import _apply_dark, _downsample

    xy = adata.obsm["spatial"]
    idx = _downsample(adata.n_obs, 400_000)
    x, y = xy[idx, 0], xy[idx, 1]
    m = np.asarray(mask, dtype=bool)[idx]
    ps = float(np.clip(height / (max(1, len(idx)) ** 0.5) * 0.42, 0.6, 2.4))
    fig = go.Figure()
    fig.add_trace(go.Scattergl(x=x[~m], y=y[~m], mode="markers", name="other",
                               marker=dict(size=ps, color="#33404f", opacity=0.5, line=dict(width=0))))
    fig.add_trace(go.Scattergl(x=x[m], y=y[m], mode="markers", name=f"match ({int(m.sum()):,})",
                               marker=dict(size=ps * 1.7, color="#e879f9", opacity=0.95, line=dict(width=0))))
    fig.update_layout(height=height, showlegend=True, margin=dict(l=0, r=0, t=26, b=0),
                      title=dict(text=title, font=dict(color="#e8eef7", size=13), x=0.01, y=0.99),
                      yaxis=dict(scaleanchor="x", scaleratio=1),
                      legend=dict(itemsizing="constant", x=0.997, y=0.997, xanchor="right", yanchor="top",
                                  bgcolor="rgba(8,11,16,.66)", bordercolor="#1d2531", borderwidth=1,
                                  font=dict(size=10)))
    _apply_dark(fig, hide_axes=True)
    fig.update_yaxes(autorange="reversed")
    return fig


# --------------------------------------------------------------------------- #
# Field/gene views the copilot can request ("show me CD8A", "dotplot of T markers", ...)
# --------------------------------------------------------------------------- #
_FIELD_ALIAS = {
    "qc": "total_counts", "counts": "total_counts", "depth": "total_counts", "total counts": "total_counts",
    "genes": "n_genes", "n genes": "n_genes", "malignant": "malignant_score", "tumor": "malignant_score",
    "tumour": "malignant_score", "malignant score": "malignant_score", "cell type": "cell_type",
    "celltype": "cell_type", "type": "cell_type", "annotation": "cell_type", "niche": "niche",
    "cluster": "leiden", "leiden": "leiden", "confidence": "annotation_confidence",
}


def _gene_vector(adata, gene, layer=None):
    import numpy as np
    sub = adata[:, gene]
    x = sub.layers[layer] if (layer and layer in adata.layers) else sub.X
    return np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()


def resolve_obs_field(adata, color_by):
    """The obs column ``color_by`` resolves to (via the field aliases), or ``None`` if it is a gene /
    unknown. Lets the copilot's ``show_spatial`` drive the MAIN specimen canvas only for real obs
    fields (a gene stays an inline figure since the big canvas colours by obs columns)."""
    kl = str(color_by).strip().lower()
    resolved = _FIELD_ALIAS.get(kl, str(color_by).strip())
    return resolved if resolved in adata.obs else None


def color_field_for(adata, color_by, ctx=None):
    """Resolve a copilot colour request to a PERSISTENT obs field the MAIN canvas can colour by.

    An obs field (via the field aliases) is returned as-is; an on-panel GENE has its expression
    written to ``obs[<gene>]`` (recorded in ``uns['gene_color_fields']`` so the app lists it in the
    colour dropdown) and its name returned; anything else returns ``(None, note)``. This lets
    'colour the map by CD3E' drive the MAIN spatial plot rather than dump an inline figure. Same
    field aliases + on-demand materialisation (malignant_score / QC) as :func:`spatial_view`.
    Returns ``(field|None, note)``.
    """
    key = str(color_by).strip()
    kl = key.lower()
    resolved = _FIELD_ALIAS.get(kl, key)
    if resolved == "malignant_score" and "malignant_score" not in adata.obs:
        from . import cnv
        cnv.malignant_score(adata, tissue=getattr(ctx, "tissue", "melanoma") if ctx else "melanoma")
    if resolved == "n_genes" and "n_genes" not in adata.obs and "total_counts" not in adata.obs:
        from . import qc as _qc
        _qc.compute_qc(adata)
    if resolved in adata.obs:
        return resolved, f"Coloured the map by obs['{resolved}']."
    var = {str(v).lower(): str(v) for v in adata.var_names}
    if kl in var:
        g = var[kl]
        adata.obs[g] = _gene_vector(adata, g)                 # persist so /api/points can colour by it
        gcols = list(adata.uns.get("gene_color_fields", []))
        if g not in gcols:
            gcols.append(g)
            adata.uns["gene_color_fields"] = gcols
        return g, f"Coloured the map by {g} expression (on-panel)."
    return None, (f"'{color_by}' is neither an obs field nor an on-panel gene. Try a gene symbol, "
                  f"'cell type', 'total counts', 'malignant', or 'niche'.")


def spatial_view(adata, color_by: str, ctx=None):
    """Spatial map colored by an obs field OR a gene's expression. Returns ``(fig|None, resolved, note)``.

    Grounded: only a real obs column or an on-panel gene is accepted (fuzzy field aliases like
    'qc'->total_counts, 'tumor'->malignant_score). An unknown key returns ``(None, key, note)``.
    """
    from .plots import spatial_figure

    key = str(color_by).strip()
    kl = key.lower()
    resolved = _FIELD_ALIAS.get(kl, key)

    if resolved == "malignant_score" and "malignant_score" not in adata.obs:
        from . import cnv
        cnv.malignant_score(adata, tissue=getattr(ctx, "tissue", "melanoma") if ctx else "melanoma")
    if resolved == "n_genes" and "n_genes" not in adata.obs and "total_counts" not in adata.obs:
        from . import qc as _qc
        _qc.compute_qc(adata)

    if resolved in adata.obs:
        fig, _ = spatial_figure(adata, resolved)
        return fig, f"spatial · {resolved}", f"Spatial map colored by obs['{resolved}']."

    var = {str(v).lower(): str(v) for v in adata.var_names}
    if kl in var:
        g = var[kl]
        adata.obs["_ss_expr"] = _gene_vector(adata, g)
        fig, _ = spatial_figure(adata, "_ss_expr")
        del adata.obs["_ss_expr"]
        return fig, f"spatial · {g}", f"Spatial map colored by {g} expression (on-panel)."
    return None, key, (f"'{color_by}' is neither an obs field nor an on-panel gene. Try a gene symbol, "
                       f"'cell type', 'total counts', 'malignant', or 'niche'.")


# Abstention labels (annotate.apply_confidence) are not cell types, so they are never keyword-matched:
# "Unknown" contains "nk" and would otherwise be filed under T/NK lymphoid.
_ABSTENTION_TOKENS = ("unassigned", "ambiguous", "unresolvable", "uncertain", "novel", "unknown")


# Deterministic keyword -> category map for the OFFLINE marker-overview grouping (no LLM). First
# category whose ANY keyword matches the (lower-cased) cell-type name wins; unmatched types fall to
# "Other". Total by construction (every type lands in exactly one category). Order here IS the row
# order of the offline dot-plot, so keep biologically-adjacent lineages adjacent.
_FALLBACK_CATEGORIES: dict[str, list[str]] = {
    "T/NK lymphoid": ["t cell", "t-cell", "cd4", "cd8", "treg", "tcm", "trm", "nk", "nkt", "ilc",
                      "lymphoid", "cytotox", "exhaust", "gd t", "gamma"],
    "B/plasma": ["b cell", "b-cell", "b/plasma", "plasma", "mzb", "germinal", "bcell"],
    "Myeloid": ["myeloid", "macrophage", "mono", "dendritic", " dc", "dc ", "tam", "microglia",
                "mast", "neutrophil", "granulocyte", "kupffer", "osteoclast", "langerhans"],
    "Epithelial/tumour": ["epithel", "tumor", "tumour", "malignant", "melanocyte", "melanoma",
                          "carcinoma", "luminal", "basal", "keratinocyte", "ductal", "acinar",
                          "hepatocyte", "alveolar", "club", "goblet", "secretory"],
    "Stromal": ["stroma", "caf", "fibroblast", "myofibro", "pericyte", "mural", "adipocyte",
                "chondro", "osteo", "smooth muscle", "myoepithel", "mesenchym"],
    "Vascular": ["endothel", "vascular", "lymphatic", "vein", "artery", "capillary"],
}


def _kw_hit(lc: str, kw: str) -> bool:
    """True when ``kw`` occurs in ``lc`` at the start of a token.

    The keywords are deliberately PREFIXES ("epithel" must match "epithelial"), so only the LEFT edge
    is anchored, and only for keywords that begin with an alphanumeric - the space-prefixed ones
    (" dc") carry their own boundary and would never match if anchored. Without the left anchor "nk"
    matches "u-nk-nown" and files an untyped cell under T/NK lymphoid.
    """
    import re

    if not kw[:1].isalnum():
        return kw in lc
    return re.search(rf"(?<![a-z0-9]){re.escape(kw)}", lc) is not None


def _fallback_group(cell_types) -> dict[str, list[str]]:
    """Deterministic offline grouping of cell-type NAMES into ``_FALLBACK_CATEGORIES`` (+ "Other").

    Total: every input lands in exactly one category, in ``_FALLBACK_CATEGORIES`` order (with "Other"
    last). Pure keyword matching, so it is stable and testable with no LLM. Abstention labels go
    straight to "Other" - they name the absence of a cell type, not a lineage to group."""
    groups: dict[str, list[str]] = {c: [] for c in _FALLBACK_CATEGORIES}
    other: list[str] = []
    for ct in cell_types:
        lc = str(ct).lower()
        placed = False
        if not lc.startswith(_ABSTENTION_TOKENS):
            for cat, kws in _FALLBACK_CATEGORIES.items():
                if any(_kw_hit(lc, k) for k in kws):
                    groups[cat].append(str(ct))
                    placed = True
                    break
        if not placed:
            other.append(str(ct))
    out = {c: v for c, v in groups.items() if v}
    if other:
        out["Other"] = other
    return out


def _category_groups(adata, ctx):
    """(groups, source): group the section's cell types into biologically related categories.

    ``ctx.use_llm`` -> ``llm.group_cell_types`` (falling back to the keyword grouping if the LLM is
    unavailable / returns nothing); else the deterministic keyword grouping. ``source`` is ``"LLM"`` or
    ``"keyword fallback"`` so the note can say which produced the categories."""
    types = list(map(str, adata.obs["cell_type"].astype("category").cat.categories)) \
        if "cell_type" in adata.obs else []
    if getattr(ctx, "use_llm", False):
        from . import llm as _llm
        g = _llm.group_cell_types(types, tissue=getattr(ctx, "tissue", ""))
        if g:
            return g, "LLM"
    return _fallback_group(types), "keyword fallback"


def category_overview(adata, ctx):
    """Marker dot-plot overview: rows = the section's cell types grouped so related types are adjacent,
    columns = the top-2 markers PER CATEGORY, ranked by expression (never by the LLM).

    Returns ``(genes, group_order, note)``:
      * ``genes`` - concatenation, in category order, of each category's top-2 discriminating markers.
        Candidates are the union of the category's member types' curated on-panel markers; they are
        ranked by one-vs-rest mean-expression difference (in-category cells vs the rest) COMPUTED from
        the section, and de-duplicated across categories so no column repeats.
      * ``group_order`` - the flat list of cell types in category order (the dot-plot rows).
      * ``note`` - kept quiet (the redundant "rows grouped by ..." caption was dropped as UI noise);
        the grouping SOURCE is recorded structurally in ``uns['marker_overview_source']`` instead.
    Also stashes ``adata.uns['marker_overview_groups'] = {category: [types]}`` so the capability wrapper
    can report the grouping. Degrades to ``([], [], note)`` when there are no cell types / markers.
    """
    import numpy as np

    groups, source = _category_groups(adata, ctx)
    adata.uns["marker_overview_groups"] = {c: list(t) for c, t in groups.items()}
    adata.uns["marker_overview_source"] = source     # 'LLM' | 'keyword fallback' (offline-honest, structural)
    panel = {str(v).lower(): str(v) for v in adata.var_names}
    marks = ctx.markers(adata) if ctx is not None else {}
    ctcol = adata.obs["cell_type"].astype(str).to_numpy() if "cell_type" in adata.obs else np.array([])

    genes_out: list[str] = []
    order_out: list[str] = []
    used: set[str] = set()
    # Top-2 markers PER CELL TYPE, not per category: a coarse type would otherwise dominate its
    # category's top-2 and leave its siblings with NO marker column (Basal vs Epithelial/Tumor,
    # NK/B/Plasma/Mast vs T cell). The category grouping is kept only for ROW ORDER (related types
    # adjacent). Each type's own on-panel markers are ranked by how much IT over-expresses them vs
    # the rest of the section; the top-2 new (not-yet-used) ones become its columns.
    for cat, types in groups.items():
        for t in types:
            order_out.append(t)
            if not ctcol.size:
                continue
            cand: list[str] = []
            for g in marks.get(t, []):
                pg = panel.get(str(g).lower())
                if pg and pg not in cand:
                    cand.append(pg)
            if not cand:
                continue
            intype = ctcol == t
            scored = []
            for g in cand:
                v = _gene_vector(adata, g)
                diff = (float(v[intype].mean() - v[~intype].mean())
                        if intype.any() and (~intype).any() else 0.0)
                scored.append((diff, g))
            scored.sort(key=lambda x: -x[0])
            added = 0                               # at most the top-2 NEW markers for this type
            for _, g in scored:
                if g in used:
                    continue
                genes_out.append(g)
                used.add(g)
                added += 1
                if added >= 2:
                    break
    note = ""   # grouping still applies; the "rows grouped ... (LLM)" caption was dropped as UI noise
    return genes_out, order_out, note


def _resolve_genes(adata, genes, cell_type, ctx):
    """Resolve a gene list: explicit ``genes`` (on-panel only) else the ``cell_type``'s canonical
    markers on the panel, else the grouped marker overview. Returns ``(genes_on_panel, note)``.

    The no-argument default delegates to :func:`category_overview` (dropping its row order); the one
    caller that needs the row order (``_cap_marker_dotplot``) calls ``category_overview`` directly."""
    panel = {str(v).lower(): str(v) for v in adata.var_names}
    if genes:
        gl = [panel[g.lower()] for g in genes if g.lower() in panel]
        miss = [g for g in genes if g.lower() not in panel]
        return gl, (f"{len(miss)} requested gene(s) not on panel: {miss[:6]}" if miss else "")
    marks = ctx.markers() if ctx is not None else {}
    if cell_type:
        for ct, gs in marks.items():
            if ct.lower() == str(cell_type).lower():
                gl = [panel[g.lower()] for g in gs if g.lower() in panel]
                return gl, (f"{ct} markers on panel: {gl}" if gl else f"no {ct} markers on the panel")
    genes_out, _order, note = category_overview(adata, ctx)
    return genes_out, note


# Null-like tokens `astype(str)` manufactures from a NaN category. `cluster.cluster` deliberately
# leaves low-signal cells with a NaN `leiden`, so a dot-plot grouped by leiden would otherwise draw a
# phantom "nan" row for them (10,588 cells on the bundled breast demo).
_NULL_GROUPS = ("nan", "none", "na", "<na>")


def dotplot_view(adata, genes, groupby: str = "cell_type", height: int = 460,
                 layer=None, title_suffix: str = "", group_order: list | None = None,
                 group_noun: str = "types", width: int | None = None):
    """Plotly dot-plot: mean expression (color) + % expressing (size) per (group, gene).

    ``layer`` reads expression from ``adata.layers[layer]`` instead of ``X`` (e.g. the
    SPLIT-purified counts, to see marker specificity after spillover removal). ``group_order``, when
    given, draws the rows in that order (biologically related types adjacent): types absent from the
    section are skipped, types present but not listed are appended (in frequency order). ``None`` keeps
    the default frequency ordering. ``group_noun`` names the rows in the title ("types"/"clusters").
    ``width`` is stamped on the figure alongside ``height`` so the PNG renderer can size a wide plot
    instead of squashing it into its default box. Null-like groups (a NaN category) are dropped."""
    import numpy as np
    import plotly.graph_objects as go

    from .plots import _apply_dark

    freq = adata.obs[groupby].astype(str).value_counts()
    eligible = [g for g in freq.index
                if str(g).lower() not in _NULL_GROUPS
                and (adata.obs[groupby].astype(str) == g).sum() >= 5]
    if group_order:
        seen = set(group_order)
        groups = [g for g in group_order if g in eligible] + [g for g in eligible if g not in seen]
    else:
        groups = eligible
    M = np.zeros((len(groups), len(genes))); P = np.zeros_like(M)
    for j, gene in enumerate(genes):
        v = _gene_vector(adata, gene, layer=layer)
        for i, grp in enumerate(groups):
            m = (adata.obs[groupby].astype(str) == grp).to_numpy()
            if m.sum():
                M[i, j] = float(v[m].mean()); P[i, j] = float((v[m] > 0).mean())
    Ms = M / (M.max(0, keepdims=True) + 1e-9)   # scale each gene 0..1 across groups
    # Marker `size` is in PIXELS (zoom-invariant), so a fixed max bubble (was 26 px) bled onto the
    # neighbouring category when rows pack tightly (22 types in ~460 px = ~19 px/row). Derive the max
    # dot from the cell pitch so a 100%-expressing bubble fits inside its own row/column.
    row_pitch = max(1.0, height - 40.0) / max(len(groups), 1)
    col_pitch = (max(1.0, width - 40.0) / max(len(genes), 1)) if width else row_pitch
    max_dot = max(5.0, min(row_pitch, col_pitch, 26.0) * 0.82)
    min_dot = min(3.0, max_dot * 0.4)
    xs, ys, sizes, colors = [], [], [], []
    for i in range(len(groups)):
        for j in range(len(genes)):
            xs.append(j); ys.append(i); sizes.append(min_dot + P[i, j] * (max_dot - min_dot)); colors.append(Ms[i, j])
    fig = go.Figure(go.Scatter(x=xs, y=ys, mode="markers",
                    marker=dict(size=sizes, color=colors, colorscale="Turbo", showscale=True,
                                colorbar=dict(title=dict(text="relative<br>(per gene)", font=dict(color="#c7d0dc")),
                                              tickfont=dict(color="#7d8a9c"), outlinewidth=0, thickness=12),
                                line=dict(width=0))))
    fig.update_layout(height=height, margin=dict(l=0, r=0, t=28, b=0),
                      title=dict(text=f"marker dot-plot{title_suffix} ({len(genes)} genes × {len(groups)} {group_noun}) · size = % expressing",
                                 font=dict(color="#e8eef7", size=12), x=0.01, y=0.99),
                      xaxis=dict(tickmode="array", tickvals=list(range(len(genes))), ticktext=genes, tickangle=-60),
                      # reversed so the FIRST group is at the TOP: rows read top-to-bottom and the
                      # per-type marker diagonal runs from the top-left, not the bottom-left.
                      yaxis=dict(tickmode="array", tickvals=list(range(len(groups))), ticktext=groups, autorange="reversed"))
    if width:
        fig.update_layout(width=int(width))
    _apply_dark(fig)
    return fig


def violin_view(adata, gene: str, groupby: str = "cell_type", height: int = 420):
    """Plotly violin of a gene's expression across ``groupby`` categories. Returns ``(fig|None, note)``."""
    import plotly.graph_objects as go

    from .plots import _PALETTE, _apply_dark

    panel = {str(v).lower(): str(v) for v in adata.var_names}
    if str(gene).lower() not in panel:
        return None, f"'{gene}' is not on the panel."
    g = panel[str(gene).lower()]
    v = _gene_vector(adata, g)
    groups = list(adata.obs[groupby].astype(str).value_counts().index)
    fig = go.Figure()
    for i, grp in enumerate(groups):
        m = (adata.obs[groupby].astype(str) == grp).to_numpy()
        fig.add_trace(go.Violin(y=v[m], name=str(grp), line=dict(color=_PALETTE[i % len(_PALETTE)]),
                                meanline_visible=True, points=False, opacity=0.85))
    fig.update_layout(height=height, showlegend=False, margin=dict(l=0, r=0, t=28, b=0),
                      title=dict(text=f"{g} expression by {groupby}", font=dict(color="#e8eef7", size=12), x=0.01, y=0.99),
                      xaxis=dict(tickangle=-40))
    _apply_dark(fig)
    return fig, f"{g} across {len(groups)} {groupby} groups."


def load_boundaries(path):
    """Load Xenium cell-segmentation polygons. ``path`` = a Xenium run dir (finds
    ``cell_boundaries.parquet``) or a direct parquet. Returns a DataFrame ``[cell_id, vertex_x,
    vertex_y]`` (one row per polygon vertex), or ``None`` if not found."""
    import pathlib
    import pandas as pd

    p = pathlib.Path(path)
    if p.is_dir():
        cand = list(p.glob("cell_boundaries.parquet")) or list(p.glob("*cell_boundaries*.parquet"))
        if not cand:
            return None
        p = cand[0]
    if p.suffix != ".parquet" or not p.exists():
        return None
    df = pd.read_parquet(p, columns=["cell_id", "vertex_x", "vertex_y"])
    df["cell_id"] = df["cell_id"].astype(str)
    return df


def segmentation_figure(adata, boundaries, color_by: str = "cell_type", max_cells: int = 4000, height: int = 620):
    """Draw actual cell-segmentation POLYGONS (not centroids) colored by ``color_by``. Subsamples to
    ``max_cells`` (polygon rendering is heavy). Returns ``(fig, n_shown, note)``."""
    from collections import defaultdict

    import numpy as np
    import plotly.graph_objects as go

    from .plots import _PALETTE, _apply_dark

    ids = adata.obs_names.astype(str)
    have = set(boundaries["cell_id"].unique())
    keep = [c for c in ids if c in have]
    note = ""
    if not keep:
        return None, 0, "The boundary file shares no cell ids with this section."
    if len(keep) > max_cells:
        keep = list(np.random.default_rng(0).choice(np.array(keep), max_cells, replace=False))
        note = f"showing {max_cells:,} of {len(have):,} cells (polygon rendering is heavy)"
    keepset = set(keep)
    grp = adata.obs[color_by].astype(str) if color_by in adata.obs else None
    b = boundaries[boundaries["cell_id"].isin(keepset)]
    gx, gy = defaultdict(list), defaultdict(list)
    for cid, sub in b.groupby("cell_id", sort=False):
        g = str(grp[cid]) if grp is not None and cid in grp.index else "cell"
        xs, ys = sub["vertex_x"].to_numpy(), sub["vertex_y"].to_numpy()
        gx[g] += [*xs.tolist(), float(xs[0]), None]      # close the ring + a gap between cells
        gy[g] += [*ys.tolist(), float(ys[0]), None]
    from .plots import _ordered_group_labels
    groups = _ordered_group_labels(adata.obs[color_by], gx) if color_by in adata.obs else sorted(gx)
    cmap = {g: _PALETTE[i % len(_PALETTE)] for i, g in enumerate(groups)}
    fig = go.Figure()
    for g in groups:
        fig.add_trace(go.Scatter(x=gx[g], y=gy[g], mode="lines", name=str(g), fill="toself",
                                 fillcolor=cmap[g], opacity=0.4, line=dict(color=cmap[g], width=1)))
    fig.update_layout(height=height, showlegend=True, margin=dict(l=0, r=0, t=26, b=0),
                      title=dict(text=f"cell segmentation ({len(keep):,} cells)", font=dict(color="#e8eef7", size=12), x=0.01, y=0.99),
                      yaxis=dict(scaleanchor="x", scaleratio=1),
                      legend=dict(itemsizing="constant", x=0.997, y=0.997, xanchor="right", yanchor="top",
                                  bgcolor="rgba(8,11,16,.66)", bordercolor="#1d2531", borderwidth=1, font=dict(size=10)))
    _apply_dark(fig, hide_axes=True)
    fig.update_yaxes(autorange="reversed")
    return fig, len(keep), note


def composition_view(adata, groupby: str | None = None, height: int = 380):
    """Plotly bar chart of ``groupby`` composition (cell counts per category).

    Shares its denominator with :func:`export.composition_table` and the delivered PNG: abstained
    cells are excluded (they are not cell types) and their share is stated in the title. The
    interactive chart and the report must never disagree about what fraction of a section is "T cell".
    """
    import plotly.graph_objects as go

    from . import annotate as _an
    from . import export as _ex
    from .plots import _PALETTE, _apply_dark

    groupby = groupby or _an.annotation_key(adata)
    vc, abstained = _ex.composition_table(adata, groupby)
    fig = go.Figure(go.Bar(x=list(map(str, vc.index)), y=vc.to_numpy(),
                           marker=dict(color=[_PALETTE[i % len(_PALETTE)] for i in range(len(vc))])))
    fig.update_layout(height=height, margin=dict(l=0, r=0, t=28, b=0),
                      title=dict(text=f"{groupby} composition ({int(vc.sum()):,} typed cells; "
                                      f"{abstained:.1%} abstained)",
                                 font=dict(color="#e8eef7", size=12), x=0.01, y=0.99),
                      xaxis=dict(tickangle=-40), yaxis=dict(title="cells"))
    _apply_dark(fig)
    return fig
