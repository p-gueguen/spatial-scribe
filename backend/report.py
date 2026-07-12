"""The shareable HTML report - a self-contained, on-brand summary of a SpatialScribe session.

This is the "hand it to your PI" deliverable: one `.html` file, no server, no CDN, no external
fonts. Every figure is matplotlib rasterised to a base64 PNG data URI (NO kaleido, NO plotly
``to_image`` - neither is on the deploy). The report is dark-field to match how imaging spatial
data is actually viewed (docs/BRAND.md report surface); the surrounding chrome is quiet so the
figures are the hero.

Entry point: ``report_html(session, *, summary, colours, obs_fields, device) -> str``. The four
keyword helpers come from ``backend/app.py`` (``_summary`` / ``_colours`` / ``_obs_fields`` /
``_device``) and are passed in so this module never imports ``app`` (which would be circular).
They may be passed as callables (the functions) or, for ``summary``/``device``, as pre-computed
values; both are accepted.

Every figure degrades: the builder for a section guards on the obs column(s) it needs and yields
nothing when they are absent, so a section that has only been loaded and QC'd renders cleanly
without the annotation figures - no traceback, no empty panels. Nothing here invents a number:
every value shown is computed from the session AnnData or the summary dict.
"""

from __future__ import annotations

import html as _html
from typing import Callable

import numpy as np

from spatialscribe.analysis import export as _ex
from spatialscribe.analysis import plots

_esc = _html.escape

# The report brand tokens (docs/BRAND.md report surface). Kept in sync with export._R_* (the
# matplotlib figures) so the page and its figures read as one dark-field surface.
_REPORT_CSS = (
    ":root{color-scheme:dark;--bg:#0b0c10;--panel:#14161c;--panel2:#171a21;--ink:#e7e9ee;"
    "--muted:#8892a6;--faint:#5b6274;--line:#20242e;--violet:#a896f2;--pass:#46e39b;"
    "--warn:#f2b24c;--fail:#f7746e}"
    "*{box-sizing:border-box}"
    "body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.6 system-ui,-apple-system,"
    "'Segoe UI',Roboto,sans-serif;-webkit-font-smoothing:antialiased}"
    "header{display:flex;gap:15px;align-items:center;padding:26px 32px;border-bottom:1px solid "
    "var(--line);background:linear-gradient(180deg,#14161c,#0b0c10)}"
    ".mark{width:44px;height:44px;border-radius:12px;display:grid;place-items:center;flex:0 0 auto;"
    "background:#14161c;border:1px solid rgba(168,150,242,.55);color:var(--violet);font-size:22px}"
    "h1{margin:0;font-size:20px;letter-spacing:.2px;font-weight:650}"
    ".sub{color:var(--muted);font-size:12.5px;margin-top:4px;font-family:ui-monospace,"
    "SFMono-Regular,Menlo,monospace;letter-spacing:.02em}"
    "main{max-width:960px;margin:0 auto;padding:8px 24px 56px}"
    "section{margin-top:40px}"
    ".kick{font-size:11px;text-transform:uppercase;letter-spacing:.16em;color:var(--violet);"
    "font-family:ui-monospace,SFMono-Regular,Menlo,monospace;margin:0 0 4px}"
    "h2{font-size:19px;font-weight:640;letter-spacing:-.01em;margin:0 0 4px;color:var(--ink)}"
    ".lede{color:var(--muted);font-size:13px;margin:0 0 16px;max-width:70ch}"
    ".hr{height:1px;background:var(--line);border:0;margin:0 0 18px}"
    ".grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(120px,1fr));gap:12px}"
    ".card{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:15px 16px}"
    ".card .v{font-size:22px;font-weight:650;color:#fff;font-family:ui-monospace,SFMono-Regular,"
    "Menlo,monospace;line-height:1.1}"
    ".card .l{font-size:11.5px;color:var(--muted);margin-top:5px;text-transform:uppercase;"
    "letter-spacing:.05em}"
    ".fig{background:var(--panel);border:1px solid var(--line);border-radius:14px;padding:14px 16px 6px;"
    "margin:14px 0}"
    ".fig img{width:100%;height:auto;border-radius:8px;display:block}"
    ".cap{color:var(--muted);font-size:12px;margin:10px 2px 8px;line-height:1.45}"
    ".cap b{color:var(--ink);font-weight:600}"
    ".maps{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:14px;margin:14px 0}"
    ".maps .fig{margin:0}"
    ".legend{display:flex;flex-wrap:wrap;gap:6px 12px;margin:8px 2px 6px}"
    ".lg{display:inline-flex;align-items:center;gap:6px;font-size:11.5px;color:#c7cbd6;"
    "font-family:ui-monospace,SFMono-Regular,Menlo,monospace}"
    ".lg i{width:10px;height:10px;border-radius:3px;display:inline-block;flex:0 0 auto}"
    ".chips{display:flex;flex-wrap:wrap;gap:8px;margin:6px 0 4px}"
    ".chip{font-size:12px;padding:5px 11px;border-radius:999px;background:var(--panel2);"
    "border:1px solid var(--line);color:#c7cbd6;display:inline-flex;gap:7px;align-items:center}"
    ".chip .d{width:8px;height:8px;border-radius:99px;flex:0 0 auto}"
    ".chip.ok{border-color:rgba(70,227,155,.4)}.chip.ok .d{background:var(--pass)}"
    ".chip.warn{border-color:rgba(242,178,76,.4)}.chip.warn .d{background:var(--warn)}"
    ".chip.bad{border-color:rgba(247,116,110,.42)}.chip.bad .d{background:var(--fail)}"
    ".chip .m{color:var(--muted)}"
    ".subhead{font-size:12px;font-weight:600;color:var(--ink);margin:20px 0 2px;"
    "text-transform:uppercase;letter-spacing:.05em}"
    ".muted{color:var(--muted)}"
    "table{border-collapse:collapse;width:100%;margin-top:8px;font-size:13px}"
    "th,td{text-align:left;padding:8px 10px;border-bottom:1px solid var(--line)}"
    "td:last-child,th:last-child{text-align:right;font-family:ui-monospace,SFMono-Regular,Menlo,"
    "monospace}"
    "th{color:var(--muted);font-weight:600;font-size:10.5px;text-transform:uppercase;letter-spacing:.08em}"
    ".code{background:#0f1116;border:1px solid var(--line);border-radius:12px;padding:15px 16px;"
    "overflow-x:auto;font:12.5px/1.55 ui-monospace,SFMono-Regular,Menlo,monospace;color:#c7cbd6;"
    "white-space:pre;margin-top:6px}"
    "footer{max-width:960px;margin:0 auto;padding:26px 24px;color:var(--faint);font-size:12px;"
    "border-top:1px solid var(--line)}"
)


def _img(fig, alt: str) -> str:
    """Rasterise a matplotlib figure to an inline ``<img>`` (base64 PNG) and close it. '' if None."""
    if fig is None:
        return ""
    src = _ex.fig_to_base64(fig, dpi=150)
    import matplotlib.pyplot as plt
    plt.close(fig)
    if not src:
        return ""
    return f'<img src="{src}" alt="{_esc(alt)}"/>'


def _fig_block(fig, alt: str, caption: str = "") -> str:
    """A framed figure card: the image plus an optional muted caption. '' when the figure is absent."""
    tag = _img(fig, alt)
    if not tag:
        return ""
    cap = f'<div class="cap">{caption}</div>' if caption else ""
    return f'<div class="fig">{tag}{cap}</div>'


def _map_png(a, color_by: str, colours, max_points: int = 60_000):
    """Inline base64 spatial-map PNG (dark-field) + its categorical legend, coloured by ``colours``
    (app._colours - the same palette the interactive canvas uses). Returns (data_uri|None, legend)."""
    try:
        import base64 as _b64
        import io as _io
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return None, []
    if "spatial" not in a.obsm or color_by not in a.obs:
        return None, []
    idx = plots._downsample(a.n_obs, max_points)
    xy = np.asarray(a.obsm["spatial"])[idx]
    rgb, legend = colours(a, color_by, idx)
    cols = np.array(rgb, dtype=float).reshape(-1, 3) / 255.0
    fig, ax = plt.subplots(figsize=(6.4, 6.4))
    fig.patch.set_facecolor("#0b0c10"); ax.set_facecolor("#0b0c10")
    ax.scatter(xy[:, 0], xy[:, 1], s=1.6, c=cols, linewidths=0, marker=".")
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.axis("off")
    buf = _io.BytesIO()
    fig.savefig(buf, format="png", dpi=120, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig); buf.seek(0)
    return "data:image/png;base64," + _b64.b64encode(buf.read()).decode(), legend


def _legend_html(legend) -> str:
    """A colour-swatch legend from ``_colours``' ``[{label, color:[r,g,b]}]``. '' when empty."""
    if not legend:
        return ""
    return '<div class="legend">' + "".join(
        f'<span class="lg"><i style="background:rgb({c[0]},{c[1]},{c[2]})"></i>{_esc(str(l["label"]))}</span>'
        for l in legend for c in [l["color"]]) + "</div>"


def _map_fig(a, color_by: str, colours, title: str) -> tuple[str, dict]:
    """A framed spatial map (title, image, legend). Returns (html, {label: rgb-float}) - the colour
    map is reused so the composition bar / confidence violins match the map exactly. ('' when absent)."""
    img, legend = _map_png(a, color_by, colours)
    if not img:
        return "", {}
    colour_map = {str(l["label"]): tuple(v / 255.0 for v in l["color"]) for l in legend}
    html = (f'<div class="fig"><div class="cap" style="margin:2px 2px 8px"><b>{_esc(title)}</b></div>'
            f'<img src="{img}" alt="{_esc(title)}"/>{_legend_html(legend)}</div>')
    return html, colour_map


def _cards(items) -> str:
    return '<div class="grid">' + "".join(
        f'<div class="card"><div class="v">{_esc(str(v))}</div><div class="l">{_esc(str(l))}</div></div>'
        for v, l in items) + "</div>"


def _chips(items, cls: str = "") -> str:
    if not items:
        return '<p class="muted">none</p>'
    return '<div class="chips">' + "".join(
        f'<span class="chip {cls}">{_esc(str(x))}</span>' for x in items) + "</div>"


# Section-flag metric -> human label + whether it reads as a percentage (qc.qc_summary metrics).
_FLAG_META = {
    "median_genes_per_cell": ("median genes / cell", False),
    "median_transcripts_per_cell": ("median counts / cell", False),
    "fraction_empty_cells": ("empty cells", True),
    "pct_counts_control": ("control transcripts", True),
    "annotation_usability": ("annotation usability", None),
}
_FLAG_CLASS = {"ok": "ok", "warn": "warn", "error": "bad", "fail": "bad"}


def _section_head(kicker: str, title: str, lede: str = "") -> str:
    ld = f'<p class="lede">{_esc(lede)}</p>' if lede else ""
    return (f'<p class="kick">{_esc(kicker)}</p><h2>{_esc(title)}</h2>'
            f'<hr class="hr">{ld}')


def report_html(session: dict, *, summary: "Callable[[dict], dict] | dict", colours: "Callable",
                obs_fields: "Callable", device: "Callable[[], str] | str") -> str:
    """Assemble the self-contained HTML report for a session.

    ``summary`` / ``device`` may be the app helper (called here) or a pre-computed value;
    ``colours`` (app._colours) and ``obs_fields`` (app._obs_fields) are always callables. Taking them
    as arguments rather than importing them is what keeps ``backend.app`` -> ``backend.report`` a
    one-way edge. Every section is optional and self-guarding, so the report renders for any session
    state - from a freshly loaded + QC'd section (QC figures only) through a fully annotated one.
    """
    import datetime as _dt
    from spatialscribe.analysis import qc as _qc

    a = session["adata"]
    # Probe the section BEFORE calling `summary`: app._summary -> qc.qc_summary back-fills
    # obs['total_counts'] when QC has not run, so a flag read afterwards always says "QC ran" and the
    # report would show a Quality-control section for a section the user only loaded.
    has_ct = "cell_type" in a.obs
    has_verdict = "annotation_verdict" in a.obs
    has_qc = "total_counts" in a.obs

    sm = summary(session) if callable(summary) else summary
    dev = device() if callable(device) else device
    tissue = session.get("tissue", "")
    platform = str(a.uns.get("platform", "spatial"))
    when = _dt.datetime.now().strftime("%Y-%m-%d %H:%M")

    sec: list[str] = []

    # --- Overview ----------------------------------------------------------------------------- #
    ov = [(f"{a.n_obs:,}", "cells"), (f"{a.n_vars:,}", "genes"), (platform, "platform"),
          (tissue or "-", "tissue")]
    if sm.get("cluster"):
        ov.append((sm["cluster"]["n_clusters"], "clusters"))
    if sm.get("spatial"):
        ov.append((sm["spatial"]["n_niches"], "niches"))
    if sm.get("annotate", {}).get("confident_pct") is not None:
        ov.append((f'{sm["annotate"]["confident_pct"]}%', "confident"))
    sec.append("<section>" + _section_head("overview", "Section at a glance") + _cards(ov) + "</section>")

    # --- Quality control ---------------------------------------------------------------------- #
    if has_qc:
        qcfrag = []
        floor_info = _qc.suggest_count_floor(a)
        floor = float(floor_info.get("floor", 0.0))
        removed = float(floor_info.get("pct_removed", 0.0)) * 100
        if floor > 0 and removed > 0:
            hist_cap = ("Transcripts and genes detected per cell (log scale). The dashed rule is the "
                        f"panel-indexed count floor (<b>{floor:g}</b>); <b>{removed:.1f}%</b> of cells "
                        "fall below it.")
        else:
            hist_cap = ("Transcripts and genes detected per cell (log scale). The panel-indexed count "
                        "floor removes no further cells here - this section is already filtered.")
        qcfrag.append(_fig_block(_ex.fig_qc_histograms(a), "counts and genes per cell", hist_cap))
        funnel = _fig_block(
            _ex.fig_qc_funnel(a), "qc funnel",
            "Cells surviving each QC gate, applied cumulatively in pipeline order. A gate is drawn "
            "only when the section carries the column it needs.")
        if funnel:
            qcfrag.append(funnel)
        try:
            m = _qc.qc_summary(a)
            flags = m.get("flags", {})
            chips = []
            for key, status in flags.items():
                label, is_pct = _FLAG_META.get(key, (key.replace("_", " "), None))
                cls = _FLAG_CLASS.get(str(status), "")
                val = m.get(key)
                if isinstance(val, (int, float)) and is_pct is not None:
                    vtxt = f"{val * 100:.1f}%" if is_pct else f"{val:g}"
                    chips.append(f'<span class="chip {cls}"><span class="d"></span>{_esc(label)}'
                                 f'<span class="m">{_esc(vtxt)}</span></span>')
                else:
                    chips.append(f'<span class="chip {cls}"><span class="d"></span>{_esc(label)}'
                                 f'<span class="m">{_esc(str(status))}</span></span>')
            if chips:
                qcfrag.append('<div class="subhead">section flags</div>'
                              '<div class="chips">' + "".join(chips) + "</div>")
        except Exception:
            pass
        if qcfrag:
            sec.append("<section>" + _section_head(
                "quality control", "Quality control",
                "Panel-aware per-cell QC and the survivor funnel. Flags are scored against the "
                "platform's thresholds.") + "".join(qcfrag) + "</section>")

    # --- Annotation --------------------------------------------------------------------------- #
    ann = sm.get("annotate") or {}
    ct_colour: dict = {}
    if has_ct or has_verdict:
        afrag = []
        # Two spatial maps side by side: the labels, and how much to trust them.
        maps = []
        if has_ct:
            html, ct_colour = _map_fig(a, "cell_type", colours, "cell type")
            if html:
                maps.append(html)
        if has_verdict:
            html, _ = _map_fig(a, "annotation_verdict", colours, "annotation verdict")
            if html:
                maps.append(html)
        if maps:
            afrag.append('<div class="maps">' + "".join(maps) + "</div>")

        if has_verdict:
            conf = ann.get("confident_pct")
            cap = (f"<b>{conf}%</b> confident, <b>{ann.get('tentative_pct')}%</b> tentative, "
                   f"<b>{ann.get('abstained_pct')}%</b> abstained across {a.n_obs:,} cells."
                   ) if conf is not None else "Confident / tentative / abstained split per cell."
            afrag.append(_fig_block(_ex.fig_confidence_meter(a), "annotatability meter", cap))

        if has_ct:
            afrag.append(_fig_block(
                _ex.fig_composition_bar(a, colour_map=ct_colour), "cell-type composition",
                "Cell-type composition as a share of all cells (tail merged into 'other')."))
            ridges = _fig_block(
                _ex.fig_confidence_ridges(a, colour_map=ct_colour), "confidence by cell type",
                "Per-cell annotation confidence by cell type; the light tick marks the median.")
            if ridges:
                afrag.append(ridges)

        rej = _fig_block(_ex.fig_rejection_bar(a), "rejection reasons",
                         "Why abstained cells were left unlabelled.")
        if rej:
            afrag.append(rej)

        if afrag:
            sec.append("<section>" + _section_head(
                "annotation", "Cell-type annotation",
                "Consensus labels with a per-cell confidence verdict. The maps place the labels and "
                "their trust in space.") + "".join(afrag) + "</section>")

    # --- Panel check -------------------------------------------------------------------------- #
    p = sm.get("panel")
    if p:
        h = ('<div class="subhead">resolvable</div>' + _chips(p["resolvable"], "ok")
             + '<div class="subhead">weakly resolved</div>' + _chips(p["weak"], "warn")
             + '<div class="subhead">cannot resolve</div>' + _chips(p["cannot"], "bad"))
        if p.get("confusable"):
            h += ('<div class="subhead">confusable pairs</div>' + _chips(
                [" ~ ".join(pr) if isinstance(pr, (list, tuple)) else str(pr) for pr in p["confusable"]]))
        sec.append("<section>" + _section_head(
            "panel check", "Panel adequacy",
            "Which cell types this gene panel can resolve, weakly resolve, or cannot separate.")
            + h + "</section>")

    # --- Spatial niches ----------------------------------------------------------------------- #
    spn = sm.get("spatial")
    if spn and spn.get("niches"):
        rows = "".join(
            f"<tr><td>{_esc(str(n['id']) if str(n['id']).lower().startswith('niche') else 'niche ' + str(n['id']))}</td>"
            f"<td>{n['count']:,}</td></tr>" for n in spn["niches"])
        sec.append("<section>" + _section_head(
            "spatial niches", "Spatial niches",
            "Recurrent neighbourhoods of co-located cell types.")
            + f"<table><thead><tr><th>niche</th><th>cells</th></tr></thead><tbody>{rows}</tbody></table>"
            + "</section>")

    # --- Re-runnable analysis ----------------------------------------------------------------- #
    body = _ex.build_runnable_script(session.get("log") or [], adata=a, tissue=tissue,
                                     source_path=session.get("source_path"))
    sec.append("<section>" + _section_head(
        "reproducibility", "Re-runnable analysis",
        "The exact pipeline that produced this section, as a runnable Python script.")
        + f'<pre class="code">{_esc(body)}</pre>' + "</section>")

    body_html = "".join(sec)
    return (
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width,initial-scale=1'>"
        "<title>SpatialScribe report</title><style>" + _REPORT_CSS + "</style></head>"
        "<body><header><div class='mark'>&#9672;</div><div><h1>SpatialScribe report</h1>"
        f"<div class='sub'>{_esc(f'{a.n_obs:,}')} cells &middot; {_esc(str(tissue))} &middot; "
        f"{_esc(str(dev))} &middot; {_esc(when)}</div></div></header>"
        f"<main>{body_html}</main>"
        "<footer>Generated by SpatialScribe &middot; spatial-transcriptomics annotation copilot. "
        "Self-contained; open in any browser.</footer>"
        "</body></html>"
    )
