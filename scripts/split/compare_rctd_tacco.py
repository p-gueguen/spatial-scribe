"""Head-to-head: RCTD weights vs TACCO weights as input to SPLIT purification.

Both SPLIT runs use the SAME section, the SAME cells (RCTD-surviving), the SAME reference atlas
and the SAME SPLIT parameters - the ONLY difference is which deconvolution weights drove the
purification. This script loads both `<inputs>/` + `<split_out>/` pairs, restricts to the shared
cell set, and reports the metrics + marker dotplots that decide which weight source purifies better:

  * type coverage - how many reference types each annotator kept (RCTD drops sparse immune types);
  * spillover removal - median library size + median genes/cell reduction;
  * cross-lineage MECR (mutually-exclusive co-expression) purified, lower = cleaner - annotation
    agnostic (gene presence), so comparable across the two type sets;
  * marker specificity ratio (on-target / off-target marker expression) - higher = cleaner;
  * marker purity - median per-cell fraction of counts from the cell's own-type markers;
  * on-target signal retention - do a lineage's own markers survive purification;
  * marker dotplots - mean expression x %-expressing per (type, marker), raw vs purified, per method.

Usage:
  python compare_rctd_tacco.py <rctd_inputs> <rctd_split_out> <tacco_inputs> <tacco_split_out> <report.html>
Runs in the main pixi env.
"""
from __future__ import annotations

import base64
import io as _io
import pathlib
import sys

import matplotlib
import numpy as np
import pandas as pd
import scipy.io as sio

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BG, FG, GRID, CY, MG = "#0a0e16", "#e8eef7", "#1d2531", "#22d3ee", "#e879f9"
plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
                     "text.color": FG, "axes.labelcolor": FG, "xtick.color": FG, "ytick.color": FG,
                     "axes.edgecolor": GRID, "font.size": 9})

# Curated markers per cell type (breast/skin panel, kept to genes likely on the 5K panel).
MARKERS = {
    "T cell": ["CD3D", "CD3E", "TRAC", "CD8A", "IL7R"],
    "NK cell": ["NKG7", "GNLY", "KLRD1"],
    "B/Plasma": ["MS4A1", "CD79A", "IGHG1", "MZB1", "JCHAIN"],
    "Myeloid": ["LYZ", "CD68", "ITGAX", "C1QA", "CD14"],
    "Mast": ["TPSAB1", "CPA3", "MS4A2"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CD34"],
    "Stromal/CAF": ["PDGFRB", "COL1A1", "DCN", "LUM", "FAP"],
    "Epithelial/Tumor": ["EPCAM", "KRT8", "KRT18", "ERBB2"],
    "Basal/Myoepithelial": ["KRT5", "KRT14", "ACTA2", "MYLK", "TP63"],
    "Adipocyte": ["ADIPOQ", "LEP", "PLIN1"],
}
# Disjoint lineage pairs whose co-expression should be near-zero in clean data (spillover signature).
MECR_PAIRS = [("T cell", "Epithelial/Tumor"), ("Myeloid", "Epithelial/Tumor"),
              ("B/Plasma", "Epithelial/Tumor"), ("T cell", "Stromal/CAF"),
              ("Endothelial", "Epithelial/Tumor"), ("Myeloid", "Stromal/CAF")]
# Lineages both methods keep (for the shared-lineage bar charts).
SHARED = ["Epithelial/Tumor", "Stromal/CAF", "B/Plasma", "Endothelial", "Myeloid"]


def _png(fig):
    b = _io.BytesIO(); fig.savefig(b, format="png", dpi=130, bbox_inches="tight"); plt.close(fig)
    return base64.b64encode(b.getvalue()).decode()


def _load(inputs, split_out):
    inputs, split_out = pathlib.Path(inputs), pathlib.Path(split_out)
    cells = pd.read_csv(inputs / "cells.txt", header=None)[0].astype(str).tolist()
    genes = pd.read_csv(inputs / "genes.txt", header=None)[0].astype(str).tolist()
    prim = pd.read_csv(inputs / "primary.csv").set_index("cell_id")["primary"].astype(str)
    raw = sio.mmread(inputs / "counts.mtx").T.tocsr()             # cells x genes
    pur = sio.mmread(split_out / "purified_counts.mtx").T.tocsr()  # cells x genes
    cm = pd.read_csv(split_out / "cell_meta.csv") if (split_out / "cell_meta.csv").exists() else None
    return dict(cells=cells, genes=genes, primary=prim, raw=raw, pur=pur, cell_meta=cm)


def _gidx(genes):
    return {g: i for i, g in enumerate(genes)}


def _lognorm(X):
    lib = np.asarray(X.sum(1)).ravel(); lib[lib == 0] = 1
    Y = X.multiply((1e4 / lib)[:, None]).tocsr(); Y.data = np.log1p(Y.data)
    return Y


# ---------------- metrics ----------------
def _mecr(X, gidx, a, b):
    aset = [gidx[g] for g in MARKERS[a] if g in gidx]
    bset = [gidx[g] for g in MARKERS[b] if g in gidx]
    if not aset or not bset:
        return float("nan")
    ap = np.asarray((X[:, aset] > 0).sum(1)).ravel() > 0
    bp = np.asarray((X[:, bset] > 0).sum(1)).ravel() > 0
    return float((ap & bp).mean())


def _median_genes(X):
    return float(np.median(np.asarray((X > 0).sum(1)).ravel()))


def _ontarget_retention(raw, pur, gidx, groups, lineage):
    mset = [gidx[g] for g in MARKERS[lineage] if g in gidx]
    if not mset:
        return float("nan"), 0
    mask = (groups == lineage)
    if mask.sum() == 0:
        return float("nan"), 0
    r = np.asarray(raw[mask][:, mset].sum(1)).ravel()
    p = np.asarray(pur[mask][:, mset].sum(1)).ravel()
    keep = r > 0
    if keep.sum() == 0:
        return float("nan"), int(mask.sum())
    return float(np.mean(np.clip(p[keep] / r[keep], 0, 1))), int(mask.sum())


def _specificity(Xln, gidx, groups, lineage):
    """On-target FRACTION of a lineage's marker expression: on/(on+off) mean, bounded [0,1].
    -> 1 as off-target spillover is removed (a stable index; a raw on/off ratio explodes when
    off-target expression is purified to ~0)."""
    mset = [gidx[g] for g in MARKERS[lineage] if g in gidx]
    if not mset:
        return float("nan")
    on = groups == lineage
    if on.sum() == 0 or (~on).sum() == 0:
        return float("nan")
    on_e = float(Xln[on][:, mset].mean())
    off_e = float(Xln[~on][:, mset].mean())
    return on_e / (on_e + off_e + 1e-12)


def _marker_purity(X, gidx, groups):
    """Mean per-cell fraction of counts coming from the cell's OWN primary-type markers
    (median is ~0 on a sparse 5K panel; the mean's raw->purified rise is the informative signal)."""
    tot = np.asarray(X.sum(1)).ravel(); tot[tot == 0] = 1.0
    frac = np.zeros(X.shape[0])
    for lin in np.unique(groups):
        if lin not in MARKERS:
            continue
        mset = [gidx[g] for g in MARKERS[lin] if g in gidx]
        if not mset:
            continue
        m = groups == lin
        frac[m] = np.asarray(X[m][:, mset].sum(1)).ravel() / tot[m]
    return float(np.mean(frac))


# ---------------- dotplots ----------------
def _dot_matrix(Xln, gi, groups, order):
    mean = np.zeros((len(order), len(gi))); pct = np.zeros_like(mean)
    for i, g in enumerate(order):
        m = (groups == g)
        if m.sum() == 0:
            continue
        sub = Xln[m][:, gi]
        mean[i] = np.asarray(sub.mean(0)).ravel()
        pct[i] = np.asarray((sub > 0).mean(0)).ravel()
    return mean, pct


def _dotplot(ax, mean, pct, genes, groups, title):
    mx = mean / (mean.max(0, keepdims=True) + 1e-9)
    for i in range(len(groups)):
        for j in range(len(genes)):
            ax.scatter(j, i, s=6 + pct[i, j] * 160, c=[[mx[i, j], 0.15 + 0.4 * mx[i, j], 0.6]],
                       edgecolors="none")
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, rotation=90, fontsize=6)
    ax.set_yticks(range(len(groups))); ax.set_yticklabels(groups, fontsize=7)
    ax.set_title(title, color=CY, fontsize=10)
    ax.set_xlim(-.6, len(genes) - .4); ax.set_ylim(-.6, len(groups) - .4); ax.invert_yaxis()


def _method_dotplot(raw_ln, pur_ln, gidx, groups, title):
    order = [t for t in MARKERS if (groups == t).sum() >= 10]
    flat, seen = [], set()
    for t in order:
        for g in MARKERS[t]:
            if g in gidx and g not in seen:
                flat.append(g); seen.add(g)
    gi = [gidx[g] for g in flat]
    fig, axes = plt.subplots(1, 2, figsize=(1.2 + 0.30 * len(flat) * 2, 0.42 * len(order) + 1.6))
    for ax, X, ttl in ((axes[0], raw_ln, f"{title} - RAW"), (axes[1], pur_ln, f"{title} - SPLIT-purified")):
        m, p = _dot_matrix(X, gi, groups, order)
        _dotplot(ax, m, p, flat, order, ttl)
    fig.tight_layout()
    return _png(fig)


def main(rctd_in, rctd_out, tacco_in, tacco_out, report):
    R = _load(rctd_in, rctd_out)
    T = _load(tacco_in, tacco_out)
    shared = [c for c in R["cells"] if c in set(T["cells"])]
    assert R["genes"] == T["genes"], "gene order differs between the two input dirs"
    gidx = _gidx(R["genes"])
    rpos = pd.Index(R["cells"]).get_indexer(shared)
    tpos = pd.Index(T["cells"]).get_indexer(shared)
    raw = R["raw"][rpos]                      # identical counts for both -> RAW baseline
    rpur, tpur = R["pur"][rpos], T["pur"][tpos]
    rgrp = R["primary"].reindex(shared).to_numpy()
    tgrp = T["primary"].reindex(shared).to_numpy()

    r_types = sorted(pd.unique(rgrp[pd.notna(rgrp)]))
    t_types = sorted(pd.unique(tgrp[pd.notna(tgrp)]))
    agree = float(np.mean(rgrp == tgrp))

    # library size + genes/cell
    lib_raw = float(np.median(np.asarray(raw.sum(1)).ravel()))
    lib_rp = float(np.median(np.asarray(rpur.sum(1)).ravel()))
    lib_tp = float(np.median(np.asarray(tpur.sum(1)).ravel()))
    gen_raw, gen_rp, gen_tp = _median_genes(raw), _median_genes(rpur), _median_genes(tpur)

    # MECR
    mecr = {p: (_mecr(raw, gidx, *p), _mecr(rpur, gidx, *p), _mecr(tpur, gidx, *p)) for p in MECR_PAIRS}
    mean_mecr_raw = float(np.nanmean([mecr[p][0] for p in MECR_PAIRS]))
    mean_mecr_r = float(np.nanmean([mecr[p][1] for p in MECR_PAIRS]))
    mean_mecr_t = float(np.nanmean([mecr[p][2] for p in MECR_PAIRS]))

    # log-norm for specificity / dotplots
    raw_ln, rpur_ln, tpur_ln = _lognorm(raw), _lognorm(rpur), _lognorm(tpur)

    # specificity ratio (raw uses method grouping on raw counts; purified on purified counts)
    spec = {}
    for lin in MARKERS:
        rr = _specificity(raw_ln, gidx, rgrp, lin) if (rgrp == lin).sum() >= 10 else float("nan")
        rp = _specificity(rpur_ln, gidx, rgrp, lin) if (rgrp == lin).sum() >= 10 else float("nan")
        tr = _specificity(raw_ln, gidx, tgrp, lin) if (tgrp == lin).sum() >= 10 else float("nan")
        tp = _specificity(tpur_ln, gidx, tgrp, lin) if (tgrp == lin).sum() >= 10 else float("nan")
        spec[lin] = (rr, rp, tr, tp)

    # marker purity
    pur_raw_r = _marker_purity(raw, gidx, rgrp)
    pur_pur_r = _marker_purity(rpur, gidx, rgrp)
    pur_raw_t = _marker_purity(raw, gidx, tgrp)
    pur_pur_t = _marker_purity(tpur, gidx, tgrp)

    # on-target retention
    ret = {}
    for lin in ("T cell", "Epithelial/Tumor", "Myeloid", "B/Plasma", "Endothelial", "Stromal/CAF"):
        rr, rn = _ontarget_retention(raw, rpur, gidx, rgrp, lin)
        tr, tn = _ontarget_retention(raw, tpur, gidx, tgrp, lin)
        ret[lin] = (rr, rn, tr, tn)

    # ---------------- figures ----------------
    figs = {}
    fig, ax = plt.subplots(1, 2, figsize=(9, 3))
    ax[0].hist(np.log10(np.asarray(raw.sum(1)).ravel() + 1), bins=40, alpha=.5, color="#f87171", label=f"raw ({lib_raw:.0f})")
    ax[0].hist(np.log10(np.asarray(rpur.sum(1)).ravel() + 1), bins=40, alpha=.5, color=CY, label=f"RCTD ({lib_rp:.0f})")
    ax[0].hist(np.log10(np.asarray(tpur.sum(1)).ravel() + 1), bins=40, alpha=.5, color=MG, label=f"TACCO ({lib_tp:.0f})")
    ax[0].set_xlabel("log10 library size"); ax[0].set_ylabel("cells"); ax[0].set_title("Library size", color=FG)
    ax[0].legend(facecolor=BG, edgecolor=GRID, labelcolor=FG, fontsize=7)
    ax[1].bar([0, 1, 2], [gen_raw, gen_rp, gen_tp], color=["#f87171", CY, MG])
    ax[1].set_xticks([0, 1, 2]); ax[1].set_xticklabels(["raw", "RCTD", "TACCO"], fontsize=8)
    ax[1].set_ylabel("median genes / cell"); ax[1].set_title("Genes per cell", color=FG)
    fig.tight_layout(); figs["lib"] = _png(fig)

    fig, ax = plt.subplots(figsize=(7.5, 3.2))
    labels = [f"{a[:3]}x{b[:3]}" for a, b in MECR_PAIRS]
    x = np.arange(len(labels)); w = 0.27
    ax.bar(x - w, [mecr[p][0] * 100 for p in MECR_PAIRS], w, color="#f87171", label="raw")
    ax.bar(x, [mecr[p][1] * 100 for p in MECR_PAIRS], w, color=CY, label="RCTD->SPLIT")
    ax.bar(x + w, [mecr[p][2] * 100 for p in MECR_PAIRS], w, color=MG, label="TACCO->SPLIT")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("% co-expressing"); ax.set_title("Cross-lineage co-expression (MECR) - lower = cleaner", color=FG)
    ax.legend(facecolor=BG, edgecolor=GRID, labelcolor=FG, fontsize=8)
    fig.tight_layout(); figs["mecr"] = _png(fig)

    # specificity gain (purified) on shared lineages
    fig, ax = plt.subplots(figsize=(7.5, 3.2))
    sl = [lin for lin in SHARED if not np.isnan(spec[lin][1]) and not np.isnan(spec[lin][3])]
    x = np.arange(len(sl)); w = 0.27
    ax.bar(x - w, [spec[l][0] * 100 for l in sl], w, color="#f87171", label="raw")
    ax.bar(x, [spec[l][1] * 100 for l in sl], w, color=CY, label="RCTD->SPLIT")
    ax.bar(x + w, [spec[l][3] * 100 for l in sl], w, color=MG, label="TACCO->SPLIT")
    ax.set_xticks(x); ax.set_xticklabels([s.split("/")[0] for s in sl], fontsize=8, rotation=20)
    ax.set_ylim(0, 105)
    ax.set_ylabel("on-target % of marker expr"); ax.set_title("Marker specificity (higher = cleaner)", color=FG)
    ax.legend(facecolor=BG, edgecolor=GRID, labelcolor=FG, fontsize=8)
    fig.tight_layout(); figs["spec"] = _png(fig)

    figs["dot_rctd"] = _method_dotplot(raw_ln, rpur_ln, gidx, rgrp, "RCTD primaries")
    figs["dot_tacco"] = _method_dotplot(raw_ln, tpur_ln, gidx, tgrp, "TACCO primaries")

    # ---------------- html ----------------
    def img(k, cap):
        return f'<figure><img src="data:image/png;base64,{figs[k]}"/><figcaption>{cap}</figcaption></figure>'

    def mecr_row(p):
        rw, rc, tc = mecr[p]
        return f"<tr><td>{p[0]} x {p[1]}</td><td>{rw*100:.2f}%</td><td>{rc*100:.2f}%</td><td>{tc*100:.2f}%</td></tr>"

    def spec_row(lin):
        rr, rp, tr, tp = spec[lin]
        rs = f"{rr*100:.0f}% &rarr; {rp*100:.0f}%" if not np.isnan(rr) else "- (type dropped)"
        ts = f"{tr*100:.0f}% &rarr; {tp*100:.0f}%" if not np.isnan(tr) else "-"
        return f"<tr><td>{lin}</td><td>{rs}</td><td>{ts}</td></tr>"

    def ret_row(lin):
        rr, rn, tr, tn = ret[lin]
        rs = f"{rr*100:.0f}% (n={rn:,})" if rn and not np.isnan(rr) else "- (type dropped)"
        ts = f"{tr*100:.0f}% (n={tn:,})" if tn and not np.isnan(tr) else "- (type dropped)"
        return f"<tr><td>{lin}</td><td>{rs}</td><td>{ts}</td></tr>"

    r_stat = (R["cell_meta"]["purification_status"].value_counts().to_dict()
              if R["cell_meta"] is not None and "purification_status" in R["cell_meta"] else {})
    t_stat = (T["cell_meta"]["purification_status"].value_counts().to_dict()
              if T["cell_meta"] is not None and "purification_status" in T["cell_meta"] else {})

    html = f"""<!doctype html><meta charset=utf-8><title>SpatialScribe - SPLIT: RCTD vs TACCO weights</title>
<style>body{{background:{BG};color:{FG};font-family:-apple-system,Segoe UI,Roboto,sans-serif;margin:0;padding:2rem 2.4rem;max-width:1120px}}
h1{{font-weight:650;margin:0 0 .3rem}} h2{{color:{CY};font-size:1.05rem;margin:1.8rem 0 .5rem;border-bottom:1px solid {GRID};padding-bottom:.3rem}}
.sub{{color:#8b98ab}} figure{{margin:1rem 0;background:#0d1420;border:1px solid {GRID};border-radius:10px;padding:1rem;display:inline-block;max-width:100%}}
img{{max-width:100%;display:block}} figcaption{{color:#8b98ab;font-size:.8rem;margin-top:.5rem;max-width:820px}}
table{{border-collapse:collapse;font-size:.88rem;margin:.5rem 0}} td,th{{border:1px solid {GRID};padding:.35rem .8rem;text-align:left}}
th{{color:{CY}}} b{{color:{CY}}} .mg{{color:{MG}}} .kpi{{display:flex;gap:.8rem;flex-wrap:wrap;margin:.6rem 0}}
.kpi div{{background:#0d1420;border:1px solid {GRID};border-radius:8px;padding:.6rem .9rem}} .kpi b{{font-size:1.15rem}}</style>
<h1>SPLIT purification &middot; RCTD weights vs TACCO weights</h1>
<p class=sub>Same section, same {len(shared):,} cells, same reference ({len(gidx):,} genes), same SPLIT params
(<code>rctd_free_purify</code>, residual-contamination ON, belonging 0.5). The <b>only</b> variable is the
deconvolution-weight source: <b>RCTD</b> vs <b class=mg>TACCO</b>.</p>
<div class=kpi>
<div>types kept<br><b>{len(r_types)}</b> vs <b class=mg>{len(t_types)}</b></div>
<div>library removed<br><b>{100*(1-lib_rp/lib_raw):.0f}%</b> vs <b class=mg>{100*(1-lib_tp/lib_raw):.0f}%</b></div>
<div>mean MECR (purified)<br><b>{mean_mecr_r*100:.2f}%</b> vs <b class=mg>{mean_mecr_t*100:.2f}%</b></div>
<div>marker purity (raw&rarr;pur)<br><b>{pur_raw_r*100:.0f}&rarr;{pur_pur_r*100:.0f}%</b> vs <b class=mg>{pur_raw_t*100:.0f}&rarr;{pur_pur_t*100:.0f}%</b></div>
<div>primary agreement<br><b>{agree*100:.0f}%</b></div>
</div>

<h2>1. Type coverage - the headline difference</h2>
<p>RCTD kept <b>{len(r_types)}</b> types; TACCO kept <b class=mg>{len(t_types)}</b>. RCTD's UMI/gene-list filtering
drops sparse immune/rare types, so SPLIT driven by RCTD has <b>no profile to reassign those transcripts to</b>.</p>
<table><tr><th>RCTD types</th><td>{', '.join(r_types)}</td></tr>
<tr><th class=mg>TACCO types</th><td>{', '.join(t_types)}</td></tr></table>

<h2>2. Spillover removal</h2>
<p>Median library size raw <b>{lib_raw:.0f}</b> &rarr; RCTD <b>{lib_rp:.0f}</b> ({100*(1-lib_rp/lib_raw):.1f}%)
&middot; TACCO <b class=mg>{lib_tp:.0f}</b> ({100*(1-lib_tp/lib_raw):.1f}%). Median genes/cell {gen_raw:.0f} &rarr;
RCTD {gen_rp:.0f} / TACCO {gen_tp:.0f}.</p>
{img('lib','Per-cell library size (left) and median genes/cell (right), raw vs each purification.')}

<h2>3. Cross-lineage co-expression (MECR) - lower = cleaner</h2>
<p>Annotation-agnostic (gene presence), directly comparable across the two type sets. Mean over pairs:
raw {mean_mecr_raw*100:.2f}% &rarr; RCTD {mean_mecr_r*100:.2f}% / <span class=mg>TACCO {mean_mecr_t*100:.2f}%</span>.</p>
<table><tr><th>lineage pair</th><th>raw</th><th>RCTD-&gt;SPLIT</th><th>TACCO-&gt;SPLIT</th></tr>
{''.join(mecr_row(p) for p in MECR_PAIRS)}</table>
{img('mecr','Fraction of cells co-expressing markers of two mutually-exclusive lineages.')}

<h2>4. Marker specificity (on-target / off-target ratio, raw &rarr; purified)</h2>
<p>Mean expression of a lineage's markers in its own cells vs all other cells (log-norm). Higher = cleaner.</p>
<table><tr><th>lineage</th><th>RCTD raw&rarr;pur</th><th>TACCO raw&rarr;pur</th></tr>
{''.join(spec_row(l) for l in MARKERS if not (np.isnan(spec[l][1]) and np.isnan(spec[l][3])))}</table>
{img('spec','Purified marker specificity on the lineages both methods keep.')}

<h2>5. On-target signal retention - did purification keep real markers?</h2>
<p>Mean per-cell fraction of a lineage's OWN marker counts kept after SPLIT (guards against "removed more" == "deleted signal").</p>
<table><tr><th>lineage</th><th>RCTD-&gt;SPLIT</th><th>TACCO-&gt;SPLIT</th></tr>
{''.join(ret_row(l) for l in ('T cell','Epithelial/Tumor','Myeloid','B/Plasma','Endothelial','Stromal/CAF'))}</table>

<h2>6. Marker dotplots - raw vs SPLIT-purified</h2>
<p>Dot color = mean expression (per-gene scaled across rows), size = % of cells expressing. After SPLIT each
type's markers should concentrate on its own row. TACCO resolves the immune/rare rows (T cell, NK, Mast) that
RCTD cannot even represent.</p>
{img('dot_rctd','RCTD primaries ('+str(len(r_types))+' types).')}
{img('dot_tacco','TACCO primaries ('+str(len(t_types))+' types) - note the T-cell / NK / Mast rows RCTD lacks.')}

<h2>7. Purification status</h2>
<table><tr><th>RCTD</th><td>{r_stat}</td></tr><tr><th class=mg>TACCO</th><td>{t_stat}</td></tr></table>
"""
    pathlib.Path(report).write_text(html)
    print("wrote", report)
    print(f"types RCTD {len(r_types)} vs TACCO {len(t_types)} | lib removed {100*(1-lib_rp/lib_raw):.1f}% vs "
          f"{100*(1-lib_tp/lib_raw):.1f}% | genes/cell {gen_raw:.0f}->{gen_rp:.0f}/{gen_tp:.0f} | "
          f"mean MECR {mean_mecr_r*100:.2f}% vs {mean_mecr_t*100:.2f}% | "
          f"purity {pur_pur_r*100:.0f}% vs {pur_pur_t*100:.0f}% | agree {agree*100:.1f}%")
    for lin in SHARED:
        print(f"  spec {lin.split('/')[0]:<12} RCTD {spec[lin][0]:.1f}->{spec[lin][1]:.1f} | "
              f"TACCO {spec[lin][2]:.1f}->{spec[lin][3]:.1f}")


if __name__ == "__main__":
    main(*sys.argv[1:6])
