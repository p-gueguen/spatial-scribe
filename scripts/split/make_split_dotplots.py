"""Before/after-SPLIT marker dotplots + spillover metrics for the SPLIT track.

Loads the RCTD-surviving section counts (raw) and the SPLIT-purified counts, annotates by the
RCTD primary type, and shows that canonical markers become MORE cell-type-specific after SPLIT
(spillover removed). Also computes library-size + MECR (mutually-exclusive co-expression rate)
raw-vs-purified. Writes a self-contained dark HTML. Runs with the main pixi env.

Usage: python make_split_dotplots.py <section.h5ad> <split_inputs_dir> <split_out_dir> <report.html>
"""
from __future__ import annotations

import base64
import io as _io
import sys
import pathlib

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BG, FG, GRID = "#0a0e16", "#e8eef7", "#1d2531"
plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
                     "text.color": FG, "axes.labelcolor": FG, "xtick.color": FG, "ytick.color": FG,
                     "axes.edgecolor": GRID, "font.size": 9})

# Curated breast-panel markers per demo cell type (kept to genes likely on the 5K panel).
MARKERS = {
    "T cell": ["CD3D", "CD3E", "TRAC", "CD8A", "IL7R"],
    "Myeloid": ["LYZ", "CD68", "ITGAX", "C1QA", "CD14"],
    "B/Plasma": ["MS4A1", "CD79A", "IGHG1", "MZB1", "JCHAIN"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CD34"],
    "Stromal/CAF": ["PDGFRB", "COL1A1", "DCN", "LUM", "FAP"],
    "Epithelial/Tumor": ["EPCAM", "KRT8", "KRT18", "ERBB2"],
    "Basal/Myoepithelial": ["KRT5", "KRT14", "ACTA2", "MYLK", "TP63"],
    "Mast": ["TPSAB1", "CPA3", "MS4A2"],
    "Adipocyte": ["ADIPOQ", "LEP", "PLIN1"],
    "NK cell": ["NKG7", "GNLY", "KLRD1"],
}


def _png(fig):
    b = _io.BytesIO(); fig.savefig(b, format="png", dpi=130, bbox_inches="tight"); plt.close(fig)
    return base64.b64encode(b.getvalue()).decode()


def _lognorm(X):
    X = sp.csr_matrix(X).astype(float)
    lib = np.asarray(X.sum(1)).ravel(); lib[lib == 0] = 1
    Y = X.multiply((1e4 / lib)[:, None]).tocsr(); Y.data = np.log1p(Y.data)
    return Y


def _dot_matrix(Xln, genes_idx, groups, order):
    """mean-expression (color) + pct-expressing (size) per (group, gene)."""
    mean = np.zeros((len(order), len(genes_idx))); pct = np.zeros_like(mean)
    for i, g in enumerate(order):
        m = (groups == g)
        if m.sum() == 0:
            continue
        sub = Xln[m][:, genes_idx]
        mean[i] = np.asarray(sub.mean(0)).ravel()
        pct[i] = np.asarray((sub > 0).mean(0)).ravel()
    return mean, pct


def _dotplot(ax, mean, pct, genes, groups, title):
    mx = mean / (mean.max(0, keepdims=True) + 1e-9)  # scale each gene 0..1 across groups
    for i in range(len(groups)):
        for j in range(len(genes)):
            ax.scatter(j, i, s=8 + pct[i, j] * 180, c=[[mx[i, j], 0.15 + 0.4 * mx[i, j], 0.6]],
                       edgecolors="none")
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, rotation=90, fontsize=6)
    ax.set_yticks(range(len(groups))); ax.set_yticklabels(groups, fontsize=7)
    ax.set_title(title, color="#22d3ee", fontsize=10); ax.set_xlim(-.6, len(genes) - .4); ax.set_ylim(-.6, len(groups) - .4)
    ax.invert_yaxis()


def main(sec_h5, indir, splitdir, report):
    indir, splitdir = pathlib.Path(indir), pathlib.Path(splitdir)
    cells = pd.read_csv(indir / "cells.txt", header=None)[0].astype(str).tolist()
    genes = pd.read_csv(indir / "genes.txt", header=None)[0].astype(str).tolist()
    prim = pd.read_csv(indir / "primary.csv").set_index("cell_id")["primary"].astype(str)
    raw = sio.mmread(indir / "counts.mtx").T.tocsr()           # cells x genes
    pur = sio.mmread(splitdir / "purified_counts.mtx").T.tocsr()  # cells x genes
    groups = prim.reindex(cells).to_numpy()

    gidx = {g: i for i, g in enumerate(genes)}
    order = [t for t in MARKERS if (groups == t).sum() >= 10]
    marker_flat, marker_types = [], []
    for t in order:
        for g in MARKERS[t]:
            if g in gidx and g not in marker_flat:
                marker_flat.append(g); marker_types.append(t)
    gi = [gidx[g] for g in marker_flat]

    rawln, purln = _lognorm(raw), _lognorm(pur)
    figs = {}
    fig, axes = plt.subplots(1, 2, figsize=(1.2 + 0.28 * len(marker_flat) * 2, 0.42 * len(order) + 1.5))
    for ax, X, ttl in ((axes[0], rawln, "RAW counts"), (axes[1], purln, "SPLIT-purified")):
        m, p = _dot_matrix(X, gi, groups, order)
        _dotplot(ax, m, p, marker_flat, order, ttl)
    fig.suptitle("Canonical markers per RCTD type - raw vs SPLIT-purified", color=FG, y=1.02)
    figs["dot"] = _png(fig)

    # spillover metrics: library size + MECR (rounded-integer presence) raw vs purified
    lib_raw = np.asarray(raw.sum(1)).ravel(); lib_pur = np.asarray(pur.sum(1)).ravel()
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.hist(np.log10(lib_raw + 1), bins=40, alpha=.6, label=f"raw (med {np.median(lib_raw):.0f})", color="#f87171")
    ax.hist(np.log10(lib_pur + 1), bins=40, alpha=.6, label=f"purified (med {np.median(lib_pur):.0f})", color="#22d3ee")
    ax.set_xlabel("log10 library size"); ax.set_ylabel("cells"); ax.legend(facecolor=BG, edgecolor=GRID, labelcolor=FG)
    ax.set_title("Library size: SPLIT removes spillover transcripts", color=FG)
    figs["lib"] = _png(fig)

    # MECR: fraction of cells co-expressing markers of two mutually-exclusive lineages (T vs Epithelial)
    def mecr(X):
        Xi = X.copy(); Xi.data = np.ones_like(Xi.data)  # presence on rounded>0
        tset = [gidx[g] for g in MARKERS["T cell"] if g in gidx]
        eset = [gidx[g] for g in MARKERS["Epithelial/Tumor"] if g in gidx]
        tpos = np.asarray(Xi[:, tset].sum(1)).ravel() > 0
        epos = np.asarray(Xi[:, eset].sum(1)).ravel() > 0
        return float((tpos & epos).mean())
    mecr_raw = mecr(sp.csr_matrix(np.rint(raw.toarray())) if raw.shape[0] < 20000 else raw)
    mecr_pur = mecr(sp.csr_matrix(np.rint(pur.toarray())) if pur.shape[0] < 20000 else pur)

    cm = pd.read_csv(splitdir / "cell_meta.csv") if (splitdir / "cell_meta.csv").exists() else None
    status = (cm["purification_status"].value_counts().to_dict()
              if cm is not None and "purification_status" in cm else {})

    def img(k, cap):
        return f'<figure><img src="data:image/png;base64,{figs[k]}"/><figcaption>{cap}</figcaption></figure>' if k in figs else ""
    html = f"""<!doctype html><meta charset=utf-8><title>SpatialScribe - SPLIT purification</title>
<style>body{{background:{BG};color:{FG};font-family:-apple-system,Segoe UI,Roboto,sans-serif;margin:0;padding:2rem 2.4rem;max-width:1100px}}
h1{{font-weight:650;margin:0 0 .2rem}} h2{{color:#22d3ee;font-size:1.05rem;margin:1.8rem 0 .5rem;border-bottom:1px solid {GRID};padding-bottom:.3rem}}
.sub{{color:#8b98ab}} figure{{margin:1rem 0;background:#0d1420;border:1px solid {GRID};border-radius:10px;padding:1rem;display:inline-block}}
img{{max-width:100%;display:block}} figcaption{{color:#8b98ab;font-size:.8rem;margin-top:.5rem}} b{{color:#22d3ee}}
table{{border-collapse:collapse;font-size:.9rem}} td,th{{border:1px solid {GRID};padding:.35rem .8rem}}</style>
<h1>SpatialScribe &middot; SPLIT residual-contamination purification</h1>
<p class=sub>RCTD (rctd-py) full-mode weights &rarr; <b>SPLIT::rctd_free_purify</b> with <b>DO_remove_residual_contamination=TRUE</b> on {len(cells):,} RCTD-surviving cells. SPLIT reassigns spillover transcripts to their true source cell using the deconvolution weights + per-type reference profile.</p>
<p><b>Spillover removed:</b> median library size {np.median(lib_raw):.0f} &rarr; {np.median(lib_pur):.0f} ({100*(1-np.median(lib_pur)/np.median(lib_raw)):.1f}% fewer transcripts/cell). <b>T-cell &times; Epithelial co-expression (MECR):</b> {mecr_raw*100:.1f}% &rarr; {mecr_pur*100:.1f}% (lower = cleaner). purification_status: {status}</p>
<h2>Marker specificity (raw vs SPLIT-purified)</h2>{img('dot','Dot color = mean expression (per-gene scaled), size = % of cells expressing. After SPLIT, each type&#39;s markers concentrate in that type and off-target spillover shrinks.')}
<h2>Spillover removal</h2>{img('lib','Per-cell library size before/after SPLIT.')}
"""
    pathlib.Path(report).write_text(html)
    print("wrote", report, "| MECR raw->pur:", round(mecr_raw, 3), "->", round(mecr_pur, 3),
          "| lib med:", int(np.median(lib_raw)), "->", int(np.median(lib_pur)))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
