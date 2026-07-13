// QC (selling point): section metrics, the six-layer annotation-QC funnel, and region QC - box-select a
// region on the map to QC it, then exclude it or crop to it (region filtering).
import React, { useEffect, useState } from "react";
import * as api from "../api";
import { capKey, getCached, setCached } from "../capCache";
import { Markdown, Help } from "./viz";

export default function QCPanel({ summary, sid, selection, onMutate, sess, onRan }:
  { summary?: any; sid?: string; selection?: number[]; onMutate?: () => Promise<void>; sess?: api.State; onRan?: () => void }) {
  const q = summary?.qc;
  const [funnel, setFunnel] = useState<any>(null);
  const [fbusy, setFbusy] = useState(false);
  const [region, setRegion] = useState<any>(null);
  const [rbusy, setRbusy] = useState(false);
  const [verdict, setVerdict] = useState<string | null>(null);
  const [vbusy, setVbusy] = useState(false);

  const runFunnel = async (withVerdict = true) => {
    if (!sid) return; setFbusy(true);
    try { const r = await api.runCap(sid, "qc_funnel"); if (r.ok) { setFunnel(r.value); setCached(capKey(sid, "qc_funnel"), r.value); onRan?.(); } } catch { /* keep */ }
    setFbusy(false);
    // Plain-language QC verdict (LLM). Auto-fetch only in a PROD build; dev/test skips it to save the
    // API call - the "refresh" button still fetches it on demand.
    if (withVerdict) {
      setVbusy(true);
      try { const v = await api.qcVerdict(sid!); const vv = v.verdict || v.note || null; setVerdict(vv); setCached(capKey(sid!, "qc_verdict"), vv); } catch { /* keep */ }
      setVbusy(false);
    }
  };
  // Auto-run the funnel on FIRST open, then serve it from the session cache so re-entering the QC tab
  // does not recompute. capCache.clearSession() drops these on any section mutation (region filter,
  // annotate, ...), so a cached hit is never stale.
  useEffect(() => {
    if (!sid) return;
    const cached = getCached(capKey(sid, "qc_funnel"));
    if (cached) {
      setFunnel(cached);
      setVerdict(getCached<string>(capKey(sid, "qc_verdict")) ?? null);
      setRegion(null);
      return;
    }
    setFunnel(null); setRegion(null); setVerdict(null); runFunnel(import.meta.env.PROD);
    /* eslint-disable-next-line */
  }, [sid]);

  const doRegionQc = async () => {
    if (!sid || !selection?.length) return; setRbusy(true);
    try { setRegion(await api.regionQc(sid, selection)); } catch { /* keep */ }
    setRbusy(false);
  };
  const doFilter = async (mode: "exclude" | "keep") => {
    if (!sid || !selection?.length) return; setRbusy(true);
    try { await api.regionFilter(sid, selection, mode); setRegion(null); if (onMutate) await onMutate(); } catch { /* keep */ }
    setRbusy(false);
  };

  const flags = funnel?.section?.flags || {};
  const cf = funnel?.count_floor || {};
  const seg = funnel?.segmentation || {};

  // Quick-scan QC health: colour draws the eye straight to the metrics that need a second look -
  // shallow depth (amber/red) or too many empties vs. clean signal (green). Heuristic Xenium
  // thresholds; the AI verdict at the bottom is the authoritative read, this is a fast triage cue.
  const TONE: Record<string, string> = { pass: "var(--pass)", warn: "var(--warn)", fail: "var(--fail)" };
  const TL: Record<string, string> = { pass: "green", warn: "amber", fail: "red" };
  const depthTone = (v?: number, warn = 100, fail = 30) => v == null ? undefined : v < fail ? "fail" : v < warn ? "warn" : "pass";
  const emptyTone = (f?: number) => f == null ? undefined : f > 0.15 ? "fail" : f > 0.05 ? "warn" : "pass";
  const metric = (k: string, v: React.ReactNode, t?: string) => (
    <div className="metric">
      <div className="k" style={{ display: "flex", alignItems: "center", justifyContent: "space-between" }}>
        <span>{k}</span>{t && <span className={"tl " + TL[t]} style={{ marginRight: 0 }} />}
      </div>
      <div className="v" style={{ color: t === "warn" || t === "fail" ? TONE[t] : undefined }}>{v}</div>
    </div>
  );

  return (
    <div className="panel">
      <h4>Signal you can trust</h4>
      {q ? (
        <div className="metricgrid">
          {metric("cells", Number(q.n_cells).toLocaleString())}
          {metric("median genes", q.median_genes, depthTone(q.median_genes, 50, 20))}
          {metric("median counts", q.median_counts, depthTone(q.median_counts, 100, 30))}
          {metric("empty cells", `${(q.empty_frac * 100).toFixed(1)}%`, emptyTone(q.empty_frac))}
          {/* neg-control %: the driver behind the "pct counts control" flag chip. Shown so a red chip
              is never an unexplained flag while the four other metrics read green (2% warn / 5% fail). */}
          {q.pct_counts_control != null && metric("neg. control %", `${q.pct_counts_control.toFixed(1)}%`,
            q.pct_counts_control > 5 ? "fail" : q.pct_counts_control > 2 ? "warn" : "pass")}
        </div>
      ) : <div className="pmuted">Run QC to compute per-cell metrics.</div>}

      <div className="subhead">six-layer QC funnel {fbusy && <span className="running"><i />running</span>}</div>
      {!funnel ? <div className="pmuted">{fbusy ? "scoring segmentation, count floor, purity, coherence …" : "not run"}</div> : (
        <>
          {/* Only surface metrics that need attention. The four green "all ok" pills just duplicated
              the metric cards above; showing nothing on a clean section keeps the tab uncluttered. */}
          {Object.entries(flags).filter(([, v]) => v !== "ok").length > 0 && (
            <div className="flagrow">
              {Object.entries(flags).filter(([, v]) => v !== "ok").map(([k, v]) => (
                <span className="flag" key={k}><span className={"tl " + (v === "warn" ? "amber" : "red")} />{String(k).replace(/_/g, " ")}</span>
              ))}
            </div>
          )}
          <div className="cmp">
            <div className="cmpcell"><div className="k">count floor</div><div className="v">{cf.floor != null ? Math.round(cf.floor) : "-"}</div>
              <div className="d">{cf.mode || ""}{cf.pct_removed != null && <> · <span style={{ color: cf.pct_removed > 0.1 ? "var(--warn)" : undefined }}>{(cf.pct_removed * 100).toFixed(1)}% removed</span></>}</div></div>
            {(() => {
              const st = seg.status ? String(seg.status) : "";
              const t = st === "ok" ? "pass" : (!st || st.startsWith("skip")) ? null : "warn";
              return (
                <div className="cmpcell"><div className="k">segmentation</div>
                  <div className="v" style={{ fontSize: ".95rem", color: t ? TONE[t] : "var(--muted)" }}>
                    {t && <span className={"tl " + TL[t]} />}{seg.status || "-"}</div>
                  <div className="d">{seg.pct_area_outlier != null ? `${(seg.pct_area_outlier * 100).toFixed(1)}% area outliers` : ""}</div></div>
              );
            })()}
          </div>
          {/* The label-dependent layers (annotation quality, spatial coherence, per-cell confidence)
              score the CURRENT cell-type labels, so they live in the Annotate tab now - not here in QC,
              which runs on the raw section BEFORE any annotation exists. (They used to appear here only
              because the bundled demos ship pre-annotated.) */}
        </>
      )}

      {/* "why weren't some cells confidently typed?" lives in the Annotate tab, next to the labels it
          describes - it is annotation-dependent (rejection reasons only exist once cells are typed), so
          showing it here in QC (which reads the RAW signal before annotation) just duplicated it. */}

      <div className="subhead">
        AI verdict {vbusy && <span className="running"><i />asking Claude</span>}
        <button className="pbtn sm" style={{ marginLeft: "auto" }} disabled={vbusy || !sid}
                onClick={async () => { setVbusy(true); try { const v = await api.qcVerdict(sid!); setVerdict(v.verdict || v.note || null); } catch { /* keep */ } setVbusy(false); }}>{vbusy ? "…" : "refresh"}</button>
      </div>
      {vbusy && !verdict ? <div className="pmuted">asking Claude whether this section's signal is trustworthy …</div>
        : verdict ? <Markdown text={verdict} /> : <div className="pmuted">No AI verdict yet - click refresh.</div>}

      <div className="subhead">region QC (H1)</div>
      {!selection?.length ? (
        <div className="pmuted">Click <b>select region</b> on the map and drag a box, then QC or filter it here.</div>
      ) : (
        <>
          <div className="pmuted" style={{ marginBottom: 4 }}>{selection.length.toLocaleString()} cells selected.</div>
          <div className="pills" style={{ marginBottom: 6 }}>
            <button className="pbtn pri sm" disabled={rbusy} onClick={doRegionQc}>QC region</button>
            <button className="pbtn sm" disabled={rbusy} onClick={() => doFilter("exclude")}>Exclude</button>
            <button className="pbtn sm" disabled={rbusy} onClick={() => doFilter("keep")}>Keep only</button>
          </div>
          {region && (
            <div className="cmp">
              <div className="cmpcell"><div className="k">region genes</div><div className="v">{region.region?.median_genes_per_cell ?? "-"}</div><div className="d">section {region.section?.median_genes_per_cell ?? "-"}</div></div>
              <div className="cmpcell"><div className="k">region counts</div><div className="v">{region.region?.median_transcripts_per_cell ?? "-"}</div><div className="d">section {region.section?.median_transcripts_per_cell ?? "-"}</div></div>
            </div>
          )}
        </>
      )}
    </div>
  );
}
