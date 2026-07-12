// Panel check (selling point): how well the probe panel + reference markers resolve each cell type.
// Global concordance metrics, per-type coverage + one-vs-rest identifiability AUC, confusable pairs,
// a merge nudge, and an automatic plain-language LLM verdict (Haiku - cheap, prompt-cached).
import React, { useEffect, useState } from "react";
import * as api from "../api";
import { capKey, getCached, setCached } from "../capCache";
import { Help, Markdown } from "./viz";
import { ReferenceSection } from "./LoadPanel";

export default function PanelCheckPanel({ sid, onRan, sess, onMutate, onColor }:
  { summary?: any; sid?: string; onRan?: () => void; sess?: api.Session | null; onMutate?: () => Promise<void>; onColor?: (f: string) => void }) {
  const [rep, setRep] = useState<any>(null);
  const [verdict, setVerdict] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);
  const [vbusy, setVbusy] = useState(false);

  const load = async (withVerdict = true) => {
    if (!sid) return;
    setBusy(true);
    try { const r = await api.panelReport(sid); setRep(r); setCached(capKey(sid, "panel_report"), r); onRan?.(); } catch { /* keep */ }
    setBusy(false);
    if (withVerdict) {
      setVbusy(true);
      try { const v = await api.panelVerdict(sid); const vv = v.verdict || v.note || null; setVerdict(vv); setCached(capKey(sid, "panel_verdict"), vv); } catch { /* keep */ }
      setVbusy(false);
    }
  };
  // Serve the report + verdict from the session cache so switching tabs and coming back does not
  // recompute; capCache.clearSession() drops these on any section mutation (annotate/rename/...), so
  // a cached hit is never stale. Only compute when nothing is cached yet.
  useEffect(() => {
    if (!sid) return;
    const cached = getCached(capKey(sid, "panel_report"));
    if (cached) { setRep(cached); setVerdict(getCached<string>(capKey(sid, "panel_verdict")) ?? null); return; }
    // Auto-fetch the LLM verdict only in a PROD build (import.meta.env.PROD); dev/test skips it to save
    // the API call - the "refresh" button (load(true)) still fetches it on demand.
    setRep(null); setVerdict(null); load(import.meta.env.PROD);
    /* eslint-disable-next-line */
  }, [sid]);

  if (!rep) return <div className="panel"><h4>Panel resolvability {busy && <span className="running"><i />checking</span>}</h4><div className="pmuted">{busy ? "scoring marker coverage, identifiability + confusable pairs …" : "no panel report"}</div></div>;

  const g = rep.global || {};
  const cov = Math.round((g.overall_marker_coverage || 0) * 100);
  const auc = g.mean_identifiability_auc;
  const basisLabel = (b?: string) => b === "depth_matched_f1" ? "depth-matched F1 (reference)"
    : b === "identifiability_auc" ? "identifiability AUC" : "marker coverage (weak)";
  // Hover a marker-coverage chip to see WHICH cell types it counts (saves the vertical space of a list).
  // undefined (older/cached report without the lists) -> no tooltip; [] -> "none"; else the type names.
  const covTip = (types?: string[]) => Array.isArray(types) ? (types.length ? types.join(", ") : "none") : undefined;

  return (
    <div className="panel">
      <h4>Panel ↔ reference resolvability <Help text="How well this probe panel + its reference markers can resolve each cell type. Marker presence is necessary, not sufficient." />{busy && <span className="running"><i />checking</span>}</h4>

      <div className="metricgrid">
        <div className="metric" title="One-vs-rest AUC of each type's on-panel markers vs the section's own cell-type labels - computed once cell types exist (run Annotate). Depth-matched F1 supersedes it when a reference is loaded.">
          <div className="k">mean identifiability</div>
          <div className="v" style={auc != null ? { color: auc >= 0.75 ? "var(--pass)" : auc >= 0.6 ? "var(--warn)" : "var(--fail)" } : undefined}>
            {auc != null ? auc.toFixed(2)
            : <span className="pmuted" style={{ fontSize: ".62rem", fontWeight: 400 }}>after annotate</span>}</div>
        </div>
        <div className="metric"><div className="k">marker coverage</div><div className="v">{cov}%</div></div>
        {g.n_typability_assessed > 0
          ? <div className="metric" title={"decided from " + basisLabel(g.typability_basis)}><div className="k">confidently typable</div><div className="v">{g.n_confidently_typable}/{g.n_typability_assessed}</div></div>
          : <div className="metric"><div className="k">types resolved</div><div className="v">{g.n_resolvable}/{g.n_types}</div></div>}
      </div>
      {g.n_typability_assessed > 0 && (
        <div className="pmuted" style={{ margin: "2px 0 6px" }}>
          Typability decided from <b>{basisLabel(g.typability_basis)}</b> - a computed discriminability score, not marker presence.
        </div>
      )}
      {/* Marker COVERAGE (are the type's markers on the panel), a different axis from typability
          (can the panel actually discriminate the type). Labelled explicitly so the two are not
          conflated - a type can be well-covered yet not typable (few markers, low identifiability). */}
      <div className="flagrow">
        <span className="pmuted" style={{ fontSize: 11 }}>marker coverage:</span>
        <span className="flag" title={covTip(g.resolvable_types)} style={covTip(g.resolvable_types) ? { cursor: "help" } : undefined}><span className="tl green" />{g.n_resolvable} good</span>
        <span className="flag" title={covTip(g.weak_types)} style={covTip(g.weak_types) ? { cursor: "help" } : undefined}><span className="tl amber" />{g.n_weak} weak</span>
        <span className="flag" title={covTip(g.cannot_types)} style={covTip(g.cannot_types) ? { cursor: "help" } : undefined}><span className="tl red" />{g.n_cannot} none</span>
        {g.n_confusable_pairs > 0 && <span className="flag">{g.n_confusable_pairs} confusable pair(s)</span>}
      </div>

      {rep.merge_nudge && (
        <div className={"callout" + (String(rep.merge_nudge).startsWith("No panel") ? "" : " warn")}>{rep.merge_nudge}</div>
      )}

      <div className="subhead">per cell type</div>
      <table className="dtable">
        <thead><tr><th>type</th><th>markers</th><th>identifiability (AUC)</th><th>typable?</th></tr></thead>
        <tbody>
          {rep.per_type.map((r: any) => (
            <tr key={r.cell_type} title={r.missing?.length ? "missing: " + r.missing.join(", ") : "all canonical markers on panel"}>
              {/* No status dot on the name: the markers count IS the coverage, and the only colour-coded
                  verdict in the row is "typable?". */}
              <td>{r.cell_type}</td>
              <td>{r.n_present}/{r.n_markers}</td>
              <td>
                {r.auc == null ? <span className="pmuted">-</span> : (
                  <span style={{ display: "flex", alignItems: "center", gap: 6 }}>
                    {/* Neutral magnitude bar (not a red->green traffic light) so the AUC reads as a value,
                        not a second pass/fail signal competing with the typable dot. */}
                    <span className="aucbar neutral" style={{ flex: 1 }}><i style={{ width: `${Math.round(r.auc * 100)}%` }} /></span>
                    <span style={{ minWidth: 28 }}>{r.auc.toFixed(2)}</span>
                  </span>
                )}
              </td>
              <td title={(r.not_typable_reason ? r.not_typable_reason + " · " : "") + (r.typability_basis ? `${basisLabel(r.typability_basis)}${r.typability_score != null ? " = " + r.typability_score : ""}` : "")}>
                {r.confidently_typable == null
                  ? <span className="pmuted">-</span>
                  : (() => {
                      // Small-n ("too few cells") is NOT a hard fail: the cells may well express their
                      // markers (they do for the atera NK/Adipocyte in the dot-plot) - we just cannot
                      // CERTIFY typability on <50 cells, where the one-vs-rest AUC is noise in either
                      // direction. Show it AMBER "unsure", not red "no" (a genuine low-AUC / no-marker
                      // fail stays red). The headline "confidently typable" count is unchanged.
                      const smallN = r.confidently_typable === false
                        && String(r.not_typable_reason || "").startsWith("too few cells");
                      return <><span className={"tl " + (r.confidently_typable ? "green" : smallN ? "amber" : "red")} />
                        {r.confidently_typable ? "yes" : smallN ? "unsure" : "no"}</>;
                    })()}
                {/* Surface the reason inline ONLY for the COUNTERINTUITIVE 'no's - a good AUC / marker
                    count that is still 'no' (a small-n discount like "too few cells (n=11 < 50)", or an
                    F1 confusion). The low-AUC case ("identifiability below X") is self-evident from the
                    AUC column, so it stays hover-only - that is the caption we removed to declutter. */}
                {r.confidently_typable === false && r.not_typable_reason &&
                  !String(r.not_typable_reason).startsWith("identifiability below") &&
                  <span className="pmuted" style={{ fontSize: 10, display: "block", lineHeight: 1.15 }}>{r.not_typable_reason}</span>}
              </td>
            </tr>
          ))}
        </tbody>
      </table>

      {rep.confusable?.length > 0 && (
        <>
          <div className="subhead">confusable (no private marker)</div>
          <div className="pills">{rep.confusable.map((p: any, i: number) =>
            <span className="pill2" key={i}>{Array.isArray(p) ? p.join(" ~ ") : String(p)}</span>)}</div>
        </>
      )}

      <div className="subhead">reference (sharpens resolvability)
        <Help text="Add a single-cell reference to upgrade the typability verdict from marker-coverage to depth-matched per-class F1, and to enable supervised label transfer in the Annotate step. The reference is shared across the whole session." />
      </div>
      <ReferenceSection sid={sid} sess={sess} onMutate={onMutate} onColor={onColor} showTransfer={false} />

      <div className="subhead">
        AI verdict
        <button className="pbtn sm" style={{ marginLeft: "auto" }} disabled={vbusy} onClick={() => load(true)}>{vbusy ? "…" : "refresh"}</button>
      </div>
      {vbusy && !verdict ? <div className="pmuted">asking Claude for a plain-language verdict …</div>
        : verdict ? <Markdown text={verdict} /> : <div className="pmuted">No AI verdict yet - click refresh.</div>}
    </div>
  );
}
