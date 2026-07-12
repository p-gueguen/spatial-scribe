// Specimen overview shown on the Load step. Also keeps the panel column reserved on every
// step so the canvas never shifts sideways. Includes a "load your own section" control that
// ingests a Xenium / CosMx / MERSCOPE output folder or a .h5ad directly from disk (server-side),
// and a "custom single-cell reference" control (upload a .h5ad or paste a server path) that scores
// how well the panel can transfer that reference's cell types and can transfer labels onto the section.
import React, { useEffect, useRef, useState } from "react";
import * as api from "../api";
import type { Session } from "../api";
import { capKey, setCached } from "../capCache";
import { getRecent, clearRecent, shortPath, type RecentSection } from "../recentSections";
import { Help } from "./viz";

const TISSUES = ["melanoma", "breast", "lung", "skin", "colon", "generic"];

// Friendly labels for the pipeline stages surfaced during "Run full analysis". The progress `label` is
// the raw capability name; `stage` is the constant "pipeline" (the job name), so it must NOT be shown.
const STAGE_LABEL: Record<string, string> = {
  starting: "starting", compute_qc: "quality control", panel_check: "panel adequacy",
  reference_match: "reference match", annotation_strategy: "annotation strategy",
  cluster: "clustering", annotate: "annotating", reference_transfer: "annotating",
  annotate_rctd: "annotating (RCTD)", niches: "spatial niches",
  malignant_concordance: "malignant calling", split_purify: "spillover purification",
  self_heal: "self-verify",
};

// After a full run, warm the diagnostic panels so opening Panel-check / QC is INSTANT: the pipeline ran
// the base caps, but the panels display the richer reports + LLM verdicts, and onMutate's invalidate
// clears any cache. Best-effort, non-blocking; the panels read these same cache keys on mount. The
// entry each panel CHECKS first (panel_report / qc_funnel) is cached LAST, after its dependent verdict
// and rejection entries, so a visit mid-prefetch either finds the whole set or falls back to a full
// compute - never a report with a missing verdict.
async function prefetchReports(sid: string) {
  // The two plain-language verdicts cost an LLM call each. Prewarm them ONLY in a production build
  // (import.meta.env.PROD: true for `vite build` = the prod bundle, false for `vite dev` =
  // the dev/test deploy), so testing does not burn the metered vLLM. On dev the panels fetch the
  // verdict on demand (their "refresh" button). The reports/funnel below are compute-only, always warm.
  const prod = import.meta.env.PROD;
  // panel: verdict first (prod only), then the report (the entry the panel checks) LAST - so a visit
  // mid-prefetch either finds the whole set or falls back to a full compute, never a half-warmed one.
  if (prod) { try { const pv = await api.panelVerdict(sid); setCached(capKey(sid, "panel_verdict"), pv.verdict || pv.note || null); } catch { /* best-effort */ } }
  try { const r = await api.panelReport(sid); setCached(capKey(sid, "panel_report"), r); } catch { /* best-effort */ }
  // qc: rejection + verdict first, then the funnel (the entry the panel checks) LAST.
  try { const rr = await api.runCap(sid, "rejection_reasons"); if (rr.ok) setCached(capKey(sid, "rejection_reasons"), rr.value?.breakdown || []); } catch { /* best-effort */ }
  if (prod) { try { const qv = await api.qcVerdict(sid); setCached(capKey(sid, "qc_verdict"), qv.verdict || qv.note || null); } catch { /* best-effort */ } }
  try { const f = await api.runCap(sid, "qc_funnel"); if (f.ok) setCached(capKey(sid, "qc_funnel"), f.value); } catch { /* best-effort */ }
}

// The "does this sample contain malignant cells?" flag - a DATA-step config (it gates the malignant
// callers), so it stays on the Load tab.
function TumourFlag({ sid, sess, onMutate }:
  { sid?: string; sess?: Session | null; onMutate?: () => Promise<void> }) {
  const isTumour = sess?.is_tumour ?? null;
  const toggleTumour = async (checked: boolean) => {
    if (!sid) return;
    await api.setIsTumour(sid, checked);
    if (onMutate) await onMutate();
  };
  if (!sid) return null;
  return (
    <>
      <label className="checkline" style={{ marginTop: 10, marginBottom: 6 }}>
        <input type="checkbox" checked={isTumour === true} onChange={(e) => toggleTumour(e.target.checked)} />
        <span>This sample contains malignant cells
          <Help text="Runs the malignant marker-score + Cancer-Finder callers at the Annotate step (InSituCNV is not run interactively - its infercnvpy pass is minutes-long). Leave unticked for normal tissue. When you do not answer, the app guesses from the tissue name, which wrongly fires on 'normal breast' and wrongly skips 'glioblastoma'." /></span>
      </label>
      {isTumour === null && (
        <div className="pmuted" style={{ fontSize: 11, marginBottom: 6 }}>
          Not answered: malignant calling will be guessed from the tissue name ({sess?.tissue || "unset"}).
        </div>
      )}
    </>
  );
}

// The SUPERVISED annotation route (reference transfer). Exported so it renders at the top of the
// Annotate tab - transferring a reference's labels IS annotation, so it belongs with the other
// annotation methods, not on the Load step. Load only attaches the reference as an input is gone;
// the picker + match + transfer all live here now.
export function ReferenceSection({ sid, sess, onMutate, onColor, showTransfer = true }:
  { sid?: string; sess?: Session | null; onMutate?: () => Promise<void>; onColor?: (f: string) => void; showTransfer?: boolean }) {
  const [path, setPath] = useState("");
  const [busy, setBusy] = useState(false);
  const [res, setRes] = useState<api.RefResult | null>(null);
  const [msg, setMsg] = useState<string | null>(null);
  const [xfer, setXfer] = useState<any>(null);
  const [autoBusy, setAutoBusy] = useState(false);
  const [autoRes, setAutoRes] = useState<api.AutoRefResult | null>(null);
  const [allowFetch, setAllowFetch] = useState(false);
  const [refTissue, setRefTissue] = useState("");   // editable tissue for auto-select (defaults to the section's)
  const fileRef = useRef<HTMLInputElement>(null);

  // Free-text tissue -> auto-choose + load the best-matched pre-computed reference (or a live gget
  // fetch). Reflects the returned panel match into `res` so the match readout below renders as usual.
  const doAuto = async () => {
    if (!sid) return;
    setAutoBusy(true); setMsg(null); setXfer(null);
    try {
      const r = await api.autoReference(sid, (refTissue || sess?.tissue || "").trim(), allowFetch);
      setAutoRes(r);
      if (r.ok && r.match) {
        setRes({ ok: true, n_ref_cells: r.n_ref_cells || 0, label_key: r.label_key || null,
                 n_labels: null, match: r.match, obs_fields: r.obs_fields || [], done: r.done || {} });
        if (onMutate) await onMutate();
      }
    } catch (e: any) {
      setAutoRes({ ok: false, auto: { status: "error", chosen: null, message: String(e?.message || e) } });
    }
    setAutoBusy(false);
  };

  const register = async (fn: () => Promise<api.RefResult>) => {
    if (!sid) return;
    setBusy(true); setMsg(null); setXfer(null);
    try { const r = await fn(); setRes(r); if (onMutate) await onMutate(); }
    catch (e: any) { setMsg("could not load reference: " + (e?.message || e)); }
    setBusy(false);
  };
  const onUpload = (f?: File) => f && register(() => api.uploadReference(sid!, f));
  const onPath = () => path.trim() && register(() => api.setReference(sid!, path.trim()));
  const transfer = async (setAsPrimary: boolean) => {
    if (!sid) return;
    setBusy(true); setXfer(null);
    try {
      const r = await api.runCap(sid, "reference_transfer", { set_as_primary: setAsPrimary });
      if (r.ok) { setXfer(r.value); if (onMutate) await onMutate(); if (onColor) onColor(setAsPrimary ? "cell_type" : "celltypist_label"); }
      else setXfer({ error: r.error });
    } catch (e: any) { setXfer({ error: String(e) }); }
    setBusy(false);
  };

  const m = res?.match;
  const g = m?.global || {};
  const notRes: string[] = m?.not_resolvable || [];
  const verdictClass = g.verdict === "good" ? "tl green" : g.verdict === "fair" ? "tl amber" : "tl red";

  return (
    <>
      <h4 style={{ marginBottom: 6 }}>{showTransfer ? "Annotate from a single-cell reference" : "Add a single-cell reference"}
        <Help text={showTransfer
          ? "The SUPERVISED annotation route: bring your own scRNA / snRNA reference .h5ad (label column auto-detected), or auto-select one for the tissue. It scores how well this panel can transfer the reference's cell types, then transfers those labels onto the section. The alternative is to annotate from the clusters (below)."
          : "Bring your own scRNA / snRNA reference .h5ad (or auto-select one for the tissue) to SHARPEN panel resolvability - the typability verdict upgrades from marker-coverage to depth-matched per-class F1 - and to enable supervised label transfer in the Annotate step."} />
      </h4>
      {!sid && <div className="pmuted">Load a section first.</div>}
      {sid && (
        <>
          {/* Auto-select from a free-text tissue - editable, defaults to the loaded section's tissue. */}
          <div className="field" style={{ marginBottom: 6 }}>
            <label>tissue for auto-select <Help text="Edit to pick a reference for a different tissue than the section was loaded as (e.g. type 'liver'). Auto-select finds the best-matched pre-computed reference for this tissue." /></label>
            <input value={refTissue} placeholder={sess?.tissue || "e.g. liver"}
                   onChange={(e) => setRefTissue(e.target.value)}
                   onKeyDown={(e) => { if (e.key === "Enter" && !autoBusy && !busy) doAuto(); }} />
          </div>
          <label className="checkline" style={{ marginBottom: 6 }}>
            <input type="checkbox" checked={allowFetch} onChange={(e) => setAllowFetch(e.target.checked)} />
            <span>Fetch a CELLxGENE reference if nothing local fits
              <Help text="If no local reference matches the tissue, fetch a small CELLxGENE reference live via gget (needs a gget env + network)." /></span>
          </label>
          <button className="pbtn" style={{ width: "100%", marginBottom: autoBusy ? 4 : 8 }} disabled={autoBusy || busy}
                  onClick={doAuto}>{autoBusy ? "selecting reference …" : "Auto-select reference for this tissue"}</button>
          {/* A large atlas (e.g. the 621k-cell breast reference) is minutes to LOAD off disk the first
              time - the button alone reads as frozen, so say so honestly. */}
          {autoBusy && (
            <div className="running" style={{ width: "100%", justifyContent: "center", marginBottom: 8 }}>
              <i />loading the reference - a large atlas can take 1-2 min
            </div>
          )}
          {autoRes && (
            <div className={"callout" + (autoRes.ok ? "" : " warn")} style={{ marginBottom: 8 }}>
              {autoRes.ok
                ? <><b>Selected the '{autoRes.auto.chosen}' reference</b> ({(autoRes.n_ref_cells || 0).toLocaleString()} cells).
                    {autoRes.recommended_mode && <> Recommended: <b>{autoRes.recommended_mode === "cluster"
                      ? "cluster (reference too weak for this panel)" : "supervised transfer"}</b>.</>}</>
                : <>{autoRes.auto.message || "No suitable local reference; upload one or cluster instead."}</>}
            </div>
          )}
          <input ref={fileRef} type="file" accept=".h5ad" style={{ display: "none" }}
                 onChange={(e) => onUpload(e.target.files?.[0])} />
          <button className="pbtn pri" style={{ width: "100%" }} disabled={busy}
                  onClick={() => fileRef.current?.click()}>{busy ? "loading reference …" : "Upload reference .h5ad"}</button>
          <div className="field" style={{ marginTop: 8 }}><label>…or a server-side path</label>
            <input value={path} placeholder="/path/to/reference.h5ad" onChange={(e) => setPath(e.target.value)}
                   onKeyDown={(e) => { if (e.key === "Enter") onPath(); }} /></div>
          {path.trim() && <button className="pbtn" style={{ width: "100%" }} disabled={busy} onClick={onPath}>Register path</button>}
          {msg && <div className="err">{msg}</div>}
        </>
      )}

      {res && (
        <>
          <div className="callout" style={{ marginTop: 10 }}>
            Reference loaded: <b>{res.n_ref_cells.toLocaleString()}</b> cells · label
            column <b>{res.label_key}</b>{res.n_labels ? <> · <b>{res.n_labels}</b> types</> : null}.
          </div>
          {m && m.status === "ok" && (
            <>
              <div className="subhead">reference &harr; panel match</div>
              <div className="cmp">
                <div className="cmpcell"><div className="k">match score</div>
                  <div className="v"><span className={verdictClass} />{m.match_score != null ? m.match_score.toFixed(2) : "-"}</div>
                  <div className="d">{g.verdict} - single 0-1 reference-to-panel fit</div></div>
                <div className="cmpcell"><div className="k">panel-gene overlap</div>
                  <div className="v">{Math.round((g.gene_overlap_frac || 0) * 100)}%</div>
                  <div className="d">{g.n_shared_genes}/{g.n_panel_genes} genes</div></div>
                <div className="cmpcell"><div className="k">resolvable types</div>
                  <div className="v">{g.n_resolvable}/{g.n_types}</div>
                  <div className="d">mean per-type F1 {g.mean_f1 ?? "-"}</div></div>
              </div>
              {notRes.length > 0 && (
                <div className="pmuted" style={{ marginTop: 6 }}>Panel cannot resolve: {notRes.slice(0, 8).join(", ")}{notRes.length > 8 ? " …" : ""}</div>
              )}
              {m.recommendation?.headline && (
                <div className={"callout" + (m.recommendation.action === "supervised_transfer" ? "" : " warn")}>
                  <b>Recommendation:</b> {m.recommendation.headline}
                  <div className="pmuted" style={{ marginTop: 3 }}>{m.recommendation.why}</div>
                </div>
              )}
              {(m.alternatives || []).length > 0 && (
                <div className="pmuted" style={{ marginTop: 6 }}>
                  Better-matched references for <b>{sess?.tissue}</b>: {(m.alternatives || []).map((r: any) => r.tissue_key).join(", ")}.
                </div>
              )}
            </>
          )}
          {m && m.status !== "ok" && (
            <div className="callout warn" style={{ marginTop: 8 }}>Match not computed: {m.reference_source || m.status}</div>
          )}

          {showTransfer && (<>
          <div className="subhead">transfer labels from the reference <Help text={'Trains a CellTypist model on the reference (in-env) and predicts per-cell labels; also runs TACCO OT when installed. "Adopt" overwrites cell_type with the transferred labels.'} /></div>
          <div style={{ display: "flex", gap: 8 }}>
            <button className="pbtn" style={{ flex: 1 }} disabled={busy} onClick={() => transfer(false)}>{busy ? "transferring …" : "Transfer (advisory)"}</button>
            <button className="pbtn pri" style={{ flex: 1 }} disabled={busy} onClick={() => transfer(true)}>Transfer &amp; adopt</button>
          </div>
          {xfer?.error && <div className="err">{String(xfer.error?.hint || xfer.error?.message || xfer.error)}</div>}
          {xfer && !xfer.error && (
            <table className="dtable" style={{ marginTop: 6 }}>
              <thead><tr><th>arm</th><th>status</th></tr></thead>
              <tbody>{(xfer.arms || []).map((a: any, i: number) => (
                <tr key={i}><td>{a.method}</td><td>{a.status === "ok" ? `ok${a.n_labels ? ` · ${a.n_labels} types` : ""}` : a.status}</td></tr>
              ))}
              {xfer.agreement_with_cell_type != null && (
                <tr><td>agreement vs cell_type</td><td>{Math.round(xfer.agreement_with_cell_type * 100)}%</td></tr>
              )}</tbody>
            </table>
          )}
          </>)}
        </>
      )}
    </>
  );
}

export default function LoadPanel({
  sess, sid, onMutate, onColor, onLoad, onLoadSynthetic, onLoadDemo, busy,
}: {
  sess: Session | null; sid?: string; onMutate?: () => Promise<void>; onColor?: (f: string) => void;
  onLoad?: (path: string, tissue: string) => void; onLoadSynthetic?: () => void;
  onLoadDemo?: (name: string) => void; busy?: boolean;
}) {
  const [path, setPath] = useState("");
  const [tissue, setTissue] = useState("melanoma");
  // "recently loaded" history (localStorage). Re-read whenever a section loads (sess changes), since
  // App.loadSection pushes the entry on success.
  const [recent, setRecent] = useState<RecentSection[]>(() => getRecent());
  useEffect(() => { setRecent(getRecent()); }, [sess?.session_id]);

  const metrics = sess
    ? [
        { k: "cells", v: sess.n_obs.toLocaleString() },
        { k: "genes", v: sess.n_vars.toLocaleString() },
        { k: "tissue", v: sess.tissue },
        { k: "device", v: sess.device || "-" },
      ]
    : [];

  // "Run full analysis": kick the whole spine as a background job, poll its status, and reload the
  // section when it finishes (so the app becomes a viewer of precomputed results). One writer at a
  // time - the backend 409s per-click steps while it runs, and the button disables itself.
  const [plRun, setPlRun] = useState(false);
  const [plStat, setPlStat] = useState<api.PipelineStatus | null>(null);
  const runFull = async () => {
    if (!sid || plRun) return;
    setPlRun(true); setPlStat(null);
    try {
      await api.runPipeline(sid, { tumour: sess?.is_tumour ?? null });
      for (;;) {
        await new Promise((r) => setTimeout(r, 1000));
        let st: api.PipelineStatus;
        try { st = await api.pipelineStatus(sid); } catch { continue; }
        setPlStat(st);
        if (st.done) {
          if (!st.error && onMutate) await onMutate();
          if (!st.error && sid) void prefetchReports(sid);   // warm Panel-check / QC so opening them is instant
          break;
        }
      }
    } catch (e: any) {
      setPlStat({ running: false, frac: 0, label: "", stage: null, done: true, error: String(e),
                  route: null, stages: [], summary: {}, state: null });
    }
    setPlRun(false);
  };

  return (
    <div className="panel">
      <h4>Your specimen</h4>
      {sess && (
        <div className="metricgrid">
          {metrics.map((m) => (
            <div className="metric" key={m.k}>
              <div className="k">{m.k}</div>
              <div className="v" style={{ fontSize: String(m.v).length > 8 ? "1.05rem" : "1.5rem" }}>{m.v}</div>
            </div>
          ))}
        </div>
      )}
      {sess?.panel_name && (
        <div className="callout" style={{ marginTop: 10 }}>Panel &middot; <b>{sess.panel_name}</b></div>
      )}

      {sid && (
        <div style={{ marginTop: 12 }}>
          <button className="pbtn pri" style={{ width: "100%" }} disabled={busy || plRun} onClick={runFull}
                  title="Run the whole pipeline (QC, panel check, cluster, annotate, niches, malignant calling, self-heal) in one background job, then load the results.">
            {plRun ? "Running full analysis …" : "Run full analysis"}
          </button>
          {plRun && (
            <div style={{ marginTop: 8 }}>
              <div className="meter"><div className="seg" style={{ width: `${Math.round((plStat?.frac || 0) * 100)}%`, background: "var(--accent, var(--pass))" }} /></div>
              <div className="pmuted" style={{ marginTop: 4, fontSize: ".8rem" }}>{STAGE_LABEL[plStat?.label || "starting"] || plStat?.label || "starting"} &hellip;</div>
            </div>
          )}
          {!plRun && plStat?.done && (
            plStat.error
              ? <div className="callout warn" style={{ marginTop: 8 }}>Full analysis failed: {plStat.error}</div>
              : <div className="callout" style={{ marginTop: 8 }}>
                  Full analysis done &middot; route <b>{plStat.route}</b> &middot; {(plStat.summary?.ok?.length ?? 0)} ran
                  {(plStat.summary?.skipped?.length ?? 0) > 0 ? `, ${plStat.summary!.skipped!.length} skipped` : ""}
                </div>
          )}
        </div>
      )}

      <div className="pmuted" style={{ marginTop: 10 }}>
        {sess?.demo?.source || "Demo - 10x Prime 5K breast Xenium (public, ~20k-cell subset)."}{" "}
        <Help text="Walk the steps on the left, or just ask the copilot below." />
      </div>
      {sess?.demo?.role && (
        <div className="callout" style={{ marginTop: 8 }}>{sess.demo.role}</div>
      )}

      {(onLoadDemo || onLoadSynthetic) && (
        <>
          <hr style={{ border: "none", borderTop: "1px solid var(--line)", margin: "14px 0 10px" }} />
          <h4 style={{ marginBottom: 6 }}>
            Bundled sections{" "}
            <Help text="breast = a public Xenium Prime 5K breast section at ~50 counts/cell: only 1/10 cell types are confidently typable, and that is a depth verdict, not a panel defect. synthetic = an instant simulated melanoma section (nothing to fetch), for a quick tour of the workflow." />
          </h4>
          {onLoadDemo && (
            <button className="pbtn" style={{ width: "100%" }} disabled={busy} onClick={() => onLoadDemo("breast")}>
              Xenium Prime 5K breast (public)
            </button>
          )}
          {onLoadSynthetic && (
            <button className="pbtn" style={{ width: "100%", marginTop: onLoadDemo ? 8 : 0 }} disabled={busy} onClick={onLoadSynthetic}>
              Load synthetic melanoma (instant demo)
            </button>
          )}
        </>
      )}

      {onLoad && (
        <>
          <hr style={{ border: "none", borderTop: "1px solid var(--line)", margin: "14px 0 10px" }} />
          <h4 style={{ marginBottom: 6 }}>Load your own section <Help text="Server-side path to a Xenium / CosMx / MERSCOPE output folder, or a .h5ad file." /></h4>
          <input
            className="pill" style={{ width: "100%", marginBottom: 8 }}
            value={path} onChange={(e) => setPath(e.target.value)}
            placeholder="/path/to/xenium_output"
            onKeyDown={(e) => { if (e.key === "Enter" && path.trim()) onLoad(path.trim(), tissue); }}
          />
          <div style={{ display: "flex", gap: 8 }}>
            <input className="pill" list="tissue-suggestions" value={tissue} placeholder="tissue (free text)"
                   onChange={(e) => setTissue(e.target.value)} style={{ flex: "none", width: 150 }}
                   title="Free-text tissue/tumour context (e.g. 'uveal melanoma', 'mouse brain'). Drives the auto-selected reference and marker set." />
            <datalist id="tissue-suggestions">{TISSUES.map((t) => <option key={t} value={t} />)}</datalist>
            <button className="btn primary" style={{ flex: 1, minHeight: "auto", padding: ".5rem" }}
                    disabled={busy || !path.trim()} onClick={() => onLoad(path.trim(), tissue.trim() || "generic")}>
              {busy ? "loading …" : "Load section"}
            </button>
          </div>

          {recent.length > 0 && (
            <div style={{ marginTop: 10 }}>
              <div className="pmuted" style={{ fontSize: ".72rem", marginBottom: 4, display: "flex", justifyContent: "space-between" }}>
                <span>recently loaded</span>
                <span style={{ cursor: "pointer" }} title="clear this browser's history" onClick={() => { clearRecent(); setRecent([]); }}>clear</span>
              </div>
              {recent.map((r) => (
                <button key={r.path} title={`${r.path}  (${r.tissue || "generic"})`} disabled={busy}
                        onClick={() => { setPath(r.path); setTissue(r.tissue || "generic"); onLoad(r.path, r.tissue || "generic"); }}
                        style={{ display: "flex", justifyContent: "space-between", gap: 8, width: "100%", textAlign: "left",
                                 background: "transparent", border: "1px solid var(--line)", borderRadius: 6,
                                 padding: "5px 8px", marginBottom: 4, cursor: busy ? "default" : "pointer",
                                 color: "inherit", fontSize: ".76rem" }}>
                  <span style={{ overflow: "hidden", textOverflow: "ellipsis", whiteSpace: "nowrap" }}>{shortPath(r.path)}</span>
                  <span className="pmuted" style={{ flex: "none" }}>{r.tissue || "generic"}</span>
                </button>
              ))}
            </div>
          )}
        </>
      )}

      {/* Custom reference: available on every step (the reference feeds panel-match and transfer)
          - but only once a section is loaded. */}
      {sess && <TumourFlag sid={sid} sess={sess} onMutate={onMutate} />}
    </div>
  );
}
