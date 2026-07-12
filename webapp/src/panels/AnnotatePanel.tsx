// Annotate: per-cell confidence + abstention, self-verify (marker-agreement audit + remediation),
// the granular "why weren't cells typed" breakdown, curation (rename/merge + subcluster), cell
// states, and cross-checked malignant calling (marker + Cancer-Finder + InSituCNV).
import React, { useEffect, useState } from "react";
import * as api from "../api";
import { capKey, getCached, setCached } from "../capCache";
import { Heat, Fig, Help } from "./viz";

type Tab = "confidence" | "verify" | "curate" | "states" | "markers";
// Order tells a review story: how sure are the labels (confidence), are they right (verify ->
// markers evidence), what states, then curate (EDIT last, after you've reviewed).
const TABS: Tab[] = ["confidence", "verify", "markers", "states", "curate"];
const TAB_LABELS: Record<Tab, string> = {
  confidence: "confidence", verify: "verify",
  curate: "curate", states: "states", markers: "markers",
};

// App remounts this panel (through its React key) whenever the section mutates, so a stale result held
// in local state cannot outlive a rename / subcluster / re-run. The chosen subtab is UI state, not a
// result, so it is kept module-side and survives the remount - otherwise curating a type would bounce
// the user back to "confidence" every time.
let _lastTab: Tab = "confidence";

// The reference lives session-wide (set in the Panel-check step, its first consumer). Annotate shows
// only its STATUS + the supervised transfer action - not the uploader - so this tab stays a review
// surface. "change" jumps back to Panel check; with no reference the user is pointed there.
function ReferenceTransfer({ sid, sess, onMutate, onColor, onGoStep }:
  { sid?: string; sess?: api.State | null; onMutate?: () => Promise<void>; onColor?: (f: string) => void; onGoStep?: (k: string) => void }) {
  const [busy, setBusy] = useState(false);
  const [xfer, setXfer] = useState<any>(null);
  const ref = sess?.reference || null;
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
  return (
    <div style={{ margin: "2px 0 6px" }}>
      <div className="subhead">supervised: transfer from a reference
        <Help text="Transfer a single-cell reference's cell-type labels onto this section (CellTypist in-env; TACCO OT when installed). The reference is set in the Panel-check step - use 'change' to set or swap it." />
      </div>
      {ref ? (
        <>
          <div className="covrow">
            <span className="dot" style={{ background: "var(--pass)" }} />
            reference: <b>{ref.name || ref.label_key || "custom"}</b>
            {ref.n_types != null ? <> · {ref.n_types} types</> : null} · {ref.n_cells.toLocaleString()} cells
            <button className="pbtn sm" style={{ marginLeft: "auto" }} onClick={() => onGoStep?.("panel")}>change</button>
          </div>
          <div style={{ display: "flex", gap: 8, marginTop: 6 }}>
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
        </>
      ) : (
        <div className="pmuted">No reference set. <button className="pbtn sm" onClick={() => onGoStep?.("panel")}>Add one in Panel check</button> for supervised transfer - or annotate from the clusters below.</div>
      )}
    </div>
  );
}

export default function AnnotatePanel({ summary, sid, onMutate, onColor, onPreviewColor, onEmphasize, sess, onRan, onGoStep }:
  { summary?: any; sid?: string; onMutate?: () => Promise<void>; onColor?: (f: string) => void;
    onPreviewColor?: (f: string | null, onlyType?: string, posMin?: number) => void; onEmphasize?: (labels: string[] | null) => void;
    sess?: api.State; onRan?: () => void; onGoStep?: (k: string) => void }) {
  const d = summary?.annotate;
  // Opening this tab when annotation is already done (via Run-full-analysis or a prior run) greens the
  // rail step - it reflects work done this session, and the panel does not itself re-run the step.
  useEffect(() => { if (sid && sess?.done?.annotate) onRan?.(); /* eslint-disable-next-line */ }, [sid, sess?.done?.annotate]);
  const [tab, setTabState] = useState<Tab>(_lastTab);
  const setTab = (t: Tab) => { _lastTab = t; setTabState(t); };
  const types: any[] = d?.cell_types ?? [];
  const typeNames = types.map((t) => t.name);

  // confidence: rejection breakdown
  const [rej, setRej] = useState<any>(null);
  useEffect(() => {
    setRej(null);
    if (sid && d?.confident_pct != null && tab === "confidence") {
      api.runCap(sid, "rejection_reasons").then((r) => { if (r.ok) setRej(r.value); }).catch(() => {});
    } /* eslint-disable-next-line */
  }, [sid, tab, d?.abstained_pct]);

  // Annotation-quality battery (label-dependent: it scores the CURRENT labels, so it belongs here in
  // Annotate, not in the QC tab which runs on the raw section before annotation). Reuses the qc_funnel
  // capability - now that cell_type exists, the funnel computes its label-dependent layers - cached
  // per session so revisiting the tab does not recompute.
  const [aq, setAq] = useState<any>(null);
  useEffect(() => {
    setAq(null);
    if (sid && d?.confident_pct != null && tab === "confidence") {
      const key = capKey(sid, "qc_funnel");
      const cached = getCached<any>(key);
      if (cached) { setAq(cached); return; }
      api.runCap(sid, "qc_funnel").then((r) => { if (r.ok) { setAq(r.value); setCached(key, r.value); } }).catch(() => {});
    } /* eslint-disable-next-line */
  }, [sid, tab, d?.abstained_pct]);

  // verify (dedicated tab): marker-agreement audit + remediation
  const [ver, setVer] = useState<any>(null); const [vbusy, setVbusy] = useState(false);
  const [led, setLed] = useState<any>(null);   // trust ledger (resolvable x coherent x agreed per type)
  const [heal, setHeal] = useState<any>(null); const [hbusy, setHbusy] = useState(false);
  useEffect(() => {
    if (sid && tab === "verify" && !ver) {
      const key = capKey(sid, "verify_report");
      const cached = getCached(key);
      if (cached) { setVer(cached); return; }
      setVbusy(true);
      api.verifyReport(sid).then((r) => { setVer(r); if (r?.status === "ok") setCached(key, r); })
        .catch((e) => setVer({ status: "error", error: String(e) })).finally(() => setVbusy(false));
    }
    if (sid && tab === "verify" && !led) {                    // trust ledger runs alongside the audit
      const lk = capKey(sid, "trust_ledger");
      const lc = getCached<any>(lk);
      if (lc) setLed(lc);
      else api.runCap(sid, "trust_ledger").then((r) => { if (r.ok) { setLed(r.value); setCached(lk, r.value); } }).catch(() => {});
    } /* eslint-disable-next-line */
  }, [sid, tab]);
  // "Optimize annotation": run the safe self-heal pass (subcluster the flagged heterogeneous types,
  // then re-verify), then refresh the audit + recolour by the new subtypes. self_heal only ADDS a
  // `subtype` column (it never changes cell_type / renames / merges), so no full remount is needed.
  const runHeal = async () => {
    if (!sid) return; setHbusy(true); setHeal(null);
    try {
      const r = await api.runCap(sid, "self_heal");
      if (r.ok) {
        setHeal(r.value);
        const vr = await api.verifyReport(sid); setVer(vr);        // show the improved agreement
        if (vr?.status === "ok") setCached(capKey(sid, "verify_report"), vr);
        if (onColor) onColor("subtype");
      } else setHeal({ error: r.error });
    } catch (e: any) { setHeal({ error: String(e) }); }
    setHbusy(false);
  };

  // curate
  const [old, setOld] = useState(""); const [neu, setNeu] = useState(""); const [cbusy, setCbusy] = useState(false);
  const [subType, setSubType] = useState(""); const [sub, setSub] = useState<any>(null); const [sbusy, setSbusy] = useState(false);
  const [mergePrev, setMergePrev] = useState<any>(null); const [mgbusy, setMgbusy] = useState(false);
  // Reset a selection to a CURRENT type when the old one is gone (a rename/merge changes typeNames);
  // `s || typeNames[0]` kept the stale name, so after renaming Oligo->Coucou "Subcluster Oligo" hit 0
  // cells. Keep the pick if still valid, else fall back to the first type.
  useEffect(() => {
    if (typeNames.length) {
      setOld((o) => (typeNames.includes(o) ? o : typeNames[0]));
      setSubType((s) => (typeNames.includes(s) ? s : typeNames[0]));
    } /* eslint-disable-next-line */
  }, [typeNames.join(",")]);

  const rename = async () => {
    if (!sid || !old || !neu.trim()) return; setCbusy(true);
    try { await api.renameCelltype(sid, old, neu.trim()); setNeu(""); if (onMutate) await onMutate(); } catch { /* keep */ }
    setCbusy(false);
  };
  const doSubcluster = async () => {
    if (!sid || !subType) return; setSbusy(true); setSub(null);
    try { const r = await api.runCap(sid, "subcluster", { cell_type: subType }); if (r.ok) { setSub(r.value?.rows || []); if (onMutate) await onMutate(); if (onColor) onColor("subtype"); } else setSub({ error: r.error }); } catch (e: any) { setSub({ error: String(e) }); }
    setSbusy(false);
  };
  // Suggest merges (panel + ontology): PREVIEW the confusable types the panel cannot separate, grouped
  // only within a shared Cell-Ontology lineage (dry_run), then Apply collapses each group into one label.
  const suggestMerges = async () => {
    if (!sid) return; setMgbusy(true); setMergePrev(null);
    try { const r = await api.runCap(sid, "merge_types", { dry_run: true }); setMergePrev(r.ok ? r.value : { error: r.error }); }
    catch (e: any) { setMergePrev({ error: String(e) }); }
    setMgbusy(false);
  };
  const applyMerges = async () => {
    if (!sid) return; setMgbusy(true);
    try { const r = await api.runCap(sid, "merge_types", {}); if (r.ok) { setMergePrev(null); if (onColor) onColor("cell_type"); if (onMutate) await onMutate(); } else setMergePrev({ error: r.error }); }
    catch (e: any) { setMergePrev({ error: String(e) }); }
    setMgbusy(false);
  };

  // states + malignant concordance
  const [states, setStates] = useState<any>(null); const [mal, setMal] = useState<any>(null); const [stbusy, setStbusy] = useState(false); const [mbusy, setMbusy] = useState(false);
  const [activeState, setActiveState] = useState<string | null>(null);   // heatmap column currently painted on the map
  // On (re)entering the states tab, serve state_by_celltype from the session cache, and rehydrate the
  // button-driven results (assign_cell_states, malignant_concordance) so leaving and coming back does
  // not recompute. score_fields maps each program to its per-cell obs column for the hover recolour.
  const scoreField = (st: string) => states?.score_fields?.[st] ?? ("state_" + st.replace(/[/ ]/g, "_")); // fallback for older servers lacking score_fields
  useEffect(() => {
    if (sid && tab === "states") {
      if (!states) {
        const key = capKey(sid, "state_by_celltype");
        const cached = getCached(key);
        if (cached) setStates(cached);
        else { setStbusy(true); api.runCap(sid, "state_by_celltype").then((r) => { if (r.ok) { setStates(r.value); setCached(key, r.value); } }).catch(() => {}).finally(() => setStbusy(false)); }
      }
      if (!cstate) { const c = getCached(capKey(sid, "assign_cell_states")); if (c) setCstate(c); }
      if (!mal) { const c = getCached(capKey(sid, "malignant_concordance")); if (c) setMal(c); }
    } /* eslint-disable-next-line */
  }, [sid, tab]);
  const runMal = async () => {
    if (!sid) return; setMbusy(true);
    try { const r = await api.runCap(sid, "malignant_concordance"); if (r.ok) { setMal(r.value); setCached(capKey(sid, "malignant_concordance"), r.value); if (onColor) onColor("malignant_score"); } else setMal({ error: r.error }); }
    catch (e: any) { setMal({ error: String(e) }); }
    setMbusy(false);
  };
  // CyteType-style per-cell state typing: assign each cell a dominant program label (cell_state).
  const [cstate, setCstate] = useState<any>(null); const [csbusy, setCsbusy] = useState(false);
  const runAssignStates = async () => {
    if (!sid) return; setCsbusy(true);
    try { const r = await api.runCap(sid, "assign_cell_states"); if (r.ok) { setCstate(r.value); setCached(capKey(sid, "assign_cell_states"), r.value); if (onColor) onColor("cell_state"); if (onMutate) await onMutate(); } else setCstate({ error: r.error }); }
    catch (e: any) { setCstate({ error: String(e) }); }
    setCsbusy(false);
  };

  // markers: dot-plot (optionally SPLIT-purified before/after) + SPLIT spillover purification
  const [mkGenes, setMkGenes] = useState(""); const [mkType, setMkType] = useState("");
  const [figs, setFigs] = useState<api.Artifact[]>([]); const [dpbusy, setDpbusy] = useState(false); const [dpnote, setDpnote] = useState<string | null>(null);
  const [split, setSplit] = useState<any>(null); const [spbusy, setSpbusy] = useState(false); const [corrected, setCorrected] = useState(false);
  useEffect(() => { if (typeNames.length) setMkType((s) => (typeNames.includes(s) ? s : typeNames[0])); /* eslint-disable-line */ }, [typeNames.join(",")]);
  const drawDotplot = async () => {
    if (!sid) return; setDpbusy(true); setDpnote(null);
    const params: any = { corrected };
    const gl = mkGenes.split(/[,\s]+/).map((s) => s.trim()).filter(Boolean);
    if (gl.length) params.genes = gl; else if (mkType) params.cell_type = mkType;
    try { const r = await api.runCap(sid, "marker_dotplot", params); if (r.ok) { setFigs(r.artifacts || []); setDpnote(r.value?.note || null); } else setDpnote(String(r.error?.message || r.error)); }
    catch (e: any) { setDpnote(String(e)); }
    setDpbusy(false);
  };
  const runSplit = async () => {
    if (!sid) return; setSpbusy(true);
    try { const r = await api.runCap(sid, "split_purify"); const v = r.ok ? r.value : { status: String(r.error?.message || r.error) }; setSplit(v); if (r.ok && r.value?.status === "ok") { setCorrected(true); setCached(capKey(sid, "split_purify"), v); } }
    catch (e: any) { setSplit({ status: String(e) }); }
    setSpbusy(false);
  };
  // Rehydrate the button-driven SPLIT result on (re)entering the markers tab so it is not recomputed.
  useEffect(() => {
    if (sid && tab === "markers" && !split) { const c = getCached(capKey(sid, "split_purify")); if (c) { setSplit(c); if (c?.status === "ok") setCorrected(true); } }
    /* eslint-disable-next-line */
  }, [sid, tab]);
  // Auto-draw the category overview on tab open (rows = cell types ordered by category, cols = top
  // markers per category), so a figure is visible immediately - matching the verify/states
  // tabs. Empty {} params hit the backend's category-overview default; the result is memoized per
  // session so leaving and returning does not recompute. Manual controls below refine it.
  useEffect(() => {
    if (sid && tab === "markers" && !figs.length && !dpbusy) {
      const key = capKey(sid, "marker_dotplot", "overview");
      const cached = getCached<{ figs: api.Artifact[]; note: string | null }>(key);
      if (cached) { setFigs(cached.figs); setDpnote(cached.note); return; }
      setDpbusy(true);
      api.runCap(sid, "marker_dotplot", {})
        .then((r) => { if (r.ok) { const v = { figs: r.artifacts || [], note: r.value?.note || null }; setFigs(v.figs); setDpnote(v.note); setCached(key, v); } })
        .catch(() => {})
        .finally(() => setDpbusy(false));
    } /* eslint-disable-next-line */
  }, [sid, tab]);

  // Annotation is CREATED here (not just reviewed): the supervised route (ReferenceSection = transfer a
  // reference's labels) and the unsupervised route (annotate from cluster markers) both live in this tab,
  // so the panel renders the method choice when nothing is annotated yet, then the review subtabs once
  // cell_type exists. (`runStep('annotate')` returns {ok} in a 200 body - surface a failure, don't hide it.)
  const [annBusy, setAnnBusy] = useState(false);
  const [annErr, setAnnErr] = useState<string | null>(null);
  const runAnnotate = async () => {
    if (!sid) return; setAnnBusy(true); setAnnErr(null);
    try {
      const r = await api.runStep(sid, "annotate");
      if (r && r.ok === false) setAnnErr(String(r.error?.hint || r.error?.message || r.error));
      else if (onMutate) await onMutate();
    } catch (e: any) { setAnnErr(String(e)); }
    setAnnBusy(false);
  };

  const maxCount = types.reduce((m: number, t: any) => Math.max(m, t.count), 0) || 1;

  return (
    <div className="panel">
      <h4>{d ? "Named, with confidence" : "Name the cells"}</h4>
      {/* The reference is set in the Panel-check step (its first consumer - it sharpens resolvability
          there). Here we only show its status + the supervised label-transfer action, so Annotate is a
          REVIEW surface, not a second uploader (which cramped this tab). */}
      <ReferenceTransfer sid={sid} sess={sess} onMutate={onMutate} onColor={onColor} onGoStep={onGoStep} />
      {!d ? (
        <>
          <div className="subhead">or annotate from the clusters</div>
          {sess?.done?.cluster
            ? <button className="pbtn pri" style={{ width: "100%" }} disabled={annBusy} onClick={runAnnotate}>{annBusy ? "annotating …" : "Annotate from cluster markers"}</button>
            : <div className="pmuted">Cluster the section first (Cluster tab), then annotate from the cluster markers.</div>}
          {annErr && <div className="err">{annErr}</div>}
        </>
      ) : (
        <>
          <div className="subtabs">
            {TABS.map((t) => <span key={t} className={"subtab" + (tab === t ? " on" : "")} onClick={() => setTab(t)}>{TAB_LABELS[t]}</span>)}
          </div>

      {tab === "confidence" && (
        <>
          {d.confident_pct != null && (
            <>
              <div className="meter">
                <div className="seg" style={{ width: d.confident_pct + "%", background: "var(--pass)" }} />
                <div className="seg" style={{ width: d.tentative_pct + "%", background: "var(--warn)" }} />
                <div className="seg" style={{ width: d.abstained_pct + "%", background: "var(--fail)" }} />
              </div>
              {/* PASS / WARN / FAIL: the same three words obs['annotation_verdict'] paints on the map. */}
              <div className="meterlbl">pass <b>{d.confident_pct}%</b> · warn <b>{d.tentative_pct}%</b> · fail <b>{d.abstained_pct}%</b>
                <Help text="The per-cell annotation verdict. PASS = confidently typed, WARN = tentative, FAIL = abstained. Colour the map by annotation_verdict to see where each falls." />
              </div>
            </>
          )}
          {(() => {
            const q = aq?.annotation_quality;
            const aqi = q?.aqi;
            const coh = aq?.spatial_coherence?.mean_coherence;
            // Preferred headline: the AQI (the one robust, validated 0-1 index). Fall back to the raw
            // internal-validity signal only when the AQI could not be computed (guarded skip).
            if (aqi && typeof aqi.aqi === "number") {
              const score: number = aqi.aqi;
              const comp = aqi.components || {};
              const abst = aqi.abstention || {};
              const tone = score >= 0.6 ? "var(--pass)" : score >= 0.45 ? "var(--warn)" : "var(--fail)";
              const pct = (v: any) => (typeof v === "number" ? v.toFixed(2) : "-");
              return (
                <div className="metric" style={{ textAlign: "left", padding: "10px 12px", margin: "8px 0", borderLeft: `3px solid ${tone}` }}>
                  <div className="k">annotation quality index (AQI) <Help text="One robust 0-1 index of whether the cell-type labels hold up: the soft-min of marker purity (C) and marker-program fidelity (M), capped by the panel/depth resolvability ceiling (A) when a reference is loaded. Validated to order sections by true accuracy (Spearman 1.0 on 3 expert-ground-truth sections). It is an INDEX, not a correctness rate, and uses NO ground truth at runtime. Cross-method agreement is reported separately as the abstention signal, since its level does not transfer across sections." /></div>
                  <div className="v" style={{ fontSize: "1.55rem", color: tone }}>{score.toFixed(2)}<span className="pmuted" style={{ fontSize: ".6rem", fontWeight: 400 }}> / 1.00 index</span></div>
                  <div className="meter" style={{ marginTop: 4 }}><div className="seg" style={{ width: `${Math.round(score * 100)}%`, background: tone }} /></div>
                  {aqi.argmin && <div className="d" style={{ marginTop: 4 }}>limited by <b>{aqi.argmin}</b></div>}
                  {/* All sub-scores on ONE line; per-metric definitions + remedies live in the Help
                      tooltips, not the card body, so it stays lean (coherence is NOT an AQI term). */}
                  <div className="d" style={{ marginTop: 4 }}>
                    purity <b>{pct(comp.C)}</b> · fidelity <b>{pct(comp.M)}</b>
                    {typeof comp.A === "number" && <> · ceiling <b>{pct(comp.A)}</b></>}
                    {typeof coh === "number" && <> · coherence <b>{coh.toFixed(2)}</b></>}
                    {" "}<Help text={`purity (C) = marker contamination; fidelity (M) = marker-program AUC; ${typeof comp.A === "number" ? "ceiling (A) = panel/depth resolvability cap; " : ""}coherence = fraction of same-type neighbours (not an AQI term). AQI = min(A, soft-min(C, M)).`} />
                  </div>
                  <div className="d" style={{ marginTop: 4 }}>
                    per-cell signal <b>{abst.available ? `ensemble · ${abst.n_voters} methods` : "advisory"}</b>
                    {" "}<Help text={abst.available
                      ? "Cross-method agreement backs the per-cell abstention (trustworthy)."
                      : "Heuristic only. Run ≥3 reference methods in Panel-check for a trustworthy per-cell abstention signal."} />
                  </div>
                  {aqi.flags?.coverage_limited && (
                    <div className="callout warn" style={{ marginTop: 6 }}>
                      Purity coverage-limited: <b>{Math.round((aqi.no_dict_frac || 0) * 100)}%</b> of cells have a label with no on-panel markers, so the AQI can read low.
                      {" "}<Help text="C is scored over the covered types only. Load a reference in Panel-check to type the rest (reference-derived markers)." />
                    </div>
                  )}
                  {Array.isArray(aqi.coverage?.missing) && aqi.coverage.missing.length > 0 && (
                    <div className="callout warn" style={{ marginTop: 6 }}>
                      Markers suggest <b>{aqi.coverage.missing.join(", ")}</b>{typeof aqi.coverage.missing_frac === "number" && <> (~{Math.round(aqi.coverage.missing_frac * 100)}% of cells)</>} but no cells carry the label.
                      {" "}<Help text="The AQI scores only the labels present, so it can't penalize a dropped lineage - advisory. Recover by re-clustering at a higher resolution or running reference methods." />
                    </div>
                  )}
                  {Array.isArray(aqi.coverage?.missing_vs_reference) && aqi.coverage.missing_vs_reference.length > 0 && (
                    <div className="callout warn" style={{ marginTop: 6 }}>
                      Reference lineages with no matching label: <b>{aqi.coverage.missing_vs_reference.slice(0, 8).join(", ")}{aqi.coverage.missing_vs_reference.length > 8 ? ` +${aqi.coverage.missing_vs_reference.length - 8} more` : ""}</b>.
                      {" "}<Help text="A dropped compartment or a naming difference - catches compartments even when their markers are off-panel. If real, use supervised reference transfer." />
                    </div>
                  )}
                  {score >= 0.6 && !abst.available && !(aqi.coverage?.missing?.length > 0) && (
                    <div className="callout warn" style={{ marginTop: 6 }}>
                      High coherence, but no reference ensemble - a whole-cluster mislabel would go unseen.
                      {" "}<Help text="An internally coherent, marker-consistent cluster can still be the wrong label. Add a reference in Panel-check, run the methods, then check the verify tab." />
                    </div>
                  )}
                </div>
              );
            }
            const score: number | null = typeof q?.internal_validity?.integrated === "number" ? q.internal_validity.integrated : null;
            if (score == null) return null;
            const auc: number | null = typeof q?.marker_fidelity?.mean_auc === "number" ? q.marker_fidelity.mean_auc : null;
            const tone = score >= 0.65 ? "var(--pass)" : score >= 0.5 ? "var(--warn)" : "var(--fail)";
            return (
              <div className="metric" style={{ textAlign: "left", padding: "10px 12px", margin: "8px 0", borderLeft: `3px solid ${tone}` }}>
                <div className="k">cell-type annotation quality (internal validity) <Help text="Reference-free internal-validity signal (spatial-anno-metrics: silhouette, neighborhood purity, Ward proportion-match). Shown when the full AQI could not be computed. 0-1, higher is better - NOT a correctness rate, and no ground truth is used." /></div>
                <div className="v" style={{ fontSize: "1.55rem", color: tone }}>{score.toFixed(2)}<span className="pmuted" style={{ fontSize: ".6rem", fontWeight: 400 }}> / 1.00</span></div>
                <div className="meter" style={{ marginTop: 4 }}><div className="seg" style={{ width: `${Math.round(score * 100)}%`, background: tone }} /></div>
                {auc != null && <div className="d" style={{ marginTop: 4 }}>marker-program fidelity <b>{auc.toFixed(2)}</b> AUC · internal validity from silhouette, purity + Ward match</div>}
                {typeof coh === "number" && <div className="d" style={{ marginTop: 4 }}>spatial coherence <b>{coh.toFixed(2)}</b> · same-type neighbors</div>}
              </div>
            );
          })()}
          {types.slice(0, 35).map((t: any) => (
            <div className="barrow" key={t.name}>
              <div className="lab">{t.name}</div>
              <div className="track"><div className="fill" style={{ width: (t.count / maxCount) * 100 + "%" }} /></div>
              <div className="n">{t.count.toLocaleString()}</div>
            </div>
          ))}
          <div className="subhead">why weren't some cells typed?</div>
          {!rej ? <div className="pmuted">checking …</div> : (rej.breakdown?.length ? (
            <table className="dtable">
              <thead><tr><th>reason</th><th>cells</th><th title="Each reason's share of the UNTYPED cells (these sum to ~100%) - NOT a fraction of all cells in the section.">% of untyped</th></tr></thead>
              <tbody>{(() => {
                // Discreet in-cell data bar: width scaled to the LARGEST reason so the descending trend
                // reads at a glance; the printed number is still the true % of untyped.
                const maxPct = Math.max(...rej.breakdown.map((r: any) => r.pct_of_untyped || 0), 1);
                return rej.breakdown.map((r: any) => {
                  const p = Math.round(r.pct_of_untyped || 0);
                  const w = Math.round(100 * (r.pct_of_untyped || 0) / maxPct);
                  return (
                    <tr key={r.reason} title={r.description}>
                      <td>{r.label}</td>
                      <td>{Number(r.n_cells).toLocaleString()}</td>
                      <td style={{ background: `linear-gradient(to right, rgba(124,108,196,0.25) ${w}%, transparent ${w}%)` }}>{p}%</td>
                    </tr>
                  );
                });
              })()}</tbody>
            </table>
          ) : <div className="pmuted">Every cell was confidently typed.</div>)}
          {rej?.panel_warnings?.map((w: any, i: number) => (<div className="callout warn" key={i}>{w.message}</div>))}
        </>
      )}


      {tab === "verify" && (
        <>
          {led && led.per_type && (
            <>
              <div className="subhead">trust ledger <Help text="Three INDEPENDENT per-type verdicts and their disagreements: resolvable (can the panel separate this type?) x coherent (do its cells express its markers?) x agreed (do the reference methods back the label, when >=3 voted). The rows that matter are the CONTRADICTIONS - especially 'coherent but disputed', a whole-cluster mislabel the AQI index alone cannot see." /></div>
              {led.n_flagged > 0
                ? led.per_type.filter((r: any) => r.flags.length).slice(0, 6).map((r: any) => (
                    <div className="callout warn" key={r.cell_type}><b>{r.cell_type}</b> ({r.n_cells.toLocaleString()}) — {r.flags.join("; ")}</div>))
                : <div className="pmuted">No contradictions: resolvable, coherent{led.has_ensemble ? " and reference-agreed" : ""} line up for every type{led.has_ensemble ? "" : " (no ensemble - agreement not checked)"}.</div>}
              <table className="dtable" style={{ marginTop: 6 }}>
                <thead><tr><th>cell type</th><th>resolvable</th><th>coherent</th><th>agreed</th></tr></thead>
                <tbody>{led.per_type.slice(0, 16).map((r: any) => (
                  <tr key={r.cell_type} title={r.flags.length ? r.flags.join("; ") : ""} style={r.flags.length ? { color: "var(--warn)" } : undefined}>
                    <td>{r.cell_type}</td>
                    <td>{r.resolvable == null ? "-" : r.resolvable ? "yes" : "no"}</td>
                    <td>{r.coherent == null ? "-" : r.coherent.toFixed(2)}</td>
                    <td>{r.agreed == null ? "-" : r.agreed.toFixed(2)}</td>
                  </tr>))}</tbody>
              </table>
              {!led.has_ensemble && <div className="pmuted" style={{ marginTop: 4 }}>{led.note}</div>}
            </>
          )}
          <div className="subhead">label audit {vbusy && <span className="running"><i />auditing</span>}<Help text="Do the labels agree with their canonical markers? Marker-argmax agreement + one-vs-rest AUC per type, with grounded fixes for the ones that fail. Advisory - it never changes labels." /></div>
          {vbusy && !ver ? <div className="pmuted">auditing labels …</div>
            : !ver || ver.status !== "ok" ? <div className="pmuted">{ver?.note || ver?.error || "not available"}</div> : (
            <>
              <div className="cmp">
                <div className="cmpcell"><div className="k">section agreement</div><div className="v">{ver.section_agreement != null ? Math.round(ver.section_agreement * 100) + "%" : "-"}</div><div className="d">cells topping their own markers</div></div>
                <div className="cmpcell"><div className="k">mean AUC</div><div className="v">{ver.mean_auc != null ? ver.mean_auc.toFixed(2) : "-"}</div><div className="d">{ver.failed?.length || 0} type(s) flagged</div></div>
              </div>
              <table className="dtable" style={{ marginTop: 6 }}>
                <thead><tr><th>cell type</th><th>agree</th><th>AUC</th><th></th></tr></thead>
                <tbody>{(ver.per_type || []).slice(0, 16).map((r: any) => {
                  // Why a row is NOT green - so an amber dot with blank agree/AUC (a rare or
                  // marker-less type) explains itself instead of reading as an unexplained flag,
                  // exactly like the panel-check typability reason.
                  const reason = r.status === "skipped_small" ? `too few cells (n=${r.n_cells})`
                    : r.status === "unscoreable" ? "no on-panel markers"
                    : r.status === "fail" ? (r.confused_with ? `confused with ${r.confused_with}` : (r.cause || "low marker agreement"))
                    : null;
                  return (
                  <tr key={r.cell_type} title={reason || ""}>
                    <td>{r.cell_type}</td>
                    <td>{r.argmax_agreement != null ? Math.round(r.argmax_agreement * 100) + "%" : "-"}</td>
                    <td>{r.auc != null ? r.auc.toFixed(2) : "-"}</td>
                    <td>
                      <span className={"tl " + (r.status === "fail" ? "red" : r.status === "pass" ? "green" : "amber")} />
                      {reason && <span className="pmuted" style={{ fontSize: 10, marginLeft: 4 }}>{reason}</span>}
                    </td>
                  </tr>);
                })}</tbody>
              </table>
              {(ver.suggestions || []).length > 0 && (
                <>
                  <div className="subhead">suggested fixes</div>
                  {ver.suggestions.slice(0, 6).map((s: any, i: number) => (
                    <div className={"callout" + (s.advisory ? "" : " warn")} key={i}><b>{s.cell_type}</b> · {s.action}: {s.reason}</div>
                  ))}
                </>
              )}
              <div className="subhead" style={{ marginTop: 12 }}>optimize {hbusy && <span className="running"><i />re-running</span>}<Help text="Runs the safe self-heal pass: subclusters the flagged heterogeneous types to resolve their hidden substructure, then re-verifies. It only adds a subtype column - it never renames or merges (those stay in Curate)." /></div>
              <button className="pbtn pri" disabled={hbusy || !sid} onClick={runHeal}>{hbusy ? "re-running …" : "Optimize annotation (self-heal)"}</button>
              {heal && (heal.error
                ? <div className="err">{String(heal.error?.message || heal.error)}</div>
                : <div className="callout" style={{ marginTop: 8 }}>Fixed <b>{heal.n_fixed ?? 0}</b> of {heal.initial_failed?.length ?? 0} flagged type(s){heal.n_auto_actions ? ` in ${heal.n_auto_actions} action(s)` : ""}. {heal.final_failed?.length ? <>Still flagged: {heal.final_failed.join(", ")}.</> : <>None remain.</>}</div>)}
            </>
          )}
        </>
      )}

      {tab === "curate" && (
        <>
          <div className="subhead">rename / merge <Help text="Reuse an existing name as the target to merge two types." /></div>
          <div className="field"><label>rename</label>
            <select value={old} onChange={(e) => setOld(e.target.value)}>{typeNames.map((n) => <option key={n} value={n}>{n}</option>)}</select></div>
          <div className="field"><label>to (reuse a name to merge)</label>
            <input value={neu} placeholder="new label" onChange={(e) => setNeu(e.target.value)} onKeyDown={(e) => { if (e.key === "Enter") rename(); }} /></div>
          <button className="pbtn pri" style={{ width: "100%" }} disabled={cbusy || !neu.trim()} onClick={rename}>{cbusy ? "applying …" : "Apply rename"}</button>

          <div className="subhead" style={{ marginTop: 12 }}>suggest merges (panel + ontology) {mgbusy && <span className="running"><i />working</span>}<Help text="Preview the cell types the PANEL cannot separate (no private on-panel marker), grouped ONLY within a shared Cell-Ontology lineage (never across compartments - immune stays immune). Review the groups, then Apply collapses each into one coarser label." /></div>
          <button className="pbtn" style={{ width: "100%" }} disabled={mgbusy || !sid} onClick={suggestMerges}>Suggest merges</button>
          {mergePrev?.error && <div className="err">{String(mergePrev.error?.hint || mergePrev.error?.message || mergePrev.error)}</div>}
          {mergePrev && !mergePrev.error && (mergePrev.would_merge?.length ? (
            <>
              <table className="dtable" style={{ marginTop: 6 }}>
                <thead><tr><th>merge into</th><th>cells</th></tr></thead>
                <tbody>{mergePrev.would_merge.map((m: any, i: number) => (
                  <tr key={i} title={"members: " + (m.group || []).join(", ")}>
                    <td>{m.label}</td><td>{Number(m.n_cells).toLocaleString()}</td>
                  </tr>
                ))}</tbody>
              </table>
              <button className="pbtn pri" style={{ width: "100%", marginTop: 6 }} disabled={mgbusy} onClick={applyMerges}>{`Apply ${mergePrev.would_merge.length} merge${mergePrev.would_merge.length > 1 ? "s" : ""}`}</button>
            </>
          ) : <div className="pmuted" style={{ marginTop: 6 }}>{mergePrev.note || "The panel resolves every type - no confusable groups to merge."}</div>)}

          <div className="subhead">subcluster a type</div>
          <div className="field"><label>cell type</label>
            <select value={subType} onChange={(e) => setSubType(e.target.value)}>{typeNames.map((n) => <option key={n} value={n}>{n}</option>)}</select></div>
          <button className="pbtn pri" style={{ width: "100%" }} disabled={sbusy} onClick={doSubcluster}>{sbusy ? "subclustering (markers + LLM) …" : `Subcluster ${subType}`}</button>
          {sub?.error && <div className="err">{String(sub.error?.hint || sub.error?.message || sub.error)}</div>}
          {Array.isArray(sub) && (
            <table className="dtable" style={{ marginTop: 6 }}>
              <thead><tr><th>subtype</th><th>cells</th><th>conf</th></tr></thead>
              <tbody>{sub.map((r: any, i: number) => (
                <tr key={i} title={(r.rationale || "") + (r.top_markers ? "\nmarkers: " + (r.top_markers || []).join(", ") : "")}>
                  <td>{r.label || r.subcluster}</td><td>{Number(r.n_cells).toLocaleString()}</td><td>{r.confidence || "-"}</td>
                </tr>
              ))}</tbody>
            </table>
          )}
        </>
      )}

      {tab === "states" && (
        <>
          <div className="subhead">cell-type × state signatures {stbusy && <span className="running"><i />scoring</span>}<Help text="KNOWN cell-state signatures (e.g. exhaustion, proliferation) scored per cell type - distinct from the spatial tab's de-novo NMF gene programs, which are unsupervised. Mean signature score per cell type; hover a row to colour the map by its per-cell score." /></div>
          {stbusy && !states ? <div className="pmuted">scoring cell states …</div>
            : states?.cell_types?.length ? <Heat rows={states.cell_types} cols={states.states} matrix={states.matrix}
                activeCol={activeState}
                onHover={(row, col) => { setActiveState(col); onPreviewColor?.(scoreField(col), row); }}
                onHoverCol={(c) => { setActiveState(c); onPreviewColor?.(scoreField(c)); }}
                onLeave={() => { setActiveState(null); onPreviewColor?.(null); }} />
              : <div className="pmuted">No state signatures had enough on-panel genes.</div>}

          <div className="subhead">type cells by state <Help text="Label each cell with its dominant program (obs['cell_state']) and colour the map by it. State is what a cell is DOING, orthogonal to what it IS." /></div>
          {!cstate ? <button className="pbtn pri" style={{ width: "100%" }} disabled={csbusy} onClick={runAssignStates}>{csbusy ? "typing cell states …" : "Type cell states"}</button>
            : cstate.error ? <div className="err">{String(cstate.error?.hint || cstate.error?.message || cstate.error)}</div>
            : (
              <>
                {(Object.entries(cstate.distribution || {}) as [string, any][])
                  .sort((a, b) => Number(b[1]) - Number(a[1])).map(([st, nCells]) => {
                  const tot = (Object.values(cstate.distribution || {}) as any[]).reduce((s: number, v) => s + Number(v), 0) || 1;
                  const pct = (Number(nCells) / tot) * 100;
                  return (
                    <div className="barrow" key={st}>
                      <div className="lab">{st}</div>
                      <div className="track"><div className="fill" style={{ width: pct + "%" }} /></div>
                      <div className="n">{Number(nCells).toLocaleString()}</div>
                    </div>
                  );
                })}
                {cstate.per_celltype && Object.keys(cstate.per_celltype).length > 0 && (
                  <table className="dtable" style={{ marginTop: 6 }}>
                    <thead><tr><th>cell type</th><th>dominant state</th><th>% stateful</th></tr></thead>
                    <tbody>{Object.entries(cstate.per_celltype).slice(0, 14).map(([ct, d]: any) => (
                      <tr key={ct} title={cstate.llm_labels?.[ct]?.rationale || ""}>
                        <td>{ct}</td>
                        <td>{cstate.llm_labels?.[ct]?.state ? cstate.llm_labels[ct].state : d.top_state}</td>
                        <td>{Math.round((d.pct_stateful || 0) * 100)}%</td>
                      </tr>
                    ))}</tbody>
                  </table>
                )}
                <button className="pbtn sm" style={{ marginTop: 6 }} onClick={() => onColor && onColor("cell_state")}>colour map by state</button>
              </>
            )}

          <div className="subhead">malignant / tumor calling (cross-checked) {mbusy && <span className="running"><i />calling</span>}<Help text="Marker score always; Cancer-Finder runs in an isolated env when configured. Reports each caller's %-malignant. (InSituCNV is not run interactively - its infercnvpy pass is minutes-long.)" /></div>
          {!mal ? <button className="pbtn" disabled={mbusy} onClick={runMal}>{mbusy ? "calling malignant cells …" : "Call malignant cells"}</button>
            : mal.error ? <div className="err">{String(mal.error?.hint || mal.error?.message || mal.error)}</div>
            : mal.status && mal.status !== "ok" ? <div className="pmuted">{mal.status}</div> : (
              <>
                <div className="pmuted" style={{ marginBottom: 4 }}>Hover a caller to light up its malignant cells on the map.</div>
                <table className="dtable">
                  <thead><tr><th>caller</th><th>% malignant</th><th>status</th></tr></thead>
                  <tbody>{(mal.callers || []).map((c: any, i: number) => {
                    // Each caller writes a per-cell obs column (magma featureplot) - hover to recolour the
                    // map by it so the malignant cells glow. Only the callers that actually ran have one.
                    const n = String(c.name || "").toLowerCase();
                    const field = c.status !== "ok" ? null
                      : n.includes("cancer") ? "cancerfinder_prob"
                      : (n.includes("cnv") || n.includes("insitu")) ? "cnv_score"
                      : n.includes("marker") ? "malignant_score" : null;
                    return (
                      <tr key={i} style={field ? { cursor: "pointer" } : undefined}
                          onMouseEnter={() => field && onPreviewColor?.(field, undefined, typeof c.threshold === "number" ? c.threshold : 0.5)}
                          onMouseLeave={() => field && onPreviewColor?.(null)}>
                        <td>{c.name}</td>
                        <td>{c.status === "ok" ? Math.round((c.pct_malignant || 0) * 100) + "%" : "-"}</td>
                        <td>{c.status === "ok" ? "ok" : "skipped"}</td>
                      </tr>
                    );
                  })}</tbody>
                </table>
                {mal.concordance != null && (
                  <div className="covrow" style={{ marginTop: 6 }}><span className="dot" style={{ background: "var(--pass)" }} />Cancer-Finder vs InSituCNV agree<b>{Math.round(mal.concordance * 100)}%</b></div>
                )}
                {(mal.notes || []).map((n: string, i: number) => <div className="pmuted" key={i} style={{ marginTop: 4 }}>{n}</div>)}
              </>
            )}
        </>
      )}

      {tab === "markers" && (
        <>
          {/* The dot-plot is the point of this tab: top-2 markers per current cell type, auto-drawn on
              open (see the tab useEffect). Keep it at the TOP; the refine controls sit below it. */}
          <div className="subhead">top markers per cell type {dpbusy && <span className="running"><i />drawing</span>}
            <Help text="Top-2 markers per cell type - dot colour = mean expression, size = % of cells expressing. Drawn automatically on open; refine to specific genes or one type below." /></div>
          {figs.length ? figs.map((f, i) => <Fig key={i} a={f} />)
            : <div className="pmuted">{dpbusy ? "drawing marker dot-plot …" : "no on-panel markers to plot"}</div>}
          {dpnote && <div className="pmuted" style={{ marginTop: 4 }}>{dpnote}</div>}

          {/* Refine (secondary): swap in specific genes, one type's canonical markers, or the purified layer. */}
          <div className="subhead">refine</div>
          <div className="field">
            <input value={mkGenes} placeholder="genes e.g. CD3D, EPCAM (optional)" onChange={(e) => setMkGenes(e.target.value)} onKeyDown={(e) => { if (e.key === "Enter") drawDotplot(); }} /></div>
          <div className="field">
            <select value={mkType} onChange={(e) => setMkType(e.target.value)}>{typeNames.map((n) => <option key={n} value={n}>{n}</option>)}</select></div>
          <label className="pmuted" style={{ display: "flex", alignItems: "center", gap: 8, margin: "2px 0 8px" }}>
            <input type="checkbox" checked={corrected} onChange={(e) => setCorrected(e.target.checked)} /> show SPLIT-purified
          </label>
          <button className="pbtn" style={{ width: "100%" }} disabled={dpbusy} onClick={drawDotplot}>{dpbusy ? "drawing …" : "Redraw"}</button>

          <div className="subhead">SPLIT spillover purification</div>
          <div className="pmuted" style={{ marginBottom: 6 }}>Decontaminate transcript spillover between neighbouring cells, then tick &ldquo;show SPLIT-purified&rdquo; above to compare raw vs decontaminated.</div>
          <button className="pbtn" disabled={spbusy} onClick={runSplit}>{spbusy ? "purifying …" : "Run SPLIT purification"}</button>
          {split && (split.status === "ok" ? (
            <>
              {split.method && <div className="pmuted" style={{ marginTop: 4 }}>method · {split.method === "marker_neighbour" ? "in-app neighbour + marker" : "RCTD -> SPLIT"}</div>}
              <div className="cmp" style={{ marginTop: 6 }}>
                <div className="cmpcell"><div className="k">spillover removed</div><div className="v">{Math.round((split.pct_removed || 0) * 100)}%</div><div className="d">median library {Math.round(split.median_lib_before)} -&gt; {Math.round(split.median_lib_after)}</div></div>
                <div className="cmpcell"><div className="k">cells purified</div><div className="v">{Number(split.n_purified).toLocaleString()}</div><div className="d">{Math.round((split.coverage || 0) * 100)}% of section</div></div>
                {split.method === "marker_neighbour" && split.coexpr_pair && (
                  <div className="cmpcell"><div className="k">co-expression {split.coexpr_pair}</div><div className="v">{Math.round((split.coexpr_before || 0) * 100)}% -&gt; {Math.round((split.coexpr_after || 0) * 100)}%</div><div className="d">cross-lineage spillover</div></div>
                )}
              </div>
              {split.note && <div className="pmuted" style={{ marginTop: 4 }}>{split.note}</div>}
            </>
          ) : <div className="pmuted" style={{ marginTop: 4 }}>{split.status}</div>)}
        </>
      )}
        </>
      )}
    </div>
  );
}
