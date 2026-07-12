// Spatial: TME niches, the neighborhood-enrichment heatmap, and de-novo gene programs (NMF).
import React, { useEffect, useRef, useState } from "react";
import * as api from "../api";
import { capKey, getCached, setCached } from "../capCache";
import { Heat, Help } from "./viz";

type Tab = "niches" | "neighbors" | "programs";

// See AnnotatePanel: App remounts this panel on any section mutation so stale local results die with
// it. The subtab is UI state, so it survives the remount.
let _lastTab: Tab = "niches";

export default function SpatialPanel({ summary, sid, sess, onColor, onPreviewColor, onEmphasize, onMutate }:
  { summary?: any; sid?: string; sess?: any; onColor?: (f: string) => void;
    onPreviewColor?: (f: string | null) => void;
    onEmphasize?: (labels: string[] | null) => void; onMutate?: () => Promise<void> | void }) {
  const d = summary?.spatial;
  const [tab, setTabState] = useState<Tab>(_lastTab);
  const setTab = (t: Tab) => { _lastTab = t; setTabState(t); };
  const typeNames: string[] = (summary?.annotate?.cell_types ?? []).map((t: any) => t.name);

  // neighbors. nhood is seeded from the per-session cache so returning to this tab does not recompute
  // (App clears the cache on any section-mutating op / new load - see capCache).
  const [nhood, setNhood] = useState<any>(() => getCached(capKey(sid, "neighborhood_enrichment")) ?? null);
  const [nbusy, setNbusy] = useState(false);
  const nhoodAuto = useRef<string | undefined>(undefined);   // session we've kicked the auto-enrichment for
  const emphColoured = useRef(false);   // have we recoloured the map to cell_type for the current heatmap
  // recolour=true on the manual recompute (colour the map by cell_type so hovering a heatmap pair can
  // highlight it); the automatic tab-open run passes false so it never fights the user's colouring.
  const runNhood = async (recolour: boolean) => {
    if (!sid) return; setNbusy(true); setNhood(null);
    if (recolour) { onColor?.("cell_type"); emphColoured.current = true; }
    try {
      const r = await api.runCap(sid, "neighborhood_enrichment");
      const v = r.ok ? r.value : { error: r.error }; setNhood(v);
      if (r.ok) setCached(capKey(sid, "neighborhood_enrichment"), v);
    } catch (e: any) { setNhood({ error: String(e) }); }
    setNbusy(false);
  };
  // Compute the enrichment automatically the first time the neighbors tab is opened this session
  // (once - the cache short-circuits a repeat). No recolour, so it does not disturb the current map.
  useEffect(() => {
    if (tab !== "neighbors" || !sid || !typeNames.length) return;
    if (getCached(capKey(sid, "neighborhood_enrichment")) || nhoodAuto.current === sid) return;
    nhoodAuto.current = sid; runNhood(false);
    /* eslint-disable-next-line */
  }, [tab, sid, typeNames.length]);
  // Hovering a heatmap pair: recolour to cell_type once (so the emphasis has cell_type colours to
  // match), then emphasise. After the first hover the map stays cell_type, so later hovers are exact.
  const onHeatHover = (r: string, c: string) => {
    if (!emphColoured.current) { onColor?.("cell_type"); emphColoured.current = true; }
    onEmphasize?.(r === c ? [r] : [r, c]);
  };

  // programs. Seeded from the per-session cache so leaving and returning keeps the NMF result.
  const [prog, setProg] = useState<any>(() => getCached(capKey(sid, "programs")) ?? null);
  const [pbusy, setPbusy] = useState(false);
  const runProg = async () => {
    if (!sid) return; setPbusy(true); setProg(null);
    try {
      const r = await api.runCap(sid, "discover_programs", { n_programs: 8 });
      const v = r.ok ? r.value : { error: r.error }; setProg(v);
      if (r.ok) setCached(capKey(sid, "programs"), v);
    } catch (e: any) { setProg({ error: String(e) }); }
    setPbusy(false);
  };
  // LLM-name each program from its top loading genes; relabels obs["program"] server-side + recolours.
  // Do NOT call onMutate here: it invalidates the caches and REMOUNTS this panel, and the remount
  // re-seeded `prog` from the just-cleared cache -> the freshly-named table vanished and the tab
  // snapped back to the "Discover programs" button. name_programs only renames labels (no cell change),
  // so a local setProg + re-cache + onColor("program") (refetches points with the new labels) suffices.
  const nameProg = async () => {
    if (!sid) return; setPbusy(true);
    try {
      const r = await api.runCap(sid, "name_programs");
      if (r.ok) { setProg(r.value); setCached(capKey(sid, "programs"), r.value); onColor?.("program"); }
      else setProg({ error: r.error });
    } catch (e: any) { setProg({ error: String(e) }); }
    setPbusy(false);
  };
  // niches. Auto-computed the first time the niches tab is opened (novae is heavy, so it runs ONCE and
  // caches; the pill shows progress). The result is held locally (like nhood/programs) so it survives a
  // remount without an onMutate that would loop the auto-run. Falls back to the summary's spatial-step
  // niches when those exist. Rows are normalised to {id, count}.
  const [nicheRes, setNicheRes] = useState<any>(() => getCached(capKey(sid, "niches")) ?? null);
  const [nichesBusy, setNichesBusy] = useState(false);
  const nicheAuto = useRef<string | undefined>(undefined);
  const runNiches = async (recolour: boolean) => {
    if (!sid) return; setNichesBusy(true); setNicheRes(null);
    try {
      const r = await api.runCap(sid, "niches");
      const v = r.ok ? (r.value || []) : { error: r.error }; setNicheRes(v);
      if (r.ok) { setCached(capKey(sid, "niches"), v); if (recolour) onColor?.("niche"); }
    } catch (e: any) { setNicheRes({ error: String(e) }); }
    setNichesBusy(false);
  };
  const nicheRows: any[] = Array.isArray(nicheRes) ? nicheRes.map((r: any) => ({ id: r.niche, count: r.n_cells }))
    : (d?.niches ?? []);
  useEffect(() => {
    if (tab !== "niches" || !sid || !typeNames.length) return;
    if (getCached(capKey(sid, "niches")) || nicheRows.length || nicheAuto.current === sid) return;
    nicheAuto.current = sid; runNiches(false);
    /* eslint-disable-next-line */
  }, [tab, sid, typeNames.length]);
  // The niches tab colours the MAIN canvas by niche - that is the tab's whole point. Recolour once
  // niches are available (freshly auto-computed OR cached/summary from a prior visit), once per tab
  // entry so a manual recolour by the user while on the tab is not clobbered. onColor("niche")
  // refetches points by obs['niche'] (present server-side once niches ran), matching the rail-run's
  // step.color behaviour. Reset on leaving so re-entering re-applies it.
  const nicheColoured = useRef(false);
  useEffect(() => { if (tab !== "niches") nicheColoured.current = false; }, [tab]);
  useEffect(() => {
    if (tab === "niches" && nicheRows.length && !nicheColoured.current) {
      nicheColoured.current = true; onColor?.("niche");
    }
    /* eslint-disable-next-line */
  }, [tab, nicheRows.length]);
  const maxN = nicheRows.reduce((m: number, n: any) => Math.max(m, n.count), 0) || 1;

  return (
    <div className="panel">
      <h4>The tissue in context</h4>
      <div className="subtabs">
        {(["niches", "neighbors", "programs"] as Tab[]).map((t) =>
          <span key={t} className={"subtab" + (tab === t ? " on" : "")} onClick={() => setTab(t)}>{t}</span>)}
      </div>

      {tab === "niches" && (
        <>
          <div className="subhead">TME niches {nichesBusy && <span className="running"><i />computing niches</span>}
            <Help text="Zero-shot spatial domains from the novae graph foundation model. Computed automatically when you open this tab (a heavier model - it can take a moment); cached after." /></div>
          {nicheRes?.error && <div className="err">{String(nicheRes.error?.hint || nicheRes.error?.message || nicheRes.error)}</div>}
          {nicheRows.length ? (
            <>
              <div className="metricgrid"><div className="metric"><div className="k">niches</div><div className="v">{nicheRows.length}</div></div></div>
              {nicheRows.map((n: any) => (
                <div className="barrow" key={n.id}>
                  <div className="lab" onClick={() => onColor && onColor("niche")} style={{ cursor: "pointer" }}>{String(n.id).startsWith("Niche") ? n.id : `Niche ${n.id}`}</div>
                  <div className="track"><div className="fill" style={{ width: (n.count / maxN) * 100 + "%" }} /></div>
                  <div className="n">{n.count.toLocaleString()}</div>
                </div>
              ))}
            </>
          ) : !nicheRes?.error && <div className="pmuted">{nichesBusy ? "computing spatial niches (novae) …" : "no niches yet"}</div>}
        </>
      )}

      {tab === "neighbors" && (
        <>
          <div className="subhead">neighborhood enrichment {nbusy && <span className="running"><i />computing</span>}<Help text="Which cell types sit next to each other more (red) or less (violet) than chance. Computed automatically when you open this tab; hover a pair to highlight both on the map." /></div>
          {!typeNames.length ? <div className="pmuted">Annotate first.</div> : (
            <>
              <button className="pbtn" onClick={() => runNhood(true)} disabled={nbusy}>{nbusy ? "computing …" : (nhood && !nhood.error ? "Recompute" : "Compute enrichment")}</button>
              {nhood?.error && <div className="err">{String(nhood.error?.hint || nhood.error?.message || nhood.error)}</div>}
              {nhood && !nhood.error && nhood.categories && (
                <>
                  <Heat rows={nhood.categories} cols={nhood.categories} matrix={nhood.zscore} diverging
                        onHover={onHeatHover}
                        onLeave={() => onEmphasize?.(null)} />
                  {/* Which map colour is which when a pair is hovered (App.onEmphasize: row->red, col->blue). */}
                  <div className="pmuted" style={{ marginTop: 4, fontSize: 11 }}>
                    hover a pair &rarr; <b style={{ color: "#EF4444" }}>row</b> cell red, <b style={{ color: "#3898FF" }}>column</b> cell blue on the map
                  </div>
                </>
              )}
            </>
          )}
        </>
      )}

      {tab === "programs" && (
        <>
          <div className="subhead">de-novo gene programs (NMF) {pbusy && <span className="running"><i />scoring signatures</span>}<Help text="Hover a program row to recolour the map by that program's per-cell signature score. The 'signature' column scores how well each program's own cells express its genes (AUROC)." /></div>
          <button className="pbtn pri" style={{ width: "100%" }} disabled={pbusy} onClick={runProg}>{pbusy ? "discovering programs + scoring signatures …" : "Discover programs"}</button>
          {prog?.error && <div className="err">{String(prog.error?.hint || prog.error?.message || prog.error)}</div>}
          {Array.isArray(prog) && (
            <>
              {/* Lean layout: cell count moved to the row tooltip (low value inline); top genes capped at
                  5 in a smaller mono size so each row is ~2 lines, not 4. */}
              <table className="dtable" style={{ marginTop: 6 }}>
                <thead><tr><th>program</th><th>signature</th><th>top genes</th></tr></thead>
                <tbody>{prog.map((p: any, i: number) => (
                  <tr key={i} title={`${Number(p.n_cells).toLocaleString()} cells${p.rationale ? " · " + p.rationale : (p.confidence ? " · " + p.confidence : "")}`}
                      onMouseEnter={() => p.score_field && onPreviewColor?.(p.score_field)}
                      onMouseLeave={() => onPreviewColor?.(null)} style={{ cursor: "pointer" }}>
                    <td style={{ verticalAlign: "middle", whiteSpace: "nowrap" }}>{p.program}</td>
                    <td style={{ verticalAlign: "middle" }}>{p.specificity != null ? <ProgSpec s={p.specificity} label /> : <span className="pmuted">-</span>}</td>
                    <td style={{ verticalAlign: "middle", whiteSpace: "normal", fontSize: 11, letterSpacing: "-0.01em", color: "var(--muted)" }}>{(p.top_genes || []).slice(0, 5).join(", ")}</td></tr>
                ))}</tbody>
              </table>
              {/* What the traffic-light next to each signature score means. */}
              <div className="pmuted" style={{ marginTop: 6, display: "flex", gap: 10, alignItems: "center", flexWrap: "wrap", fontSize: 11 }}>
                <span>signature specificity (AUROC):</span>
                <span><span className="tl green" /> ≥0.85 strong</span>
                <span><span className="tl amber" /> 0.70-0.85 moderate</span>
                <span><span className="tl red" /> &lt;0.70 weak</span>
              </div>
              {sess?.has_key && (
                <button className="pbtn" style={{ width: "100%", marginTop: 6 }} disabled={pbusy} onClick={nameProg}>{pbusy ? "labeling …" : "Label programs with AI"}</button>
              )}
            </>
          )}
        </>
      )}

    </div>
  );
}

// signature specificity = AUROC of a program's own signature score separating its cells from the rest.
const specTone = (s: number) => (s >= 0.85 ? "green" : s >= 0.7 ? "amber" : "red");

function ProgSpec({ s, label }: { s: number; label?: boolean }) {
  const pct = Math.max(0, Math.min(1, Number(s))) * 100;
  return (
    <span className="specwrap" title={`signature specificity (AUROC) ${Number(s).toFixed(2)}`}>
      <span className="aucbar" style={{ minWidth: label ? 64 : 46 }}><i style={{ width: pct + "%" }} /></span>
      <span className={"tl " + specTone(s)} />
      {label && <b className="specnum">{Number(s).toFixed(2)}</b>}
    </span>
  );
}
