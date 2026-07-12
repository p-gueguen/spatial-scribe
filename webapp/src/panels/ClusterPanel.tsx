// Cluster: Leiden communities on the shared neighbour graph, with a live resolution slider
// (higher = more, finer clusters). Re-runs on the server and recolours the map. Below it, the top-N
// differentially-expressed genes of every cluster as a dot-plot (click it to open full-screen).
import React, { useEffect, useState } from "react";
import * as api from "../api";
import { capKey, getCached, setCached } from "../capCache";
import { Fig, Help } from "./viz";

export default function ClusterPanel({ summary, sid, onMutate, sess, onRan }:
  { summary?: any; sid?: string; onMutate?: () => Promise<void>; sess?: api.State; onRan?: () => void }) {
  const d = summary?.cluster;
  // Opening this tab when clustering is already done (via Run-full-analysis or a prior run) greens the
  // rail step - it reflects work done this session, and the panel does not itself re-run the step.
  useEffect(() => { if (sid && sess?.done?.cluster) onRan?.(); /* eslint-disable-next-line */ }, [sid, sess?.done?.cluster]);
  // Seed the slider from the resolution the CURRENT clusters were made at, not a hardcoded 0.5. App
  // remounts this panel on every section mutation (a re-cluster is one), so a plain useState(0.5) reset
  // the slider back to 0.5 right after a run at, say, 0.1 - the reported "snaps back to 0.5" bug.
  const [res, setRes] = useState<number>(d?.resolution ?? 0.5);
  // useState only reads its initial value once, and the panel remounts BEFORE the fresh summary
  // arrives (App bumps dataEpoch, then awaits fetchSummary), so seeding alone snapped back to the
  // stale 0.5. Re-sync to the server's authoritative resolution whenever IT changes (a completed
  // re-cluster / full run) - a plain slider drag does not change d.resolution, so it is never clobbered.
  useEffect(() => { if (d?.resolution != null) setRes(d.resolution); }, [d?.resolution]);
  const [busy, setBusy] = useState(false);

  // DEGs. Seeded from the per-session cache, so leaving the tab and coming back does not recompute;
  // any re-cluster goes through App.invalidate(), which drops the cache and remounts this panel.
  const DEG_KEY = capKey(sid, "cluster_markers", "2");
  const cached = getCached<{ figs: api.Artifact[]; value: any }>(DEG_KEY);
  const [degs, setDegs] = useState<any>(cached?.value ?? null);
  const [degFigs, setDegFigs] = useState<api.Artifact[]>(cached?.figs ?? []);
  const [degBusy, setDegBusy] = useState(false);
  const [degErr, setDegErr] = useState<string | null>(null);

  const [clErr, setClErr] = useState<string | null>(null);
  const recluster = async () => {
    if (!sid) return; setBusy(true); setClErr(null);
    // run_step returns {ok,error} in a 200 body - a failed cluster (e.g. an all-zero panel gene) is NOT
    // an HTTP error. Surface it: silently ignoring r.ok is what made "clustering does nothing / empty
    // leiden" invisible on a Xenium section with unexpressed panel genes.
    try {
      const r = await api.runStep(sid, "cluster", { resolution: res });
      if (r && r.ok === false) setClErr(String(r.error?.hint || r.error?.message || r.error));
      else if (onMutate) await onMutate();
    } catch (e: any) { setClErr(String(e)); }
    setBusy(false);
  };

  const runDegs = async () => {
    if (!sid) return;
    setDegBusy(true); setDegErr(null);
    try {
      const r = await api.runCap(sid, "cluster_markers", { n_genes: 2 });
      if (r.ok) {
        const figs = r.artifacts || [];
        setDegs(r.value); setDegFigs(figs);
        setCached(DEG_KEY, { figs, value: r.value });
      } else setDegErr(String(r.error?.hint || r.error?.message || r.error));
    } catch (e: any) { setDegErr(String(e)); }
    setDegBusy(false);
  };

  return (
    <div className="panel">
      <h4>Cells grouped by expression</h4>
      <div className="metricgrid">
        <div className="metric"><div className="k">leiden clusters</div><div className="v">{d ? d.n_clusters : "-"}</div></div>
        {d?.resolution != null && <div className="metric"><div className="k">resolution</div><div className="v">{d.resolution}</div></div>}
      </div>
      <div className="subhead">re-cluster {busy && <span className="running"><i />clustering</span>}</div>
      <div className="rngrow">
        <span className="lab">resolution
          <Help text="Higher resolution splits the section into more, finer Leiden communities." /></span>
        <input className="rng" type="range" min={0.1} max={2} step={0.1} value={res} onChange={(e) => setRes(parseFloat(e.target.value))} />
        <b>{res.toFixed(1)}</b>
      </div>
      <button className="pbtn pri" style={{ width: "100%", marginTop: 6 }} disabled={busy} onClick={recluster}>
        {busy ? "clustering …" : `Cluster at ${res.toFixed(1)}`}
      </button>
      {clErr && <div className="err">{clErr}</div>}

      <button className="pbtn" style={{ width: "100%", marginTop: 8 }} disabled={busy || degBusy || !d} onClick={runDegs}>
        {degBusy ? "ranking genes …" : "Top DEGs per cluster"}
      </button>
      {degErr && <div className="err">{degErr}</div>}
      {degs && !degErr && (
        <div className="pmuted" style={{ marginTop: 4 }}>
          {degs.genes?.length} genes across {degs.n_clusters} clusters &middot; {degs.source} ({degs.device}).
          Click the plot to enlarge.
        </div>
      )}
      {degFigs.map((f, i) => <Fig key={i} a={f} />)}
    </div>
  );
}
