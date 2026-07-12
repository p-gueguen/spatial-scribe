import React, { useEffect, useMemo, useRef, useState } from "react";
import SpecimenCanvas from "./SpecimenCanvas";
import * as api from "./api";
import { clearSession, clearAll } from "./capCache";
import { pushRecent } from "./recentSections";
import QCPanel from "./panels/QCPanel";
import PanelCheckPanel from "./panels/PanelCheckPanel";
import ClusterPanel from "./panels/ClusterPanel";
import AnnotatePanel from "./panels/AnnotatePanel";
import SpatialPanel from "./panels/SpatialPanel";
import ReportPanel from "./panels/ReportPanel";
import LoadPanel from "./panels/LoadPanel";

const PANEL_FOR: Record<string, any> = { panel: PanelCheckPanel, qc: QCPanel, cluster: ClusterPanel, annotate: AnnotatePanel, spatial: SpatialPanel, report: ReportPanel };

// Fun, honest status lines while a step runs (annotate is the slow one - it calls Claude).
const ANNOTATE_QUIPS = [
  "Reading each cluster's top marker genes …",
  "Asking Claude to name the clusters …",
  "Weighing CD3D, MS4A1, PECAM1, EPCAM …",
  "Cross-checking markers against the panel …",
  "Scoring an honest per-cell confidence …",
  "Flagging ambiguous and low-signal cells …",
  "Deciding which cells to leave unlabeled …",
];
const GENERIC_QUIPS = ["Crunching on the GPU …", "Working through your section …", "Almost there …"];

// Determinate when the engine reports progress (frac>0 via /api/{sid}/progress), else the rotating
// quips carry the wait so nothing regresses when a step emits no ticks.
function RunProgress({ label, quips, frac, serverLabel }:
  { label: string; quips: string[]; frac: number; serverLabel: string }) {
  const [i, setI] = useState(0);
  useEffect(() => {
    const t = setInterval(() => setI((x) => (x + 1) % quips.length), 1600);
    return () => clearInterval(t);
  }, [quips.length]);
  const determinate = frac > 0;
  return (
    <div className="runoverlay">
      <div className="spinner" />
      <div className="runlabel">{label}</div>
      <div className="runquip">{determinate && serverLabel ? serverLabel : quips[i]}</div>
      {determinate && <div className="runbar"><div className="fill" style={{ width: `${Math.round(frac * 100)}%` }} /></div>}
    </div>
  );
}

// A small, quiet colour scale for a continuous featureplot (gene / QC / program score) - the
// specimen legend is categorical, so a feature colouring otherwise has no key for which end is
// high/low. The gradient renders the EXACT stops the backend painted the map with (`pts.ramp`:
// magma for a feature score, brand violet for QC/depth), so the key can never disagree with the map.
function FeatureScale({ field, ramp }: { field: string; ramp?: string[] | null }) {
  const stops = ramp && ramp.length >= 2 ? ramp : ["#23203a", "#6c5ce0", "#a896f2", "#ebe6ff"];
  return (
    <div className="fscale" title={`colour scale for ${field}`}>
      <span className="fslo">low</span>
      <span className="fsbar" style={{ background: `linear-gradient(90deg, ${stops.join(", ")})` }} />
      <span className="fshi">high</span>
    </div>
  );
}

const svg = (children: React.ReactNode) => (
  <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.8} strokeLinecap="round" strokeLinejoin="round">{children}</svg>
);
const ICONS: Record<string, React.ReactNode> = {
  load: svg(<><path d="M12 3v10" /><path d="M8 9l4 4 4-4" /><path d="M4 16v3a1 1 0 0 0 1 1h14a1 1 0 0 0 1-1v-3" /></>),
  panel: svg(<><rect x="3.5" y="3.5" width="7" height="7" rx="1.6" /><rect x="13.5" y="3.5" width="7" height="7" rx="1.6" /><rect x="3.5" y="13.5" width="7" height="7" rx="1.6" /><path d="M13.8 17.4l1.9 1.9 4-4.2" /></>),
  qc: svg(<><path d="M4 15.5a8 8 0 0 1 16 0" /><path d="M12 15.5l4.2-3.4" /><circle cx="12" cy="15.5" r="1.15" /></>),
  cluster: svg(<><circle cx="7" cy="8" r="2.3" /><circle cx="16.5" cy="7" r="2.3" /><circle cx="8.5" cy="16.5" r="2.3" /><circle cx="17" cy="16" r="2.3" /><path d="M9.1 9.3l5.4 5.1M9.3 8l4.9-.7" /></>),
  annotate: svg(<><path d="M4 11.4V5.6A1.6 1.6 0 0 1 5.6 4h5.8a1.6 1.6 0 0 1 1.13.47l7 7a1.6 1.6 0 0 1 0 2.26l-5.8 5.8a1.6 1.6 0 0 1-2.26 0l-7-7A1.6 1.6 0 0 1 4 11.4z" /><circle cx="8.2" cy="8.2" r="1.35" /></>),
  spatial: svg(<><path d="M12 3l7.5 4.3v9.4L12 21l-7.5-4.3V7.3z" /><circle cx="12" cy="12" r="2.1" /></>),
  report: svg(<><path d="M6 3h7.5L19 8.5V20a1 1 0 0 1-1 1H6a1 1 0 0 1-1-1V4a1 1 0 0 1 1-1z" /><path d="M13.3 3.2V8.2H18.6" /><path d="M8.5 13h7M8.5 16.4h5" /></>),
};
const MIC = svg(<><rect x="9" y="3" width="6" height="11" rx="3" /><path d="M5.5 11.5a6.5 6.5 0 0 0 13 0" /><path d="M12 18v3" /></>);

interface Step { key: string; n: string; label: string; title: string; sub: string; hint: string; run: string | null; color: string; }
const STEPS: Step[] = [
  { key: "load", n: "01", label: "Load", title: "Your specimen, in view", sub: "A public 10x Prime 5K breast Xenium section (~20k-cell demo subset).", hint: "Bring a section in - a public breast Xenium demo.", run: null, color: "cell_type" },
  { key: "panel", n: "02", label: "Panel check", title: "What the panel can see", sub: "Which cell types this panel can actually resolve, before you over-trust a label.", hint: "Which cell types this panel can actually resolve.", run: "panel", color: "cell_type" },
  { key: "qc", n: "03", label: "QC", title: "Signal you can trust", sub: "Section QC and the six-layer funnel - lasso a region to QC it on the map.", hint: "Metrics, the six-layer funnel, and lasso region QC.", run: "qc", color: "total_counts" },
  { key: "cluster", n: "04", label: "Cluster", title: "Cells grouped by expression", sub: "Leiden clustering on the neighbour graph.", hint: "Leiden clustering - raise resolution to subdivide.", run: "cluster", color: "leiden" },
  { key: "annotate", n: "05", label: "Annotate", title: "Named, with confidence", sub: "Marker + consensus cell types with an honest per-cell confidence.", hint: "Marker + consensus types with an honest confidence.", run: "annotate", color: "cell_type" },
  { key: "spatial", n: "06", label: "Spatial + niches", title: "The tissue in context", sub: "Immune exclusion, neighborhoods, and TME niches.", hint: "Immune exclusion, neighborhoods, and TME niches.", run: "spatial", color: "niche" },
  { key: "report", n: "07", label: "Report", title: "The whole run, packaged", sub: "A re-runnable script and a shareable summary.", hint: "Re-runnable script, annotated .h5ad, shareable HTML.", run: null, color: "cell_type" },
];
const CHIPS = ["How many cells couldn't be confidently labeled, and why?", "Color the map by CD8A expression", "Which cell types neighbor the tumour?"];

export default function App() {
  const [sess, setSess] = useState<api.Session | null>(null);
  const [stepKey, setStepKey] = useState("load");
  const [colorBy, setColorBy] = useState("cell_type");
  const [basis, setBasis] = useState<"spatial" | "umap">("spatial");   // map coords: spatial (hero) or UMAP
  const [pts, setPts] = useState<api.Points | null>(null);
  // Client cache of colourings by `${basis}:${color_by}` -> the recolour payload (rgb/legend/ramp).
  // A pure recolour reuses the current pts' x/y/ids, so a revisited colouring is a synchronous swap
  // (0 network). Same lifetime as capCache: invalidate() clears it on any section mutation/change.
  const ptsCache = useRef<Map<string, Pick<api.Points, "rgb" | "legend" | "ramp" | "color_by">>>(new Map());
  type Emph = { match: [number, number, number]; paint: [number, number, number] };
  const [emphasize, setEmphasize] = useState<Emph[] | null>(null);   // heatmap-hover highlight
  const [sessionLost, setSessionLost] = useState(false);   // server restarted / session timed out -> offer a reload
  const lastLoad = useRef<null | (() => Promise<void>)>(null);   // how to re-establish the current section
  // A hovered enrichment-heatmap pair: match each type by its on-map colour but repaint it a fixed
  // red/blue, so the two stay distinct even when their own map colours are similar (a single type on
  // the diagonal keeps its own colour). Canvas fades everything else. No-op unless the map is coloured
  // by a matching categorical field.
  const onEmphasize = (labels: string[] | null) => {
    if (!labels || !pts) { setEmphasize(null); return; }
    const PAIR: [number, number, number][] = [[239, 68, 68], [56, 152, 255]];   // red / blue
    const solo = labels.length === 1;
    const items = labels.map((l, i) => {
      const col = pts.legend.find((le) => le.label === l)?.color as [number, number, number] | undefined;
      return col ? { match: col, paint: solo ? col : (PAIR[i] ?? PAIR[0]) } : null;
    }).filter(Boolean) as Emph[];
    setEmphasize(items.length ? items : null);
  };
  const [answer, setAnswer] = useState<{ text: string; tag?: string; thinking?: boolean }>({ text: "Ask about your section in plain language - I run the real analysis and recolour the map." });
  const [busy, setBusy] = useState(false);
  const [runningStep, setRunningStep] = useState<string | null>(null);
  const [progress, setProgress] = useState<api.Progress | null>(null);
  const [note, setNote] = useState<string | null>(null);
  const [selection, setSelection] = useState<number[]>([]);
  const [figures, setFigures] = useState<api.Artifact[]>([]);
  const [zoomSrc, setZoomSrc] = useState<string | null>(null);   // a copilot figure enlarged to a modal
  useEffect(() => {
    if (!zoomSrc) return;
    const onKey = (e: KeyboardEvent) => { if (e.key === "Escape") setZoomSrc(null); };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [zoomSrc]);

  // Dictate a question. Click to record, click again to stop -> the blob is transcribed on the
  // server (faster-whisper, local) and dropped into the ask box for you to check before sending.
  // Browsers only expose a microphone on a SECURE context (https://, or localhost). On a plain-http
  // origin `navigator.mediaDevices` is undefined outright - so say why, rather than fail on click.
  const micSupported = typeof navigator !== "undefined" && !!navigator.mediaDevices?.getUserMedia;
  const [recording, setRecording] = useState(false);
  const [transcribing, setTranscribing] = useState(false);
  const recRef = useRef<MediaRecorder | null>(null);
  const startRec = async () => {
    setNote(null);
    let stream: MediaStream;
    try { stream = await navigator.mediaDevices.getUserMedia({ audio: true }); }
    catch (e: any) { setNote("microphone unavailable: " + (e?.message || e)); return; }
    const chunks: Blob[] = [];
    const rec = new MediaRecorder(stream);   // let the browser pick: webm/opus, ogg/opus or mp4 - all decode server-side
    rec.ondataavailable = (e) => { if (e.data.size) chunks.push(e.data); };
    rec.onstop = async () => {
      stream.getTracks().forEach((t) => t.stop());   // release the mic, drop the recording indicator
      setRecording(false); setTranscribing(true);
      try {
        const text = await api.transcribe(new Blob(chunks, { type: rec.mimeType }));
        if (text && inputRef.current) { inputRef.current.value = text; inputRef.current.focus(); }
        else if (!text) setNote("Didn't catch that - try again, closer to the microphone.");
      } catch (e: any) { setNote("transcription failed: " + (e?.message || e)); }
      setTranscribing(false);
    };
    recRef.current = rec; rec.start(); setRecording(true);
  };
  const toggleMic = () => { if (recording) recRef.current?.stop(); else startRec(); };
  useEffect(() => () => { if (recRef.current?.state === "recording") recRef.current.stop(); }, []);
  const inputRef = useRef<HTMLInputElement>(null);
  const step = STEPS.find((s) => s.key === stepKey)!;

  const [summary, setSummary] = useState<any>(null);
  const fetchPoints = async (sid: string, cb: string, b: string = basis) => { try { const p = await api.getPoints(sid, cb, b); setPts(p); ptsCache.current.set(`${p.basis ?? b}:${p.color_by}`, { rgb: p.rgb, legend: p.legend, ramp: p.ramp, color_by: p.color_by }); return p; } catch (e) { /* keep last */ } };
  const fetchSummary = async (sid: string) => { try { setSummary(await api.getSummary(sid)); } catch (e) { /* keep last */ } };
  // Recolour by an obs field WITHOUT refetching x/y/ids (positions don't change on a recolour, so this
  // also spares deck.gl a getPosition GPU re-upload). Dedups a no-op; serves a previously-seen colouring
  // from ptsCache synchronously (0 network); else fetches the colours-only payload (~2.3x smaller) and
  // merges. Full fetch only when there are no points to merge into (first load / after a basis change).
  const recolor = async (sid: string, cb: string, b: string = basis) => {
    if (pts && b === basis) {
      const hit = ptsCache.current.get(`${b}:${cb}`);
      if (hit && hit.rgb.length === pts.rgb.length) {                         // canonical colouring cached -> apply synchronously (0 network)
        setColorBy(cb);
        // pts.rgb === hit.rgb means it is already shown (skip a no-op re-render); a masked preview left
        // a different rgb under the same field, so restoring the cached canonical here also un-masks it.
        if (pts.rgb !== hit.rgb) setPts(prev => prev ? { ...prev, rgb: hit.rgb, color_by: hit.color_by, legend: hit.legend, ramp: hit.ramp } : prev);
        return;
      }
      const p = await api.getPoints(sid, cb, b, true);                        // colours-only fetch + merge (reuse x/y/ids, ~2.3x smaller)
      ptsCache.current.set(`${b}:${cb}`, { rgb: p.rgb, legend: p.legend, ramp: p.ramp, color_by: p.color_by });
      setColorBy(cb); setPts(prev => (prev && p.n === prev.n) ? { ...prev, rgb: p.rgb, color_by: p.color_by, legend: p.legend, ramp: p.ramp } : prev); return;
    }
    setColorBy(cb); await fetchPoints(sid, cb, b);                            // no points yet / basis change -> full fetch
  };

  useEffect(() => { (async () => {
    setBusy(true);
    // Prefer the precomputed real-data demo; if its cache isn't configured on the server
    // (SPATIALSCRIBE_DEMO_CACHE unset), fall back to the offline synthetic so the app always loads.
    let s: api.Session;
    try { s = await api.loadDemo(); }
    catch { try { s = await api.loadSynthetic(); setNote("Demo cache not configured - loaded the synthetic melanoma demo."); }
            catch (e: any) { setNote("could not load a section: " + (e?.message || e)); setBusy(false); return; } }
    setSess(s);
    const first = s.obs_fields.includes("cell_type") ? "cell_type" : (s.obs_fields[0] || "total_counts");
    setColorBy(first); await fetchPoints(s.session_id, first); await fetchSummary(s.session_id); setBusy(false);
  })(); }, []);

  // A new / resized section invalidates any region selection (indices no longer map to the same cells).
  useEffect(() => { setSelection([]); setBasis("spatial"); }, [sess?.n_obs]);

  // Server restart / session timeout -> every call 404s "unknown session". Show a reload prompt (see
  // api.onSessionLost) and re-establish the SAME section if we know how, else fall back to a page reload.
  useEffect(() => { api.onSessionLost(() => setSessionLost(true)); }, []);
  const recoverSession = async () => {
    setSessionLost(false);
    if (lastLoad.current) { try { await lastLoad.current(); return; } catch { /* fall through to reload */ } }
    window.location.reload();
  };

  // Poll the server's coarse progress while a step runs so RunProgress can show a determinate bar.
  // The run endpoints are sync defs served from FastAPI's threadpool, so this GET answers mid-run.
  useEffect(() => {
    const sid = sess?.session_id;
    if (!runningStep || !sid) return;
    let alive = true;
    const t = setInterval(async () => {
      try { const p = await api.getProgress(sid); if (alive) setProgress(p); } catch { /* keep last */ }
    }, 400);
    return () => { alive = false; clearInterval(t); };   // stop on completion and on unmount
  }, [runningStep, sess?.session_id]);

  const gotoStep = async (s: Step) => {
    setStepKey(s.key); setNote(null);
    if (sess) { const cb = sess.obs_fields.includes(s.color) ? s.color : colorBy; await recolor(sess.session_id, cb); }
  };
  const doRun = async () => {
    if (!sess || !step.run) return;
    invalidate(sess.session_id);   // a run mutates the section - drop its memoized capability results
    setBusy(true); setRunningStep(step.key); setProgress(null); setNote(null);
    const r = await api.runStep(sess.session_id, step.run);
    if (r.error) setNote(r.error.hint || r.error.message || "step failed");
    const ns = { ...sess, ...r.state }; setSess(ns);
    const cb = r.state.obs_fields.includes(step.color) ? step.color : colorBy;
    setColorBy(cb); await fetchPoints(sess.session_id, cb); await fetchSummary(sess.session_id);
    setRunningStep(null); setBusy(false);
  };
  const changeColor = async (cb: string) => { setEmphasize(null); if (sess) await recolor(sess.session_id, cb, basis); else setColorBy(cb); };
  // A rail-step panel finished its first-visit auto-compute (QC funnel, panel report): green the step.
  // Idempotent, and the backend returns fresh state so `ran` (and the rail) update.
  const markRan = async (stepKey: string) => {
    if (!sess || (sess.ran || []).includes(stepKey)) return;
    try { const st = await api.markRan(sess.session_id, stepKey); setSess((s) => (s ? { ...s, ...st } : s)); } catch { /* keep */ }
  };
  // Transient hover preview: recolour the canvas WITHOUT touching colorBy, so the dropdown/HUD keep
  // the real field and passing null restores it. A sequence guard makes the LATEST hover win
  // regardless of fetch-resolution order - without it a slow enter-fetch could land after the
  // leave-fetch and strand the map on a program colour (or vice-versa: revert while still hovering).
  const previewSeq = useRef(0);
  // The per-cell program score is precomputed server-side (obs['program_score_k']); the only cost on
  // hover is refetching the recoloured points, which lags on a big section. Show a tiny "colouring"
  // chip, but only if the fetch is actually slow (>140ms) so a fast recolour does not flicker it -
  // without any feedback the map looks frozen and the user thinks the hover is broken.
  const [previewBusy, setPreviewBusy] = useState(false);
  const previewColor = async (field: string | null, onlyType?: string, posMin?: number) => {
    if (!sess) return;
    const seq = ++previewSeq.current;
    const t = field ? setTimeout(() => { if (seq === previewSeq.current) setPreviewBusy(true); }, 140) : undefined;
    try {
      const p = await api.getPoints(sess.session_id, field ?? colorBy, basis, true, onlyType, posMin);   // colours-only recolour (optionally masked to one cell type / thresholded)
      // Reuse the existing x/y/ids (idx is deterministic) and swap only the colours - ~2x faster than a
      // full refetch. If the point count no longer matches (section changed) keep the current view; a
      // full fetch will reconcile it.
      if (seq === previewSeq.current)
        setPts(prev => (prev && p.n === prev.n) ? { ...prev, rgb: p.rgb, color_by: p.color_by, legend: p.legend, ramp: p.ramp } : prev);
    } catch { /* keep last */ }
    finally { clearTimeout(t); if (seq === previewSeq.current) setPreviewBusy(false); }
  };
  const changeBasis = async (b: "spatial" | "umap") => { setBasis(b); if (sess) await fetchPoints(sess.session_id, colorBy, b); };
  // Anything that changes the section invalidates BOTH caches of a capability result: the module-level
  // capCache AND the panels' own useState copies. Clearing capCache alone is not enough - a panel that
  // never unmounted still holds the pre-mutation value in local state and would show it as current
  // (subcluster from the curate tab, return to states, read stale programs). Bumping `dataEpoch` is
  // part of the panel's React key, so the panel remounts and its local state goes with it; the panels
  // remember their own subtab across the remount.
  const [dataEpoch, setDataEpoch] = useState(0);
  const invalidate = (sid?: string) => {
    if (sid) clearSession(sid); else clearAll();
    ptsCache.current.clear();   // colourings depend on the section's data - a mutation restains them
    setDataEpoch((e) => e + 1);
  };
  // Refresh after a mutating panel action (rename / merge / region filter / subcluster).
  const onMutate = async () => {
    if (!sess) return;
    try {
      const st = await api.getState(sess.session_id);
      setSess({ ...sess, ...st });
      // Land the FRESH summary BEFORE bumping dataEpoch (the remount). A summary-seeded control - the
      // cluster resolution slider - otherwise re-mounts on the STALE resolution and then a useEffect
      // snaps it to the new one: the visible flicker (chosen -> previous -> chosen). Invalidating after
      // the summary lands means the remount already reads the new value, so there is nothing to snap.
      await fetchSummary(sess.session_id);
      invalidate(sess.session_id);   // rename / merge / region filter / subcluster / re-cluster: drop stale caches + remount
      const cb = st.obs_fields.includes(colorBy) ? colorBy : (st.obs_fields.includes("cell_type") ? "cell_type" : (st.obs_fields[0] || colorBy));
      setColorBy(cb); await fetchPoints(sess.session_id, cb);
    } catch { /* keep */ }
  };
  const loadSynthetic = async () => {
    lastLoad.current = loadSynthetic;
    invalidate();   // a different section - no cached result from the old one may leak in
    setBusy(true); setNote(null); setSelection([]);
    try {
      const s = await api.loadSynthetic(); setSess(s);
      const first = s.obs_fields.includes("cell_type") ? "cell_type" : (s.obs_fields[0] || "total_counts");
      setColorBy(first); await fetchPoints(s.session_id, first); await fetchSummary(s.session_id);
      setNote(`Loaded synthetic melanoma (${s.n_obs.toLocaleString()} cells)`);
    } catch (e: any) { setNote("load failed: " + (e?.message || e)); }
    setBusy(false);
  };
  // Swap between the bundled demo sections ("breast" = the shallow hard case, "atera5k" = the deep
  // positive control). Same shape as loadSynthetic; the server holds the processed cache.
  const loadDemoNamed = async (name: string) => {
    lastLoad.current = () => loadDemoNamed(name);
    invalidate();   // a different section - no cached result from the old one may leak in
    setBusy(true); setNote(null); setSelection([]);
    try {
      const s = await api.loadDemo(name); setSess(s);
      const first = s.obs_fields.includes("cell_type") ? "cell_type" : (s.obs_fields[0] || "total_counts");
      setColorBy(first); await fetchPoints(s.session_id, first); await fetchSummary(s.session_id);
      setNote(`Loaded ${name} demo (${s.n_obs.toLocaleString()} cells)`);
    } catch (e: any) { setNote("load failed: " + (e?.message || e)); }
    setBusy(false);
  };
  const loadSection = async (path: string, tissue: string) => {
    lastLoad.current = () => loadSection(path, tissue);
    invalidate();   // a different section - no cached result from the old one may leak in
    setBusy(true); setNote(null); setAnswer({ text: "Loading your section - this can take a moment for a full run." });
    try {
      const s = await api.loadSection(path, tissue);
      setSess(s);
      pushRecent(path, tissue);   // remember this section for the Load tab's "recently loaded" list
      const first = s.obs_fields.includes("cell_type") ? "cell_type" : (s.obs_fields[0] || "total_counts");
      setColorBy(first); await fetchPoints(s.session_id, first); await fetchSummary(s.session_id);
      setNote(`Loaded ${s.n_obs.toLocaleString()} cells from ${path}`);
    } catch (e: any) {
      setNote("load failed: " + (e?.message || e));
    }
    setBusy(false);
  };
  const send = async (prompt: string) => {
    if (!sess || !prompt.trim()) return;
    if (!sess.has_key) { setAnswer({ text: "Copilot is disabled: ANTHROPIC_API_KEY is not set on the server." }); return; }
    setBusy(true); setAnswer({ text: "thinking …", thinking: true }); setFigures([]); setEmphasize(null);
    let viewShown = false;   // did a view tool repaint the map this turn? (a load then needn't re-repaint)
    try {
      await api.streamCopilot(sess.session_id, prompt, async (ev) => {
        // `thinking` drives the animated violet "running"-style pill (status / pre-answer phase); it
        // clears the moment real answer tokens start streaming, so the answer itself renders as text.
        if (ev.type === "status") setAnswer({ text: ev.text, thinking: true });
        else if (ev.type === "map_view") {
          // Drive the MAIN specimen canvas: recolour, and if the copilot named a category to spotlight
          // ("where are the T cells"), emphasise it (the canvas fades the rest) instead of embedding a
          // thumbnail. Use the freshly-fetched legend, not the stale `pts` closure.
          viewShown = true;
          const cb = ev.color_by; setColorBy(cb);
          const p = await fetchPoints(sess.session_id, cb);
          const hl = ev.highlight as string | undefined;
          const m = hl && p?.legend ? p.legend.find((l: any) => String(l.label) === String(hl)) : undefined;
          setEmphasize(m ? [{ match: m.color, paint: m.color }] : null);
          setAnswer((a) => ({ text: a?.text || "", tag: hl ? `showing ${hl} on the map` : `recoloured the map by ${cb}`, thinking: a?.thinking }));
        }
        else if (ev.type === "figures") setFigures(ev.figures || []);
        else if (ev.type === "token") setAnswer((a) => ({ text: ev.text, tag: a?.tag, thinking: false }));
        else if (ev.type === "done") {
          const ns = { ...sess, ...ev.state }; setSess(ns);
          if (ev.loaded) {
            // The copilot swapped in a NEW section this turn: drop the old (demo) section's cached
            // panel results, and - unless a view tool already repainted the map - repaint it from the
            // new section so the old points stop lingering on the spatial plot.
            invalidate(sess.session_id);
            if (!viewShown) {
              const cb = ns.obs_fields.includes("cell_type") ? "cell_type" : (ns.obs_fields[0] || "total_counts");
              setColorBy(cb); setEmphasize(null); setSelection([]);
              await fetchPoints(sess.session_id, cb);
            }
          } else if (ev.mutated) {
            // The copilot MUTATED the section in place (merge / relabel / subcluster / self_heal /
            // cell-states): drop cached panel results and repaint the map, else the change is applied
            // server-side but never shows (the "optimize merging did nothing" report). Mirrors the
            // wizard's onMutate. Repaint by the annotation column so a merge/relabel is visible.
            invalidate(sess.session_id);
            if (!viewShown) {
              const cb = ns.obs_fields.includes("cell_type_final") ? "cell_type_final"
                : (ns.obs_fields.includes("cell_type") ? "cell_type" : colorBy);
              setColorBy(cb); setEmphasize(null);
              await fetchPoints(sess.session_id, cb);
            }
          }
          await fetchSummary(sess.session_id);
        }
      });
    } catch (e: any) { setAnswer({ text: "copilot error: " + (e?.message || e), thinking: false }); }
    setBusy(false);
  };
  const onSelect = (idx: number[]) => { setSelection(idx); setNote(idx.length ? `Selected ${idx.length.toLocaleString()} cells - QC or filter on the QC tab` : null); };
  // Click the model tag -> flip the copilot to the OTHER configured provider (local vLLM <-> Anthropic).
  const toggleLlm = async () => {
    const provs = sess?.llm?.providers || [];
    if (!sess || provs.length < 2) return;
    const next = provs.find((p) => p !== sess.llm?.provider) || provs[0];
    try { const r = await api.setLlmProvider(next); if (r.ok && r.llm) setSess({ ...sess, llm: r.llm }); } catch { /* keep */ }
  };

  // The rail greens off the steps the USER ran this session (`ran`), not the data-derived `done` -
  // otherwise a bundled demo (shipped pre-processed) reads as already-run on first load. `load` is
  // always done once a section is open; `report` completes once every compute step has been run.
  const ran = useMemo(() => new Set<string>(sess?.ran || []), [sess]);
  const COMPUTE_STEPS = ["panel", "qc", "cluster", "annotate", "spatial"];
  const stepDone = (k: string) =>
    k === "load" ? !!sess : k === "report" ? COMPUTE_STEPS.every((c) => ran.has(c)) : ran.has(k);
  const doneCount = sess ? STEPS.filter((s) => stepDone(s.key)).length : 0;
  const hud = useMemo(() => ({ k: "specimen", v: sess ? sess.n_obs.toLocaleString() : "-", sub: `cells · by ${colorBy}` }), [sess, colorBy]);
  // Leiden badge reads the authoritative cluster count (the legend is capped at 35 for display, so
  // a >35-cluster section would otherwise under-report on the map).
  const nLeiden = summary?.cluster?.n_clusters;
  const stat = pts ? (colorBy === "leiden" ? `${nLeiden ?? (new Set(pts.legend.map((l) => l.label)).size || pts.legend.length)} clusters` : (pts.legend.length ? `${pts.legend.length} ${colorBy === "niche" ? "niches" : "categories"}` : null)) : null;

  return (
    <div className="app">
      {sessionLost && (
        <div className="runoverlay" style={{ zIndex: 60 }}>
          <div className="runlabel">Session expired</div>
          <div className="runquip" style={{ maxWidth: 340, textAlign: "center" }}>
            The server restarted (or your session timed out). Reload to bring your section back.
          </div>
          <button className="pbtn pri" style={{ marginTop: 12 }} onClick={recoverSession}>Reload section</button>
        </div>
      )}
      <aside className="rail">
        <div className="brand">SpatialScribe</div>
        <div className="brandsub">spatial-transcriptomics copilot</div>
        <div className="railmeta">{sess ? `${sess.device || "GPU"} · ${sess.n_obs.toLocaleString()} cells · ${sess.tissue}` : "loading …"}</div>
        {sess?.panel_name && <div className="railmeta" title="Xenium panel (from experiment.xenium)" style={{ marginTop: 2 }}>panel · {sess.panel_name}</div>}
        <div className="railbar"><div className="top"><span>run progress</span><span>{doneCount} / {STEPS.length} done</span></div>
          <div className="track"><div className="fill" style={{ width: `${Math.round((doneCount / STEPS.length) * 100)}%` }} /></div></div>
        <hr />
        <div className="steps">
          {STEPS.map((s) => {
            const cls = "step" + (s.key === stepKey ? " active" : (stepDone(s.key) ? " done" : ""));
            return (
              <div key={s.key} className={cls} onClick={() => gotoStep(s)}>
                <span className="sdot">{ICONS[s.key]}</span>
                <span className="slab">{s.label}</span>
                <span className="shint">{s.hint}</span>
              </div>
            );
          })}
        </div>
      </aside>

      <div className="main">
        <div className="topbar">
          <div>
            <span className="kicker">{step.n} · {step.label}</span>
            <div className="title">{step.title}</div>
            <div className="sub">{note ? <span style={{ color: "var(--violet)" }}>{note}</span> : step.sub}</div>
          </div>
          <div className="tbctl">
            {step.run && <button className="btn primary" onClick={doRun} disabled={busy}>{busy ? "running …" : `Run ${step.label}`}</button>}
            {sess?.done?.cluster && (
              <div className="seg" role="group" aria-label="map view">
                {(["spatial", "umap"] as const).map((b) => (
                  <button key={b} className={"segbtn" + (basis === b ? " on" : "")} onClick={() => changeBasis(b)}>
                    {b === "spatial" ? "Spatial" : "UMAP"}
                  </button>
                ))}
              </div>
            )}
            {sess && (
              <select className="pill" value={colorBy} onChange={(e) => changeColor(e.target.value)}>
                {sess.obs_fields.map((f) => <option key={f} value={f}>{f}</option>)}
              </select>
            )}
          </div>
        </div>
        <div className="stage">
          <div className="panelcol">
            {PANEL_FOR[step.key] && summary
              ? (() => { const P = PANEL_FOR[step.key]; return <P key={`${sess?.session_id ?? "none"}:${dataEpoch}`} summary={summary} sid={sess?.session_id} sess={sess} selection={selection} onMutate={onMutate} onColor={changeColor} onPreviewColor={previewColor} onEmphasize={onEmphasize} onRan={() => markRan(step.key)} onGoStep={(k: string) => setStepKey(k)} />; })()
              : <LoadPanel sess={sess} sid={sess?.session_id} onMutate={onMutate} onColor={changeColor}
                           onLoad={step.key === "load" ? loadSection : undefined}
                           onLoadSynthetic={step.key === "load" ? loadSynthetic : undefined}
                           onLoadDemo={step.key === "load" ? loadDemoNamed : undefined} busy={busy} />}
          </div>
          {pts ? (
            <div className="viewport">
              <SpecimenCanvas x={pts.x} y={pts.y} rgb={pts.rgb} ids={pts.ids} n={pts.n} hud={hud} legend={pts.legend} stat={stat} onSelect={onSelect} emphasize={emphasize ?? undefined} />
              {pts.color_by && (!pts.legend || pts.legend.length === 0) && <FeatureScale field={pts.color_by} ramp={pts.ramp} />}
              {previewBusy && <div className="previewchip"><i />colouring …</div>}
              {runningStep && (
                <RunProgress label={`Running ${step.label} …`}
                             quips={runningStep === "annotate" ? ANNOTATE_QUIPS : GENERIC_QUIPS}
                             frac={progress?.frac ?? 0} serverLabel={progress?.label ?? ""} />
              )}
            </div>
          ) : <div className="empty">{busy ? "loading specimen …" : "no section"}</div>}
        </div>

        <div className="drawer">
          <div className="cav">S</div>
          <div className="composer">
            {/* Fixed-height row (see .answer): the thinking pill is taller than a plain answer line,
                so both live in a min-height flex row and the pane no longer jumps when it appears/goes. */}
            <div className="answer" title={answer.thinking ? undefined : answer.text}>{answer.thinking
              ? <span className="running phrase"><i />{answer.text}</span>
              : <span className="txt">{answer.text}{answer.tag && <span className="mtag"> · {answer.tag}</span>}</span>}</div>
            {figures.length > 0 && (
              <div style={{ display: "flex", gap: 8, overflowX: "auto", padding: "2px 0" }}>
                {figures.map((f, i) => f.png && <img key={i} src={f.png} alt={f.title || "figure"} title={f.title || "click to enlarge"}
                  onClick={() => setZoomSrc(f.png!)}
                  style={{ height: 128, borderRadius: 9, border: "1px solid var(--line2)", flex: "none", cursor: "zoom-in" }} />)}
              </div>
            )}
            {zoomSrc && (
              <div className="figmodal" onClick={() => setZoomSrc(null)}>
                <button className="close" onClick={() => setZoomSrc(null)}>close</button>
                <img src={zoomSrc} alt="figure enlarged" />
              </div>
            )}
            <div className="cinput">
              <input ref={inputRef} placeholder={sess?.has_key ? "ask the copilot - it drives the map …" : "set ANTHROPIC_API_KEY to enable the copilot"}
                     onKeyDown={(e) => { if (e.key === "Enter") { send((e.target as HTMLInputElement).value); (e.target as HTMLInputElement).value = ""; } }} />
              {/* Which model backs the copilot - tucked into the bar (violet dot = local/self-hosted
                  endpoint, green = Anthropic API). Clickable to switch provider when >1 is configured. */}
              {sess?.llm?.available && sess.llm.model && (() => {
                const local = sess.llm.provider === "openai";
                const canToggle = (sess.llm.providers?.length || 0) > 1;
                return (
                  <span className={"llmtag" + (canToggle ? " sw" : "")}
                        onClick={canToggle ? toggleLlm : undefined}
                        title={canToggle ? "click to switch the copilot LLM (local vLLM <-> Anthropic)"
                                         : (local ? "local / self-hosted LLM endpoint" : "Anthropic API")}>
                    <i style={{ background: local ? "var(--violet)" : "var(--pass)" }} />
                    {sess.llm.model} · {local ? "local" : "Anthropic"}
                  </span>
                );
              })()}
              {sess?.has_stt && (
                <button className={"mic" + (recording ? " rec" : "")} onClick={toggleMic}
                        disabled={!micSupported || transcribing} aria-pressed={recording}
                        title={!micSupported
                          ? "Dictation needs a secure page - open the app over https:// or via localhost"
                          : transcribing ? "Transcribing …" : recording ? "Stop and transcribe" : "Dictate your question"}>
                  {MIC}
                </button>
              )}
              <button className="snd" onClick={() => { if (inputRef.current) { send(inputRef.current.value); inputRef.current.value = ""; } }}>↑</button>
            </div>
            <div className="chips">{CHIPS.map((c) => <span key={c} className="chip" onClick={() => send(c)}>{c}</span>)}</div>
          </div>
        </div>
      </div>
    </div>
  );
}
