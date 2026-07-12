// Thin fetch wrappers over the FastAPI backend. Same-origin /api (vite proxies to the API in dev;
// in production the container serves the SPA and reverse-proxies /api to the GPU backend).
//
// Every URL is prefixed by u() with the app's base path (vite `base`: "/" in dev,
// "/app/spatialscribe/" on a reverse proxy) so the SPA works under any mount point without CORS.
const BASE = (import.meta.env.BASE_URL || "/").replace(/\/$/, "");
const u = (p: string) => BASE + p;

// The backend keeps sessions IN MEMORY, so a restart (dev redeploys are frequent) or a timeout drops
// them and every later call 404s "unknown session". Surface that as a recoverable state (a reload
// prompt) instead of a bare "Error" on every button. App registers the handler.
let sessionLostHandler: (() => void) | null = null;
export const onSessionLost = (fn: () => void) => { sessionLostHandler = fn; };
const json = (r: Response) => {
  if (!r.ok) return r.text().then((t) => {
    if (r.status === 404 && /unknown session/i.test(t)) sessionLostHandler?.();
    return Promise.reject(new Error(t));
  });
  return r.json();
};
const post = (url: string, body?: any) =>
  fetch(u(url), { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body ?? {}) }).then(json);

export interface State {
  n_obs: number; n_vars: number; tissue: string;
  // null = unanswered, so the engine falls back to matching the tissue string against a keyword list
  // that fires on "normal breast" and misses "glioblastoma". An explicit answer always wins.
  is_tumour?: boolean | null;
  // `done` = an output exists (data-derived; gates functional UI like the UMAP toggle).
  // `ran` = the steps the USER ran this session; the rail greens off `ran` so a pre-processed demo
  // does not show as already-done on first load.
  obs_fields: string[]; done: Record<string, boolean>; ran?: string[]; has_key: boolean;
  // which model backs the copilot (provider "openai" = a local/self-hosted endpoint, "anthropic" = API).
  // `providers` = the backends configured this process; >1 means the model tag can toggle between them.
  llm?: { available: boolean; provider: string; model: string | null; providers?: string[] };
  has_stt?: boolean;   // server has faster-whisper -> the ask box can be dictated
  device?: string;
  panel_name?: string | null;   // e.g. "hMulti_100g (Human, Multi, 480 targets)" for a Xenium run
  demo?: { source: string; role: string } | null;   // set for the bundled demo sections only
  reference?: RefSummary | null;   // the session's single-cell reference (set in Panel-check); null if none
}
export interface RefSummary { name: string | null; n_cells: number; label_key: string | null; n_types: number | null; }
export interface Session extends State { session_id: string; }
export interface Points {
  x: number[]; y: number[]; rgb: number[]; ids: number[]; n: number;
  color_by: string; legend: { label: string; color: [number, number, number] }[];
  basis?: string;   // which coordinates the map is drawn in: "spatial" (default) or "umap"
  ramp?: string[] | null;   // continuous featureplot only: the hex stops (low->high) the map was painted with
}

// name: "breast" (Xenium 5K, the hard case) | "atera5k" (the deep positive control). See docs/DATASETS.md.
export const loadDemo = (name: string = "breast"): Promise<Session> =>
  post(`/api/load_demo?name=${encodeURIComponent(name)}`);
export const listDemos = (): Promise<{ demos: { name: string; available: boolean }[] }> =>
  fetch(u("/api/demos")).then(json);
export const loadSection = (path: string, tissue?: string, isTumour?: boolean | null): Promise<Session> =>
  post("/api/load_section", { path, tissue: tissue || "melanoma", is_tumour: isTumour ?? null });
// Data-step checkbox: does this section contain malignant cells? Gates Cancer-Finder + InSituCNV at
// Annotate. Pass null to clear the answer and fall back to the tissue-keyword heuristic.
export const setIsTumour = (sid: string, isTumour: boolean | null): Promise<State> =>
  post(`/api/${sid}/is_tumour`, { is_tumour: isTumour });
export const getState = (sid: string): Promise<State> => fetch(u(`/api/${sid}/state`)).then(json);
// Switch the copilot LLM backend at runtime (process-wide); returns the resulting llm descriptor.
export const setLlmProvider = (provider: string): Promise<{ ok: boolean; llm: State["llm"] }> =>
  post(`/api/llm/provider`, { provider });
// colorsOnly: a RECOLOUR - the server returns just rgb/legend/ramp (no x/y/ids), the caller reuses
// its existing point positions. ~2x faster than a full points fetch (see backend points()).
// onlyType restricts the colouring to one cell type (the rest fade) - a states-heatmap cell hover
// shows that state's score only on its cell type.
export const getPoints = (sid: string, colorBy: string, basis: string = "spatial", colorsOnly = false, onlyType?: string, posMin?: number): Promise<Points> =>
  fetch(u(`/api/${sid}/points?color_by=${encodeURIComponent(colorBy)}&basis=${encodeURIComponent(basis)}${colorsOnly ? "&colors_only=true" : ""}${onlyType ? "&only_type=" + encodeURIComponent(onlyType) : ""}${posMin != null ? "&pos_min=" + posMin : ""}`)).then(json);
export const runStep = (sid: string, step: string, params?: any): Promise<{ ok: boolean; error: any; state: State }> =>
  post(`/api/${sid}/run/${step}`, { params: params ?? null });
export const askCopilot = (sid: string, prompt: string): Promise<{ reply: string; map_view: { color_by: string } | null; state: State }> =>
  post(`/api/${sid}/copilot`, { prompt });
export const getSummary = (sid: string): Promise<any> => fetch(u(`/api/${sid}/summary`)).then(json);

// Coarse progress of the currently-running step/capability. Both run endpoints are sync `def`, so
// FastAPI serves them from the threadpool and this GET is answerable while a step runs.
export interface Progress { running: boolean; step: string | null; frac: number; label: string }
export const getProgress = (sid: string): Promise<Progress> =>
  fetch(u(`/api/${sid}/progress`)).then(json);

// "Run full analysis": run the whole spine (QC -> gate -> annotate -> niches -> malignant -> self-heal)
// as a BACKGROUND job, then poll pipelineStatus and reload the section when done.
export interface PipelineStatus {
  running: boolean; frac: number; label: string; stage: string | null;
  done: boolean; error: string | null; route: string | null;
  stages: { name: string; status: string }[];
  summary: { ok?: string[]; skipped?: string[]; failed?: string[] };
  state: State | null;
}
export const runPipeline = (
  sid: string, opts?: { tumour?: boolean | null; rctd?: boolean; split?: boolean; resolution?: number },
): Promise<{ started: boolean }> => post(`/api/${sid}/run_pipeline`, opts ?? {});
export const pipelineStatus = (sid: string): Promise<PipelineStatus> =>
  fetch(u(`/api/${sid}/pipeline_status`)).then(json);

// Figure/map artifacts a capability or the copilot emits.
export interface Artifact { kind: "figure" | "map_view"; title?: string; png?: string; color_by?: string; note?: string; }
export interface CapResult { ok: boolean; error: any; value: any; artifacts: Artifact[]; state: State; }

// Run ANY registered capability (qc_funnel, immune_exclusion, neighborhood_enrichment,
// state_by_celltype, malignant_score, discover_programs, subcluster, rejection_reasons, ...).
export const runCap = (sid: string, name: string, params?: any): Promise<CapResult> =>
  post(`/api/${sid}/run_cap/${name}`, { params: params ?? null });

// Green a rail step once its panel finished the first-visit auto-compute; returns fresh state.
export const markRan = (sid: string, step: string): Promise<Partial<State>> =>
  post(`/api/${sid}/ran/${step}`, {});
export const regionQc = (sid: string, indices: number[]): Promise<any> =>
  post(`/api/${sid}/region_qc`, { indices });
export const regionFilter = (sid: string, indices: number[], mode: "exclude" | "keep"): Promise<Session> =>
  post(`/api/${sid}/region_filter`, { indices, mode });
export const renameCelltype = (sid: string, oldName: string, newName: string): Promise<Session> =>
  post(`/api/${sid}/rename_celltype`, { old: oldName, new: newName });
export const panelReport = (sid: string): Promise<any> => fetch(u(`/api/${sid}/panel_report`)).then(json);
export const panelVerdict = (sid: string): Promise<any> => fetch(u(`/api/${sid}/panel_verdict`)).then(json);
export const qcVerdict = (sid: string): Promise<any> => fetch(u(`/api/${sid}/qc_verdict`)).then(json);
export const loadSynthetic = (): Promise<Session> => post("/api/load_synthetic");

// Custom single-cell reference: upload a .h5ad or register a server-side path. Returns the panel
// reference<->match readout so the UI can show the fit immediately.
export interface RefResult {
  ok: boolean; n_ref_cells: number; label_key: string | null; n_labels: number | null;
  match: any; obs_fields: string[]; done: Record<string, boolean>;
}
export const setReference = (sid: string, path: string, labelKey?: string): Promise<RefResult> =>
  post(`/api/${sid}/reference`, { path, label_key: labelKey ?? null });
export const uploadReference = (sid: string, file: File): Promise<RefResult> => {
  const fd = new FormData();
  fd.append("file", file);
  // No Content-Type header - the browser sets the multipart boundary.
  return fetch(u(`/api/${sid}/reference/upload`), { method: "POST", body: fd }).then(json);
};
// Free-text tissue -> auto-choose + load the best-matched pre-computed reference (registry, or a
// live CELLxGENE gget fetch when allowFetch). Returns which reference was chosen, the panel match,
// and the supervised-vs-clustering route. ok=false means nothing suitable was found (cluster instead).
export interface AutoRefResult {
  ok: boolean; n_ref_cells?: number; label_key?: string | null;
  auto: { status: string; chosen: string | null; source?: string; message?: string;
          ranked?: { tissue_key: string; score: number; available: boolean; description?: string }[] };
  match?: any; recommended_mode?: string | null; route?: any;
  obs_fields?: string[]; done?: Record<string, boolean>;
}
export const autoReference = (sid: string, tissue: string, allowFetch = false): Promise<AutoRefResult> =>
  post(`/api/${sid}/reference/auto`, { tissue, allow_fetch: allowFetch });
export const verifyReport = (sid: string, neighborhood = false): Promise<any> =>
  fetch(u(`/api/${sid}/verify_report?neighborhood=${neighborhood}`)).then(json);

// Dictated question -> text. Post the raw MediaRecorder blob; the server transcribes it locally
// (faster-whisper) and the audio never leaves the cluster. "" means the recording was silent.
export const transcribe = async (blob: Blob): Promise<string> => {
  const fd = new FormData();
  fd.append("audio", blob, "question.webm");   // no Content-Type header - the browser sets the boundary
  const r = await fetch(u("/api/stt"), { method: "POST", body: fd });
  return (await json(r)).text as string;
};

export const exportScriptUrl = (sid: string) => u(`/api/${sid}/export/script`);
export const exportH5adUrl = (sid: string) => u(`/api/${sid}/export/h5ad`);
export const exportReportUrl = (sid: string) => u(`/api/${sid}/export/report`);

// Copilot with SSE-over-POST: onEvent gets {type:"status"|"map_view"|"token"|"done", ...}.
export const streamCopilot = async (sid: string, prompt: string, onEvent: (e: any) => void): Promise<void> => {
  const res = await fetch(u(`/api/${sid}/copilot/stream`), {
    method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ prompt }),
  });
  if (!res.ok || !res.body) throw new Error(await res.text().catch(() => "stream failed"));
  const reader = res.body.getReader();
  const dec = new TextDecoder();
  let buf = "";
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    buf += dec.decode(value, { stream: true });
    let i: number;
    while ((i = buf.indexOf("\n\n")) >= 0) {
      const line = buf.slice(0, i).split("\n").find((l) => l.startsWith("data:"));
      buf = buf.slice(i + 2);
      if (line) { try { onEvent(JSON.parse(line.slice(5).trim())); } catch { /* skip */ } }
    }
  }
};
