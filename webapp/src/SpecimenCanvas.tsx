import React, { useEffect, useMemo, useRef, useState } from "react";
import { Deck, OrthographicView } from "@deck.gl/core";
import { ScatterplotLayer } from "@deck.gl/layers";

// The GPU specimen viewport: a deck.gl WebGL scatter (binary attributes -> scales past 1e6 cells),
// a floating glass HUD, native zoom/pan, and a box "select region" that calls onSelect(indices).
type Legend = { label: string; color: [number, number, number] };
export interface SpecimenProps {
  x: number[]; y: number[]; rgb: number[]; ids: number[]; n: number;
  hud?: { k?: string; v?: string; sub?: string };
  legend?: Legend[]; stat?: string | null;
  onSelect?: (indices: number[]) => void;
  emphasize?: { match: [number, number, number]; paint: [number, number, number] }[];   // external highlight (e.g. hover a heatmap pair): match a type's map colour, repaint it `paint` (fixed red/blue so the pair stays distinct even if their own colours are similar) - fade the rest
}

// Typical nearest-neighbour distance in DATA units, estimated from a sample via a hash grid (O(n)
// build + O(sample) query). Drawing dots in COMMON (data) units at ~half this makes them tile
// WITHOUT overlapping at every zoom: the dot-to-spacing ratio is then zoom-invariant, unlike a fixed
// pixel radius (which stays put while cells pack together on zoom-out, so it overlaps). A grid of
// ~sqrt(n) bins/axis puts ~1 cell per bin, so each cell's nearest neighbour is within its 3x3
// neighbourhood; this also tracks tissue shape (empty bins carry no weight) unlike a bbox-area guess.
function cellPitch(x: number[], y: number[], n: number): number {
  if (n < 2) return 1;
  let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
  for (let i = 0; i < n; i++) { const px = x[i], py = y[i]; if (px < minX) minX = px; if (px > maxX) maxX = px; if (py < minY) minY = py; if (py > maxY) maxY = py; }
  const w = (maxX - minX) || 1, h = (maxY - minY) || 1;
  const G = Math.max(1, Math.round(Math.sqrt(n)));
  const cw = w / G, ch = h / G;
  const gx = (px: number) => Math.min(G - 1, Math.max(0, ((px - minX) / cw) | 0));
  const gy = (py: number) => Math.min(G - 1, Math.max(0, ((py - minY) / ch) | 0));
  const bins = new Map<number, number[]>();
  for (let i = 0; i < n; i++) { const k = gy(y[i]) * G + gx(x[i]); const a = bins.get(k); if (a) a.push(i); else bins.set(k, [i]); }
  const S = Math.min(2000, n), step = Math.max(1, (n / S) | 0);
  const d: number[] = [];
  for (let i = 0; i < n; i += step) {
    const bx = gx(x[i]), by = gy(y[i]); let best = Infinity;
    for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
      const nx = bx + dx, ny = by + dy; if (nx < 0 || ny < 0 || nx >= G || ny >= G) continue;
      const cell = bins.get(ny * G + nx); if (!cell) continue;
      for (const j of cell) { if (j === i) continue; const ex = x[j] - x[i], ey = y[j] - y[i]; const dd = ex * ex + ey * ey; if (dd < best) best = dd; }
    }
    if (best < Infinity) d.push(Math.sqrt(best));
  }
  if (!d.length) return Math.sqrt((w * h) / n);
  d.sort((a, b) => a - b);
  return d[d.length >> 1] || Math.sqrt((w * h) / n);   // median nearest-neighbour distance
}

const glass: React.CSSProperties = {
  position: "absolute", background: "rgba(13,15,20,.72)", backdropFilter: "blur(10px)",
  border: "1px solid #2E323D", borderRadius: 13, boxShadow: "0 20px 50px -30px #000",
  color: "#C4C8D2", fontSize: 12, padding: "11px 14px", pointerEvents: "none",
  fontFamily: "'JetBrains Mono', ui-monospace, monospace",
};

export default function SpecimenCanvas({ x, y, rgb, ids, n, hud = {}, legend = [], stat, onSelect, emphasize }: SpecimenProps) {
  const wrapRef = useRef<HTMLDivElement>(null);
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const deckRef = useRef<any>(null);
  const [selectMode, setSelectMode] = useState(false);
  const [box, setBox] = useState<{ x0: number; y0: number; x1: number; y1: number } | null>(null);

  const positions = useMemo(() => {
    const p = new Float32Array(n * 2);
    for (let i = 0; i < n; i++) { p[i * 2] = x[i]; p[i * 2 + 1] = y[i]; }
    return p;
  }, [x, y, n]);
  const colors = useMemo(() => Uint8Array.from(rgb), [rgb]);
  // Two emphasis channels that both fade the rest of the section so a subset pops:
  //  - hi: transient legend-hover (match a category colour)
  //  - selIdx: a persistent box-selection (rendered-point indices) - selected cells also grow
  const [hi, setHi] = useState<[number, number, number] | null>(null);
  const [selIdx, setSelIdx] = useState<Set<number> | null>(null);
  // A new section (load / region filter changes the point count) invalidates a selection.
  useEffect(() => { setSelIdx(null); }, [n]);

  const emph = useMemo<{ mask: Uint8Array; paint: Uint8Array | null } | null>(() => {
    // External emphasis (e.g. hovering a cell-type pair in the enrichment heatmap) takes priority:
    // light up every cell whose colour matches a given type's map colour and repaint it that type's
    // fixed `paint` colour (red/blue), fading the rest. Internal hover (hi) keeps its own colour.
    const ext = emphasize && emphasize.length ? emphasize : (hi ? [{ match: hi, paint: hi }] : null);
    if (ext) {
      const mask = new Uint8Array(n), paint = new Uint8Array(n * 3);
      for (let i = 0; i < n; i++) {
        const r = colors[i * 3], g = colors[i * 3 + 1], b = colors[i * 3 + 2];
        const hit = ext.find((c) => c.match[0] === r && c.match[1] === g && c.match[2] === b);
        if (hit) { mask[i] = 1; paint[i * 3] = hit.paint[0]; paint[i * 3 + 1] = hit.paint[1]; paint[i * 3 + 2] = hit.paint[2]; }
      }
      return { mask, paint };
    }
    if (selIdx && selIdx.size) {
      const mask = new Uint8Array(n); selIdx.forEach((i) => { if (i >= 0 && i < n) mask[i] = 1; }); return { mask, paint: null };
    }
    return null;
  }, [hi, emphasize, selIdx, colors, n]);

  const displayColors = useMemo(() => {
    if (!emph) return colors;
    const { mask, paint } = emph;
    const out = new Uint8Array(colors.length);
    for (let i = 0; i < n; i++) {
      const r = colors[i * 3], g = colors[i * 3 + 1], b = colors[i * 3 + 2];
      if (mask[i]) {
        if (paint) { out[i * 3] = paint[i * 3]; out[i * 3 + 1] = paint[i * 3 + 1]; out[i * 3 + 2] = paint[i * 3 + 2]; }
        else { out[i * 3] = r; out[i * 3 + 1] = g; out[i * 3 + 2] = b; }
      } else { out[i * 3] = (r * 0.14 + 8) | 0; out[i * 3 + 1] = (g * 0.14 + 9) | 0; out[i * 3 + 2] = (b * 0.14 + 13) | 0; }
    }
    return out;
  }, [colors, emph, n]);

  // Refit the camera only when the coordinates ACTUALLY change (a basis switch or a new/filtered
  // section), not on every points re-fetch. App hands us a fresh x/y array on each fetch even when the
  // values are identical (same section, new step), so keying the view on the array identity made
  // useMemo return a new object every step -> deck.gl re-applied initialViewState and snapped the zoom
  // by a hair (px depends on the panel size, which reflows on tab change) = the visible jitter on
  // load -> panel-check. A content signature is stable across identical re-fetches, so the view holds.
  const viewSig = n ? `${n}:${x[0]},${y[0]}:${x[n - 1]},${y[n - 1]}` : "0";
  const initialViewState = useMemo(() => {
    let a = Infinity, b = -Infinity, c = Infinity, d = -Infinity;
    for (let i = 0; i < n; i++) { const px = x[i], py = y[i]; if (px < a) a = px; if (px > b) b = px; if (py < c) c = py; if (py > d) d = py; }
    const cx = (a + b) / 2, cy = (c + d) / 2, span = Math.max(b - a, d - c) || 1;
    const px = wrapRef.current ? Math.min(wrapRef.current.clientWidth, wrapRef.current.clientHeight) : 700;
    return { target: [cx, cy, 0], zoom: Math.log2(px / span) - 0.25 } as any;
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [viewSig]);

  // Dot radius in DATA units = ~0.7 x the typical cell spacing, drawn in "common" units so it scales
  // WITH the zoom (see cellPitch). Sized so dots slightly overlap - a SOLID, filled tissue look with
  // no gaps - while staying zoom-invariant, so it never degenerates into an over-packed blob on
  // zoom-out the way a fixed pixel radius did. radiusMin/MaxPixels only clamp the extremes (a
  // visibility floor zoomed all the way out, a cap so a lone cell zoomed all the way in is not a blob).
  const radius = useMemo(() => 0.7 * cellPitch(x, y, n), [x, y, n]);
  // Selected cells also grow (not just stay bright), so a region reads instantly on a dense map.
  const displayRadii = useMemo(() => {
    if (!selIdx || !selIdx.size) return null;
    const r = new Float32Array(n); r.fill(radius);
    selIdx.forEach((i) => { if (i >= 0 && i < n) r[i] = radius * 3; });
    return r;
  }, [selIdx, n, radius]);

  useEffect(() => {
    const attrs: any = { getPosition: { value: positions, size: 2 }, getFillColor: { value: displayColors, size: 3 } };
    if (displayRadii) attrs.getRadius = { value: displayRadii, size: 1 };
    const layer = new ScatterplotLayer({
      id: "cells",
      data: { length: n, attributes: attrs },
      radiusUnits: "common", getRadius: radius, radiusMinPixels: 0.6, radiusMaxPixels: 24,
      pickable: true, stroked: false, opacity: 0.92,
      updateTriggers: { getFillColor: displayColors, getRadius: displayRadii },
    });
    if (!deckRef.current) {
      deckRef.current = new Deck({
        canvas: canvasRef.current!, views: [new OrthographicView({ flipY: true })],
        initialViewState, controller: { dragPan: !selectMode, scrollZoom: true, doubleClickZoom: true },
        layers: [layer], getCursor: () => (selectMode ? "crosshair" : "grab"),
      });
    } else {
      // Pass initialViewState so a coordinate-basis switch (Spatial <-> UMAP) refits the camera to the
      // new range. deck.gl deep-compares it and only resets when it changed; its useMemo deps are
      // [x,y,n], so colour/hover/selection re-renders keep the same object -> user pan/zoom is preserved.
      deckRef.current.setProps({ initialViewState, layers: [layer], controller: { dragPan: !selectMode, scrollZoom: true, doubleClickZoom: true } });
    }
  }, [positions, displayColors, displayRadii, n, radius, selectMode, initialViewState]);

  useEffect(() => () => { deckRef.current?.finalize?.(); deckRef.current = null; }, []);

  const rel = (e: React.PointerEvent) => { const r = (e.currentTarget as HTMLElement).getBoundingClientRect(); return { x: e.clientX - r.left, y: e.clientY - r.top }; };
  const onDown = (e: React.PointerEvent) => {
    if (!selectMode) return;
    e.currentTarget.setPointerCapture?.(e.pointerId);   // keep move/up on this overlay, not deck's canvas
    const p = rel(e); setBox({ x0: p.x, y0: p.y, x1: p.x, y1: p.y });
  };
  const onMove = (e: React.PointerEvent) => { if (!selectMode || !box) return; const p = rel(e); setBox({ ...box, x1: p.x, y1: p.y }); };
  const onUp = () => {
    if (!selectMode || !box) { setBox(null); return; }
    const x0 = Math.min(box.x0, box.x1), y0 = Math.min(box.y0, box.y1), w = Math.abs(box.x1 - box.x0), h = Math.abs(box.y1 - box.y0);
    let picked: number[] = [];
    if (w > 3 && h > 3 && deckRef.current) {
      const objs = deckRef.current.pickObjects({ x: x0, y: y0, width: w, height: h, layerIds: ["cells"] });
      picked = objs.map((o: any) => (ids ? ids[o.index] : o.index));
      setSelIdx(objs.length ? new Set<number>(objs.map((o: any) => o.index)) : null);  // highlight on the map
    } else {
      setSelIdx(null);
    }
    onSelect?.(picked);
    setBox(null);
  };
  const clearSelection = () => { setSelIdx(null); onSelect?.([]); };

  const btnStyle: React.CSSProperties = {
    position: "absolute", top: 16, right: 16, pointerEvents: "auto", cursor: "pointer",
    fontFamily: "'JetBrains Mono', monospace", fontSize: 11, letterSpacing: ".08em", textTransform: "uppercase",
    padding: ".34rem .6rem", borderRadius: 9, zIndex: 3,
    color: selectMode ? "#0B0C10" : "#A896F2", background: selectMode ? "#A896F2" : "rgba(168,150,242,.12)",
    border: "1px solid " + (selectMode ? "transparent" : "rgba(168,150,242,.34)"),
  };

  return (
    <div ref={wrapRef} style={{ position: "relative", width: "100%", height: "100%", background: "#0A0C11", borderRadius: 16, overflow: "hidden" }}>
      <canvas ref={canvasRef} style={{ width: "100%", height: "100%" }} />
      {/* Drag-capture overlay: active only in select mode, so deck.gl keeps pan/zoom otherwise.
          Sitting above the canvas, it reliably receives the whole drag (deck never sees it). */}
      <div onPointerDown={onDown} onPointerMove={onMove} onPointerUp={onUp}
           style={{ position: "absolute", inset: 0, zIndex: 2, cursor: selectMode ? "crosshair" : "default",
                    pointerEvents: selectMode ? "auto" : "none" }} />
      {(hud.v || hud.k) && (
        <div style={{ ...glass, top: 18, left: 18 }}>
          {hud.k && <div style={{ fontSize: 8, letterSpacing: ".16em", textTransform: "uppercase", color: "#9096A4" }}>{hud.k}</div>}
          {hud.v && <div style={{ fontFamily: "'Space Grotesk', sans-serif", fontWeight: 700, fontSize: 24, color: "#F3F4F8", marginTop: 2 }}>{hud.v}</div>}
          {hud.sub && <div style={{ fontSize: 10, color: "#6A7080", marginTop: 2 }}>{hud.sub}</div>}
        </div>
      )}
      {stat && <div style={{ ...glass, left: 18, bottom: 16, padding: ".45rem .7rem", color: "#A896F2", fontWeight: 600 }}>{stat}</div>}
      {selIdx && selIdx.size > 0 && (
        <div style={{ ...glass, left: 18, bottom: stat ? 56 : 16, padding: ".4rem .65rem", color: "#F3F4F8",
                      pointerEvents: "auto", cursor: "pointer", display: "flex", alignItems: "center", gap: 8 }}
             onPointerDown={(e) => { e.stopPropagation(); clearSelection(); }}>
          <span style={{ color: "#A896F2", fontWeight: 700 }}>{selIdx.size.toLocaleString()} selected</span>
          <span style={{ color: "#9096A4" }}>· clear ✕</span>
        </div>
      )}
      {/* The always-on "scroll to zoom · drag to pan" hint overlapped the featureplot colour scale
          (both bottom-right) and stated the obvious - scroll/drag are standard map gestures. Keep only
          the contextual "drag to select" cue, which is non-obvious and shown just while selecting. */}
      {selectMode && <div style={{ ...glass, left: 16, bottom: 16, color: "#9096A4", fontSize: 11 }}>drag to select a region</div>}
      {legend.length > 0 && (
        <div style={{ ...glass, top: 54, right: 16, maxHeight: "70%", overflow: "auto", pointerEvents: "auto" }}
             onPointerLeave={() => setHi(null)}>
          {legend.slice(0, 35).map((l, i) => {
            const active = !!hi && hi[0] === l.color[0] && hi[1] === l.color[1] && hi[2] === l.color[2];
            const dim = !!hi && !active;
            return (
              <div key={i} onPointerEnter={() => setHi(l.color)}
                   style={{ display: "flex", alignItems: "center", gap: 7, margin: "3px 0", cursor: "pointer",
                            opacity: dim ? 0.4 : 1, fontWeight: active ? 700 : 400, color: active ? "#F3F4F8" : undefined, transition: "opacity .12s" }}>
                {/* flex:"0 0 auto" so a two-char label ("10".."29") can't shrink the swatch to a
                    sliver (default flex-shrink:1 did, hence "tiny bars, only colours 0-9"). */}
                <span style={{ flex: "0 0 auto", width: 9, height: 9, borderRadius: 99, background: `rgb(${l.color[0]},${l.color[1]},${l.color[2]})`, boxShadow: `0 0 8px rgb(${l.color[0]},${l.color[1]},${l.color[2]})` }} />
                <span style={{ whiteSpace: "nowrap" }}>{l.label}</span>
              </div>
            );
          })}
        </div>
      )}
      <div style={btnStyle} onPointerDown={(e) => { e.stopPropagation(); setSelectMode(!selectMode); }}>{selectMode ? "selecting" : "select region"}</div>
      {box && selectMode && (
        <div style={{ position: "absolute", left: Math.min(box.x0, box.x1), top: Math.min(box.y0, box.y1), width: Math.abs(box.x1 - box.x0), height: Math.abs(box.y1 - box.y0), border: "1.5px dashed #A896F2", background: "rgba(168,150,242,.12)", borderRadius: 6, pointerEvents: "none", zIndex: 3 }} />
      )}
    </div>
  );
}
