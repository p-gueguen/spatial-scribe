// Shared presentational helpers for the analysis panels: inline figures, a minimal markdown
// renderer for LLM verdicts, and a compact heatmap (diverging for z-scores, violet ramp for scores).
import React from "react";
import * as api from "../api";

// A small "?" affordance that reveals its explanation on hover (native title), so panels can carry
// their guidance without a wall of inline text. Place it next to a heading or control.
// A "?" affordance that reveals its explanation on hover. Uses a cursor-anchored position:fixed
// popover (NOT a native `title`, which is unreliable, and NOT an absolute child, which the rail's
// overflow-y:auto would clip) so the text is always visible. pointer-events:none so it never eats
// the hover. Place it next to a heading or control.
export function Help({ text }: { text: string }) {
  const [tip, setTip] = React.useState<{ x: number; y: number } | null>(null);
  const at = (e: React.MouseEvent) =>
    setTip({ x: Math.min(e.clientX + 12, window.innerWidth - 292), y: Math.min(e.clientY + 16, window.innerHeight - 150) });
  return (
    <span className="help" role="img" aria-label={text} tabIndex={0}
          onMouseEnter={at} onMouseMove={at} onMouseLeave={() => setTip(null)}>
      ?
      {tip && <span className="helptip" style={{ left: tip.x, top: tip.y }}>{text}</span>}
    </span>
  );
}

// Every figure is click-to-zoom: click opens a full-viewport modal with the same PNG at natural
// size. Escape, the close control, or a backdrop click closes it. Implemented here so every panel
// that renders a Fig gets the behaviour for free. The .figmodal / .figimg CSS lives in styles.css.
export function Fig({ a }: { a: api.Artifact }) {
  const [zoom, setZoom] = React.useState(false);
  React.useEffect(() => {
    if (!zoom) return;
    const onKey = (e: KeyboardEvent) => { if (e.key === "Escape") setZoom(false); };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [zoom]);
  if (a.kind !== "figure" || !a.png) return null;
  return (
    <div>
      {a.title && <div className="pmuted" style={{ fontSize: ".68rem", marginTop: 4 }}>{a.title}</div>}
      <img className="figimg" src={a.png} alt={a.title || "figure"} onClick={() => setZoom(true)} />
      {zoom && (
        <div className="figmodal" onClick={() => setZoom(false)}>
          <button className="close" onClick={() => setZoom(false)}>close</button>
          <img src={a.png} alt={a.title || "figure"} />
        </div>
      )}
    </div>
  );
}

function inline(s: string): React.ReactNode {
  return s.split(/(\*\*[^*]+\*\*)/g).map((p, i) =>
    p.startsWith("**") && p.endsWith("**")
      ? <strong key={i}>{p.slice(2, -2)}</strong>
      : <React.Fragment key={i}>{p}</React.Fragment>);
}

export function Markdown({ text }: { text: string }) {
  const lines = (text || "").split("\n");
  const out: React.ReactNode[] = [];
  let list: string[] = [];
  const flush = (k: string) => { if (list.length) { out.push(<ul key={k}>{list.map((l, j) => <li key={j}>{inline(l)}</li>)}</ul>); list = []; } };
  lines.forEach((ln, i) => {
    const t = ln.trim();
    if (/^#{1,6}\s/.test(t)) { flush("u" + i); const lvl = (t.match(/^#+/) as RegExpMatchArray)[0].length; const txt = t.replace(/^#+\s/, ""); out.push(lvl <= 2 ? <h2 key={i}>{inline(txt)}</h2> : <h3 key={i}>{inline(txt)}</h3>); }
    else if (/^[-*]\s/.test(t)) { list.push(t.replace(/^[-*]\s/, "")); }
    else { flush("u" + i); if (t) out.push(<p key={i} style={{ margin: "3px 0" }}>{inline(t)}</p>); }
  });
  flush("uend");
  return <div className="verdict">{out}</div>;
}

// onHoverCol fires with a COLUMN (a state/program) when the cursor enters that column's header OR
// any cell in it - the user's mental model is "hover the program", so both surfaces recolour the
// map. onHover (row,col) is unchanged for the enrichment-pair highlight; onLeave clears both.
// onHover(row,col) fires when the cursor enters a CELL; onHoverCol(col) fires only from the COLUMN
// HEADER - two distinct gestures so a states-heatmap cell can highlight its cell TYPE while hovering
// the state title recolours by that state's score. `activeCol` gets a discreet violet outline so the
// user can see which column is currently painted on the map.
export function Heat({ rows, cols, matrix, diverging = false, activeCol, onHover, onHoverCol, onLeave }:
  { rows: string[]; cols: string[]; matrix: number[][]; diverging?: boolean; activeCol?: string | null;
    onHover?: (row: string, col: string) => void; onHoverCol?: (col: string) => void; onLeave?: () => void }) {
  if (!rows?.length || !cols?.length) return <div className="pmuted">no data to plot</div>;
  const flat = matrix.flat().filter((v) => Number.isFinite(v)) as number[];
  // For a square self-pair matrix (neighborhood enrichment), the self-self DIAGONAL z-scores are
  // vastly larger than any cross-type value (e.g. Ependymal-Ependymal ~842 vs cross pairs |z|<30).
  // Including them in amax crushed every off-diagonal cell into a near-constant faint floor
  // (|z|=30 -> alpha 0.15, |z|=0 -> 0.12), so enriched (violet) and avoided (red) cross-type pairs
  // were indistinguishable dark boxes. Normalize the diverging scale on the OFF-DIAGONAL cells so
  // cross-type interactions get the full colour range; the diagonal clamps to full intensity
  // (self-self is trivially maximally colocalized).
  const offdiag = rows.length === cols.length
    ? (matrix.flatMap((row, ri) => row.filter((v, ci) => ri !== ci && Number.isFinite(v))) as number[])
    : flat;
  const amax = Math.max(0.001, ...(offdiag.length ? offdiag : flat).map((v) => Math.abs(v)));
  const lo = Math.min(...flat), hi = Math.max(...flat);
  const color = (v: number) => {
    if (!Number.isFinite(v)) return "#14161c";
    if (diverging) { const t = Math.max(-1, Math.min(1, v / amax)); return t >= 0 ? `rgba(247,116,110,${0.12 + 0.88 * t})` : `rgba(108,92,224,${0.12 + 0.88 * -t})`; }
    const t = (v - lo) / ((hi - lo) || 1); return `rgba(168,150,242,${0.1 + 0.9 * t})`;
  };
  return (
    <div style={{ overflowX: "auto" }} onPointerLeave={() => onLeave?.()}>
      <div className="heat" style={{ gridTemplateColumns: `78px repeat(${cols.length}, minmax(15px,1fr))` }}>
        <div />
        {cols.map((c, i) => <div key={i} className="heatlbl col" title={c} onPointerEnter={() => onHoverCol?.(c)}
                                  style={c === activeCol ? { color: "#C9BCFF", fontWeight: 600 } : undefined}>{c}</div>)}
        {rows.map((r, ri) => (
          <React.Fragment key={ri}>
            <div className="heatlbl" title={r} style={{ alignSelf: "center" }}>{r}</div>
            {cols.map((_, ci) => {
              const v = matrix[ri]?.[ci];
              const active = cols[ci] === activeCol;
              return <div key={ci} className="heatcell"
                          title={`${r} · ${cols[ci]}: ${Number.isFinite(v) ? (v as number).toFixed(2) : "-"}`}
                          onPointerEnter={() => onHover?.(r, cols[ci])}
                          style={{ background: color(v as number),
                                   boxShadow: active ? "inset 0 0 0 1.5px rgba(168,150,242,.85)" : undefined }} />;
            })}
          </React.Fragment>
        ))}
      </div>
    </div>
  );
}
