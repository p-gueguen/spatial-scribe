// Client-side history of successfully-loaded sections (server-side path + tissue), newest-first, in
// localStorage. Pure convenience - no server state; capped and deduped by path. Only successful loads
// are recorded (App.loadSection calls pushRecent after the section resolves).
export type RecentSection = { path: string; tissue: string; at: number };

const KEY = "spatialscribe.recentSections";
const CAP = 8;

export function getRecent(): RecentSection[] {
  try {
    const v = JSON.parse(localStorage.getItem(KEY) || "[]");
    return Array.isArray(v) ? v.filter((r) => r && typeof r.path === "string") : [];
  } catch {
    return [];
  }
}

export function pushRecent(path: string, tissue: string): void {
  const p = (path || "").trim();
  if (!p) return;
  const rest = getRecent().filter((r) => r.path !== p);          // dedupe by path, keep newest
  const next = [{ path: p, tissue: (tissue || "").trim(), at: Date.now() }, ...rest].slice(0, CAP);
  try { localStorage.setItem(KEY, JSON.stringify(next)); } catch { /* quota / private mode - ignore */ }
}

export function clearRecent(): void {
  try { localStorage.removeItem(KEY); } catch { /* ignore */ }
}

// last two path segments, e.g. "…/internal/xenium_output" - full path goes in the row's title tooltip.
export function shortPath(p: string): string {
  const parts = p.replace(/\/+$/, "").split("/").filter(Boolean);
  return parts.length <= 2 ? p : "…/" + parts.slice(-2).join("/");
}
