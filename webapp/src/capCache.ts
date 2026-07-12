// Session-scoped memo for expensive capability results, so leaving a tab and coming back does not
// recompute. Keyed by session id, so loading a new section cannot serve another section's numbers.
// Deliberately a plain module Map: results are small JSON, the page owns one session at a time, and
// nothing here needs to outlive a reload.
const _cache = new Map<string, any>();

// "::" separates the sid/name/extra fields. sids are opaque ids and the capability names passed here
// are hardcoded literals with no "::", so a session prefix is unambiguous for clearSession; undefined
// sid folds to "-".
const _sep = "::";
const _prefix = (sid: string | undefined): string => (sid ?? "-") + _sep;

// Build the cache key for a capability result. `extra` scopes a variant (a type_a|type_b pair, a
// cluster resolution, ...) so different parameterisations of the same capability do not overwrite.
export function capKey(sid: string | undefined, name: string, extra?: string): string {
  return _prefix(sid) + name + (extra != null ? _sep + extra : "");
}

// Read a cached result, or undefined if nothing is stored under `key`.
export function getCached<T = any>(key: string): T | undefined {
  return _cache.get(key) as T | undefined;
}

// Store a result under `key`. Values are the same small JSON the panels already hold in React state.
export function setCached(key: string, value: any): void {
  _cache.set(key, value);
}

// Drop every entry for one session. Call after any op that mutates that section's cells (rename,
// subcluster, region filter, re-cluster, annotate) so stale numbers can never be served.
export function clearSession(sid: string | undefined): void {
  const p = _prefix(sid);
  for (const k of _cache.keys()) if (k.startsWith(p)) _cache.delete(k);
}

// Drop everything. Call when a brand-new section is loaded.
export function clearAll(): void {
  _cache.clear();
}
