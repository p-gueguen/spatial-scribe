import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// Dev: vite serves the SPA and proxies /api to the FastAPI backend (a compute node:8000), so the browser
// talks to a single origin (no CORS in the browser). For production, build and serve dist/ from
// FastAPI (or any static host) behind the same origin as /api.
export default defineConfig({
  plugins: [react()],
  // Base path the built SPA is served under. "/" by default; set VITE_BASE=/sub/path/ when the
  // app is served under a sub-path behind a reverse proxy, so assets + the /api prefix resolve.
  base: process.env.VITE_BASE || "/",
  server: {
    host: "0.0.0.0",
    // Allow any host so the dev server is reachable by hostname / FQDN / IP (not just localhost),
    // which vite's default host allowlist would otherwise reject with "Blocked request".
    allowedHosts: true,
    port: 5173,
    proxy: {
      "/api": { target: process.env.VITE_API_TARGET || "http://localhost:8000", changeOrigin: true },
    },
  },
});
