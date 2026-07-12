# SpatialScribe - Design System bundle

Self-contained component preview cards for the SpatialScribe "control-desk" theme
(blue-black canvas, single lavender-violet accent, Space Grotesk / Manrope / JetBrains Mono).
Each `*.html` is a standalone `@dsCard` preview; `_fonts.css` carries the embedded fonts.

## Publish to claude.ai/design

From an interactive Claude Code terminal (desktop), in this repo:

    /design-login        # one-time: authorize design-system access
    /design-sync         # pick/create a project, review the plan, push these cards

The sync is incremental (one component at a time) and shows you the exact file plan
before writing. Cards are grouped as Tokens / Components / Overview.

## Source of truth

These mirror the live app theme in `webapp/src/styles.css` (React SPA + deck.gl HUD)
(+ `.streamlit/config.stakent.toml`). Refine there and regenerate, or edit cards directly.
