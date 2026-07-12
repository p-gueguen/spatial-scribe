# SpatialScribe - brand & design system

The one source of truth for how SpatialScribe looks. It applies to **everything the project
ships**: the React SPA (`webapp/`), the deck.gl specimen viewport, the plotly / matplotlib
figures, the self-contained HTML reports, and the `PROGRESS.html` dashboard. If you build a new
surface, pull from here so the whole project reads as one product.

> **The live theme is the "control-desk" look** (blue-black canvas, single lavender-violet accent,
> Space Grotesk / Manrope / JetBrains Mono). The earlier warm "field-notebook" theme (ivory paper,
> clay accent, Hanken Grotesk + Fira Code) was **retired with Streamlit on 2026-07-09** - the React
> SPA is the only front-end now, and its tokens live in `webapp/src/styles.css`. If any file still
> points at the warm theme (`src/spatialscribe/app/theme.py`, a `.streamlit/config*.toml`,
> Hanken/Fira/clay `#D97757`), it is stale.

## Design requirements

Hard constraints for every SpatialScribe surface. These are design decisions, not casual edits -
change them deliberately.

- **Palette: blue-black control desk + a single lavender-violet accent** (`--violet #A896F2`, deep
  `--violet-deep #6C5CE0`). No warm neutrals, no clay, no second accent channel. Emphasis that is
  not the violet is carried by **weight** and by the near-white ink, not by another colour.
- **Dark by default, and that is justified.** Imaging spatial data is viewed dark-field, so the
  whole instrument is a calm dark console and the data glows on it. The **specimen viewport** is the
  darkest inset (`#0A0C11`) within the already-dark app.
- **Three font families, and no more.** Space Grotesk (display), Manrope (all UI text), JetBrains
  Mono (data / numbers / gene names). Hierarchy comes from weight and size, never a fourth face.
- **One accent, used sparingly.** Violet marks only the active step, the primary action, key
  numbers, the kicker, links, and the brand wordmark. The **gradient wordmark is the one deliberate
  gradient-on-text**; nowhere else uses gradient text.
- **No emojis in product surfaces.** State is shown with type weight, colour, and layout (numbered
  steps, CSS status dots, a violet-tinted active row, a pulsing "running" pill), never emoji glyphs.
  (Emojis remain fine in internal docs / commit messages.)
- **Accessible by default.** Near-white ink on the dark ground clears 4.5:1; visible violet keyboard
  focus rings; `prefers-reduced-motion` respected (it stops animation and keeps the colour);
  transitions 140-500ms.
- **Distinctive, not templated.** The identity is the blue-black console, the lavender glow, the
  deck.gl **specimen viewport**, and the JetBrains Mono data language - not a generic SaaS dark mode.
  The violet is lavender (cool, soft), not a neon/electric accent.
- **Consistent across surfaces.** The React app, the plot figures, the HTML report, and
  `PROGRESS.html` all pull the tokens and fonts below, so the project reads as one product.

## Principle

**Control desk + specimen viewport.** SpatialScribe is a working instrument for imaging spatial
transcriptomics, so the interface reads like a calm, dark control desk - a blue-black console you
operate - built around a single **specimen viewport** where the spatial data glows (dark-field is
genuinely how you view fluorescent cells). The chrome is quiet and near-monochrome; the data and the
one lavender accent are the only bright things.

- **One restrained accent.** A single lavender **violet** (`#A896F2`, deep `#6C5CE0`) marks only what
  matters: the active step, the primary action, key numbers, the kicker, links, the brand. Everything
  else is a cool blue-black neutral.
- **Let the data breathe.** Minimal padding and chrome; the deck.gl specimen viewport is the dominant
  element (near-full-height, dark, GPU-rendered). The step panels are a narrow left column; the
  copilot is a compact bottom drawer.
- **Glow is the state language.** A soft violet radial sits behind the app; the active step dot, the
  primary button, and the "running" pill glow; traffic-light dots (mint/amber/coral) glow in tables.
  Motion (a pulsing dot) is what says "still working", so a computing section never reads as static.
- **Weight carries hierarchy.** Space Grotesk 700 titles + bold key numbers + the violet accent do
  the emphasis, so the near-monochrome page never reads flat.

## Color tokens

Blue-black control-desk palette (dark). Ink is near-white for 4.5:1; the violet stays legible on the
dark ground without deepening (unlike the retired clay, which needed a text-dark variant).

| Token | Hex / value | Use |
|-------|-------------|-----|
| `--bg` | `#0B0C10` | app / page background (blue-black) |
| `--rail` | `#171A22` | left rail (the step navigator) |
| `--surface` | `#15171E` | elevated surfaces (panels, cards) |
| `--surface2` | `#1C1F28` | hover, secondary buttons, inputs on a surface |
| `--sink` | `#111319` | recessed wells (inputs, code blocks, meter tracks) |
| `--line` / `--line2` | `#23262F` / `#2E323D` | hairline borders / stronger dividers |
| `--rail-line` | `#2E323D` | rail border + the stepper's connector line |
| `--ink` | `#F3F4F8` | primary text (near-white) |
| `--body` | `#C4C8D2` | secondary text, prose |
| `--muted` / `--faint` | `#9096A4` / `#6A7080` | labels / meta text |
| `--violet` | `#A896F2` | **the accent**: active step, kicker, key numbers, links, brand |
| `--violet-2` | `#8B76EC` | mid violet (wordmark gradient stop) |
| `--violet-deep` | `#6C5CE0` | primary-button / active-dot base, focus ring |
| `--violet-tint` | `rgba(168,150,242,.12)` | active-step fill, kicker / callout background |
| `--violet-line` | `rgba(168,150,242,.34)` | violet hairline (kicker border, focus-within) |
| `--glow` | `rgba(124,92,224,.55)` | violet glow on the active dot / primary button |
| `--pass` | `#46E39B` | PASS / confident / "resolvable" (mint) |
| `--warn` | `#F2B24C` | WARN / tentative / "weak" (amber) |
| `--fail` | `#F7746E` | FAIL / abstained / "cannot resolve" (coral) |
| *(viewport)* | `#0A0C11` | the deck.gl specimen canvas (the darkest inset) |

The app carries a barely-there **violet radial glow** behind everything (`.app::before`: a large
`rgba(124,92,224,.16)` glow at top-right + a fainter blue glow at top-left), so the dark ground has
depth and is never a flat black.

## Typography

**Three families**, loaded from Google Fonts in `webapp/index.html`. Hierarchy comes from weight and
size; one display face for headings, one sans for UI, one mono for data.

| Role | Family | Used for |
|------|--------|----------|
| Display | **Space Grotesk** (500/600/700) | brand wordmark (gradient), section titles, panel headings, metric values, run labels |
| Body / UI | **Manrope** (400/500/600/700) | prose, button labels, captions, sub-text |
| Data / mono | **JetBrains Mono** (400/500/600) | kickers, metric labels, gene names, tables, code, the rail stepper, the viewport HUD, all numeric readouts |

Scale: title ~1.7rem Space Grotesk 700 (tight `-.03em`); kicker 0.62rem JetBrains Mono (letter-spacing
.2em, uppercase, violet); metric value 1.5rem Space Grotesk 700; body 0.82-0.86rem Manrope. Use
`<b>` + `--ink` (or the violet) to make a key number or term stand out inline.

## Layout

- **`.app` = CSS grid `288px 1fr`**, full-viewport, `overflow:hidden`, with the ambient violet radial
  glow layer behind it.
- **Left rail** (288px, `--rail`): the gradient wordmark + a mono sub-label + a live status dot + a
  thin progress bar, then the numbered **stepper** (`01 Load` ... `07 Report`) as connected dots.
- **Main column**: a **topbar** (violet kicker pill + Space Grotesk title + muted sub + controls -
  demo picker, GPU/CPU segment, primary run button) -> the **stage** -> the **copilot drawer**.
- **Stage** = a narrow **panel column** (~380px, scrolls independently, `webapp/src/panels/*`) beside
  the **specimen viewport** (the deck.gl WebGL canvas hero, fills the rest).
- **Copilot** = a compact **bottom drawer** spanning the main column (avatar, a streaming answer line,
  a mono input with a mic/dictate toggle + send, example chips, and a live LLM model tag). *(This
  replaced the old right-rail copilot; the canvas now owns the full width above the drawer.)*

## Components

- **Kicker** - a violet mono pill with a glowing dot that encodes the wizard step (`03 · QUALITY`),
  sitting above a Space Grotesk title. (`.kicker` + `.title`)
- **Rail stepper** - connected numbered dots; the active row is violet-tinted with a gradient dot and
  a glow (`.step.active`), a completed step shows a mint check (`.step.done`). (`.steps` / `.step`)
- **Readout metric** - a `linear-gradient(180deg, --surface, --sink)` tile with an uppercase mono
  label and a large Space Grotesk value. Never truncate the value. (`.metric`)
- **Signal meter** - a segmented pass/warn/fail bar (mint / amber / coral) with a mono legend, for the
  annotatability headline (confident / tentative / abstained). (`.meter` + `.meterlbl`)
- **Reason pills** - small mono pills breaking down abstention reasons. (`.pills` / `.pill2`)
- **Traffic-light dots** - glowing `.tl.green` / `.tl.amber` / `.tl.red` dots for panel-check
  resolvable / weak / cannot-resolve in the per-type tables.
- **Specimen viewport** - `webapp/src/SpecimenCanvas.tsx`: a **deck.gl WebGL scatter** (binary
  attributes -> scales past 1e6 cells) on `#0A0C11`, with a floating **glass HUD**
  (`rgba(13,15,20,.72)` + `backdrop-filter:blur(10px)`, `--line2` border), native zoom/pan, a box
  **"select region"** tool, a categorical legend with glowing swatches, and a violet featureplot
  colour scale (`#23203a -> #6c5ce0 -> #a896f2 -> #ebe6ff`). Continuous overlays: the violet ramp for
  QC / depth / confidence, **magma** for feature scores (signature / program / malignancy);
  categoricals use the fluorescence palette below.
- **Primary button** - violet gradient `linear-gradient(135deg, #9A86F0, #6C5CE0)` with a glow;
  secondary = `--surface` with a violet-line hover. (`.btn.primary` / `.pbtn.pri`)
- **Running pill** - a violet pill with a pulsing dot for a section that is computing; the pulse is
  what says "still working" (reduced-motion stops the pulse, keeps the colour). (`.running`)

### Categorical (cell-type / cluster) palette

Bright fluorescence hues chosen to glow on the dark viewport and stay perceptually distinct - the
single source is **`plots._PALETTE`** (41 entries, gated by `tests/test_palette.py` on min pairwise
CIE76 >= 25 and min luminance), and it **cycles by design** past 41 categories:

`#22d3ee #e879f9 #34d399 #fbbf24 #f87171 #a78bfa #4ade80 #f472b6 #38bdf8 #facc15 #fb923c #2dd4bf
#c084fc #60a5fa #f9a8d4 #86efac #fde047 #fca5a5 #5eead4 #d8b4fe #93c5fd ...`

This palette lives inside the dark specimen viewport and the report figures; the surrounding UI chrome
stays blue-black + violet.

## Brand mark

`webapp/public/favicon.svg`: a dark rounded square (`#0B0C10`) with a lavender (`#A896F2`) **hexagon
outline + a centre dot** - a single cell in a field. The wordmark is "SpatialScribe" in Space Grotesk
700 with a lavender gradient (`linear-gradient(102deg, #F1EFFC, #B9A8F6 48%, #8B76EC)` clipped to
text).

## Where it lives

- **App theme + all component classes:** `webapp/src/styles.css` (`:root` tokens + the shared
  `.panel` / `.metric` / `.meter` / `.dtable` / ... vocab). **The single source of truth** - refine
  here.
- **Specimen viewport + glass HUD:** `webapp/src/SpecimenCanvas.tsx` (deck.gl scatter, HUD, select).
- **Fonts:** `webapp/index.html` (the Google Fonts `<link>`). Base64-embedded copies live in
  `design-system/_fonts.css` for the standalone design-system cards / offline self-hosting;
  regenerate with `python scripts/embed_fonts.py`.
- **Figures:** `src/spatialscribe/analysis/plots.py` (`_PALETTE`, the dark plot theme).
- **HTML report:** `backend/report.py` reuses the same blue-black + violet tokens (`#0b0c10` /
  `#070a10` grounds, `#a896f2` accent, JetBrains Mono), so a delivered report matches the app.
- **Design-system preview cards:** `design-system/*.html` - standalone `@dsCard` previews (Tokens /
  Components / Overview), publishable to `claude.ai/design` via `/design-sync`. They mirror
  `webapp/src/styles.css`; keep them in sync when the tokens change.
- **`PROGRESS.html`:** generated (`the progress renderer from `progress.yaml`); its CSS should
  track these tokens.

## Gotchas

- **Fonts load from the Google Fonts CDN** (the `<link>` in `webapp/index.html`), not self-hosted.
  Offline, behind a proxy, or under a strict CSP they fall back to system fonts and the UI looks
  generic. The base64-embedded copies in `design-system/_fonts.css` are the offline-safe source if you
  need to inline them (a standalone artifact must inline them - a webfont `<link>` is blocked by the
  Artifact CSP). This is also the "judge clones the repo offline" caveat.
- **The wordmark gradient is the ONE deliberate gradient-on-text.** Everywhere else, emphasis is
  weight + the flat violet - never gradient text on body or headings.
- **The featureplot colour ramp is deliberately two-track:** violet for a QC / depth / confidence
  field ("more of the same measurement"), magma for a feature score (signature / program /
  malignancy), so a sparse positive tail reads against the dark ground. Add a new score column ->
  register its prefix in `backend.app._is_feature_field` or it renders as a QC metric.
- **Stale-theme references bite.** The retired warm field-notebook theme (Hanken Grotesk + Fira Code +
  clay `#D97757`, `src/spatialscribe/app/theme.py`, `.streamlit/config*.toml`) is gone. The
  `design-system/README.md` still mentions `.streamlit/config.stakent.toml`; that line is stale - the
  React `styles.css` is the only live source.
