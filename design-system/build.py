"""Generate the SpatialScribe design-system bundle (self-contained @dsCard preview HTML).

Emits design-system/*.html cards + _fonts.css + README, ready for claude.ai/design via
/design-sync. Every card is self-contained: tokens + component CSS inlined, fonts @import'd
from the sibling _fonts.css (the bundled base64 woff2 in design-system/fonts_stakent.css).
"""
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "design-system"
OUT.mkdir(exist_ok=True)
shutil.copyfile(OUT / "fonts_stakent.css", OUT / "_fonts.css")

# Shared tokens + component CSS (plain classes mirroring the Streamlit theme's look).
CSS = """@import url('_fonts.css');
:root{
  --bg:#0B0C10; --rail:#171A22; --surface:#15171E; --surface2:#1C1F28; --sink:#111319;
  --rail-line:#2E323D; --line:#23262F; --line2:#2E323D;
  --ink:#F3F4F8; --body:#C4C8D2; --muted:#9096A4; --faint:#6A7080;
  --violet:#A896F2; --violet-2:#8B76EC; --violet-deep:#6C5CE0;
  --violet-tint:rgba(168,150,242,.12); --violet-line:rgba(168,150,242,.32); --glow:rgba(124,92,224,.55);
  --pass:#46E39B; --warn:#F2B24C; --fail:#F7746E;
  --sans:'Manrope',ui-sans-serif,system-ui,-apple-system,sans-serif;
  --display:'Space Grotesk',ui-sans-serif,system-ui,sans-serif;
  --mono:'JetBrains Mono',ui-monospace,'SF Mono',monospace;
  --shadow:0 1px 0 rgba(255,255,255,.02) inset,0 18px 40px -26px rgba(0,0,0,.75);
}
*{box-sizing:border-box}
body{margin:0;min-height:100vh;background:var(--bg);color:var(--body);font-family:var(--sans);
  -webkit-font-smoothing:antialiased;padding:34px;position:relative;overflow-x:hidden}
body::before{content:"";position:fixed;inset:0;pointer-events:none;z-index:0;
  background:radial-gradient(900px 520px at 92% -8%,rgba(124,92,224,.16),transparent 60%),
            radial-gradient(760px 480px at 8% -12%,rgba(70,120,220,.07),transparent 60%)}
.ds{position:relative;z-index:1;max-width:1160px;margin:0 auto}
.eyebrow{font-family:var(--mono);font-size:.62rem;letter-spacing:.24em;text-transform:uppercase;color:var(--faint);margin:0 0 1.3rem}
.grid{display:grid;gap:16px}
h1,h2,h3{font-family:var(--display);color:var(--ink);letter-spacing:-.02em;margin:0}

/* buttons */
.btn{font-family:var(--sans);font-weight:600;font-size:.92rem;border-radius:11px;border:1px solid var(--line2);
  background:var(--surface);color:var(--ink);min-height:2.7rem;padding:.55rem 1.25rem;cursor:pointer;
  display:inline-flex;align-items:center;gap:.5rem}
.btn.primary{background:linear-gradient(135deg,#9A86F0,var(--violet-deep));border-color:transparent;color:#fff;font-weight:700;
  box-shadow:0 10px 26px -10px var(--glow),inset 0 1px 0 rgba(255,255,255,.18)}
.btn.hover{background:linear-gradient(135deg,#A895F5,#7A69E6);box-shadow:0 14px 32px -10px var(--glow),inset 0 1px 0 rgba(255,255,255,.24)}
.btn.ghost{background:transparent}
.btn:disabled,.btn.disabled{opacity:.45;cursor:not-allowed}

/* metric tile */
.tile{background:linear-gradient(180deg,var(--surface),var(--sink));border:1px solid var(--line);border-radius:14px;
  padding:16px 18px 14px;box-shadow:var(--shadow);min-width:190px}
.tile .k{font-family:var(--mono);text-transform:uppercase;letter-spacing:.12em;font-size:.58rem;color:var(--muted)}
.tile .v{font-family:var(--display);font-weight:700;font-size:1.85rem;letter-spacing:-.03em;color:var(--ink);
  margin-top:.35rem;text-shadow:0 0 22px rgba(168,150,242,.22)}
.tile .d{font-family:var(--mono);font-size:.72rem;margin-top:.3rem}
.up{color:var(--pass)} .down{color:var(--fail)}

/* nav rail */
.rail{background:var(--rail);border:1px solid var(--rail-line);border-radius:16px;padding:18px 14px;width:250px}
.brand{font-family:var(--display);font-weight:700;font-size:1.4rem;letter-spacing:-.02em;line-height:1;width:max-content;
  background:linear-gradient(102deg,#F1EFFC 0%,#B9A8F6 48%,#8B76EC 100%);-webkit-background-clip:text;background-clip:text;
  -webkit-text-fill-color:transparent;color:transparent}
.brandsub{font-family:var(--mono);font-size:.56rem;letter-spacing:.22em;text-transform:uppercase;color:var(--muted);margin:.34rem 0 0}
.status{font-family:var(--mono);font-size:.66rem;color:var(--faint);margin-top:.6rem;display:flex;align-items:center;gap:.4rem}
.status::before{content:"";width:6px;height:6px;border-radius:99px;background:var(--pass);box-shadow:0 0 8px var(--pass)}
.rail hr{border:none;border-top:1px solid var(--rail-line);margin:.9rem 0}
.railbar{margin:.7rem 0 .1rem}
.railbar .top{display:flex;justify-content:space-between;font-family:var(--mono);font-size:.55rem;letter-spacing:.16em;text-transform:uppercase;color:var(--faint);margin-bottom:6px}
.railbar .track{height:5px;border-radius:99px;background:var(--sink);overflow:hidden;border:1px solid var(--line)}
.railbar .fill{height:100%;background:linear-gradient(90deg,var(--violet-deep),var(--violet));box-shadow:0 0 10px -1px var(--glow)}
.nav{font-family:var(--mono);font-size:.8rem;color:var(--body);padding:.5rem .7rem;border-radius:11px;border:1px solid transparent;margin:3px 0}
.nav.active{background:var(--violet-tint);color:var(--ink);border-color:var(--violet-line);font-weight:700;
  box-shadow:-3px 0 0 var(--violet) inset,0 8px 22px -14px var(--glow)}

/* section header */
.kicker{display:inline-flex;align-items:center;gap:.42rem;font-family:var(--mono);font-size:.62rem;letter-spacing:.2em;
  text-transform:uppercase;color:var(--violet);background:var(--violet-tint);border:1px solid var(--violet-line);
  border-radius:99px;padding:.24rem .6rem}
.kicker::before{content:"";width:6px;height:6px;border-radius:99px;background:var(--violet);box-shadow:0 0 8px var(--violet)}
.title{font-family:var(--display);font-weight:700;font-size:1.95rem;letter-spacing:-.03em;color:var(--ink);margin:.6rem 0 .35rem;line-height:1.05}
.sub{color:var(--muted);font-size:.9rem;max-width:70ch;line-height:1.55}

/* signal meter */
.meter .bar{display:flex;height:12px;border-radius:99px;overflow:hidden;border:1px solid var(--line2);background:var(--surface2)}
.meter .seg{height:100%;box-shadow:0 0 12px -2px currentColor}
.legend{display:flex;gap:1.3rem;margin-top:.6rem;font-family:var(--mono);font-size:.74rem;color:var(--muted);flex-wrap:wrap}
.legend b{color:var(--ink);font-weight:600}
.dot{display:inline-block;width:9px;height:9px;border-radius:99px;margin-right:.45rem;vertical-align:middle;background:currentColor;box-shadow:0 0 8px currentColor}

/* pills */
.pills{display:flex;gap:.45rem;flex-wrap:wrap;align-items:center}
.pills-lbl{font-family:var(--mono);font-size:.6rem;color:var(--faint);letter-spacing:.16em;text-transform:uppercase;margin-right:.15rem}
.pill{font-family:var(--mono);font-size:.72rem;color:var(--pc);background:color-mix(in srgb,var(--pc) 14%,transparent);
  border:1px solid color-mix(in srgb,var(--pc) 32%,transparent);border-radius:99px;padding:.16rem .6rem;white-space:nowrap}
.pill b{color:var(--ink);font-weight:600;margin-left:.15rem}

/* swatches + panels */
.panel{background:var(--surface);border:1px solid var(--line);border-radius:14px;padding:18px;box-shadow:var(--shadow)}
.swatch{border:1px solid var(--line);border-radius:12px;overflow:hidden;background:var(--sink)}
.sw-chip{height:64px} .sw-body{padding:9px 11px}
.sw-name{font-family:var(--mono);font-size:.72rem;color:var(--ink)}
.sw-hex{font-family:var(--mono);font-size:.66rem;color:var(--muted);margin-top:2px}
.row{display:flex;gap:14px;flex-wrap:wrap;align-items:flex-start}
.lbl{font-family:var(--mono);font-size:.6rem;letter-spacing:.16em;text-transform:uppercase;color:var(--faint);margin:0 0 .7rem}
"""


def card(name, group, subtitle, body, cls="ds"):
    marker = f'<!-- @dsCard group="{group}" name="{name}" subtitle="{subtitle}" -->'
    return (f'{marker}\n<meta charset="utf-8"><title>{name}</title>\n'
            f'<style>{CSS}</style>\n<div class="{cls}">\n{body}\n</div>\n')


def swatch(name, var, hexv):
    return (f'<div class="swatch"><div class="sw-chip" style="background:{var}"></div>'
            f'<div class="sw-body"><div class="sw-name">{name}</div><div class="sw-hex">{hexv}</div></div></div>')


# ---- Tokens: colors ----
groups = [
    ("Canvas & surfaces", [("bg", "#0B0C10"), ("rail", "#171A22"), ("surface", "#15171E"),
                           ("surface2", "#1C1F28"), ("sink", "#111319"), ("line", "#23262F")]),
    ("Text", [("ink", "#F3F4F8"), ("body", "#C4C8D2"), ("muted", "#9096A4"), ("faint", "#6A7080")]),
    ("Accent · violet", [("violet", "#A896F2"), ("violet-2", "#8B76EC"), ("violet-deep", "#6C5CE0")]),
    ("Signal", [("pass", "#46E39B"), ("warn", "#F2B24C"), ("fail", "#F7746E")]),
]
blocks = ""
for gname, items in groups:
    sw = "".join(swatch(n, f"var(--{n})", h) for n, h in items)
    blocks += (f'<div style="margin-bottom:22px"><div class="lbl">{gname}</div>'
               f'<div class="grid" style="grid-template-columns:repeat(auto-fill,minmax(150px,1fr))">{sw}</div></div>')
colors = card("Color tokens", "Tokens", "Canvas, text, violet accent, signal",
              '<div class="eyebrow">SpatialScribe · design system</div>'
              '<h2 style="margin-bottom:20px">Color</h2>' + blocks)

# ---- Tokens: type ----
typ = card("Type scale", "Tokens", "Space Grotesk · Manrope · JetBrains Mono",
    '<div class="eyebrow">SpatialScribe · design system</div><h2 style="margin-bottom:22px">Typography</h2>'
    '<div class="row">'
    '<div class="panel" style="flex:1;min-width:300px">'
    '<div class="lbl">Display · Space Grotesk</div>'
    '<div style="font-family:var(--display);color:var(--ink);font-weight:700;font-size:2.4rem;letter-spacing:-.03em">Quality control</div>'
    '<div style="font-family:var(--display);color:var(--ink);font-weight:700;font-size:2.4rem;letter-spacing:-.03em;text-shadow:0 0 22px rgba(168,150,242,.22)">100,000</div>'
    '<div style="color:var(--muted);font-family:var(--mono);font-size:.7rem;margin-top:.6rem">headings + big readouts · 600/700</div>'
    '</div>'
    '<div class="panel" style="flex:1;min-width:300px">'
    '<div class="lbl">Body · Manrope</div>'
    '<div style="color:var(--body);font-size:1rem;line-height:1.55">Load a spatial section. The interface stays quiet so the readouts and the accent do the talking.</div>'
    '<div style="color:var(--muted);font-family:var(--mono);font-size:.7rem;margin-top:.6rem">UI + prose · 400/600</div>'
    '</div>'
    '<div class="panel" style="flex:1;min-width:300px">'
    '<div class="lbl">Data · JetBrains Mono</div>'
    '<div style="font-family:var(--mono);color:var(--ink);font-size:1.05rem">CONFIDENCE 0.32 &middot; PAS 33.5%</div>'
    '<div style="font-family:var(--mono);color:var(--violet);font-size:.9rem;margin-top:.4rem">03 · QUALITY</div>'
    '<div style="color:var(--muted);font-family:var(--mono);font-size:.7rem;margin-top:.6rem">labels, tickers, code · 400/500</div>'
    '</div></div>')

# ---- Components: buttons ----
buttons = card("Buttons", "Components", "Primary CTA, secondary, ghost, disabled",
    '<div class="eyebrow">Components</div><h3 style="font-size:1.15rem;margin-bottom:18px">Buttons</h3>'
    '<div class="row" style="align-items:center">'
    '<button class="btn primary">Load breast example</button>'
    '<button class="btn primary hover">Load breast example</button>'
    '<button class="btn">Load synthetic</button>'
    '<button class="btn ghost">View profile</button>'
    '<button class="btn disabled" disabled>Run funnel</button>'
    '</div>'
    '<div class="row" style="margin-top:18px;color:var(--muted);font-family:var(--mono);font-size:.68rem">'
    '<span>primary&nbsp;·&nbsp;glossy violet gradient + glow</span><span>hover&nbsp;·&nbsp;brighter</span>'
    '<span>secondary</span><span>ghost</span><span>disabled</span></div>')

# ---- Components: metric tiles ----
tiles = card("Metric tiles", "Components", "Readout cards with delta",
    '<div class="eyebrow">Components</div><h3 style="font-size:1.15rem;margin-bottom:18px">Readout tiles</h3>'
    '<div class="row">'
    '<div class="tile"><div class="k">cells</div><div class="v">100,000</div><div class="d up">&#8593; loaded</div></div>'
    '<div class="tile"><div class="k">median genes / cell</div><div class="v">46</div><div class="d">Xenium 5K</div></div>'
    '<div class="tile"><div class="k">spatial coherence PAS</div><div class="v">33.5%</div><div class="d up">&#8593; mean 0.37</div></div>'
    '<div class="tile"><div class="k">mean confidence</div><div class="v">0.32</div><div class="d down">&#8595; heuristic</div></div>'
    '</div>')

# ---- Components: nav rail ----
navs = ["01 Load", "02 Panel check", "03 QC", "04 Cluster", "05 Annotate", "06 Spatial", "07 Report"]
navhtml = "".join(f'<div class="nav{" active" if n=="03 QC" else ""}">{n}</div>' for n in navs)
rail = card("Nav rail", "Components", "Sidebar with active violet pill",
    '<div class="eyebrow">Components</div><h3 style="font-size:1.15rem;margin-bottom:18px">Navigation rail</h3>'
    f'<div class="rail"><div class="brand">SpatialScribe</div><div class="brandsub">spatial copilot</div>'
    f'<div class="status">GPU &middot; 100,000 cells</div>'
    f'<div class="railbar"><div class="top"><span>run progress</span><span>3 / 7 done</span></div>'
    f'<div class="track"><div class="fill" style="width:43%"></div></div></div><hr>{navhtml}</div>', cls="ds")

# ---- Components: section header ----
header = card("Section header", "Components", "Kicker chip + title + sub",
    '<div class="eyebrow">Components</div><h3 style="font-size:1.15rem;margin-bottom:20px">Section header</h3>'
    '<div class="kicker">03 · Quality</div>'
    '<div class="title">Signal you can trust</div>'
    '<div class="sub">A section-level readout of the loaded sample before clustering. The kicker chip encodes the wizard step.</div>')

# ---- Components: signal meter ----
meter = card("Signal meter", "Components", "Confident / tentative / abstained",
    '<div class="eyebrow">Components</div><h3 style="font-size:1.15rem;margin-bottom:18px">Confidence signal bar</h3>'
    '<div class="meter" style="max-width:640px">'
    '<div class="bar"><div class="seg" style="width:62%;background:var(--pass);color:var(--pass)"></div>'
    '<div class="seg" style="width:24%;background:var(--warn);color:var(--warn)"></div>'
    '<div class="seg" style="width:14%;background:var(--fail);color:var(--fail)"></div></div>'
    '<div class="legend">'
    '<span><span class="dot" style="color:var(--pass)"></span>confident <b>62%</b></span>'
    '<span><span class="dot" style="color:var(--warn)"></span>tentative <b>24%</b></span>'
    '<span><span class="dot" style="color:var(--fail)"></span>abstained <b>14%</b></span>'
    '<span style="color:var(--faint)">&middot; never a confident wrong label</span></div></div>')

# ---- Components: reason pills ----
pills = card("Abstention pills", "Components", "Reason-coded exclusion chips",
    '<div class="eyebrow">Components</div><h3 style="font-size:1.15rem;margin-bottom:18px">Abstention pills</h3>'
    '<div class="pills"><span class="pills-lbl">excluded</span>'
    '<span class="pill" style="--pc:var(--muted)">low confidence <b>13,175</b></span>'
    '<span class="pill" style="--pc:var(--fail)">low quality <b>10,588</b></span>'
    '<span class="pill" style="--pc:var(--warn)">mixed / doublet <b>9,705</b></span>'
    '<span class="pill" style="--pc:var(--violet)">panel can\'t resolve <b>640</b></span></div>')

# ---- Overview: composed mini dashboard ----
overview = card("Dashboard overview", "Overview", "Rail + header + tiles + meter",
    '<div class="row" style="gap:22px;flex-wrap:nowrap">'
    f'<div class="rail" style="flex:0 0 230px"><div class="brand">SpatialScribe</div>'
    f'<div class="brandsub">spatial copilot</div><div class="status">CPU &middot; 100,000 cells</div><hr>{navhtml}</div>'
    '<div style="flex:1;min-width:0">'
    '<div class="kicker">03 · Quality</div><div class="title">Quality control</div>'
    '<div class="row" style="margin:18px 0 20px">'
    '<div class="tile"><div class="k">cells</div><div class="v">100,000</div><div class="d up">&#8593; loaded</div></div>'
    '<div class="tile"><div class="k">median genes / cell</div><div class="v">46</div><div class="d">Xenium 5K</div></div>'
    '<div class="tile"><div class="k">mean confidence</div><div class="v">0.32</div><div class="d down">&#8595; heuristic</div></div>'
    '</div>'
    '<div class="panel">'
    '<div class="meter"><div class="bar"><div class="seg" style="width:62%;background:var(--pass);color:var(--pass)"></div>'
    '<div class="seg" style="width:24%;background:var(--warn);color:var(--warn)"></div>'
    '<div class="seg" style="width:14%;background:var(--fail);color:var(--fail)"></div></div>'
    '<div class="legend"><span><span class="dot" style="color:var(--pass)"></span>confident <b>62%</b></span>'
    '<span><span class="dot" style="color:var(--warn)"></span>tentative <b>24%</b></span>'
    '<span><span class="dot" style="color:var(--fail)"></span>abstained <b>14%</b></span></div>'
    '<div class="pills" style="margin-top:16px"><span class="pills-lbl">excluded</span>'
    '<span class="pill" style="--pc:var(--muted)">low confidence <b>13,175</b></span>'
    '<span class="pill" style="--pc:var(--fail)">low quality <b>10,588</b></span>'
    '<span class="pill" style="--pc:var(--warn)">mixed / doublet <b>9,705</b></span></div>'
    '<div style="margin-top:18px"><button class="btn primary">Run QC funnel</button></div>'
    '</div></div></div>', cls="ds")

files = {
    "tokens-colors.html": colors, "tokens-type.html": typ, "buttons.html": buttons,
    "metric-tiles.html": tiles, "nav-rail.html": rail, "section-header.html": header,
    "signal-meter.html": meter, "reason-pills.html": pills, "overview.html": overview,
}
for fn, html in files.items():
    (OUT / fn).write_text(html)

readme = """# SpatialScribe - Design System bundle

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
"""
(OUT / "README.md").write_text(readme)
print("wrote", len(files) + 2, "files to", OUT)
for p in sorted(OUT.iterdir()):
    print(f"  {p.name}  ({p.stat().st_size} bytes)")
