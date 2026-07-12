# 3-minute demo storyboard (product-first, general audience)

Two guiding principles (per Anthropic's advice): speak to a **general audience**, not spatial-
transcriptomics specialists, and spend the time on **what the finished tool does**, not on how it
was built. Everything shown is real, reproducible output on the bundled public 10x Prime 5K breast
Xenium section - no staged results.

Target 3:00, voiceover ~430 words. One continuous screen-recording of the live tool.

---

**0:00 - The problem (20s).** *On screen:* a whole tumor slide zooming to individual cells; the
number "700,000 cells."
> "A cancer lab runs an experiment and gets back a map of a tumor at single-cell resolution -
> hundreds of thousands of cells, each one placed exactly where it sits in the tissue. It should be
> a goldmine. But today the biologist who ran it usually can't read it: the slide waits weeks for a
> specialist to analyze it by hand, in code. SpatialScribe lets the scientist do it themselves, in
> an afternoon."

**0:20 - Load, and an honest first answer (30s).** *On screen:* click *Load breast example*; open
**Panel check**.
> "One click loads a real section. Before it labels a single cell, it does something most tools
> skip: it tells you, in plain language, which cell types this experiment can and cannot reliably
> tell apart. It's honest about its own limits up front, instead of labeling everything and hoping."

**0:50 - Quality, and dropping the bad parts (25s).** *On screen:* the quality tiles, then draw a
loop around a damaged patch of tissue on the map; that region's quality updates on its own.
> "It checks the data quality for you. And if part of the slide is folded or dead tissue, you just
> circle it and drop it, so a bad patch doesn't poison the whole result. Biologists asked for this
> directly."

**1:15 - Cell types, with the confidence to say 'I don't know' (35s).** *On screen:* run
annotation; the map fills with colored cell types; show the confidence headline and the map colored
by confidence.
> "Now it identifies the cell types and paints the map. Here's the part that makes it trustworthy:
> it refuses to fake a label it can't support. The cells it's sure of get named. The cells it can't
> confidently call are set aside and flagged with the reason - the signal was too weak, or this
> experiment simply can't tell those two types apart. You never get a confident wrong answer, which
> is exactly what you need before you build a conclusion on it."

**1:50 - Just ask (the hero, 45s).** *On screen:* type in the chat: *"Are the immune cells getting
into the tumor, or locked out?"* The copilot runs the real analysis, draws the plot, and answers
with the finding on screen.
> "And you don't need to know the tools. You just ask. 'Are the immune cells getting into the tumor,
> or locked out?' It runs the real spatial analysis on your data and answers: the T cells are
> strongly shut out of the tumor - and it shows you exactly where. Ask it 'what neighborhoods are in
> this tissue?' and it finds them and names them. The important part: every answer comes from an
> analysis it actually ran. It is wired so it cannot make a number up."

**2:35 - Zoom in, and take it with you (20s).** *On screen:* click the T-cell population,
**Subcluster**; then **Report** and the download.
> "Click any population and split it into finer subtypes on the spot. Then export a full report
> where every figure traces back to code you can re-run - so the result is reproducible, not a
> black box."

**2:55 - Close (10s).** *On screen:* the tissue map, then the public GitHub page.
> "A wet-lab scientist did this alone, on real data, in an afternoon instead of waiting weeks. It's
> open source, and it runs on a laptop."

---

## What changed from the earlier storyboard (and why)

- **Cut the "how I built it" montage** (planning, worktrees, subagents, Playwright, dashboards). It
  is genuinely interesting, but the guidance is to show the product, not the process. Keep those
  points for the written submission / a follow-up post, not the 3-min video.
- **Removed specialist jargon** from the spoken track: no AUC, no z-scores, no "cross-method
  ensemble", no "abstention basis". The real numbers still appear on screen for credibility; the
  voiceover translates them into plain outcomes ("locked out of the tumor").
- **Led with the human problem and the payoff**, kept the copilot as the emotional peak, and made
  the close about the outcome (afternoon vs weeks), not the infrastructure.

## Recording tips

- Point the copilot at **Anthropic (Sonnet)** for the recording, not the internal endpoint: set
  `ANTHROPIC_API_KEY` and `ANTHROPIC_MODEL=claude-sonnet-5`. Sonnet abstains more honestly and
  fabricates far less than the Haiku default, and it avoids the internal-endpoint tool-parsing
  failure mode. Before you hit record, send one real copilot question and confirm it answers with a
  tool result, not raw text.
- `export SPATIALSCRIBE_FORCE_CPU=1` for a clean, reproducible run that any judge can repeat.
- Pre-click the demo once so the cache is warm and every step is instant on camera.
- Do one full dry run end to end before the real take - the live copilot and the first-render of
  the map are the two things most likely to surprise you.
