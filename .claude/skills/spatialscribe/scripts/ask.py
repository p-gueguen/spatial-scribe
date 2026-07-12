#!/usr/bin/env python
"""One-shot SpatialScribe copilot: ask a plain-language question about a section.

Loads a section and hands it to the grounded Claude copilot (`agent.tools.run_copilot`), a
whitelisted tool-use loop over the capability registry - it ACTS (runs real capabilities) and
answers only from their computed numbers, never free-text speculation. Needs ANTHROPIC_API_KEY.

Usage:
    python ask.py --path <dir_or.h5ad> --tissue "human breast" "are the T cells excluded from the tumour?"
    python ask.py --demo "which cell types can this panel resolve?"
"""
from __future__ import annotations

import argparse
import os
import sys


def main() -> int:
    ap = argparse.ArgumentParser(description="Ask the SpatialScribe copilot one question.")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--path", help="Section: platform output dir or .h5ad.")
    src.add_argument("--demo", action="store_true", help="Use the synthetic melanoma demo.")
    ap.add_argument("--tissue", default="melanoma", help="Free-text tissue/organism context.")
    ap.add_argument("--reference", help="Optional single-cell reference .h5ad.")
    ap.add_argument("--figures-dir", help="If set, save any figure artifacts (PNG) here.")
    ap.add_argument("question", help="The plain-language question.")
    args = ap.parse_args()

    if not os.environ.get("ANTHROPIC_API_KEY"):
        print("ANTHROPIC_API_KEY is not set - the copilot needs it. `export ANTHROPIC_API_KEY=sk-...`",
              file=sys.stderr)
        return 2

    from spatialscribe.analysis import demo, export, io, reference as _ref
    from spatialscribe.agent.tools import run_copilot

    adata = demo.load_demo().adata if args.demo else io.load(args.path).adata
    ref = ref_key = None
    if args.reference:
        ref, ref_key = _ref.load_reference(args.reference)

    artifacts: list = []
    answer = run_copilot(adata, args.question, tissue=args.tissue, use_llm=True,
                         artifacts=artifacts, reference=ref, ref_label_key=ref_key)
    print(answer)

    if args.figures_dir and artifacts:
        os.makedirs(args.figures_dir, exist_ok=True)
        for i, art in enumerate(artifacts):
            fig = art.get("fig") if isinstance(art, dict) else None
            if fig is not None:
                p = export.save_figure_png(fig, os.path.join(args.figures_dir, f"artifact_{i}.png"))
                if p:
                    print(f"[figure] {p}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
