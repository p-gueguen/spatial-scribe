"""Headless CLI for the full SpatialScribe analysis spine.

    python scripts/run.py --demo --out DIR                 # synthetic section, no data needed
    python scripts/run.py --section PATH --out DIR [--tissue TXT] [--tumour|--no-tumour] \
                          [--resolution R] [--reference REF.h5ad [--ref-label-key COL]] \
                          [--rctd] [--split]

Loads a spatial section from disk and runs ``analysis.pipeline.run_pipeline`` - the SAME
degrade-graceful, registry-driven spine the app's "Run full analysis" job and the copilot drive, so
the headless run can never diverge from the interactive one. Writes the standard outputs to DIR:
``annotated.h5ad`` + ``report.html`` + ``rerun.py`` (the export stage) plus ``run.json``, the
machine-readable per-stage status record an external agent parses.

Exit code: ``0`` clean, ``1`` iff a stage FAILED (an optional arm that reports ``skipped:`` is NOT a
failure), ``2`` if the section could not be loaded at all. So a caller asserts success at the process
level - never on mere presence of an output - matching the pipeline's honest status contract.

This closes the headless-CLI reference in ``analysis/pipeline.py`` and gives agents a shell-out entry.
"""
from __future__ import annotations

import argparse
import json
import os
import sys


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="spatialscribe-run",
                                 description="Run the full SpatialScribe analysis spine headlessly.")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--section", "--path", dest="section", default=None,
                     help="server-side path: a Xenium/CosMx/MERSCOPE output folder or a .h5ad")
    src.add_argument("--demo", action="store_true",
                     help="run on the built-in synthetic melanoma section (no data files needed)")
    ap.add_argument("--out", required=True,
                    help="output directory (annotated.h5ad + report.html + rerun.py + run.json)")
    ap.add_argument("--tissue", default="",
                    help="tissue context; blank = infer from the panel metadata (parity with the app)")
    tum = ap.add_mutually_exclusive_group()
    tum.add_argument("--tumour", dest="tumour", action="store_const", const=True, default=None,
                     help="the sample contains malignant cells (drives the malignant gate)")
    tum.add_argument("--no-tumour", dest="tumour", action="store_const", const=False,
                     help="the sample has no malignant cells. Omit BOTH to fall back to the "
                          "tissue-keyword heuristic.")
    ap.add_argument("--resolution", type=float, default=0.5,
                    help="leiden clustering resolution (default 0.5, matching the app's cluster cap)")
    ap.add_argument("--reference", default=None,
                    help="optional single-cell reference (.h5ad) for supervised annotation")
    ap.add_argument("--ref-label-key", default=None,
                    help="cell-type column in the reference (auto-detected if omitted)")
    ap.add_argument("--rctd", action="store_true",
                    help="use RCTD as the annotate engine (needs a reference + env; else skips)")
    ap.add_argument("--split", action="store_true",
                    help="run SPLIT spillover purification (needs a reference + R; else skips)")
    args = ap.parse_args(argv)

    from spatialscribe.analysis import capabilities as cap
    from spatialscribe.analysis import io as _io
    from spatialscribe.analysis import llm
    from spatialscribe.analysis import pipeline as _pl

    if args.demo:
        from spatialscribe.analysis import demo as _demo
        path = "synthetic-melanoma-demo"      # built in; no data files needed
        a = _demo.load_demo().adata
    else:
        # Resolve the path the way the backend load_section does (corrects a dropped mount prefix).
        path = cap._resolve_server_path((args.section or "").strip())
        if path is None:
            print(f"error: section not found on the server: {args.section!r}", file=sys.stderr)
            return 2
        try:
            a = _io.load(path).adata
        except Exception as exc:  # noqa: BLE001 - a clean usage error, not a stack trace
            print(f"error: could not load section: {type(exc).__name__}: {exc}", file=sys.stderr)
            return 2
    if "spatial" not in a.obsm:
        print("error: loaded section has no spatial coordinates (obsm['spatial'])", file=sys.stderr)
        return 2

    # Tissue: explicit > inferred from the panel > app default (parity with backend load_section).
    tissue = (args.tissue or "").strip() or _io.infer_tissue(a) or "melanoma"

    # Optional reference threaded exactly like the app's _ctx (in-memory AnnData + its label column).
    reference = ref_key = None
    if args.reference:
        from spatialscribe.analysis import reference as _ref
        reference, ref_key = _ref.load_reference(args.reference, label_key=args.ref_label_key)

    ctx = cap.RunContext(tissue=tissue, use_llm=llm.available(),
                         is_tumour=args.tumour,   # True / False / None (None => tissue-keyword gate)
                         reference=reference, ref_label_key=ref_key)
    # Stream coarse progress to STDERR so a headless/SLURM run shows life without polluting stdout.
    ctx.progress = lambda f, l: print(f"[{f * 100:5.1f}%] {l}", file=sys.stderr, flush=True)

    os.makedirs(args.out, exist_ok=True)
    rec = _pl.run_pipeline(a, ctx, _pl.PipelineOptions(
        export=True, out_dir=args.out, resolution=args.resolution, rctd=args.rctd, split=args.split))

    # Machine-readable status an agent parses (h5ad/report/rerun.py come from the export stage itself).
    with open(os.path.join(args.out, "run.json"), "w") as fh:
        json.dump({"section": path, "tissue": tissue, "route": rec["route"],
                   "stages": rec["stages"], "summary": rec["summary"]}, fh, indent=2)

    print(json.dumps({"route": rec["route"], "summary": rec["summary"]}, indent=2))
    failed = rec["summary"].get("failed", [])
    if failed:
        print(f"error: {len(failed)} stage(s) FAILED: {', '.join(failed)}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
