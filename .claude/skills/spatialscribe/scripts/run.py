#!/usr/bin/env python
"""Headless SpatialScribe driver: raw section -> QC'd, routed, annotated .h5ad + report + script.

A thin wrapper over ``analysis.pipeline.run_pipeline`` - the ONE spine both the CLI and the app's
"Run full analysis" job drive, so the autonomous run and the interactive app never diverge:

    load -> compute_qc -> panel_check -> reference_match -> annotation_strategy (the gate)
         -> cluster -> { reference_transfer | annotate | (clusters only) } -> niches
         -> malignant_concordance (tumour-gated) -> split_purify (opt-in) -> self_heal -> export

`annotation_strategy` runs the coarsen -> reselect ladder and recommends supervised label transfer
vs unsupervised de-novo clustering; the pipeline obeys that recommendation. Every stage is
degrade-graceful and its honest status (ok / skipped / failed) is printed at the end.

Usage:
    python run.py --path <xenium_dir_or.h5ad> --tissue "human breast" [--reference ref.h5ad] \
                  [--tumour | --no-tumour] [--rctd] [--split] --out results/
    python run.py --demo --out results/            # synthetic melanoma demo, no data needed
"""
from __future__ import annotations

import argparse
import os
import sys


def main() -> int:
    ap = argparse.ArgumentParser(description="Headless SpatialScribe run (full pipeline).")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--path", help="Xenium/CosMx/MERSCOPE/Atera output dir, or a .h5ad section.")
    src.add_argument("--demo", action="store_true", help="Use the synthetic melanoma demo section.")
    ap.add_argument("--tissue", default="melanoma", help="Free-text tissue/organism context.")
    ap.add_argument("--reference", help="Optional single-cell reference .h5ad for label transfer.")
    ap.add_argument("--label-key", help="Reference cell-type column (optional; auto-detected).")
    ap.add_argument("--resolution", type=float, default=1.0, help="Leiden resolution (default 1.0).")
    tum = ap.add_mutually_exclusive_group()
    tum.add_argument("--tumour", dest="tumour", action="store_const", const=True,
                     help="Force the malignant-calling gate ON (else inferred from the tissue name).")
    tum.add_argument("--no-tumour", dest="tumour", action="store_const", const=False,
                     help="Force the malignant-calling gate OFF.")
    ap.set_defaults(tumour=None)
    ap.add_argument("--rctd", action="store_true",
                    help="Use RCTD as the annotate engine (needs a reference + SPATIALSCRIBE_RCTD_PYTHON).")
    ap.add_argument("--split", action="store_true",
                    help="Run SPLIT spillover purification (needs a reference + R).")
    ap.add_argument("--out", default="spatialscribe_out", help="Output directory.")
    args = ap.parse_args()

    from spatialscribe.analysis import capabilities as cap
    from spatialscribe.analysis import demo, io, pipeline
    from spatialscribe.analysis import reference as _ref

    # 1. load
    if args.demo:
        adata = demo.load_demo().adata
        print(f"[load] synthetic demo: {adata.n_obs:,} cells x {adata.n_vars} genes")
    else:
        adata = io.load(args.path).adata
        print(f"[load] {args.path}: {adata.n_obs:,} cells x {adata.n_vars} genes "
              f"(platform={adata.uns.get('platform')})")

    ref = ref_key = None
    if args.reference:
        ref, ref_key = _ref.load_reference(args.reference, label_key=args.label_key)
        print(f"[reference] {args.reference}: label_key={ref_key!r}, {ref.n_obs:,} cells")

    ctx = cap.RunContext(tissue=args.tissue, use_llm=bool(os.environ.get("ANTHROPIC_API_KEY")),
                         reference=ref, ref_label_key=ref_key, is_tumour=args.tumour)
    opts = pipeline.PipelineOptions(resolution=args.resolution, rctd=args.rctd, split=args.split,
                                    export=True, out_dir=args.out,
                                    source_path=None if args.demo else args.path)

    # 2. run the full spine (obeys the gate; exports the annotated .h5ad + HTML report + rerun.py)
    res = pipeline.run_pipeline(adata, ctx, opts)

    # 3. report each stage's honest status
    print(f"[route] {res['route']}")
    for s in res["stages"]:
        kind, _, why = s["status"].partition(":")
        print(f"  [{kind:7}] {s['name']}" + (f"  ({why})" if why else ""))
    ran, skipped, failed = (res["summary"][k] for k in ("ok", "skipped", "failed"))
    print(f"[done] cell_type present: {'cell_type' in adata.obs}  |  out: {args.out}")
    print(f"       ran={len(ran)} skipped={len(skipped)} failed={len(failed)}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
