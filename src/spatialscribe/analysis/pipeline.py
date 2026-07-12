"""The full analysis spine as ONE degrade-graceful function.

`run_pipeline` runs load-time-annotated section -> QC -> panel adequacy -> the supervised-vs-clustering
gate -> cluster -> annotate -> niches -> malignant calling -> spillover purification -> self-heal ->
export, driving the SAME capability registry the app and copilot use. It is shared by the headless CLI
(`scripts/run.py`) and the app's background "Run full analysis" job, so the autonomous run and the
interactive app can never diverge.

Design contract (see docs/superpowers/specs/2026-07-11-spatialscribe-preprocessing-pipeline-design.md):
- Every stage is wrapped: a skip or failure is RECORDED and the pipeline continues (never aborts on an
  optional arm). The result lists each stage's honest status ``ok | skipped:<why> | failed:<why>`` -
  so a caller asserts ``status == "ok"``, never mere presence (the repo's silent-skip gotcha).
- The expensive / env-gated arms (RCTD engine, SPLIT, learned CNV) stay opt-in and degrade to
  ``skipped`` when unconfigured - exactly as the individual capabilities already do.
"""
from __future__ import annotations

import os
from dataclasses import dataclass

from . import capabilities as cap

# Nominal stage count, only used to drive the progress fraction (ctx.tick). The real set of stages run
# depends on the gate + which optional arms are configured; the fraction is advisory, not exact.
_NOMINAL_STAGES = 12


@dataclass
class PipelineOptions:
    """Knobs for a full run. The malignant gate is read from ``ctx.is_tumour`` (single source)."""
    resolution: float = 1.0
    rctd: bool = False           # use RCTD (annotate_rctd) as the annotate engine when a ref + env allow
    split: bool = False          # run SPLIT spillover purification (needs a reference + R)
    export: bool = True          # write annotated .h5ad + HTML report + rerun.py (needs out_dir)
    out_dir: str | None = None
    source_path: str | None = None  # the loaded section's path -> threaded into the re-runnable script
    max_heal_iter: int = 3


def run_pipeline(adata, ctx, opts: PipelineOptions | None = None) -> dict:
    """Run the full spine on ``adata`` (mutated in place). Returns a per-stage status record.

    Parameters
    ----------
    adata : AnnData
        A loaded section (``io.load(...).adata``); mutated in place (obs cell_type / niche / malignant
        scores, uns route / pipeline record).
    ctx : capabilities.RunContext
        Carries tissue, ``is_tumour`` (the malignant gate), an optional reference, and ``tick``.
    opts : PipelineOptions

    Returns
    -------
    dict : ``{"route", "stages": [{"name", "status"}], "action_log", "summary"}``.
    """
    from . import cnv as _cnv

    opts = opts or PipelineOptions()
    log: list[dict] = []
    stages: list[dict] = []

    def _record(name: str, status: str) -> None:
        stages.append({"name": name, "status": status})

    def _run(name: str, params: dict | None = None, *, optional: bool = False):
        # Map the capability's OWN local 0..1 progress into THIS stage's global band
        # [idx/N, (idx+1)/N] so the bar advances monotonically. Before this, the stage-start tick
        # (idx/N) and the capability's internal ticks (cluster: 0.05->1.0; compute_qc/annotate via
        # progress=ctx.tick) overwrote each other, so the bar jumped to ~1.0 inside a stage then
        # snapped back to the next stage's low fraction - the "wrong order" the user sees. Swapping
        # ctx.progress catches both direct ctx.tick calls and progress=ctx.tick sub-calls (tick reads
        # self.progress at call time). Clamped at 0.99 so a longer-than-nominal run never exceeds the bar.
        idx = len(stages)
        lo = min(idx / _NOMINAL_STAGES, 0.99)
        hi = min((idx + 1) / _NOMINAL_STAGES, 0.99)
        outer = ctx.progress

        def _banded(frac, label):
            if outer is not None:
                outer(lo + (hi - lo) * max(0.0, min(1.0, float(frac))), str(label))

        ctx.progress = _banded
        try:
            ctx.tick(0.0, name)                                  # stage start -> lo
            res = cap.run(adata, name, params or {}, ctx)
        finally:
            ctx.progress = outer                                 # restore the app sink for the next stage
        if res.ok:
            val = res.value
            # A capability can succeed (ok=True) yet self-report a degrade: {"status": "skipped: ..."}.
            if isinstance(val, dict) and str(val.get("status", "")).startswith("skipped"):
                _record(name, "skipped:" + str(val["status"]).split("skipped:", 1)[-1].strip())
            else:
                _record(name, "ok")
                if res.record:
                    log.append(res.record)
            return val
        # ok=False: a genuine error. Optional arms downgrade to a skip; the spine surfaces the failure.
        err = res.error or {}
        why = str(err.get("error_type") or err.get("message") or err)
        _record(name, ("skipped:" if optional else "failed:") + why)
        return None

    # 1-4. QC, panel adequacy, reference match (only when a reference is loaded).
    _run("compute_qc")
    _run("panel_check")
    if ctx.reference is not None:
        _run("reference_match", optional=True)

    # 5. THE GATE: self-verify the reference<->panel match, decide supervised vs de-novo clustering.
    route_val = _run("annotation_strategy")
    route = (route_val or {}).get("recommended_mode", "annotate") if isinstance(route_val, dict) else "annotate"

    # 6. Unsupervised backbone (embedding + leiden), always.
    _run("cluster", {"resolution": opts.resolution})

    # 7. Annotate per the gate. RCTD is the opt-in engine; it falls back to reference_transfer if it
    #    skips (no env), and the gate can route to de-novo clusters when the reference cannot resolve
    #    the panel (never forces supervised labels the data cannot support).
    if route == "reference_transfer" and ctx.reference is not None:
        if opts.rctd:
            _run("annotate_rctd", {"set_as_primary": True}, optional=True)
        if "cell_type" not in adata.obs:
            _run("reference_transfer", {"set_as_primary": True})
    elif route == "cluster":
        _record("annotate", "skipped:gate routed to de-novo clusters (reference cannot resolve the panel)")
    elif "cell_type" in adata.obs:
        # Idempotent: a section that already carries labels (a loaded demo cache, or any pre-annotated
        # upload) keeps them - re-running the de-novo cluster-marker annotator CLOBBERS higher-quality
        # reference labels (the atera RCTD cache collapses its 10 per-cell types to a 4-type per-cluster
        # call). Mirrors the reference_transfer branch above, which already guards on cell_type presence;
        # the Annotate tab's own button still forces a re-annotation when the user wants one.
        _record("annotate", "skipped:section already annotated (labels preserved)")
    else:
        _run("annotate")

    # 8. Spatial niches, malignant calling (tumour-gated), spillover purification, self-heal.
    if "cell_type" in adata.obs:
        _run("niches", optional=True)

        is_tumour, gate = _cnv.is_tumour_context(ctx.tissue, ctx.is_tumour)
        if is_tumour:
            _run("malignant_concordance", optional=True)
        else:
            _record("malignant_concordance", f"skipped:non-tumour context ({gate})")

        if opts.split:
            _run("split_purify", optional=True)

        # Enhanced self-verify + re-run loop (verify.autorerun): a failing marker-check re-annotates
        # or abstains, so a type can leave the failing set. Replaces the bare advisory self_verify.
        _run("self_heal", {"max_rounds": opts.max_heal_iter}, optional=True)
    else:
        _record("niches", "skipped:no cell_type (clusters-only route)")
        _record("self_heal", "skipped:no cell_type (clusters-only route)")

    # Six-layer QC funnel - runs on EVERY route, not only when annotation produced labels. run_funnel
    # ALWAYS does the section-level layers (Layer 0 section, Layer 1 segmentation, Layer 2 count floor)
    # and adds the label-dependent layers (purity / coherence / per-cell confidence + the section AQI)
    # only when cell_type exists. It used to live INSIDE `if cell_type`, so a clusters-only or
    # annotate-failed section got NO section QC at all - exactly the layers that never needed labels.
    # It also scores the FINAL (self-healed) labels + stamps the AQI on uns for the exported report;
    # optional, so a degenerate section never aborts the run.
    fun = _run("qc_funnel", optional=True)
    if isinstance(fun, dict) and fun.get("annotation_quality") is not None:
        adata.uns["annotation_quality"] = fun["annotation_quality"]

    # 9. Export (headless CLI only; the app displays the populated session in place).
    if opts.export and opts.out_dir:
        ctx.tick(0.99, "export")
        from . import export as _export
        try:
            os.makedirs(opts.out_dir, exist_ok=True)
            _export.export_h5ad(adata, os.path.join(opts.out_dir, "annotated.h5ad"))
            _export.render_analysis_report(adata, log, ctx.tissue,
                                           os.path.join(opts.out_dir, "report.html"),
                                           source_path=opts.source_path)
            _export.export_script(log, os.path.join(opts.out_dir, "rerun.py"), ctx.tissue, adata,
                                  source_path=opts.source_path)
            _record("export", "ok")
        except Exception as exc:  # noqa: BLE001 - export is best-effort; never abort the run on it
            _record("export", f"failed:{exc}")

    ctx.tick(1.0, "done")
    adata.uns["pipeline"] = {"route": route, "stages": stages}
    summary = {k: [s["name"] for s in stages if s["status"].split(":", 1)[0] == k]
               for k in ("ok", "skipped", "failed")}
    return {"route": route, "stages": stages, "action_log": log, "summary": summary}
