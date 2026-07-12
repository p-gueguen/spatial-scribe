"""Progress reporting (CONTRACT section 1, engine).

Every long capability, run with a recording ``ctx.progress`` sink, must emit a monotonically
non-decreasing sequence of fracs in [0,1] that ends at 1.0; and every capability must still run
when no sink is installed (the engine stays UI-agnostic). ``ctx.tick`` clamps and never raises.
"""
from __future__ import annotations


def _run_with_progress(adata, name, params):
    from spatialscribe.analysis import capabilities as cap

    events: list[tuple[float, str]] = []
    ctx = cap.RunContext(tissue="melanoma", use_llm=False,
                         progress=lambda f, l: events.append((f, l)))
    res = cap.run(adata, name, params, ctx)
    assert res.ok, (name, res.error)
    return events


def _assert_monotone_to_one(events, name):
    assert events, f"{name}: no progress emitted"
    fracs = [f for f, _ in events]
    assert all(0.0 <= f <= 1.0 for f in fracs), (name, fracs)
    assert all(b >= a for a, b in zip(fracs, fracs[1:])), (name, fracs)
    assert fracs[-1] == 1.0, (name, fracs)
    assert all(isinstance(l, str) and l for _, l in events), name


def test_compute_qc_progress(raw_adata):
    _assert_monotone_to_one(_run_with_progress(raw_adata, "compute_qc", {}), "compute_qc")


def test_cluster_progress(raw_adata):
    _assert_monotone_to_one(_run_with_progress(raw_adata, "cluster", {"resolution": 1.0}), "cluster")


def test_qc_funnel_progress(processed_adata):
    _assert_monotone_to_one(_run_with_progress(processed_adata, "qc_funnel", {}), "qc_funnel")


def test_annotate_progress(raw_adata, ctx):
    from spatialscribe.analysis import capabilities as cap

    assert cap.run(raw_adata, "cluster", {"resolution": 1.0}, ctx).ok   # leiden prereq
    _assert_monotone_to_one(_run_with_progress(raw_adata, "annotate", {}), "annotate")


def test_niches_progress(processed_adata):
    _assert_monotone_to_one(_run_with_progress(processed_adata, "niches", {}), "niches")


def test_progress_none_still_runs(raw_adata):
    """A capability with no progress sink installed runs exactly as before."""
    from spatialscribe.analysis import capabilities as cap

    ctx = cap.RunContext(tissue="melanoma", use_llm=False)
    assert ctx.progress is None
    for name, params in [("compute_qc", {}), ("cluster", {"resolution": 1.0})]:
        assert cap.run(raw_adata, name, params, ctx).ok, name


def test_tick_clamps_and_never_raises():
    from spatialscribe.analysis import capabilities as cap

    rec: list[float] = []
    ctx = cap.RunContext(progress=lambda f, l: rec.append(f))
    ctx.tick(2.0, "over")     # clamps to 1.0
    ctx.tick(-1.0, "under")   # clamps to 0.0
    assert rec == [1.0, 0.0]

    def _boom(_f, _l):
        raise RuntimeError("sink is broken")

    cap.RunContext(progress=_boom).tick(0.5, "x")   # swallowed, must not raise
