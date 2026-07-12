"""Best-effort background precompute - hide the latency of a deferred, off-critical-path step.

What it does
------------
SpatialScribe keeps the expensive-but-non-blocking work off the synchronous wizard flow: chiefly
the UMAP embedding, which nothing on the load->qc->cluster->annotate->spatial critical path consumes
(only the UMAP plot reads ``obsm['X_umap']``), yet is ~75% of the CPU pipeline. This module runs such
a step in a daemon thread so it is ready by the time the user opens the view that needs it, WITHOUT
blocking the step that kicked it off.

How to use it
-------------
>>> from . import preload
>>> preload.start(adata, "umap", ctx)   # returns immediately; idempotent per (adata, step)
>>> ...                                  # the wizard keeps going - leiden + the spatial canvas render now
>>> preload.ready(adata, "umap")         # True once obsm['X_umap'] is in place
>>> preload.result(adata, "umap")        # join the background task (blocks) if you need it synchronously

Thread-safety contract: a preloaded step may only WRITE ``adata`` keys the critical path never writes
concurrently. UMAP writes ``obsm['X_umap']`` + ``uns['umap']``; annotate/spatial write ``obs`` columns
and read ``X`` - disjoint, so the background embed cannot corrupt the synchronous path. The on-demand
``capabilities.ensure(adata, step)`` in the view path is the correctness guarantee; this only hides latency.

Depends on
----------
Standard-library ``threading`` + ``capabilities.ensure``. No heavy imports at module load.
"""
from __future__ import annotations

import threading

_TASKS: dict[tuple[int, str], threading.Thread] = {}
_LOCK = threading.Lock()
_warmed = False


def ready(adata, name: str) -> bool:
    """True once every key capability ``name`` produces is present on ``adata`` (the preload has landed)."""
    from . import capabilities as cap

    c = cap.REGISTRY.get(name)
    if c is None or not c.produces:
        return False
    return all(k.present(adata) for k in c.produces)


def start(adata, name: str, ctx=None) -> None:
    """Kick off ``capabilities.ensure(adata, name, ctx)`` in a daemon thread (idempotent).

    A no-op if the step is already produced, or a task for this ``(adata, name)`` is already
    running - so calling it on every Streamlit rerun starts at most one thread.
    """
    from . import capabilities as cap

    if ready(adata, name):
        return
    key = (id(adata), name)
    with _LOCK:
        t = _TASKS.get(key)
        if t is not None and t.is_alive():
            return
        t = threading.Thread(target=cap.ensure, args=(adata, name, ctx),
                             name=f"preload:{name}", daemon=True)
        _TASKS[key] = t
        t.start()


def result(adata, name: str, timeout: float | None = None) -> bool:
    """Block until the preload thread for ``(adata, name)`` finishes (or ``timeout``); return :func:`ready`."""
    with _LOCK:
        t = _TASKS.get((id(adata), name))
    if t is not None:
        t.join(timeout)
    return ready(adata, name)


def warm_backend() -> None:
    """Warm the numba JIT (pynndescent neighbors + UMAP) in a background daemon thread. Idempotent.

    scanpy switches to APPROXIMATE neighbors above ~4k cells, and pynndescent's first call JIT-compiles
    for ~30 s - a one-time, process-wide cost that otherwise lands on the user's FIRST real clustering.
    Running neighbors + UMAP once on a throwaway synthetic section at app startup compiles those kernels
    off the interactive path, so the first Cluster click is fast. Never touches user data; best-effort
    (a warm-up failure is swallowed - it only ever costs a little startup CPU)."""
    global _warmed
    with _LOCK:
        if _warmed:
            return
        _warmed = True
    threading.Thread(target=_warm_backend_impl, name="preload:warm", daemon=True).start()


def _warm_backend_impl() -> None:
    try:
        import anndata as ad
        import numpy as np

        from .backend import get_backend

        rng = np.random.default_rng(0)
        a = ad.AnnData(rng.poisson(1.0, size=(6000, 60)).astype("float32"))  # >4k -> approximate path
        be = get_backend()
        be.pca(a, n_comps=30)
        be.neighbors(a)
        be.umap(a)
    except Exception:
        pass   # warm-up is best-effort; never surface an error onto the app
