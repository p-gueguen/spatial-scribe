"""Cancer-Finder malignant-cell probability via an isolated torch subprocess (a second opinion).

Cancer-Finder (Patchouli-M/SequencingCancerFinder) is a domain-adaptation (VREx) classifier that
labels each cell malignant/normal from its transcriptome - an INDEPENDENT malignant caller
alongside the CNV path (`cnv.call_malignant_cnv`). It needs torch + the CF repo + a pretrained
checkpoint, none of which belong in the main env, so it runs in an isolated env via
`subprocesses/cancerfinder/run_cancerfinder.py` and this joins `obs['cancerfinder_prob']` +
`obs['cancerfinder_malignant']` back. Same isolation pattern as the CNV / ovrlpy subprocesses.

Config via env vars (the cluster defaults), so the committed code carries no hard dependency:
    CANCERFINDER_PYTHON  - a python with torch + scanpy + anndata
    CANCERFINDER_REPO    - the SequencingCancerFinder checkout (models/ + utils/)
    CANCERFINDER_CKPT    - the pretrained checkpoint (e.g. sc_pretrain_article.pkl)

Caveat (internal Atera benchmark): Cancer-Finder is over-sensitive on single-cell Xenium (it agrees
with the CNV caller at ~0.98 AUROC but over-calls at a 0.5 threshold) - prefer the probability
ranking, or raise the threshold.
"""
from __future__ import annotations

import os

_CF_PY = os.environ.get("CANCERFINDER_PYTHON", "")   # a python with torch + scanpy + anndata
_CF_REPO = os.environ.get("CANCERFINDER_REPO", "")   # the SequencingCancerFinder checkout
# derive the conventional checkpoint under the repo when only CANCERFINDER_REPO is set:
_CF_CKPT = os.environ.get("CANCERFINDER_CKPT") or (
    os.path.join(_CF_REPO, "checkpoints", "sc_pretrain_article.pkl") if _CF_REPO else "")


def _join_cf(adata, parquet_path, threshold: float) -> int:
    """Join a ``{cell_id, cancerfinder_prob}`` parquet onto obs (reindexed). Returns n covered."""
    import numpy as np
    import pandas as pd

    df = pd.read_parquet(parquet_path)
    df["cell_id"] = df["cell_id"].astype(str)
    aligned = df.set_index("cell_id").reindex(adata.obs_names.astype(str))
    prob = aligned["cancerfinder_prob"].to_numpy(dtype=float)
    adata.obs["cancerfinder_prob"] = prob
    adata.obs["cancerfinder_malignant"] = np.where(np.isnan(prob), False, prob > threshold)
    return int(aligned["cancerfinder_prob"].notna().sum())


# Background runs started by prefetch(), keyed by _fingerprint of their exact input; call_cancerfinder
# joins a matching one instead of starting its own. A run is only reused for byte-identical counts,
# cells, genes and parameters (max_cells subsamples by cell count, so a different cell set is a
# different answer); anything else falls back to a synchronous run.
_PENDING: dict = {}


def _config(env_python=None, repo=None, ckpt=None):
    """Resolved (python, repo, ckpt), or a skip message when any is unset/missing."""
    from pathlib import Path

    env_python = env_python or _CF_PY
    repo = repo or _CF_REPO
    # Derive the conventional checkpoint under the EFFECTIVE repo when only repo was supplied
    # (so passing repo= alone works, not just the module-level default).
    ckpt = ckpt or (os.path.join(repo, "checkpoints", "sc_pretrain_article.pkl") if repo else _CF_CKPT)
    for label, envname, path in (("cancerfinder python", "CANCERFINDER_PYTHON", env_python),
                                 ("cancerfinder repo", "CANCERFINDER_REPO", repo),
                                 ("checkpoint", "CANCERFINDER_CKPT", ckpt)):
        if not path:
            return None, f"{label} not configured; set {envname}"
        if not Path(path).exists():
            return None, f"{label} not found ({path}); set {envname}"
    return (env_python, repo, ckpt), None


def _counts(adata):
    return adata.layers["counts"] if "counts" in adata.layers else adata.X


def _fingerprint(adata, threshold: float, max_cells: int) -> str:
    import hashlib

    import numpy as np
    import scipy.sparse as sp

    h = hashlib.blake2b(digest_size=16)
    h.update(f"{threshold!r}|{int(max_cells)}|{adata.shape}".encode())
    h.update("\0".join(adata.obs_names.astype(str)).encode())
    h.update("\0".join(adata.var_names.astype(str)).encode())
    c = _counts(adata)
    if sp.issparse(c):
        c = c.tocsr()
        for arr in (c.data, c.indices, c.indptr):
            h.update(np.ascontiguousarray(arr).view(np.uint8))
    else:
        h.update(np.ascontiguousarray(np.asarray(c)).view(np.uint8))
    return h.hexdigest()


def _launch(adata, threshold: float, max_cells: int, cfg) -> dict:
    """Export the counts and start the subprocess; returns a job dict (owns a TemporaryDirectory)."""
    import subprocess
    import tempfile
    from pathlib import Path

    import anndata as ad

    env_python, repo, ckpt = cfg
    tmp = tempfile.TemporaryDirectory(prefix="sscf_")   # removed in _join, whatever happens
    d = Path(tmp.name)
    h5, out, err = d / "section.h5ad", d / "cf.parquet", d / "stderr.log"
    a_exp = ad.AnnData(X=_counts(adata).copy())
    a_exp.obs_names = adata.obs_names.astype(str)
    a_exp.var_names = adata.var_names.astype(str)
    a_exp.write_h5ad(h5)

    script = str(Path(__file__).resolve().parents[3] / "subprocesses" / "cancerfinder" / "run_cancerfinder.py")
    cmd = [env_python, script, "--h5ad", str(h5), "--out", str(out), "--repo", repo, "--ckpt", ckpt,
           "--threshold", str(threshold), "--max-cells", str(int(max_cells))]
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}   # don't leak the app's src path
    # The app caps BLAS/OMP at 8 threads for itself; the torch forward + dense normalize here scale
    # past that (measured 12.9 s -> 8.4 s at 32 threads, identical probabilities).
    n = str(min(32, os.cpu_count() or 8))
    env.update(OMP_NUM_THREADS=n, MKL_NUM_THREADS=n, OPENBLAS_NUM_THREADS=n)
    # stderr to a file, not a pipe: a background run nobody reads yet would block on a full pipe.
    with open(err, "w") as fh:
        proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=fh, text=True, env=env)
    return {"proc": proc, "tmp": tmp, "out": out, "err": err}


def _discard(job: dict) -> None:
    try:
        if job["proc"].poll() is None:
            job["proc"].kill()
            job["proc"].wait()
    finally:
        job["tmp"].cleanup()


def prefetch(adata, threshold: float = 0.5, max_cells: int = 0) -> str | None:
    """Start Cancer-Finder in the background on the CURRENT counts so a later call_cancerfinder with
    the same input and parameters only waits for it. Its ~10 s (subprocess imports + torch forward)
    then overlaps clustering/annotation instead of adding to them. Returns a handle for
    :func:`discard`, or None when nothing was started. Never raises."""
    cfg, _ = _config()
    if cfg is None:
        return None
    try:
        fp = _fingerprint(adata, threshold, max_cells)
        if fp not in _PENDING:
            _PENDING[fp] = _launch(adata, threshold, max_cells, cfg)
        return fp
    except Exception:  # noqa: BLE001 - an optimisation; the synchronous path still works
        return None


def discard(handle: str | None) -> None:
    """Kill + clean up a prefetched run nobody joined (e.g. the malignant stage was skipped)."""
    job = _PENDING.pop(handle, None) if handle else None
    if job is not None:
        _discard(job)


def call_cancerfinder(adata, threshold: float = 0.5, env_python: str | None = None,
                      repo: str | None = None, ckpt: str | None = None, max_cells: int = 0) -> dict:
    """Per-cell Cancer-Finder malignant probability, run in an isolated torch env.

    Writes ``obs['cancerfinder_prob']`` (0-1) + ``obs['cancerfinder_malignant']`` (prob > threshold).
    ``max_cells`` uniformly subsamples for tractability (0 = all). Returns a summary; on ANY failure
    (env/repo/checkpoint missing, subprocess error) returns ``{'status': 'skipped: ...'}`` and never
    raises, so the pipeline degrades gracefully when Cancer-Finder is not configured. Joins a matching
    :func:`prefetch` run when there is one.
    """
    import numpy as np

    def _skip(msg: str) -> dict:
        return {"status": f"skipped: {msg}", "pct_malignant": 0.0, "threshold": threshold}

    cfg, why = _config(env_python, repo, ckpt)
    if cfg is None:
        return _skip(why)
    job = None
    if (env_python, repo, ckpt) == (None, None, None) and _PENDING:
        try:
            job = _PENDING.pop(_fingerprint(adata, threshold, max_cells), None)
        except Exception:  # noqa: BLE001
            job = None
    try:
        if job is None:
            job = _launch(adata, threshold, max_cells, cfg)
        try:
            rc = job["proc"].wait(timeout=7200)
        except Exception as exc:
            return _skip(f"subprocess error ({exc})")
        if rc != 0 or not job["out"].exists():
            tail = job["err"].read_text(errors="replace").strip().splitlines() if job["err"].exists() else []
            return _skip(f"cancerfinder subprocess failed ({tail[-1] if tail else 'no output'})")
        covered = _join_cf(adata, job["out"], threshold)
    except Exception as exc:  # noqa: BLE001 - export/launch failure degrades like a subprocess failure
        return _skip(f"subprocess error ({exc})")
    finally:
        if job is not None:
            _discard(job)

    prob = np.asarray(adata.obs["cancerfinder_prob"], dtype=float)
    valid = ~np.isnan(prob)
    pct = float((prob[valid] > threshold).mean()) if covered else 0.0
    return {"status": "ok", "pct_malignant": pct, "threshold": threshold,
            "mean_prob": float(np.nanmean(prob)) if covered else 0.0,
            "coverage": covered / max(1, adata.n_obs)}
