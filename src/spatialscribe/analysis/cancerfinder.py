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

Caveat (an internal benchmark Atera benchmark): Cancer-Finder is over-sensitive on single-cell Xenium (it agrees
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


def call_cancerfinder(adata, threshold: float = 0.5, env_python: str | None = None,
                      repo: str | None = None, ckpt: str | None = None, max_cells: int = 0) -> dict:
    """Per-cell Cancer-Finder malignant probability, run in an isolated torch env.

    Writes ``obs['cancerfinder_prob']`` (0-1) + ``obs['cancerfinder_malignant']`` (prob > threshold).
    ``max_cells`` uniformly subsamples for tractability (0 = all). Returns a summary; on ANY failure
    (env/repo/checkpoint missing, subprocess error) returns ``{'status': 'skipped: ...'}`` and never
    raises, so the pipeline degrades gracefully when Cancer-Finder is not configured.
    """
    import os as _os
    import subprocess
    import tempfile
    from pathlib import Path

    import numpy as np

    def _skip(msg: str) -> dict:
        return {"status": f"skipped: {msg}", "pct_malignant": 0.0, "threshold": threshold}

    env_python = env_python or _CF_PY
    repo = repo or _CF_REPO
    # Derive the conventional checkpoint under the EFFECTIVE repo when only repo was supplied
    # (so passing repo= alone works, not just the module-level default).
    ckpt = ckpt or (os.path.join(repo, "checkpoints", "sc_pretrain_article.pkl") if repo else _CF_CKPT)
    for label, envname, path in (("cancerfinder python", "CANCERFINDER_PYTHON", env_python),
                                 ("cancerfinder repo", "CANCERFINDER_REPO", repo),
                                 ("checkpoint", "CANCERFINDER_CKPT", ckpt)):
        if not path:
            return _skip(f"{label} not configured; set {envname}")
        if not Path(path).exists():
            return _skip(f"{label} not found ({path}); set {envname}")

    import anndata as ad

    tmp = Path(tempfile.mkdtemp(prefix="sscf_"))
    h5, out = tmp / "section.h5ad", tmp / "cf.parquet"
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    a_exp = ad.AnnData(X=counts.copy())
    a_exp.obs_names = adata.obs_names.astype(str)
    a_exp.var_names = adata.var_names.astype(str)
    a_exp.write_h5ad(h5)

    script = str(Path(__file__).resolve().parents[3] / "subprocesses" / "cancerfinder" / "run_cancerfinder.py")
    cmd = [env_python, script, "--h5ad", str(h5), "--out", str(out), "--repo", repo, "--ckpt", ckpt,
           "--threshold", str(threshold), "--max-cells", str(int(max_cells))]
    env = {k: v for k, v in _os.environ.items() if k != "PYTHONPATH"}   # don't leak the app's src path
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=7200, env=env)
    except Exception as exc:
        return _skip(f"subprocess error ({exc})")
    if proc.returncode != 0 or not out.exists():
        tail = (proc.stderr or "").strip().splitlines()
        return _skip(f"cancerfinder subprocess failed ({tail[-1] if tail else 'no output'})")

    covered = _join_cf(adata, out, threshold)
    prob = np.asarray(adata.obs["cancerfinder_prob"], dtype=float)
    valid = ~np.isnan(prob)
    pct = float((prob[valid] > threshold).mean()) if covered else 0.0
    return {"status": "ok", "pct_malignant": pct, "threshold": threshold,
            "mean_prob": float(np.nanmean(prob)) if covered else 0.0,
            "coverage": covered / max(1, adata.n_obs)}
