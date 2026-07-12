"""Orchestrate subprocess-based annotation methods (RCTD / SingleR / scANVI / panhumanpy).

Each method runs in its OWN env via a subprocess script under ``subprocesses/annotation/`` that
writes a ``{cell_id, ...}`` parquet; this module joins it onto ``adata.obs`` (keyed on cell_id,
reindexed to obs order) and returns a coverage summary. Same pattern as ``qc.apply_ovrlpy_vsi``.
Nothing here imports torch / tensorflow / R - those live behind the subprocess boundary.
"""
from __future__ import annotations

import subprocess
from pathlib import Path


def join_parquet(adata, parquet_path, columns: list[str], prefix: str) -> dict:
    """Join ``{cell_id, <columns>}`` parquet onto ``obs`` as ``{prefix}_{col}`` (reindexed)."""
    import pandas as pd

    df = pd.read_parquet(parquet_path)
    df["cell_id"] = df["cell_id"].astype(str)
    aligned = df.set_index("cell_id").reindex(adata.obs_names.astype(str))
    for col in columns:
        adata.obs[f"{prefix}_{col}"] = aligned[col].to_numpy() if col in aligned else None
    present = [c for c in columns if c in aligned]
    covered = int(aligned[present[0]].notna().sum()) if present else 0
    return {"prefix": prefix, "coverage": float(covered / max(1, adata.n_obs)), "n_covered": covered}


def _run_subprocess(env_python: str, script: str, args: list[str]) -> Path | None:
    """Run a method subprocess (its env's python + the script). Returns None (not raises) on failure."""
    cmd = [env_python, script, *args]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    except Exception:
        return None
    if proc.returncode != 0:
        return None
    # Convention: the script prints the parquet path as its last stdout line.
    out = proc.stdout.strip().splitlines()
    return Path(out[-1]) if out and Path(out[-1]).exists() else None


def run_panhumanpy(sample, env_python: str | None = None) -> "Path | None":
    """Export the section's counts, run panhumanpy in its env, return the parquet path (or None).

    Returns None when ``env_python`` is not configured (the reference-free method still needs its
    own TF env) - the caller logs 'not run'. Never imports panhumanpy into this interpreter.
    """
    import tempfile
    from pathlib import Path

    if env_python is None:
        return None
    tmp = Path(tempfile.mkdtemp())
    h5 = tmp / "counts.h5ad"
    sample.adata.write_h5ad(h5)
    script = str(Path(__file__).resolve().parents[3] / "subprocesses" / "annotation" / "run_panhumanpy.py")
    return _run_subprocess(env_python, script, [str(h5), str(tmp / "ph.parquet")])


def join_panhumanpy(adata, parquet) -> dict:
    """Join a panhumanpy parquet onto obs as ph_broad/ph_medium/ph_fine/ph_confidence."""
    return join_parquet(adata, parquet, columns=["broad", "medium", "fine", "confidence"], prefix="ph")


def run_rctd(sample, reference_path=None, ref_label_key="cell_type", env_python: str | None = None):
    """Run RCTD doublet-mode in the rctd-py env; return the parquet Path (or None to skip).

    Skips (returns None) when no reference or no env_python is configured - the reference path
    only lights up when a reference .h5ad is supplied. Never imports rctd into this interpreter.
    """
    import tempfile
    from pathlib import Path

    if reference_path is None or env_python is None:
        return None
    tmp = Path(tempfile.mkdtemp())
    h5 = tmp / "spatial.h5ad"
    sample.adata.write_h5ad(h5)
    script = str(Path(__file__).resolve().parents[3] / "subprocesses" / "annotation" / "run_rctd.py")
    return _run_subprocess(env_python, script, [str(h5), str(reference_path), ref_label_key, str(tmp / "rctd.parquet")])


def join_rctd(adata, parquet) -> dict:
    """Join an RCTD parquet onto obs as rctd_first_type/second_type/spot_class/weight/singlet_score."""
    return join_parquet(adata, parquet,
                        columns=["first_type", "second_type", "spot_class", "weight", "singlet_score"],
                        prefix="rctd")


def run_singler(sample, reference_path=None, ref_label_key="cell_type", env_python: str | None = None):
    """Run SingleR (BiocPy) in the singler env; return the parquet Path (or None to skip).

    Skips (returns None) when no reference or env_python is configured. CPU-only method; never
    imports singler into this interpreter.
    """
    import tempfile
    from pathlib import Path

    if reference_path is None or env_python is None:
        return None
    tmp = Path(tempfile.mkdtemp())
    h5 = tmp / "spatial.h5ad"
    sample.adata.write_h5ad(h5)
    script = str(Path(__file__).resolve().parents[3] / "subprocesses" / "annotation" / "run_singler.py")
    return _run_subprocess(env_python, script, [str(h5), str(reference_path), ref_label_key, str(tmp / "singler.parquet")])


def join_singler(adata, parquet) -> dict:
    """Join a SingleR parquet onto obs as singler_label / singler_delta."""
    return join_parquet(adata, parquet, columns=["label", "delta"], prefix="singler")


def run_scanvi(sample, reference_path=None, ref_label_key="cell_type", env_python: str | None = None):
    """Run scANVI (scArches surgery) in the scanvi env; return the parquet Path (or None to skip).

    Skips (returns None) when no reference or env_python is configured. GPU-preferred; never imports
    scvi-tools/torch into this interpreter.
    """
    import tempfile
    from pathlib import Path

    if reference_path is None or env_python is None:
        return None
    tmp = Path(tempfile.mkdtemp())
    h5 = tmp / "spatial.h5ad"
    sample.adata.write_h5ad(h5)
    script = str(Path(__file__).resolve().parents[3] / "subprocesses" / "annotation" / "run_scanvi.py")
    return _run_subprocess(env_python, script, [str(h5), str(reference_path), ref_label_key, str(tmp / "scanvi.parquet")])


def join_scanvi(adata, parquet) -> dict:
    """Join a scANVI parquet onto obs as scanvi_label / scanvi_confidence / scanvi_entropy."""
    return join_parquet(adata, parquet, columns=["label", "confidence", "entropy"], prefix="scanvi")
