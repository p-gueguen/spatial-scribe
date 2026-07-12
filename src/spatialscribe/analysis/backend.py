"""GPU/CPU compute switch - the single source of GPU-native routing for the whole engine.

What it does
------------
`get_backend()` returns ONE object whose methods run on the GPU via ``rapids_singlecell`` (rsc)
when a CUDA device is present (the cluster L40S / Blackwell) and otherwise fall back to plain
``scanpy`` / ``squidpy`` / ``sklearn``. EVERY heavy ``analysis`` step calls through here, so
"is this GPU or CPU?" is a single switch and the CPU path stays reproducible for judges/CI.

Covered operations (GPU rsc  <->  CPU scanpy/squidpy):
  preprocessing : normalize_log, hvg, scale, pca, neighbors, calculate_qc_metrics
  clustering    : leiden, umap
  annotation    : score_genes, rank_genes_groups
  integration   : harmony
  spatial       : spatial_neighbors, nhood_enrichment, spatial_autocorr, co_occurrence

Operations with NO GPU equivalent stay on CPU by design and are NOT routed here: CellTypist,
TACCO, sklearn NMF / permutation tests / CV classifiers, and the isolated-env subprocess callers
(infercnvpy, Cancer-Finder). See the internal docs (the engine audit).

rapids-singlecell 0.14.x notes (the installed gpu-env build):
  * ``rsc.gr`` ships ONLY spatial_autocorr / co_occurrence / ligrec - spatial_neighbors and
    nhood_enrichment are NOT GPU-ported, so they run on squidpy on BOTH legs.
  * ``tl.rank_genes_groups`` HAS a GPU Wilcoxon in 0.14 (custom kernel, Dask-aware).
  * Out-of-core / multi-GPU via Dask: ``SPATIALSCRIBE_GPU_DASK=1`` routes the chunkable steps
    (pca covariance_eigh, normalize, hvg, score_genes, rank_genes_groups) through dask-cuda for
    datasets that exceed a single GPU's VRAM. See ``dask_enabled``. Harmony + Lanczos/Randomized
    PCA are single-GPU only.

How to use it
-------------
>>> be = get_backend()                         # "gpu" if rsc + a device, else "cpu"
>>> be.normalize_log(adata); be.pca(adata); be.neighbors(adata)
>>> be.leiden(adata, resolution=0.5); be.umap(adata)
>>> be.score_genes(adata, ["CD3D", "CD3E"], score_name="score_Tcell")
>>> be.spatial_neighbors(adata); be.nhood_enrichment(adata, "cell_type")

Depends on
----------
scanpy + squidpy (always). rapids_singlecell + cupy + cuml + cugraph (optional, GPU path).
"""

from __future__ import annotations

import functools
import os


@functools.lru_cache(maxsize=1)
def gpu_available() -> bool:
    """True iff rapids_singlecell + a usable CUDA device are importable.

    Honors ``SPATIALSCRIBE_FORCE_CPU=1`` to force the CPU path (used in CI / the
    judge-reproducibility gate). Cached - call ``gpu_available.cache_clear()`` +
    ``get_backend.cache_clear()`` to re-detect after toggling the env (the GPU/CPU bench does this).
    """
    if os.environ.get("SPATIALSCRIBE_FORCE_CPU") == "1":
        return False
    try:
        import cupy as cp  # noqa: F401
        import rapids_singlecell  # noqa: F401

        return cp.cuda.runtime.getDeviceCount() > 0
    except Exception:
        return False


@functools.lru_cache(maxsize=1)
def dask_enabled() -> bool:
    """True iff out-of-core / multi-GPU Dask routing is requested (``SPATIALSCRIBE_GPU_DASK=1``)
    AND a GPU + dask-cuda are available.

    When on, the GPU backend chunks the largest steps so a dataset can exceed a single GPU's VRAM
    (RMM host spill / multi-GPU). Off by default: single-GPU is faster and simpler below the ~10^7
    cells where onboard VRAM becomes the bottleneck. Fully realizing multi-GPU still needs a
    ``dask_cuda.LocalCUDACluster`` at the process entry point (see the internal docs); this flag gates the
    per-op chunked paths.
    """
    if os.environ.get("SPATIALSCRIBE_GPU_DASK") != "1" or not gpu_available():
        return False
    try:
        import dask_cuda  # noqa: F401

        return True
    except Exception:
        return False


class _CpuBackend:
    name = "cpu"

    # ----------------------------- preprocessing ----------------------------- #
    def normalize_log(self, adata, target_sum: float | None = None):
        """log1p total-count normalization. ``target_sum=None`` scales each cell to the MEDIAN
        transcript count - this is deliberate and is 10x's recommended Xenium normalization
        (``Seurat::NormalizeData(scale.factor = median(nCount_Xenium))``). Do NOT hardcode a fixed
        1e4 target: on sparse/targeted panels the median scale-factor beat 1e4-CPM on every
        clustering-recovery metric in the normalization benchmark (an internal benchmark harness log1p-
        median was #1 on both the 5K breast and the 405 renal sections), because 1e4 over-inflates
        low-count cells before the log. On deep WTA the two are close, so median stays a safe default."""
        import scanpy as sc

        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)

    def hvg(self, adata, n_top_genes: int | None = 2000, flavor: str = "seurat"):
        import scanpy as sc

        # For small targeted panels, HVG selection is optional; guard on gene count.
        if n_top_genes is None or adata.n_vars <= n_top_genes:
            return
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor=flavor)

    def filter_genes(self, adata, min_cells: int = 1):
        import scanpy as sc

        sc.pp.filter_genes(adata, min_cells=min_cells)

    def scale(self, adata, max_value: float | None = 10.0):
        import scanpy as sc

        sc.pp.scale(adata, max_value=max_value)

    def pca(self, adata, n_comps: int = 50):
        import scanpy as sc

        sc.pp.pca(adata, n_comps=min(n_comps, adata.n_vars - 1))

    def neighbors(self, adata, n_neighbors: int = 15):
        import scanpy as sc

        sc.pp.neighbors(adata, n_neighbors=n_neighbors)

    def calculate_qc_metrics(self, adata, qc_vars=("control",)):
        """Per-cell/per-gene QC columns (total_counts, n_genes_by_counts, pct_counts_<qc_var>)."""
        import scanpy as sc

        sc.pp.calculate_qc_metrics(
            adata, qc_vars=list(qc_vars), percent_top=None, inplace=True, log1p=False
        )

    # ------------------------------- clustering ------------------------------ #
    def leiden(self, adata, resolution: float = 0.5, key_added: str = "leiden"):
        import scanpy as sc

        # flavor="igraph" is the modern, faster CPU path.
        sc.tl.leiden(
            adata, resolution=resolution, key_added=key_added,
            flavor="igraph", n_iterations=2, directed=False,
        )

    def umap(self, adata):
        import scanpy as sc

        sc.tl.umap(adata)

    # ------------------------------ annotation ------------------------------- #
    def score_genes(self, adata, gene_list, score_name: str = "score", ctrl_size: int = 50, **kwargs):
        """Per-cell signature score (obs[score_name]). ctrl_size=50 keeps the control pool small -
        the default (large) pool dilutes to zero-count noise on <500-gene Xenium cells. Extra kwargs
        (e.g. n_bins, random_state) pass straight through to scanpy/rsc score_genes."""
        import scanpy as sc

        sc.tl.score_genes(adata, list(gene_list), score_name=score_name, ctrl_size=ctrl_size, **kwargs)

    def rank_genes_groups(self, adata, groupby: str, method: str = "wilcoxon",
                          key_added: str = "rank_genes_groups"):
        import scanpy as sc

        sc.tl.rank_genes_groups(adata, groupby=groupby, method=method,
                                key_added=key_added, use_raw=False)

    # ------------------------------ integration ------------------------------ #
    def harmony(self, adata, key: str, basis: str = "X_pca", adjusted_basis: str = "X_pca_harmony"):
        import scanpy.external as sce

        sce.pp.harmony_integrate(adata, key, basis=basis, adjusted_basis=adjusted_basis)

    # -------------------------------- spatial -------------------------------- #
    # spatial_neighbors + nhood_enrichment have no GPU port (rsc 0.14) -> squidpy on both legs.
    def spatial_neighbors(self, adata, **kwargs):
        import squidpy as sq

        sq.gr.spatial_neighbors(adata, **kwargs)

    def nhood_enrichment(self, adata, cluster_key: str, **kwargs):
        import squidpy as sq

        sq.gr.nhood_enrichment(adata, cluster_key=cluster_key, **kwargs)

    def spatial_autocorr(self, adata, mode: str = "moran", **kwargs):
        import squidpy as sq

        sq.gr.spatial_autocorr(adata, mode=mode, **kwargs)

    def co_occurrence(self, adata, cluster_key: str, **kwargs):
        import squidpy as sq

        sq.gr.co_occurrence(adata, cluster_key=cluster_key, **kwargs)


class _GpuBackend:
    name = "gpu"

    @staticmethod
    def _on_gpu(adata) -> bool:
        """True iff ``adata.X`` already lives on the device (cupy / cupyx sparse)."""
        return "cupy" in type(adata.X).__module__

    def _to_gpu(self, adata):
        # Idempotent: move the section onto the device only if it is not already there. Every heavy
        # step calls this, so the GPU path works even when normalize_log was skipped (e.g. the app
        # re-clusters a section that was loaded pre-normalized) - previously that left the data on
        # host and rapids raised "input is not a CuPy ndarray" in pca/neighbors.
        import rapids_singlecell as rsc

        if self._on_gpu(adata):
            return
        # rapids/cupy sparse kernels require FLOAT input; some exports (e.g. Atera WTA) ship raw
        # int32/int64 counts, which make the GPU qc/normalize kernels raise. Cast to float32 on host
        # before the device transfer (scanpy tolerates int, so this only bites the GPU leg).
        dt = getattr(getattr(adata, "X", None), "dtype", None)
        if dt is not None and dt.kind in "iub":
            adata.X = adata.X.astype("float32")
        rsc.get.anndata_to_GPU(adata)

    def _to_cpu(self, adata):
        """Bring the section back to host (idempotent, best-effort) so numpy/pandas downstream
        consumers (rank_genes tables, plotting, pandas obs math) get host arrays."""
        try:
            import rapids_singlecell as rsc

            rsc.get.anndata_to_CPU(adata)
        except Exception:
            pass

    # ----------------------------- preprocessing ----------------------------- #
    def normalize_log(self, adata, target_sum: float | None = None):
        # target_sum=None => median-transcript scale factor (10x's recommended Xenium normalization
        # + the benchmark-best on targeted panels). See _CpuBackend.normalize_log; keep the default.
        import rapids_singlecell as rsc

        self._to_gpu(adata)
        rsc.pp.normalize_total(adata, target_sum=target_sum)
        rsc.pp.log1p(adata)

    def hvg(self, adata, n_top_genes: int | None = 2000, flavor: str = "seurat"):
        import rapids_singlecell as rsc

        if n_top_genes is None or adata.n_vars <= n_top_genes:
            return
        self._to_gpu(adata)
        rsc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor=flavor)

    def filter_genes(self, adata, min_cells: int = 1):
        import rapids_singlecell as rsc

        self._to_gpu(adata)
        rsc.pp.filter_genes(adata, min_cells=min_cells)

    def scale(self, adata, max_value: float | None = 10.0):
        import rapids_singlecell as rsc

        self._to_gpu(adata)
        rsc.pp.scale(adata, max_value=max_value)

    def pca(self, adata, n_comps: int = 50):
        # rsc 0.14 sparse PCA: covariance_eigh (Dask-chunkable) for WTA-scale, else the fast default.
        # Lanczos/Randomized solvers exist but do NOT accept Dask, so out-of-core uses covariance_eigh.
        import rapids_singlecell as rsc

        self._to_gpu(adata)
        n = min(n_comps, adata.n_vars - 1)
        if dask_enabled():
            rsc.pp.pca(adata, n_comps=n, svd_solver="covariance_eigh")
        else:
            rsc.pp.pca(adata, n_comps=n)

    def neighbors(self, adata, n_neighbors: int = 15):
        import rapids_singlecell as rsc

        self._to_gpu(adata)
        rsc.pp.neighbors(adata, n_neighbors=n_neighbors)

    def calculate_qc_metrics(self, adata, qc_vars=("control",)):
        import rapids_singlecell as rsc

        self._to_gpu(adata)
        # log1p=False to match the CPU path (qc.py reads the raw metrics, not the log1p'd ones).
        rsc.pp.calculate_qc_metrics(adata, qc_vars=list(qc_vars), log1p=False)
        self._to_cpu(adata)   # obs metrics land host; QC math downstream is pandas/numpy

    # ------------------------------- clustering ------------------------------ #
    def leiden(self, adata, resolution: float = 0.5, key_added: str = "leiden"):
        # Leiden stays on CPU (scanpy igraph) on purpose. cugraph Leiden is verified crash-free on
        # this build but its resolution behaves differently from igraph's - on the demo's precomputed
        # neighbor graph it over-splits into singleton clusters that then break rank_genes_groups.
        # Keeping Leiden on CPU makes clustering identical to the reproducibility (CPU) path; the GPU
        # wins we want are neighbors/UMAP. Bring the section back to host first so igraph gets a
        # numpy/scipy graph.
        import scanpy as sc

        self._to_cpu(adata)
        sc.tl.leiden(adata, resolution=resolution, key_added=key_added,
                     flavor="igraph", n_iterations=2, directed=False)

    def umap(self, adata):
        # GPU UMAP via cuml (verified stable; ~560x faster than CPU on 100k cells). Last GPU step of
        # the cluster pipeline, so afterwards bring the section back to host: every downstream consumer
        # (rank_genes_groups, plotting, annotation) expects a numpy-backed AnnData. Falls back to CPU
        # scanpy UMAP if the GPU path raises.
        import rapids_singlecell as rsc

        try:
            self._to_gpu(adata)
            rsc.tl.umap(adata)
        except Exception:
            import scanpy as sc

            self._to_cpu(adata)
            sc.tl.umap(adata)
            return
        self._to_cpu(adata)

    # ------------------------------ annotation ------------------------------- #
    def score_genes(self, adata, gene_list, score_name: str = "score", ctrl_size: int = 50, **kwargs):
        import rapids_singlecell as rsc

        self._to_gpu(adata)
        rsc.tl.score_genes(adata, list(gene_list), score_name=score_name, ctrl_size=ctrl_size, **kwargs)
        self._to_cpu(adata)   # obs[score_name] + X land host for the surrounding CPU marker logic

    def rank_genes_groups(self, adata, groupby: str, method: str = "wilcoxon",
                          key_added: str = "rank_genes_groups"):
        # rsc 0.14 has a GPU Wilcoxon kernel. Falls back to scanpy on any parity/kernel issue so the
        # uns['rank_genes_groups'] structure downstream consumers read (names/scores recarrays) is safe.
        import rapids_singlecell as rsc

        try:
            self._to_gpu(adata)
            rsc.tl.rank_genes_groups(adata, groupby=groupby, method=method, key_added=key_added,
                                     use_raw=False)
        except Exception:
            import scanpy as sc

            self._to_cpu(adata)
            sc.tl.rank_genes_groups(adata, groupby=groupby, method=method, key_added=key_added,
                                    use_raw=False)
            return
        self._to_cpu(adata)

    # ------------------------------ integration ------------------------------ #
    def harmony(self, adata, key: str, basis: str = "X_pca", adjusted_basis: str = "X_pca_harmony"):
        # GPU Harmony (custom CuPy kernels, ~250-350x). Single-GPU only (no Dask). Operates on the PCA
        # embedding in obsm - no X transfer needed. Falls back to scanpy.external on failure.
        import rapids_singlecell as rsc

        try:
            rsc.pp.harmony_integrate(adata, key, basis=basis, adjusted_basis=adjusted_basis)
        except Exception:
            import scanpy.external as sce

            sce.pp.harmony_integrate(adata, key, basis=basis, adjusted_basis=adjusted_basis)

    # -------------------------------- spatial -------------------------------- #
    def spatial_neighbors(self, adata, **kwargs):
        # No GPU port in rsc 0.14 - operates on obsm['spatial'] + graph, not X, so squidpy is fine.
        import squidpy as sq

        self._to_cpu(adata)
        sq.gr.spatial_neighbors(adata, **kwargs)

    def nhood_enrichment(self, adata, cluster_key: str, **kwargs):
        import squidpy as sq

        self._to_cpu(adata)
        sq.gr.nhood_enrichment(adata, cluster_key=cluster_key, **kwargs)

    def spatial_autocorr(self, adata, mode: str = "moran", **kwargs):
        # GPU Moran's I / Geary's C (rsc.gr.spatial_autocorr) over the prebuilt spatial graph.
        import rapids_singlecell as rsc

        try:
            self._to_gpu(adata)
            rsc.gr.spatial_autocorr(adata, mode=mode, **kwargs)
            self._to_cpu(adata)
        except Exception:
            import squidpy as sq

            self._to_cpu(adata)
            sq.gr.spatial_autocorr(adata, mode=mode, **kwargs)

    def co_occurrence(self, adata, cluster_key: str, **kwargs):
        import rapids_singlecell as rsc

        try:
            self._to_gpu(adata)
            rsc.gr.co_occurrence(adata, cluster_key=cluster_key, **kwargs)
            self._to_cpu(adata)
        except Exception:
            import squidpy as sq

            self._to_cpu(adata)
            sq.gr.co_occurrence(adata, cluster_key=cluster_key, **kwargs)


@functools.lru_cache(maxsize=1)
def get_backend():
    """Return the active backend singleton ("gpu" if available, else "cpu")."""
    return _GpuBackend() if gpu_available() else _CpuBackend()
