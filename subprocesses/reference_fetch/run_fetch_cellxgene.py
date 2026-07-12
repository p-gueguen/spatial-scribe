"""Fetch a small CELLxGENE reference for a free-text tissue - isolated, metadata-first subprocess.

cellxgene-census + tiledbsoma clash with the scanpy/spatialdata main env, so this runs OUT of
process (the same isolation rule as ovrlpy and the annotators). The main app calls it through
``reference.fetch_cellxgene_reference``; it writes a downsampled ``.h5ad`` and prints ONE JSON
summary line (last stdout line) that the caller parses.

**Metadata-first routing** (the reason this exists rather than a blind ``gget.cellxgene`` pull):
it reads only the census OBS table first - cheap, even for a 30M-cell tissue - and picks which
cells to materialise via :func:`spatialscribe.analysis.reference.select_reference_cells`
(cell-type stratified + recent-chemistry preferred), so a huge tissue is never densified whole and
rare cell types survive instead of being drowned by a random downsample. Only the chosen soma
joinids are then pulled with ``get_anndata``.

Run env (resolved by the caller): a python that has cellxgene-census + tiledbsoma, e.g.
    uv run --with cellxgene-census --with tiledbsoma python \
        subprocesses/reference_fetch/run_fetch_cellxgene.py --out ref.h5ad --tissue kidney

Requires network access to the CZI CELLxGENE census. Emits {"status":"error",...} (never a bare
traceback to stderr) so the caller can degrade cleanly.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

# The cell-selection routing is the single source of truth in the main package; load it by FILE
# (not `import spatialscribe...`) so we do NOT trigger the package __init__ (which pulls scanpy and
# would fail in this stripped census env). reference.py has only stdlib at module level.
_REF_PY = Path(__file__).resolve().parents[2] / "src" / "spatialscribe" / "analysis" / "reference.py"
_spec = importlib.util.spec_from_file_location("_ss_reference_routing", _REF_PY)
_ref = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_ref)                                        # type: ignore[union-attr]
select_reference_cells = _ref.select_reference_cells
resolve_census_tissue = _ref.resolve_census_tissue

_ORGANISM_LABEL = {"homo_sapiens": "Homo sapiens", "mus_musculus": "Mus musculus"}


def _distinct_tissue_general(census, organism: str) -> list[str]:
    """The census ``tissue_general`` vocabulary for an organism (from the summary table)."""
    try:
        df = census["census_info"]["summary_cell_counts"].read().concat().to_pandas()
        lab = _ORGANISM_LABEL.get(organism, "")
        sub = df[(df["organism"] == lab) & (df["category"] == "tissue_general")]
        return sorted(sub["label"].astype(str).unique().tolist())
    except Exception:
        return []


def _build_filter(tissue: str, disease: str | None, development_stage: str | None) -> str:
    parts = ["is_primary_data == True", f"tissue_general == '{tissue}'"]
    if disease:
        parts.append(f"disease == '{disease}'")
    if development_stage:
        parts.append(f"development_stage == '{development_stage}'")
    return " and ".join(parts)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--tissue", required=True)
    ap.add_argument("--species", default="homo_sapiens")
    ap.add_argument("--disease", default=None)
    ap.add_argument("--development-stage", default=None)
    ap.add_argument("--census-version", default="2025-01-30")
    ap.add_argument("--max-cells", type=int, default=20000)
    ap.add_argument("--min-cells-per-type", type=int, default=25)
    ap.add_argument("--no-prefer-recent-chemistry", action="store_true")
    ap.add_argument("--keep-spatial", action="store_true")
    a = ap.parse_args()

    try:
        import cellxgene_census

        cols = ["soma_joinid", "cell_type", "assay", "dataset_id", "development_stage", "tissue_general"]
        with cellxgene_census.open_soma(census_version=a.census_version) as census:
            # 1) resolve the free-text tissue to a real census tissue_general label.
            tissue = a.tissue.strip()
            vocab = _distinct_tissue_general(census, a.species)
            if vocab and tissue not in vocab:
                resolved = resolve_census_tissue(tissue, vocab)
                if resolved:
                    tissue = resolved
            # 2) METADATA ONLY - the whole tissue's obs, cheap even at tens of millions of cells.
            obs = cellxgene_census.get_obs(
                census, a.species,
                value_filter=_build_filter(tissue, a.disease, a.development_stage),
                column_names=cols)
            if obs is None or len(obs) == 0:
                print(json.dumps({"status": "error",
                                  "message": f"no primary cells for tissue={tissue!r} "
                                             f"(query {a.tissue!r}), species={a.species}"}))
                return 1
            # 3) ROUTE: stratify by cell type + prefer recent chemistry, on metadata only.
            picked = select_reference_cells(
                obs, target_cells=a.max_cells, min_cells_per_type=a.min_cells_per_type,
                prefer_recent_chemistry=not a.no_prefer_recent_chemistry,
                exclude_spatial=not a.keep_spatial)
            joinids = [int(x) for x in picked["soma_joinid"].tolist()]
            # 4) materialise ONLY the chosen cells.
            adata = cellxgene_census.get_anndata(
                census, organism=a.species, obs_coords=joinids,
                obs_column_names=["cell_type", "assay", "dataset_id", "development_stage", "tissue_general"],
                var_column_names=["feature_id", "feature_name"])

        # de-Ensembl: use the human/mouse symbol as var_names (dedup), keep the id in var.
        if "feature_name" in adata.var.columns:
            names = adata.var["feature_name"].astype(str)
            adata.var["ensembl_id"] = adata.var_names
            adata.var_names = names.values
            adata.var_names_make_unique()
        adata.obs["tissue"] = tissue                                 # so the tissue-consistency guard CONFIRMS
        n_assays = int(adata.obs["assay"].nunique()) if "assay" in adata.obs else 0
        top_assay = (adata.obs["assay"].value_counts().index[0]
                     if "assay" in adata.obs and adata.n_obs else None)
        adata.write_h5ad(a.out)
        print(json.dumps({
            "status": "ok", "path": a.out, "n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars),
            "label_key": "cell_type", "tissue": tissue, "query": a.tissue,
            "n_cell_types": int(adata.obs["cell_type"].nunique()) if "cell_type" in adata.obs else None,
            "n_assays": n_assays, "top_assay": str(top_assay) if top_assay is not None else None,
        }))
        return 0
    except Exception as exc:                                          # network / census failure
        print(json.dumps({"status": "error", "message": f"{type(exc).__name__}: {exc}"}))
        return 1


if __name__ == "__main__":
    sys.exit(main())
