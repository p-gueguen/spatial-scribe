"""Per-cell Vertical Signal Integrity (VSI) via ovrlpy - runs in an ISOLATED env.

ovrlpy pins polars~=1.0 / umap-learn~=0.5 that clash with the main scanpy/spatialdata
stack, so this runs as a SEPARATE process and writes a per-cell VSI parquet that the main
app joins back onto ``adata.obs`` (see ``qc.apply_ovrlpy_vsi``). Never import ovrlpy in the
main interpreter.

Isolated env:
    uv venv /data/spatial-scribe/ovrlpy_env
    uv pip install --python .../ovrlpy_env/bin/python ovrlpy
Run:
    .../ovrlpy_env/bin/python subprocesses/ovrlpy/run_ovrlpy.py \
        --transcripts /path/transcripts.parquet --out /path/vsi.parquet

Verified constraints (docs/cell-annotation-qc.md, research pass): needs the raw molecule
table with a z-coordinate; ``unassigned="UNASSIGNED"`` (a string) on current Xenium, NOT
-1; ``fit_umap=False`` (viz only); flag low-VSI cells, do not auto-drop them.
"""

from __future__ import annotations

import argparse


def main() -> int:
    ap = argparse.ArgumentParser(description="Compute per-cell VSI with ovrlpy.")
    ap.add_argument("--transcripts", required=True, help="path to transcripts.parquet (needs x,y,z,cell_id)")
    ap.add_argument("--out", required=True, help="output per-cell VSI parquet")
    ap.add_argument("--n-workers", type=int, default=8)
    ap.add_argument("--n-components", type=int, default=20)
    ap.add_argument("--min-qv", type=int, default=20)
    ap.add_argument("--unassigned", default="UNASSIGNED", help="sentinel for unassigned transcripts")
    ap.add_argument("--vsi-threshold", type=float, default=0.7)
    args = ap.parse_args()

    import ovrlpy
    import polars as pl

    tx = ovrlpy.io.read_Xenium(args.transcripts, min_qv=args.min_qv, additional_columns=["cell_id"])
    ov = ovrlpy.Ovrlp(tx, n_components=args.n_components, n_workers=args.n_workers,
                      gene_key="gene", coordinate_keys=("x", "y", "z"), random_state=42)
    ov.analyse(min_transcripts=20, fit_umap=False)

    px = ovrlpy.cell_integrity_from_transcripts(ov, cell_id="cell_id", unassigned=args.unassigned)
    px = px.filter(pl.col("signal") > 1.5)
    per_cell = px.group_by("cell_id").agg(
        vsi=pl.col("vsi").mean(),
        frac_low_vsi=(pl.col("vsi") < args.vsi_threshold).mean(),
        vsi_n_px=pl.col("vsi").len(),
    )
    per_cell = per_cell.filter(pl.col("cell_id") != args.unassigned)
    per_cell.write_parquet(args.out)
    n_low = per_cell.filter(pl.col("vsi") < args.vsi_threshold).height
    print(f"[ovrlpy] wrote {args.out}: {per_cell.height:,} cells, "
          f"{n_low:,} ({100*n_low/max(1, per_cell.height):.1f}%) low VSI (<{args.vsi_threshold})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
