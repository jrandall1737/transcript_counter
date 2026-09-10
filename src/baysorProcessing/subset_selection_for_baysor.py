#!/usr/bin/env python3
"""Subset a Xenium transcripts.parquet to a Xenium Explorer polygon selection and
write a Baysor-ready CSV (with an integer prior-segmentation column).

Xenium Explorer's "export selection coordinates" CSV stores the polygon vertices
in microns, in the same coordinate space as transcripts.parquet x_location /
y_location, so the selection can be applied directly as a point-in-polygon test.

Example:
    python subset_selection_for_baysor.py \
        --region-dir "20250827__205437__CSU_Bouchet_v1_2025-08-27/output-XETG00230__0063687__Region-1A__20250827__205535" \
        --selection EmmaRegion1LeftSelection.csv \
        --out EmmaRegion1Left_baysor_transcripts.csv
"""

import argparse
from pyexpat.errors import codes
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from matplotlib.path import Path as MplPath

# Columns pulled from transcripts.parquet.
COLUMNS = [
    "transcript_id",
    "cell_id",
    "overlaps_nucleus",
    "feature_name",
    "x_location",
    "y_location",
    "z_location",
    "qv",
    "nucleus_distance",
    "is_gene",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--region-dir", required=True,
                   help="Xenium output folder that contains transcripts.parquet")
    p.add_argument("--selection", required=True,
                   help="Xenium Explorer selection-coordinates CSV (polygon vertices, microns)")
    p.add_argument("--out", required=True, help="Output CSV path for Baysor")
    p.add_argument("--min-qv", type=float, default=20.0,
                   help="Minimum decoding Q-score to keep (default: 20.0, matches 10x)")
    return p.parse_args()


def load_polygon(selection_csv: Path) -> MplPath:
    sel = pd.read_csv(selection_csv, comment="#")
    missing = {"X", "Y"} - set(sel.columns)
    if missing:
        sys.exit(f"Selection CSV is missing column(s): {missing}. Columns found: {list(sel.columns)}")
    names = sel["Selection"].unique() if "Selection" in sel.columns else ["(unnamed)"]
    if len(names) > 1:
        sys.exit(f"Selection CSV contains multiple selections {list(names)}; keep only one.")
    verts = sel[["X", "Y"]].to_numpy(dtype=float)
    print(f"Polygon '{names[0]}': {len(verts)} vertices, "
          f"x {verts[:,0].min():.1f}-{verts[:,0].max():.1f}, "
          f"y {verts[:,1].min():.1f}-{verts[:,1].max():.1f} um")
    return MplPath(verts)


def main():
    args = parse_args()
    region_dir = Path(args.region_dir)
    transcripts = region_dir / "transcripts.parquet"
    if not transcripts.is_file():
        sys.exit(f"Not found: {transcripts}")

    poly = load_polygon(Path(args.selection))

    cols = [c for c in COLUMNS if c in pq.ParquetFile(transcripts).schema_arrow.names]
    df = pq.read_table(transcripts, columns=cols).to_pandas()
    print(f"Loaded {len(df):,} transcripts")

    # Quality + real-gene filtering (is_gene==False covers every negative control /
    # blank / unassigned / deprecated codeword).
    keep = df["qv"] >= args.min_qv
    if "is_gene" in df.columns:
        keep &= df["is_gene"].astype(bool)
    else:
        fn = df["feature_name"].astype(str)
        for pref in ("NegControlProbe_", "NegControlCodeword_", "antisense_",
                     "BLANK_", "UnassignedCodeword_", "DeprecatedCodeword_"):
            keep &= ~fn.str.startswith(pref)
    df = df[keep]
    print(f"{len(df):,} after qv>={args.min_qv} and real-gene filter")

    # Spatial filter: point-in-polygon.
    inside = poly.contains_points(df[["x_location", "y_location"]].to_numpy())
    df = df[inside].copy()
    print(f"{len(df):,} inside the selection polygon")
    if df.empty:
        sys.exit("No transcripts fall inside the polygon - check that the selection "
                 "CSV matches this region.")

    # Prior segmentation -> integer column for Baysor (0 = unassigned / cell-free).
    # raw = df["cell_id"].astype(str)
    # unassigned = raw.isin(["UNASSIGNED", "", "-1", "nan", "None"])
    # # codes, _ = pd.factorize(raw)  # each distinct cell-id string -> 0, 1, 2, ...
    # # codes = codes + 1  # reserve 0 for "unassigned" below
    # raw[unassigned] = 0
    # df["cell_id"] = raw.astype(np.int64)
    # n_cells = int(df.loc[df["cell_id"] > 0, "cell_id"].nunique())
    # print(f"Prior segmentation: {n_cells:,} cells, "
    #       f"{int((df['cell_id'] == 0).sum()):,} unassigned transcripts")


    

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False, encoding="utf-8")
    print(f"Wrote {len(df):,} rows -> {out}  ({out.stat().st_size/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
