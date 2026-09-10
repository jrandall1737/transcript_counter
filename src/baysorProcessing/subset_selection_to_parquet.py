#!/usr/bin/env python3
"""Subset a Xenium transcripts.parquet to a Xenium Explorer polygon selection,
keeping the native Xenium parquet schema (all columns, original string cell_id,
no quality filtering).

Unlike subset_selection_for_baysor.py, this does not reshape the data for
Baysor - it just crops transcripts.parquet down to the selection polygon.

Example:
    python subset_selection_to_parquet.py \
        --region-dir "20250827__205437__CSU_Bouchet_v1_2025-08-27/output-XETG00230__0063687__Region-1A__20250827__205535" \
        --selection EmmaRegion1LeftSelection.csv \
        --out EmmaRegion1Left_subset.parquet
"""

import argparse
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from subset_selection_for_baysor import load_polygon


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--region-dir", required=True,
                   help="Xenium output folder that contains transcripts.parquet")
    p.add_argument("--selection", required=True,
                   help="Xenium Explorer selection-coordinates CSV (polygon vertices, microns)")
    p.add_argument("--out", required=True, help="Output parquet path")
    return p.parse_args()


def main():
    args = parse_args()
    region_dir = Path(args.region_dir)
    transcripts = region_dir / "transcripts.parquet"
    if not transcripts.is_file():
        raise SystemExit(f"Not found: {transcripts}")

    poly = load_polygon(Path(args.selection))

    table = pq.read_table(transcripts)
    print(f"Loaded {table.num_rows:,} transcripts")

    points = np.column_stack([
        table.column("x_location").to_numpy(),
        table.column("y_location").to_numpy(),
    ])
    inside = poly.contains_points(points)
    table = table.filter(pa.array(inside))
    print(f"{table.num_rows:,} inside the selection polygon")
    if table.num_rows == 0:
        raise SystemExit("No transcripts fall inside the polygon - check that the selection "
                          "CSV matches this region.")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, out)
    print(f"Wrote {table.num_rows:,} rows -> {out}  ({out.stat().st_size/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
