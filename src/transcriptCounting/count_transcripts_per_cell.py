#!/usr/bin/env python3
"""Count transcripts per cell -> a cell x gene matrix CSV.

Script version of the ``calculate_transcript_counts`` logic in
``transcript_counter.ipynb``, generalised to run against either:

* a Xenium ``transcripts.parquet`` (uses the 10x ``cell_id`` assignment), or
* a Baysor ``<prefix>_segmentation.csv`` (uses Baysor's re-segmented ``cell``
  assignment) -- this is the new-segmentation case.

The source is auto-detected from the file, or forced with ``--source``.

Output: one row per cell, one column per gene, integer counts. The first
column holds the cell id, so the matrix loads directly with
``scanpy.read_csv(path, first_column_names=True)``.

Examples
--------
Baysor output for the EmmaRegion1Left selection::

    .venv/Scripts/python.exe src/transcriptCounting/count_transcripts_per_cell.py \
        --transcripts emmaOutput/EmmaRegion1Left/baysor-output_segmentation.csv \
        --out emmaOutput/EmmaRegion1Left/EmmaRegion1Left_baysor_cell_counts.csv \
        --summary

Original-style run straight off a Xenium region::

    .venv/Scripts/python.exe src/transcriptCounting/count_transcripts_per_cell.py \
        --transcripts "20250827__.../output-XETG..._Region-1A_.../transcripts.parquet" \
        --out output/region1A_cell_counts.csv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

# Feature-name prefixes that are controls, not real genes.
CONTROL_PREFIXES = (
    "NegControlProbe_", "NegControlCodeword_", "antisense_",
    "BLANK_", "BLANK-", "UnassignedCodeword_", "DeprecatedCodeword_",
)
# Values in a cell-id column that mean "not in a cell".
UNASSIGNED_VALUES = {"", "0", "-1", "nan", "none", "unassigned", "na", "<na>"}


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--transcripts", required=True, type=Path,
                   help="transcripts.parquet (Xenium) or *_segmentation.csv (Baysor)")
    p.add_argument("--out", required=True, type=Path,
                   help="output cell x gene matrix CSV")
    p.add_argument("--source", choices=("auto", "xenium", "baysor"), default="auto",
                   help="input type (default: auto-detect)")
    p.add_argument("--min-qv", type=float, default=20.0,
                   help="drop transcripts with qv below this (default: 20.0, matches 10x)")
    p.add_argument("--cell-column", default=None,
                   help="name of the cell-id column "
                        "(default: 'cell' for Baysor, 'cell_id' for Xenium)")
    p.add_argument("--gene-column", default=None,
                   help="name of the gene column "
                        "(default: 'gene' for Baysor, 'feature_name' for Xenium)")
    p.add_argument("--keep-controls", action="store_true",
                   help="keep NegControl / BLANK / codeword features (default: drop)")
    p.add_argument("--keep-noise", action="store_true",
                   help="Baysor: keep transcripts flagged is_noise / with no cell "
                        "(default: drop)")
    p.add_argument("--include-unassigned", action="store_true",
                   help="keep transcripts not assigned to any cell as a row named "
                        "'UNASSIGNED' (default: drop)")
    p.add_argument("--cells", type=Path, default=None,
                   help="optional CSV listing cell ids to keep (e.g. a Xenium Explorer "
                        "region export); '#'-comment lines are ignored")
    p.add_argument("--cell-list-column", default="Cell ID",
                   help="column in --cells holding the ids (default: 'Cell ID')")
    p.add_argument("--min-counts-per-cell", type=int, default=0,
                   help="drop cells with fewer than N transcripts total (default: 0)")
    p.add_argument("--require", action="append", default=[], metavar="GENE:MIN",
                   help="keep only cells with at least MIN counts of GENE; repeatable "
                        "(generalises the notebook's 'Th > 5' filter, e.g. --require Th:5)")
    p.add_argument("--long", action="store_true",
                   help="also write <out stem>.long.csv with (cell_id, gene, count) rows")
    p.add_argument("--summary", action="store_true",
                   help="also write <out stem>.cell_summary.csv (per-cell n_transcripts, "
                        "n_genes; plus Baysor cell_stats columns when available)")
    return p.parse_args(argv)


def detect_source(path: Path) -> str:
    if path.suffix.lower() in (".parquet", ".pq"):
        return "xenium"
    if path.suffix.lower() in (".csv", ".tsv", ".gz"):
        # Baysor segmentation files have a 'cell' and 'is_noise' column.
        head = pd.read_csv(path, nrows=0)
        return "baysor" if {"cell", "gene"}.issubset(head.columns) else "xenium"
    sys.exit(f"Cannot infer source type from {path!r}; pass --source explicitly.")


def load_transcripts(path: Path, cell_col: str, gene_col: str) -> pd.DataFrame:
    """Return a frame with columns: cell, gene, qv, [is_noise].

    Reads Parquet or CSV; the source type only decides which column names to
    look for, not the file format (a Xenium export can be either).
    """
    want = [cell_col, gene_col, "qv"]
    optional = ["is_noise"]
    if path.suffix.lower() in (".parquet", ".pq"):
        import pyarrow.parquet as pq
        available = set(pq.ParquetFile(path).schema_arrow.names)
        missing = [c for c in want if c not in available]
        if missing:
            sys.exit(f"{path.name}: missing expected column(s) {missing}. "
                     f"Found: {sorted(available)}")
        cols = want + [c for c in optional if c in available]
        df = pq.read_table(path, columns=cols).to_pandas()
    else:
        import pyarrow.csv as pv
        cols = want + optional
        try:
            df = pv.read_csv(
                path,
                convert_options=pv.ConvertOptions(include_columns=cols,
                                                  include_missing_columns=True),
            ).to_pandas()
            df = df.dropna(axis=1, how="all")  # drop absent optional columns
        except Exception:  # pragma: no cover - fallback for odd dialects
            df = pd.read_csv(path, usecols=lambda c: c in cols)
        missing = [c for c in want if c not in df.columns]
        if missing:
            sys.exit(f"{path.name}: missing expected column(s) {missing}.")
    return df.rename(columns={cell_col: "cell", gene_col: "gene"})


def main(argv=None):
    args = parse_args(argv)
    if not args.transcripts.is_file():
        sys.exit(f"Not found: {args.transcripts}")

    source = args.source if args.source != "auto" else detect_source(args.transcripts)
    cell_col = args.cell_column or ("cell" if source == "baysor" else "cell_id")
    gene_col = args.gene_column or ("gene" if source == "baysor" else "feature_name")
    print(f"Source: {source}  (cell column '{cell_col}', gene column '{gene_col}')")

    df = load_transcripts(args.transcripts, cell_col, gene_col)
    print(f"Loaded {len(df):,} transcripts")

    # Quality filter.
    n = len(df)
    df = df[df["qv"].astype(float) >= args.min_qv]
    print(f"{len(df):,} after qv >= {args.min_qv}  (dropped {n - len(df):,})")

    # Drop control features.
    if not args.keep_controls:
        gene = df["gene"].astype("string")
        is_control = gene.str.startswith(tuple(CONTROL_PREFIXES), na=False)
        df = df[~is_control]
        if is_control.any():
            print(f"{len(df):,} after dropping {int(is_control.sum()):,} control-feature transcripts")

    # Baysor noise transcripts.
    if source == "baysor" and "is_noise" in df.columns and not args.keep_noise:
        noise = df["is_noise"].astype("string").str.lower().isin({"true", "1"})
        df = df[~noise.fillna(False)]
        if noise.any():
            print(f"{len(df):,} after dropping {int(noise.sum()):,} Baysor noise transcripts")

    # Normalise the cell id and handle unassigned transcripts.
    cell = df["cell"].astype("string").str.strip()
    unassigned = cell.str.lower().isin(UNASSIGNED_VALUES) | cell.isna()
    if args.include_unassigned:
        cell = cell.mask(unassigned, "UNASSIGNED")
    else:
        df = df[~unassigned]
        cell = cell[~unassigned]
        print(f"{len(df):,} after dropping unassigned/cell-free transcripts")
    df = df.assign(cell=cell)

    if df.empty:
        sys.exit("No transcripts left after filtering.")

    # Optional: restrict to a supplied list of cell ids.
    if args.cells:
        wanted = pd.read_csv(args.cells, comment="#")
        if args.cell_list_column not in wanted.columns:
            sys.exit(f"{args.cells.name}: no column {args.cell_list_column!r}. "
                     f"Found: {list(wanted.columns)}")
        wanted_ids = set(wanted[args.cell_list_column].astype("string").str.strip())
        present = set(df["cell"].unique())
        missing = wanted_ids - present
        df = df[df["cell"].isin(wanted_ids)]
        print(f"{len(df):,} transcripts in the {len(wanted_ids):,} requested cells "
              f"({len(missing):,} requested cells had no transcripts)")

    # cell x gene count matrix.
    matrix = (df.groupby(["cell", "gene"]).size()
                .unstack(fill_value=0)
                .sort_index())
    matrix.index.name = "cell_id"
    matrix = matrix.reindex(sorted(matrix.columns), axis=1)
    print(f"Matrix: {matrix.shape[0]:,} cells x {matrix.shape[1]:,} genes")

    # Per-cell filters.
    totals = matrix.sum(axis=1)
    if args.min_counts_per_cell > 0:
        keep = totals >= args.min_counts_per_cell
        print(f"{int(keep.sum()):,} cells with >= {args.min_counts_per_cell} transcripts "
              f"(dropped {int((~keep).sum()):,})")
        matrix = matrix[keep]
    for spec in args.require:
        try:
            gene_name, min_str = spec.rsplit(":", 1)
            min_val = float(min_str)
        except ValueError:
            sys.exit(f"--require expects GENE:MIN, got {spec!r}")
        if gene_name not in matrix.columns:
            sys.exit(f"--require: gene {gene_name!r} not in the matrix "
                     f"({', '.join(matrix.columns[:8])} ...)")
        before = len(matrix)
        matrix = matrix[matrix[gene_name] >= min_val]
        print(f"{len(matrix):,} cells with {gene_name} >= {min_val:g}  (dropped {before - len(matrix):,})")

    if matrix.empty:
        sys.exit("No cells left after per-cell filtering.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(args.out)
    print(f"Wrote {matrix.shape[0]:,} x {matrix.shape[1]:,} matrix -> {args.out}")

    if args.long:
        long_path = args.out.with_suffix(".long.csv")
        (matrix.rename_axis("gene", axis=1).stack().rename("count")
               .reset_index().query("count > 0")
               .to_csv(long_path, index=False))
        print(f"Wrote long format -> {long_path}")

    if args.summary:
        summ = pd.DataFrame({
            "n_transcripts": matrix.sum(axis=1),
            "n_genes": (matrix > 0).sum(axis=1),
        })
        stats_path = _sibling_cell_stats(args.transcripts) if source == "baysor" else None
        if stats_path and stats_path.is_file():
            stats = pd.read_csv(stats_path).set_index("cell")
            stats.index.name = "cell_id"
            summ = summ.join(stats, how="left", rsuffix="_baysor")
            print(f"Joined Baysor cell stats from {stats_path.name}")
        summary_path = args.out.with_suffix(".cell_summary.csv")
        summ.sort_index().to_csv(summary_path)
        print(f"Wrote per-cell summary -> {summary_path}")


def _sibling_cell_stats(seg_path: Path) -> Path | None:
    """emmaOutput/.../baysor-output_segmentation.csv -> .../baysor-output_cell_stats.csv"""
    name = seg_path.name
    if "segmentation" in name:
        return seg_path.with_name(name.replace("segmentation", "cell_stats"))
    return None


if __name__ == "__main__":
    main()
