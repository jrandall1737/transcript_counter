#!/usr/bin/env python3
"""Rerunnable end-to-end pipeline: Xenium polygon subset -> legacy Baysor
preprocessor -> Baysor -> per-cell counts for both segmentations.

Exercises the *legacy* 10x_to_baysor_dataPreProcessor.py path (as opposed to
subset_selection_for_baysor.py) so its Baysor-ready CSV and the resulting
Baysor segmentation can be compared against 10x's own segmentation over the
same polygon selection.

Prior segmentation: 10x_to_baysor_dataPreProcessor.py only zeroes cell_id when
it equals integer -1, an older Xenium convention. This dataset's cell_id is a
string ("UNASSIGNED" / "abcd1234-1"), so that check never fires. By default this
script adds a follow-up step that rewrites cell_id to an integer prior
(0 = unassigned / cell-free, 1..N = one id per 10x cell) before handing it to
Baysor. Pass --raw-prior to skip that and feed the raw strings straight through
(Baysor then treats "UNASSIGNED" as one ordinary cell and logs "No unassigned
molecules found in the prior segmentation" - see README-baysor.md).

`run_pipeline()` is importable; run_pipeline_batch.py loops it over an ROI
workbook.

Example:
    python run_pipeline_comparison.py \
        --region-dir "20250827__205437__CSU_Bouchet_v1_2025-08-27/output-XETG00230__0063687__Region-1A__20250827__205535" \
        --selection 20260903Run/Input/EmmaRegion1LeftSelection.csv \
        --out-dir 20260903Run/Output/pipeline_test
"""

import argparse
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

BAYSOR_PROCESSING_DIR = Path(__file__).resolve().parent
TRANSCRIPT_COUNTING_DIR = BAYSOR_PROCESSING_DIR.parent / "transcriptCounting"

DEFAULT_SELECTION = "20260903Run/Input/EmmaRegion1LeftSelection.csv"
DEFAULT_BAYSOR_CMD = r"C:\Users\jrand\.julia\bin\baysor.cmd"

# cell_id values (case-insensitive) that mean "not assigned to a cell".
_UNASSIGNED = {"unassigned", "", "-1", "nan", "none", "na", "<na>", "0"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--region-dir", required=True,
                   help="Xenium output folder that contains transcripts.parquet")
    p.add_argument("--selection", default=DEFAULT_SELECTION,
                   help=f"Xenium Explorer selection-coordinates CSV (default: {DEFAULT_SELECTION})")
    p.add_argument("--out-dir", required=True,
                   help="Destination folder for every intermediate and final file")
    p.add_argument("--baysor-cmd", default=DEFAULT_BAYSOR_CMD,
                   help=f"Path to baysor.cmd (default: {DEFAULT_BAYSOR_CMD})")
    p.add_argument("--raw-prior", action="store_true",
                   help="feed the legacy CSV's string cell_id straight to Baysor "
                        "(default: rewrite it to an integer 0=unassigned prior)")
    p.add_argument("--skip-existing", action="store_true",
                   help="if baysor_cell_counts.csv already exists in --out-dir, "
                        "return its stats without re-running")
    return p.parse_args()


def step(title: str):
    print(f"\n=== {title} ===")


def run(cmd, **kwargs):
    print("+ " + " ".join(str(c) for c in cmd))
    subprocess.run(cmd, check=True, **kwargs)


def _rewrite_prior_to_int(csv_path: Path) -> None:
    """Replace the string cell_id column with an integer prior for Baysor:
    unassigned -> 0, every other id -> a stable positive integer."""
    df = pd.read_csv(csv_path)
    raw = df["cell_id"].astype("string").str.strip()
    unassigned = raw.isna() | raw.str.lower().isin(_UNASSIGNED)
    codes, _ = pd.factorize(raw.where(~unassigned))  # NaN -> -1, others -> 0,1,2,...
    df["cell_id"] = np.where(unassigned, 0, codes + 1).astype("int64")
    df.to_csv(csv_path, index=False)
    n_cells = int(df.loc[df["cell_id"] > 0, "cell_id"].nunique())
    print(f"  prior segmentation -> int: {n_cells:,} cells, "
          f"{int((df['cell_id'] == 0).sum()):,} unassigned transcripts")


def _matrix_stats(counts_csv: Path):
    """(n_cells, n_transcripts_assigned) from a cell x gene count matrix CSV."""
    if not counts_csv.is_file():
        return (None, None)
    m = pd.read_csv(counts_csv, index_col=0)
    return (int(m.shape[0]), int(m.to_numpy().sum()))


def _gather_stats(xenium_subset: Path, legacy_csv: Path, xenium_counts: Path,
                  baysor_counts: Path, baysor_run_log: Path, seconds: float) -> dict:
    n_subset = pq.ParquetFile(xenium_subset).metadata.num_rows if xenium_subset.is_file() else None
    n_legacy = None
    if legacy_csv.is_file():
        with open(legacy_csv, encoding="utf-8", errors="replace") as fh:
            n_legacy = max(sum(1 for _ in fh) - 1, 0)  # minus header
    xc, xa = _matrix_stats(xenium_counts)
    bc, ba = _matrix_stats(baysor_counts)
    warn = (baysor_run_log.is_file()
            and "No unassigned molecules found" in baysor_run_log.read_text("utf-8", "replace"))
    return {
        "n_subset_transcripts": n_subset,
        "n_legacy_transcripts": n_legacy,
        "xenium_cells": xc, "xenium_assigned": xa,
        "baysor_cells": bc, "baysor_assigned": ba,
        "baysor_prior_warning": bool(warn),
        "seconds": round(seconds, 1),
    }


def run_pipeline(region_dir, selection, out_dir, *, baysor_cmd=DEFAULT_BAYSOR_CMD,
                 fix_prior=True, skip_existing=False) -> dict:
    """Run the 5-step subset -> Baysor -> counts pipeline into ``out_dir``.

    Returns a stats dict (see ``_gather_stats``). Raises on any step failure.
    """
    region_dir = Path(region_dir)
    selection = Path(selection)
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    xenium_subset = out_dir / "xenium_subset.parquet"
    legacy_csv = out_dir / "legacy_baysor_transcripts.csv"
    baysor_prefix = "baysor-output"
    baysor_segmentation = out_dir / f"{baysor_prefix}_segmentation.csv"
    baysor_run_log = out_dir / "baysor_run.log"
    xenium_counts = out_dir / "xenium_subset_cell_counts.csv"
    baysor_counts = out_dir / "baysor_cell_counts.csv"

    if skip_existing and baysor_counts.is_file() and xenium_counts.is_file():
        print(f"skip-existing: {baysor_counts.name} already present in {out_dir}")
        return _gather_stats(xenium_subset, legacy_csv, xenium_counts, baysor_counts, baysor_run_log, 0.0)

    t0 = time.time()

    step("Step 1: polygon subset -> native Xenium parquet")
    run([
        sys.executable, str(BAYSOR_PROCESSING_DIR / "subset_selection_to_parquet.py"),
        "--region-dir", str(region_dir),
        "--selection", str(selection),
        "--out", str(xenium_subset),
    ])

    step("Step 2: legacy preprocessor -> Baysor-ready CSV")
    run([
        sys.executable, str(BAYSOR_PROCESSING_DIR / "10x_to_baysor_dataPreProcessor.py"),
        "-transcript", str(xenium_subset)
    ], cwd=out_dir)
    generated_path = next(out_dir.glob("*_filtered_transcripts.csv"), None)
    if generated_path is None:
        raise RuntimeError(f"legacy preprocessor produced no *_filtered_transcripts.csv in {out_dir}")
    generated_path.replace(legacy_csv)
    print(f"{generated_path.name} -> {legacy_csv.name}")

    if fix_prior:
        step("Step 2b: rewrite cell_id -> integer prior (0 = unassigned)")
        _rewrite_prior_to_int(legacy_csv)

    step("Step 3: Baysor")
    baysor_cmd_str = (
        f'"{baysor_cmd}" run '
        f"-x x_location -y y_location -z z_location -g feature_name "
        f"--min-molecules-per-cell 1 -p --prior-segmentation-confidence 0.5 "
        f"-o ./{baysor_prefix} "
        f'"{legacy_csv.name}" :cell_id'
    )
    print("+ " + baysor_cmd_str)
    with open(baysor_run_log, "w", encoding="utf-8") as log:
        subprocess.run(baysor_cmd_str, shell=True, check=True, cwd=out_dir,
                       stdout=log, stderr=subprocess.STDOUT)
    print(f"Baysor log -> {baysor_run_log}")
    if not baysor_segmentation.is_file():
        raise RuntimeError(f"Baysor produced no {baysor_segmentation.name} (see {baysor_run_log})")

    step("Step 4: per-cell counts, Xenium (10x) segmentation")
    run([
        sys.executable, str(TRANSCRIPT_COUNTING_DIR / "count_transcripts_per_cell.py"),
        "--transcripts", str(xenium_subset),
        "--out", str(xenium_counts),
        "--summary",
    ])

    step("Step 5: per-cell counts, Baysor segmentation")
    run([
        sys.executable, str(TRANSCRIPT_COUNTING_DIR / "count_transcripts_per_cell.py"),
        "--transcripts", str(baysor_segmentation),
        "--out", str(baysor_counts),
        "--summary",
    ])

    stats = _gather_stats(xenium_subset, legacy_csv, xenium_counts, baysor_counts,
                          baysor_run_log, time.time() - t0)

    step("Done")
    print(f"Xenium (10x) counts : {xenium_counts}")
    print(f"                      {xenium_counts.with_suffix('.cell_summary.csv')}")
    print(f"Baysor counts       : {baysor_counts}")
    print(f"                      {baysor_counts.with_suffix('.cell_summary.csv')}")
    print(f"Baysor run log      : {baysor_run_log}")
    print(f"Stats               : {stats}")
    return stats


def main():
    args = parse_args()
    run_pipeline(
        args.region_dir, args.selection, args.out_dir,
        baysor_cmd=args.baysor_cmd,
        fix_prior=not args.raw_prior,
        skip_existing=args.skip_existing,
    )


if __name__ == "__main__":
    main()
