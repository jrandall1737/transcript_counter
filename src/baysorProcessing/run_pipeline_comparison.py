#!/usr/bin/env python3
"""Rerunnable end-to-end pipeline: Xenium polygon subset -> legacy Baysor
preprocessor -> Baysor -> per-cell counts for both segmentations.

Exercises the *legacy* 10x_to_baysor_dataPreProcessor.py path (as opposed to
subset_selection_for_baysor.py) so its Baysor-ready CSV and the resulting
Baysor segmentation can be compared against 10x's own segmentation over the
same polygon selection.

Note: 10x_to_baysor_dataPreProcessor.py only zeroes cell_id when it equals
integer -1, an older Xenium convention. This dataset's cell_id is a string
("UNASSIGNED" / "abcd1234-1"), so that check never fires and "UNASSIGNED"
passes straight through as the Baysor prior-segmentation label. Baysor will
treat it as one ordinary cell and print "No unassigned molecules found in the
prior segmentation" - that warning in baysor_run.log is expected here, not a
pipeline bug (see README-baysor.md).

Example:
    python run_pipeline_comparison.py \
        --region-dir "20250827__205437__CSU_Bouchet_v1_2025-08-27/output-XETG00230__0063687__Region-1A__20250827__205535" \
        --selection 20260903Run/Input/EmmaRegion1LeftSelection.csv \
        --out-dir 20260903Run/Output/pipeline_test
"""

import argparse
import subprocess
import sys
from pathlib import Path

BAYSOR_PROCESSING_DIR = Path(__file__).resolve().parent
TRANSCRIPT_COUNTING_DIR = BAYSOR_PROCESSING_DIR.parent / "transcriptCounting"

DEFAULT_SELECTION = "20260903Run/Input/EmmaRegion1LeftSelection.csv"
DEFAULT_BAYSOR_CMD = r"C:\Users\jrand\.julia\bin\baysor.cmd"

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
    return p.parse_args()


def step(title: str):
    print(f"\n=== {title} ===")


def run(cmd, **kwargs):
    print("+ " + " ".join(str(c) for c in cmd))
    subprocess.run(cmd, check=True, **kwargs)


def main():
    args = parse_args()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    xenium_subset = out_dir / "xenium_subset.parquet"
    legacy_csv = out_dir / "legacy_baysor_transcripts.csv"
    baysor_prefix = "baysor-output"
    baysor_segmentation = out_dir / f"{baysor_prefix}_segmentation.csv"
    baysor_run_log = out_dir / "baysor_run.log"
    xenium_counts = out_dir / "xenium_subset_cell_counts.csv"
    baysor_counts = out_dir / "baysor_cell_counts.csv"

    step("Step 1: polygon subset -> native Xenium parquet")
    run([
        sys.executable, str(BAYSOR_PROCESSING_DIR / "subset_selection_to_parquet.py"),
        "--region-dir", args.region_dir,
        "--selection", args.selection,
        "--out", str(xenium_subset),
    ])

    step("Step 2: legacy preprocessor -> Baysor-ready CSV")
    run([
        sys.executable, str(BAYSOR_PROCESSING_DIR / "10x_to_baysor_dataPreProcessor.py"),
        "-transcript", str(xenium_subset)
    ], cwd=out_dir)
    generated_name = (
        f"X{LEGACY_MIN_X}-{LEGACY_MAX_X}_Y{LEGACY_MIN_Y}-{LEGACY_MAX_Y}_filtered_transcripts.csv"
    )
    generated_path = out_dir / generated_name
    if not generated_path.is_file():
        sys.exit(f"Expected legacy preprocessor output not found: {generated_path}")
    generated_path.replace(legacy_csv)
    print(f"{generated_path.name} -> {legacy_csv.name}")

    step("Step 3: Baysor")
    baysor_cmd = (
        f'"{args.baysor_cmd}" run '
        f"-x x_location -y y_location -z z_location -g feature_name "
        f"--min-molecules-per-cell 1 -p --prior-segmentation-confidence 0.5 "
        f"-o ./{baysor_prefix} "
        f'"{legacy_csv.name}" :cell_id'
    )
    print("+ " + baysor_cmd)
    with open(baysor_run_log, "w", encoding="utf-8") as log:
        subprocess.run(baysor_cmd, shell=True, check=True, cwd=out_dir,
                        stdout=log, stderr=subprocess.STDOUT)
    print(f"Baysor log -> {baysor_run_log}")
    if not baysor_segmentation.is_file():
        sys.exit(f"Expected Baysor output not found: {baysor_segmentation}")

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

    step("Done")
    print(f"Xenium (10x) counts : {xenium_counts}")
    print(f"                      {xenium_counts.with_suffix('.cell_summary.csv')}")
    print(f"Baysor counts       : {baysor_counts}")
    print(f"                      {baysor_counts.with_suffix('.cell_summary.csv')}")
    print(f"Baysor run log      : {baysor_run_log}")


if __name__ == "__main__":
    main()
