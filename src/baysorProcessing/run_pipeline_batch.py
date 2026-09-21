#!/usr/bin/env python3
r"""Loop run_pipeline_comparison over every ROI tab in an Excel workbook.

Each worksheet is expected to be a verbatim Xenium Explorer "export selection
coordinates" dump in column A (a free-text label row, a blank row, then
``#Selection names: ...`` / ``#Areas ...`` / ``Selection,X,Y,Class,Color`` /
``Selection N,<x>,<y>,...`` lines). The sheet name must contain ``Region <N>``
(e.g. ``Region 1A Left``, ``Region 2D R``) so it can be matched to a Xenium
region output folder under --data-root.

For each sheet this writes ``<out-root>/<safe_sheet_name>/`` containing the
extracted ``*_selection.csv`` and the full pipeline output (xenium subset,
legacy Baysor CSV, Baysor segmentation + plots, and per-cell counts for both
the 10x and the Baysor segmentation), then writes a one-row-per-sheet
``<out-root>/<workbook stem>_run_summary.csv``.

Examples:
    # see what would run, resolve every region folder, no work done
    python run_pipeline_batch.py --dry-run

    # one sheet
    python run_pipeline_batch.py --sheets "Region 1A Left"

    # everything, resumable
    python run_pipeline_batch.py --skip-existing
"""

import argparse
import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

import pandas as pd

from run_pipeline_comparison import DEFAULT_BAYSOR_CMD, run_pipeline

DEFAULT_XLSX = "emmaInputData/Modified 10x Pilot -A11 ROI Coordinates.xlsx"
DEFAULT_DATA_ROOT = "20250827__205437__CSU_Bouchet_v1_2025-08-27"
DEFAULT_OUT_ROOT = "emmaOutput"

_MAIN = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_REGION_RE = re.compile(r"Region\s*[-_ ]?(\d[A-Za-z])", re.IGNORECASE)


# ---------------------------------------------------------------- xlsx reading

def _shared_strings(z: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in z.namelist():
        return []
    xml = z.read("xl/sharedStrings.xml").decode("utf-8", "replace")
    return ["".join(re.findall(r"<t[^>]*>(.*?)</t>", si, re.S))
            for si in re.findall(r"<si>(.*?)</si>", xml, re.S)]


def iter_sheet_selections(xlsx_path: Path):
    """Yield (sheet_name, [csv_line, ...]) for every worksheet, keeping only the
    column-A lines that start with '#' or 'Selection'."""
    z = zipfile.ZipFile(xlsx_path)
    strings = _shared_strings(z)
    workbook = z.read("xl/workbook.xml").decode("utf-8", "replace")
    rels = z.read("xl/_rels/workbook.xml.rels").decode("utf-8", "replace")
    rid_to_target = dict(re.findall(r'Id="(rId\d+)"[^>]*?Target="([^"]+)"', rels))

    for m in re.finditer(r'<sheet\b[^>]*?name="([^"]+)"[^>]*?r:id="(rId\d+)"', workbook):
        name, rid = m.group(1), m.group(2)
        target = rid_to_target.get(rid, "")
        if "worksheets/" not in target:
            continue
        root = ET.fromstring(z.read("xl/" + target.lstrip("/")))
        lines = []
        for c in root.iter(_MAIN + "c"):
            ref = c.get("r", "")
            if not (ref[:1] == "A" and ref[1:].isdigit()):
                continue
            t, v = c.get("t"), c.find(_MAIN + "v")
            if t == "s" and v is not None:
                val = strings[int(v.text)]
            elif t == "inlineStr":
                is_el = c.find(_MAIN + "is")
                val = "".join(x.text or "" for x in is_el.iter(_MAIN + "t")) if is_el is not None else ""
            elif v is not None:
                val = v.text or ""
            else:
                continue
            val = val.strip()
            if val.startswith("#") or val.startswith("Selection"):
                lines.append(val)
        yield name, lines


# ---------------------------------------------------------------- helpers

def safe_name(sheet_name: str) -> str:
    return re.sub(r"\W+", "_", sheet_name).strip("_")


def resolve_region_dir(data_root: Path, sheet_name: str) -> Path:
    m = _REGION_RE.search(sheet_name)
    if not m:
        raise ValueError(f"sheet name {sheet_name!r} has no 'Region <N>'")
    region = m.group(1).upper()
    matches = sorted(p for p in data_root.glob(f"*Region-{region}*") if p.is_dir())
    if not matches:
        raise FileNotFoundError(f"no folder matching *Region-{region}* under {data_root}")
    if len(matches) > 1:
        raise RuntimeError(f"region {region} is ambiguous: {[p.name for p in matches]}")
    return matches[0]


def n_vertices(lines: list[str]) -> int:
    return sum(1 for ln in lines if re.match(r"Selection \d", ln))


# ---------------------------------------------------------------- driver

def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--xlsx", default=DEFAULT_XLSX, type=Path,
                   help=f"ROI workbook (default: {DEFAULT_XLSX})")
    p.add_argument("--data-root", default=DEFAULT_DATA_ROOT, type=Path,
                   help=f"folder holding the output-XETG..._Region-* dirs (default: {DEFAULT_DATA_ROOT})")
    p.add_argument("--out-root", default=DEFAULT_OUT_ROOT, type=Path,
                   help=f"per-sheet output folders go here (default: {DEFAULT_OUT_ROOT})")
    p.add_argument("--baysor-cmd", default=DEFAULT_BAYSOR_CMD)
    p.add_argument("--sheets", default=None,
                   help="comma-separated sheet names to run (default: all)")
    p.add_argument("--skip-existing", action="store_true",
                   help="skip a sheet whose baysor_cell_counts.csv already exists")
    p.add_argument("--raw-prior", action="store_true",
                   help="feed 10x string cell_id straight to Baysor (no integer prior remap)")
    p.add_argument("--dry-run", action="store_true",
                   help="list resolved sheet -> region -> out-dir and exit")
    return p.parse_args()


def main():
    args = parse_args()
    if not args.xlsx.is_file():
        sys.exit(f"Not found: {args.xlsx}")
    wanted = {s.strip() for s in args.sheets.split(",")} if args.sheets else None

    rows = []
    for name, lines in iter_sheet_selections(args.xlsx):
        if wanted is not None and name not in wanted:
            continue
        safe = safe_name(name)
        out_dir = args.out_root / safe
        nv = n_vertices(lines)
        rec = {
            "sheet": name, "safe_name": safe, "region": "", "region_dir": "",
            "n_vertices": nv, "status": "", "error": "",
            "n_subset_transcripts": None, "n_legacy_transcripts": None,
            "xenium_cells": None, "xenium_assigned": None, "xenium_assigned_frac": None,
            "baysor_cells": None, "baysor_assigned": None, "baysor_assigned_frac": None,
            "baysor_prior_warning": None, "seconds": None,
        }
        try:
            region_dir = resolve_region_dir(args.data_root, name)
            rec["region"] = _REGION_RE.search(name).group(1).upper()
            rec["region_dir"] = region_dir.name
            if nv < 4:
                raise ValueError(f"only {nv} polygon vertices")

            if args.dry_run:
                rec["status"] = "dry-run"
                print(f"{name:20s} -> {region_dir.name:52s} {nv:4d} verts -> {out_dir}")
            else:
                out_dir.mkdir(parents=True, exist_ok=True)
                sel_csv = out_dir / f"{safe}_selection.csv"
                sel_csv.write_text("\n".join(lines) + "\n", encoding="utf-8")
                print(f"\n{'#' * 70}\n# {name}  ->  {region_dir.name}\n{'#' * 70}")
                stats = run_pipeline(
                    region_dir, sel_csv, out_dir,
                    baysor_cmd=args.baysor_cmd,
                    fix_prior=not args.raw_prior,
                    skip_existing=args.skip_existing,
                )
                rec.update(stats)
                rec["status"] = "ok"
        except Exception as e:  # noqa: BLE001 - isolate one sheet's failure
            rec["status"] = "failed"
            rec["error"] = f"{type(e).__name__}: {e}"
            print(f"!!! {name} FAILED: {rec['error']}", file=sys.stderr)

        denom = rec["n_subset_transcripts"]
        if denom:
            if rec["xenium_assigned"] is not None:
                rec["xenium_assigned_frac"] = round(rec["xenium_assigned"] / denom, 4)
            if rec["baysor_assigned"] is not None:
                rec["baysor_assigned_frac"] = round(rec["baysor_assigned"] / denom, 4)
        rows.append(rec)

    summary = pd.DataFrame(rows)
    print("\n" + "=" * 70 + "\nRUN SUMMARY\n" + "=" * 70)
    with pd.option_context("display.max_columns", None, "display.width", 240):
        print(summary.to_string(index=False))

    if not args.dry_run and not summary.empty:
        args.out_root.mkdir(parents=True, exist_ok=True)
        summary_path = args.out_root / f"{args.xlsx.stem}_run_summary.csv"
        summary.to_csv(summary_path, index=False)
        print(f"\nSummary -> {summary_path}")
        n_fail = int((summary["status"] == "failed").sum())
        if n_fail:
            print(f"{n_fail} sheet(s) failed - see the 'error' column.", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
