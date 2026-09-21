# Re-segmenting a Xenium selection with Baysor

Step-by-step reproduction of the **EmmaRegion1Left** run: take a region drawn in
Xenium Explorer, cut the matching transcripts out of the 10x output, and
re-segment them with [Baysor](https://www.10xgenomics.com/analysis-guides/using-baysor-to-perform-xenium-cell-segmentation).

Run everything from the repo root (`C:\Games\projects\TranscriptCounting`) in
**PowerShell**.

---

## Inputs

| What | Path used in this run |
|------|-----------------------|
| Xenium region output (contains `transcripts.parquet`) | `20250827__205437__CSU_Bouchet_v1_2025-08-27/output-XETG00230__0063687__Region-1A__20250827__205535/` |
| Xenium Explorer selection (polygon vertex export, µm) | `emmaInputData/EmmaRegion1LeftSelection.csv` |
| Scripts | `src/baysorProcessing/` |
| Output written to | `emmaOutput/EmmaRegion1Left/` |

The selection CSV is what Xenium Explorer produces from
**Selection tool → export coordinates**: a `#`-commented header, then
`Selection,X,Y,Class,Color` rows. `X`/`Y` are in microns, in the *same*
coordinate frame as `transcripts.parquet` `x_location` / `y_location`, so the
polygon is applied directly with no transform.

> If you move the source files, only the `--region-dir` / `--selection` / `--out`
> paths in Step 2 and the `cd` / input filename in Step 4 change.

---

## Step 1 — Python environment

The repo `.venv` already has what the subset script needs
(`pandas`, `numpy`, `pyarrow`, `matplotlib`). To recreate it:

```bash
python -m venv .venv
.venv/Scripts/python.exe -m pip install pandas numpy pyarrow matplotlib
```

## Step 2 — Cut the selection out of the transcripts

[`subset_selection_for_baysor.py`](subset_selection_for_baysor.py):

1. reads `transcripts.parquet` from `--region-dir`
2. keeps `qv >= 20` and real genes only — uses the `is_gene` column, which
   excludes every `NegControlProbe_` / `NegControlCodeword_` / `BLANK_` /
   `antisense_` / `UnassignedCodeword_` / `DeprecatedCodeword_` feature
3. keeps only transcripts **inside the selection polygon**
   (`matplotlib.path.Path.contains_points`)
4. converts the string `cell_id` (`"UNASSIGNED"` or `"abcd1234-1"`) to an
   **integer** prior-segmentation column: `0` = unassigned / cell-free, `1..N` =
   one id per Xenium cell. Baysor's `--prior-segmentation-confidence` needs an
   integer column with `0` as the unassigned label — feeding the raw strings
   makes Baysor treat `"UNASSIGNED"` as a normal cell and warn
   *"No unassigned molecules found in the prior segmentation"*.

```powershell
.venv\Scripts\python.exe src\baysorProcessing\subset_selection_for_baysor.py `
  --region-dir "20250827__205437__CSU_Bouchet_v1_2025-08-27\output-XETG00230__0063687__Region-1A__20250827__205535" `
  --selection emmaInputData\EmmaRegion1LeftSelection.csv `
  --out emmaOutput\EmmaRegion1Left\EmmaRegion1Left_baysor_transcripts.csv
```

This run: `6,850,638 → 6,565,219` after quality filter `→ 109,816` inside the
polygon; prior segmentation = 987 cells + 27,752 unassigned. Output CSV ≈ 8 MB.

*(For a whole region rather than a drawn selection, use the older
[`10x_to_baysor_dataPreProcessor.py`](10x_to_baysor_dataPreProcessor.py) bounding-box crop, or draw a
selection polygon that encloses the region.)*

## Step 3 — Baysor (native Windows, installed via Julia)

Baysor v0.7.1 is installed as a Julia CLI app (Comonicon), Julia 1.10.12. The
launcher is:

```
C:\Users\jrand\.julia\bin\baysor.cmd
```

It is **not on PATH** by default, so the commands below call it by full path.
Optionally add it once:

```powershell
[Environment]::SetEnvironmentVariable('Path', $env:Path + ';' + "$env:USERPROFILE\.julia\bin", 'User')
```

No WSL, no `libopenlibm` executable-stack fix, no file copying — it reads and
writes Windows paths directly. To reinstall/update Baysor, follow the "Binary
download" / Julia install in the [Baysor
README](https://github.com/kharchenkolab/Baysor) (`Pkg.add` the package, then
`using Baysor; Baysor.install_cli()`), which regenerates
`.julia\bin\baysor.cmd`.

## Step 4 — Run Baysor

`-o ./baysor-output` is a **filename prefix, not a directory** — outputs are
written as `baysor-output_*` in the working directory.

```powershell
cd emmaOutput\EmmaRegion1Left
& "C:\Users\jrand\.julia\bin\baysor.cmd" run `
    -x x_location -y y_location -z z_location -g feature_name `
    --min-molecules-per-cell 1 -p --prior-segmentation-confidence 0.5 `
    -o ./baysor-output `
    EmmaRegion1Left_baysor_transcripts.csv :cell_id 2>&1 | Tee-Object run-native.log
cd ..\..
```

Runtime: ~3 min for 110k transcripts (first run of a session also pays Julia
sysimage load). A full region (~6.5M transcripts) is ~1 h, almost all of it in
"Initializing algorithm".

Check the log for `Using the following additional information about molecules:
[:confidence, :cluster, :prior_segmentation]` — that confirms the prior
segmentation was picked up. `Using 3D coordinates` is expected (drop `-z` for a
2D run).

## Step 5 — Outputs

Written to `emmaOutput/EmmaRegion1Left/`:

| File | Contents |
|------|----------|
| `baysor-output_segmentation.csv` | per-transcript: new `cell` id, `is_noise`, `assignment_confidence`, `cluster`, … |
| `baysor-output_cell_stats.csv` | per-cell: centroid, `area`, `n_transcripts`, `elongation`, … |
| `baysor-output_counts.loom` | cell × gene count matrix |
| `baysor-output_polygons_2d.json` / `_3d.json` | cell boundary polygons |
| `baysor-output_diagnostics.html` / `_borders.html` | QC plots |
| `baysor-output_log.log` | run log |

---

## Batch — every ROI in a workbook

`run_pipeline_batch.py` loops the whole subset → Baysor → per-cell-counts flow
over every tab of an Excel workbook of ROI selections (e.g.
`emmaInputData/Modified 10x Pilot -A11 ROI Coordinates.xlsx`, 22 tabs).

**Workbook format** — each worksheet must be a raw Xenium Explorer *export
selection coordinates* dump pasted into column A (an optional label row, a blank
row, then the `#Selection names:` / `#Areas` / `Selection,X,Y,Class,Color` /
`Selection N,<x>,<y>,…` lines). The sheet name must contain `Region <N>`
(`Region 1A Left`, `Region 2D R`, …) so it can be matched to a
`…__Region-<N>__…` output folder under `--data-root`. Parsing is stdlib-only —
no `openpyxl` needed.

```powershell
# list what would run — resolves every sheet -> region folder, does no work
.venv\Scripts\python.exe src\baysorProcessing\run_pipeline_batch.py --dry-run

# one sheet
.venv\Scripts\python.exe src\baysorProcessing\run_pipeline_batch.py --sheets "Region 1A Left"

# everything (resumable); ~3-5 min per ROI, ~1.5 h for 22
.venv\Scripts\python.exe src\baysorProcessing\run_pipeline_batch.py --skip-existing
```

| flag | effect |
|------|--------|
| `--xlsx` | workbook (default: the A11 ROI file) |
| `--data-root` | folder holding the `output-XETG…__Region-*` dirs |
| `--out-root` | where per-sheet folders go (default: `emmaOutput`) |
| `--sheets "A,B"` | run only these tabs |
| `--skip-existing` | skip a tab whose `baysor_cell_counts.csv` already exists |
| `--raw-prior` | feed 10x's string `cell_id` straight to Baysor instead of remapping it to an integer `0 = unassigned` prior (Baysor then logs *"No unassigned molecules found…"*) |
| `--dry-run` | print the plan and exit |

Each tab writes `emmaOutput/<Safe_Sheet_Name>/` with `*_selection.csv`,
`xenium_subset.parquet`, `legacy_baysor_transcripts.csv`, `baysor-output_*`
(segmentation / stats / polygons / plots / log), `baysor_run.log`, and
`*_cell_counts.csv` (+ `.cell_summary.csv`) for **both** the 10x and the Baysor
segmentation.

After the loop it writes `emmaOutput/<workbook stem>_run_summary.csv`, one row
per tab: region, `n_vertices`, `status`/`error`, `n_subset_transcripts`
(all molecules in the polygon), `n_legacy_transcripts` (after qv≥20 + control
drop, i.e. what Baysor saw), cell counts and assigned-fraction for each
segmentation, `baysor_prior_warning`, and `seconds`. A failed tab is recorded
and the loop continues.

Single-ROI equivalent (what the loop calls per sheet):

```powershell
.venv\Scripts\python.exe src\baysorProcessing\run_pipeline_comparison.py `
  --region-dir "20250827__205437__CSU_Bouchet_v1_2025-08-27\output-XETG00230__0063687__Region-1A__20250827__205535" `
  --selection emmaInputData\EmmaRegion1LeftSelection.csv `
  --out-dir emmaOutput\Region_1A_Left
```

---

## This run's result

- Input: 109,816 transcripts, 50 genes, prior segmentation 987 cells.
- Output: **2,111 Baysor cells**, 248 noise transcripts, converged at 500 iters.
- Prior segmentation consumed correctly (no unassigned-value warning).

---

## Appendix — WSL fallback

A Linux Baysor build also exists at `~/baysor/bin/baysor/bin/baysor` in the
`Ubuntu` WSL distro. It only matters if the native install breaks. Two gotchas
that the native path avoids:

- Recent WSL2 kernels refuse Baysor's bundled `libopenlibm.so` (executable
  stack). Clear the flag with
  [`wsl_clear_execstack.py`](wsl_clear_execstack.py):
  ```bash
  wsl -d Ubuntu -e bash -lc 'cat > ~/clear_execstack.py' < src/baysorProcessing/wsl_clear_execstack.py
  wsl -d Ubuntu -e bash -lc 'python3 ~/clear_execstack.py \
    ~/baysor/bin/baysor/lib/julia/libopenlibm.so{,.4,.4.0}'
  ```
- `/mnt/c` throws intermittent `Input/output error` on this machine — copy the
  CSV into the distro (`cat file > ~/...`) and `tar` results back out; and keep
  the `wsl` call in the foreground (the distro stops when `wsl.exe` exits).
