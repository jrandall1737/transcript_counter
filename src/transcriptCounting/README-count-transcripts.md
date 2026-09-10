# count_transcripts_per_cell.py

Command-line version of the `calculate_transcript_counts` logic from
[`transcript_counter.ipynb`](transcript_counter.ipynb). Produces a **cell × gene
count matrix** (one row per cell, one column per gene, integer counts; first
column is the cell id, so `scanpy.read_csv(path, first_column_names=True)` loads
it directly).

It reads either input:

| Input | Cell assignment used | Auto-detected by |
|-------|----------------------|------------------|
| Xenium `transcripts.parquet` | 10x `cell_id` | `.parquet` extension |
| Baysor `<prefix>_segmentation.csv` | Baysor's re-segmented `cell` | `cell` + `gene` columns present |

Force with `--source {xenium,baysor}` if detection is wrong. Parquet or CSV are
both accepted for either source; `--cell-column` / `--gene-column` override the
column names (e.g. counting the Baysor input CSV by its `cell_id` /
`feature_name` gives the pre-Baysor 10x segmentation over the same molecules).

## Filtering (all adjustable)

1. `qv >= --min-qv` (default 20, matches 10x).
2. Drop control features — `NegControlProbe_`, `NegControlCodeword_`,
   `UnassignedCodeword_`, `DeprecatedCodeword_`, `BLANK_`, `antisense_`
   (`--keep-controls` to disable). Baysor input is usually already gene-only.
3. Baysor only: drop `is_noise` transcripts (`--keep-noise` to disable).
4. Drop transcripts not in any cell — `UNASSIGNED` / `0` / `-1` / blank
   (`--include-unassigned` keeps them as one row named `UNASSIGNED`).
5. Optional `--cells region.csv` — restrict to a list of cell ids from a column
   (`--cell-list-column`, default `Cell ID`), e.g. a Xenium Explorer region
   export. Reports how many requested cells had no transcripts.
6. Optional per-cell filters: `--min-counts-per-cell N`, and repeatable
   `--require GENE:MIN` (keeps cells with `count(GENE) >= MIN`). The notebook's
   `Th > 5` filter is `--require Th:6`.

## Outputs

- `--out` : the cell × gene matrix.
- `--long` : also `<out>.long.csv` with `cell_id,gene,count` rows (count > 0).
- `--summary` : also `<out>.cell_summary.csv` — per-cell `n_transcripts`,
  `n_genes`; for Baysor input it also joins `<prefix>_cell_stats.csv`
  (centroid, area, density, cluster, confidence, …; Baysor's own
  `n_transcripts` comes through as `n_transcripts_baysor`).

## Run it on the Baysor output

After the Baysor run (see
[`../baysorProcessing/README-baysor.md`](../baysorProcessing/README-baysor.md)):

```bash
.venv/Scripts/python.exe src/transcriptCounting/count_transcripts_per_cell.py \
  --transcripts emmaOutput/EmmaRegion1Left/baysor-output_segmentation.csv \
  --out emmaOutput/EmmaRegion1Left/EmmaRegion1Left_baysor_cell_counts.csv \
  --summary --long
```

EmmaRegion1Left result: 109,816 transcripts → 248 noise dropped →
**2,111 cells × 50 genes**, 109,568 transcripts counted (matches Baysor's own
per-cell `n_transcripts` exactly).

## Run it the original way (10x segmentation)

```bash
.venv/Scripts/python.exe src/transcriptCounting/count_transcripts_per_cell.py \
  --transcripts "20250827__205437__CSU_Bouchet_v1_2025-08-27/output-XETG00230__0063687__Region-1A__20250827__205535/transcripts.parquet" \
  --out output/region1A_cell_counts.csv
```

To reproduce the notebook's region-by-region batch, loop the command over the
per-region selection CSVs (`bouchetinputData/*.csv`) passing each as `--cells`,
then concatenate the resulting matrices.

## 10x vs Baysor on the same molecules

Count the pre-Baysor 10x segmentation over the exact transcript set that went
into Baysor (`EmmaRegion1Left_baysor_transcripts.csv`, whose `cell_id` column is
the integer-coded 10x assignment, `0` = unassigned):

```bash
.venv/Scripts/python.exe src/transcriptCounting/count_transcripts_per_cell.py \
  --transcripts emmaOutput/EmmaRegion1Left/EmmaRegion1Left_baysor_transcripts.csv \
  --source xenium --cell-column cell_id --gene-column feature_name \
  --out emmaOutput/EmmaRegion1Left/EmmaRegion1Left_xenium_cell_counts.csv --summary
```

| | 10x (Xenium) | Baysor |
|---|---|---|
| cells | 987 | 2,111 |
| transcripts assigned | 82,064 (74.7%) | 109,568 (99.8%) |
| median transcripts/cell | 40 | 20 |
| median genes/cell | 15 | 9 |
| cells < 10 transcripts | 127 | 442 |

Baysor keeps nearly every molecule (only 248 dropped as noise vs 27,752 left
cell-free by 10x) and splits the tissue into ~2× as many, smaller cells.
Per-gene assigned totals stay well correlated (Pearson r = 0.95); the biggest
relative gain is `Slc1a2` (astrocyte marker, ×2.2) — diffuse glial transcripts
that 10x's nucleus-expansion segmentation tends to leave unassigned.
See `EmmaRegion1Left_xenium_vs_baysor{_comparison,_gene_totals}.csv` and
`…_xenium_vs_baysor.png` in `emmaOutput/EmmaRegion1Left/`.
