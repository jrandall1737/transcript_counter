import os
import pandas as pd
import numpy as np

# --- Windows default paths ---
BASE_DIR = os.getcwd()
INPUT_DIR = os.path.join(BASE_DIR, "input")
OUTPUT_DIR = os.path.join(BASE_DIR, "output")

INPUT_CSV = os.path.join(INPUT_DIR, "combined_th_filtered_transcript_counts.csv")
VECTORS_TSV = os.path.join(OUTPUT_DIR, "vectors.tsv")
METADATA_TSV = os.path.join(OUTPUT_DIR, "metadata.tsv")

# --- Configuration ---
FORCE_METADATA_COLS = {"cell_id", "transcript", "cluster", "sample", "batch"}
INCLUDE_TRANSCRIPT_TEXT = True
DROP_ALL_ZERO_VECTOR_COLS = True
STANDARDIZE_VECTORS = False
# ---------------------


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    if not os.path.exists(INPUT_CSV):
        raise FileNotFoundError(
            f"Could not find input file at: {INPUT_CSV}\n"
            f"Make sure it exists inside the 'input' folder."
        )

    df = pd.read_csv(INPUT_CSV)

    if df.empty:
        raise ValueError("Input CSV appears to be empty.")

    # Identify metadata columns
    forced_present = [c for c in df.columns if c in FORCE_METADATA_COLS]
    non_numeric_cols = df.select_dtypes(exclude=[np.number]).columns.tolist()
    metadata_cols = list(dict.fromkeys(forced_present + non_numeric_cols))

    if not INCLUDE_TRANSCRIPT_TEXT and "transcript" in metadata_cols:
        metadata_cols = [c for c in metadata_cols if c != "transcript"]

    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    vector_cols = [c for c in numeric_cols if c not in metadata_cols]

    if not vector_cols:
        raise ValueError("No numeric vector columns detected.")

    vectors = df[vector_cols].replace([np.inf, -np.inf], np.nan).fillna(0)

    if DROP_ALL_ZERO_VECTOR_COLS:
        nonzero_mask = (vectors != 0).any(axis=0)
        vectors = vectors.loc[:, nonzero_mask]
        vector_cols = vectors.columns.tolist()

    if vectors.shape[1] == 0:
        raise ValueError("No vector dimensions remain after filtering.")

    if STANDARDIZE_VECTORS:
        means = vectors.mean(axis=0)
        stds = vectors.std(axis=0).replace(0, 1)
        vectors = (vectors - means) / stds

    if metadata_cols:
        metadata = df[metadata_cols].copy()
    else:
        metadata = pd.DataFrame({"row_index": np.arange(len(df))})

    # Write files for TensorFlow Projector
    vectors.to_csv(VECTORS_TSV, sep="\t", header=False, index=False)
    metadata.to_csv(METADATA_TSV, sep="\t", header=True, index=False)

    print("Conversion complete.")
    print(f"Vectors written to:  {VECTORS_TSV}")
    print(f"Metadata written to: {METADATA_TSV}")
    print(f"Rows: {vectors.shape[0]}, Vector dimensions: {vectors.shape[1]}")


if __name__ == "__main__":
    main()