# Baysor guide

- Following instructions from [this guide](https://www.10xgenomics.com/analysis-guides/using-baysor-to-perform-xenium-cell-segmentation)
- Required installing and using WSL environment on Windows
- Command used:
```
./baysor-x86_x64-linux-v0.7.1_build/bin/baysor/bin/baysor run   -x x_location   -y y_location   -z z_location   -g feature_name   --min-molecules-per-cell 1   -p --prior-segmentation-confidence 0.5   -o ./baysor-output   X0.0-24000.0_Y0.0-24000.0_filtered_transcripts.csv :cell_id
``` 
