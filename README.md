# Receptor Probe Filtering

A reproducible pipeline to rank and select padlock probes for odorant receptor targets, given a pool of candidate probes.

## Requirements
- Python >= 3.9
- Standard scientific stack: `numpy`, `pandas`, `matplotlib`, `scipy`

> ⚠️ Some inputs are large (e.g., a 7.5 GB heterodimer matrix pickle file). Make sure your machine has at least 16 GB RAM (>=32 GB recommended). Close other apps before loading large pickles.

## Data
All data live in the project Dropbox:
- **Dropbox folder:** `Odorant_Receptors/data_monohan_probe_filtering_2025/`

To run locally:
1. Copy the entire contents of `data_monohan_probe_filtering_2025/` into the repo as a folder named **`data/`**.
2. Expect large files.


## Filter algorithm 
![Flow chart](flow.png)

See <a href="algo.pdf" target="_blank">algo.pdf</a> for details.

## Filtering results

#### Global:
<img src="output/fig/round1_sor.png" alt="Flow chart" width="500">

#### Pergene feature filtering:
<img src="output/fig/round2_sor.png" alt="Flow chart" width="800">


Number of kept probes in worst 10% for rank_off_target_combined: 0.77%
Number of kept probes in worst 10% for rank_tm_melting_diff: 0.25%
Number of kept probes in worst 10% for rank_binding_fraction: 0.30%