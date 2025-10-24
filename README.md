# receptor_probe_filtering

A reproducible pipeline to rank and select padlock probes for odorant receptor targets.

## Requirements
- Python ≥ 3.9
- Standard scientific stack (stock Anaconda is fine): `numpy`, `pandas`, `matplotlib`

> ⚠️ Some inputs are large (e.g., a 7.5 GB heterodimer matrix pickle file). Make sure your machine has at least 16 GB RAM (>=32 GB recommended). Close other apps before loading large pickles.

## Data
All data live in the project Dropbox:

- **Dropbox folder:** Odorant_Receptors: `data_monohan_probe_filtering_2025`

To run locally:
1. Copy the entire contents of `data_monohan_probe_filtering_2025/` into the repo as a folder named **`data/`**.
2. Expect large files.

<!-- ## Quick Start
- Open the main notebook and run all cells. The notebook will:
  1) Pre-filter globally offensive probes using the heterodimer fraction matrix,
  2) Score remaining probes with three z-scored criteria,
  3) Combine scores via exponential scaling,
  4) Keep at most the **10 best** probes per gene. -->
<!-- 
## Workflow Overview

### 1) Global offender pre-filter (heterodimer fractions)
We compute, for each probe \(i\), the summed heterodimer **fraction bound** across all partners:
\[
s_i = \sum_j P_{ij},\quad \text{(diagonal zero)}.
\]
Greedily remove the worst \(1\%\) by iteratively dropping the probe with the largest \(s_i\) and updating the sums (subtract that probe’s column). This trims extreme cross-binders before per-gene ranking.

### 2) Feature computation (per probe)
We use three per-probe features, oriented so **smaller = worse**:
1. **`arm_tm_diff_score`**: absolute Tm difference between left/right arms (°C).  
2. **`binding_score`**: scaled on-target binding probability to its own target.  
3. **`off_target_score`**: sum of Tm to other genes (°C), using a cutoff of 30 °C to include only low-affinity binders.

Each feature is **z-scored** across the (pre-filtered) probe set:
\[
z_{ik}=\frac{x_{ik}-\mu_k}{\sigma_k}.
\]

### 3) Exponential (Boltzmann-style) combination
For each metric \(k\) with weight \(w_k\) and slope \(\beta_k\):
\[
r_{ik}=w_k\,e^{-\beta_k\,z_{ik}},\qquad
S_i=\sum_k r_{ik}.
\]
- Smaller \(z_{ik}\) (worse) ⇒ larger \(r_{ik}\) ⇒ larger total score \(S_i\) (worse overall).
- Default \(w_k=1\), \(\beta_k=1\). (You can tune \(\beta_k\) to change per-SD sensitivity.)

### 4) Per-gene selection (cap at 10)
For each gene, sort by \(S_i\) **ascending** (smallest = best) and keep the **top 10**.  
If a gene has ≤10 probes, keep them all.

## Outputs
- A table with:
  - The three raw features and their z-scores,
  - Per-metric exponential terms \(r_{ik}\),
  - Total score \(S_i\),
  - Overall rank/percentile (optional),
  - Per-gene rank and a boolean **`keep`** flag.
- Plots: histograms of features (before/after pre-filter), and QC visuals showing kept vs. dropped tails.

## Notes & Tips
- **Memory:** Pickle loads fully into RAM. If you only need aggregates, compute and persist summaries (e.g., per-probe sums) then `del` the big arrays.
- **Determinism:** Ties are broken by `padlock_name` to make rankings reproducible.
- **Tuning:** If one feature dominates, cap extreme z-scores (e.g., \(|z|\le 6\)) or adjust \(\beta_k\)/\(w_k\).

## Reproducing the results
1. Ensure `data/` is populated as described.
2. Run the notebook end-to-end.
3. Review the QC tables/plots and the final per-gene keep list. -->
