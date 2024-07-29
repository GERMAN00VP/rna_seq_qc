# rna_seq_qc

This repository contains a script that allows for depth quality control by rarefaction in RNA sequencing (RNA-seq) data analysis. It includes functions to compute unique counts in parallel and track progress efficiently.

## Functions

### `compute_sample_sizes_and_probabilities`

This function processes the input DataFrame (`data`) to eliminate very low count values, determine the maximum count values, calculate the desired sample sizes, and compute the probability of gene appearance.

#### Parameters:
- `data` (pd.DataFrame): The input DataFrame with gene counts.

#### Returns:
- `sample_sizes` (pd.DataFrame): A DataFrame with the desired values per sample at each percentage.
- `df_prob` (pd.DataFrame): A DataFrame with the probability of gene appearance per sample.

#### Example Usage:

```python
import pandas as pd

# Load your data
data = pd.read_csv('your_data.csv')

# Compute sample sizes and probabilities
sample_sizes, df_prob = compute_sample_sizes_and_probabilities(data)

### `compute_unique_counts`

Calculates the number of unique indices that appear in a sample given the sample size and the probabilities associated with the indices.

#### Parameters:
- `sample_size` (int): The size of the sample to be drawn.
- `prob_series` (pd.Series): A series containing the probabilities associated with each index.
- `index` (tuple): The index of the task (row, column).

#### Returns:
- `tuple`: The index of the task and the number of unique indices in the sample.

### `parallel_compute`

Executes parallel computations to count unique indices in samples and show progress using `tqdm`.

#### Parameters:
- `df_prob` (pd.DataFrame): DataFrame containing the probabilities for each index.
- `sample_sizes` (pd.DataFrame): DataFrame containing the sample sizes for each column.

#### Returns:
- `pd.DataFrame`: DataFrame containing the count of unique indices for each sample size.
