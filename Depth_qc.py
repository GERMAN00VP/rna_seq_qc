import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm


def compute_sample_sizes_and_probabilities(data):
    """
    Eliminates rows with very low counts, determines the maximum count values,
    calculates desired sample sizes, and computes the probability of gene appearance.

    Parameters:
    data (pd.DataFrame): The input DataFrame with gene counts.

    Returns:
    tuple: Two DataFrames - sample_sizes and df_prob.
    """
    # Eliminate rows with very low count values
    data = data.loc[data.sum(axis=1) > 1000]

    # Determine the maximum count value per sample
    max_value = data.sum(axis=0).max()

    # Calculate the decreasing sequence of percentages
    percentages = np.arange(100, 0, -0.5)  # Percentages from 100% to 0%

    # Calculate desired values for counts per sample at each percentage
    desired_values = np.round(max_value * (percentages / 100))

    # Calculate sample sizes DataFrame
    # DataFrame with desired values per sample at each percentage
    # If the sample is smaller than the desired value, use the total count
    sample_sizes = data.apply(lambda x: [min(x.sum(), desired_values[i]) for i in range(len(desired_values))])

    # Calculate probability DataFrame
    # Probability of gene appearance per sample, based on observed frequency
    df_prob = pd.DataFrame([(data[col] / data[col].sum()).values for col in data.columns], 
                           index=data.columns, 
                           columns=data.index).T

    return sample_sizes, df_prob


def compute_unique_counts(sample_size, prob_series, index):
    """
    Calculate the number of unique indices that appear in a sample given
    the sample size and the probabilities associated with the indices.

    Parameters:
    sample_size (int): The size of the sample to be drawn.
    prob_series (pd.Series): A series containing the probabilities associated with each index.
    index (tuple): The index of the task (row, column).

    Returns:
    tuple: The index of the task and the number of unique indices in the sample.
    """
    indices = prob_series.index
    probabilities = prob_series.values
    
    # Generate samples and count the unique indices
    sampled_indices = np.random.choice(indices, size=sample_size, p=probabilities, replace=True)
    unique_counts = np.sum(np.unique(sampled_indices, return_counts=True)[1] > 0)
    
    return index, unique_counts


def parallel_compute(df_prob, sample_sizes):
    """
    Execute parallel computations to count unique indices in samples
    and show progress using tqdm.

    Parameters:
    df_prob (pd.DataFrame): DataFrame containing the probabilities for each index.
    sample_sizes (pd.DataFrame): DataFrame containing the sample sizes for each column.

    Returns:
    pd.DataFrame: DataFrame containing the count of unique indices for each sample size.
    """
    # Create a list to store the results
    results = []
    
    # Prepare the tasks to be executed in parallel
    tasks = []
    for i in range(sample_sizes.shape[0]):
        for col_index, col in enumerate(df_prob.columns):
            tasks.append((sample_sizes[col][i].astype(int), df_prob[col], (i, col_index)))
    
    # Use tqdm to show progress
    with ProcessPoolExecutor(max_workers=31) as executor:
        futures = [executor.submit(compute_unique_counts, *task) for task in tasks]
        
        # Monitor progress
        for future in tqdm(as_completed(futures), total=len(futures), desc="Progress"):
            results.append(future.result())
    
    # Sort results according to the original index
    results.sort(key=lambda x: x[0])
    sorted_results = [result for _, result in results]
    
    # Convert the results to a DataFrame
    sorted_results = np.array(sorted_results).reshape(sample_sizes.shape[0], len(df_prob.columns))
    df_non_zeros = pd.DataFrame(sorted_results, columns=df_prob.columns)
    
    return df_non_zeros

