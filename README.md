# Entropy-Based Knot Selection Algorithm

This repository contains R code for an entropy-based knot selection algorithm. The purpose of this algorithm is to select knots from a spatial dataframe containing an additional continuous value column, which allows for faster and more accurate results than other methods that only factor in the spatial part. The algorithm's code is located in the `ka_v6_functions.R` file and includes the main function: the base algorithm.

## Data

["8hour_44201_2021.csv"](https://aqs.epa.gov/aqsweb/airdata/download_files.html#eighthour)

## Usage

### Importing the Algorithm
```
source("ka_v6_functions.R")
```

**Note:**

A spatial dataframe, in this context, is a dataframe that contains the columns `x`, `y`, and `z`.

### Function: Entropy Maximization

**Parameters:**
- `df.points`: A spatial dataframe containing the column `signal`.
- `n.neighbors`: The number of neighbors to consider when calculating metrics for each knot.
- `radius.mult`: The value to multiply each point's radius by when removing other points from the list of knot candidates.
- `max.knots`: The maximum number of iterations to pick knots (not guaranteed to always reach the max).

**Return:**
- The final list of knots: A dataframe containing the `x`, `y`, `z`, `signal` (continuous value), and any other metrics used to determine the knots.

**Example:**
```
knots.entropy = entropy_max(df.points=test.data, n.neighbors=5, radius.mult=3, max.knots=25)
```


### Function: Generate Knot Metrics in Parallel

**Parameters:**
- `df.points`: A spatial dataframe containing the column `signal`.
- `n.neighbors`: The number of neighbors to consider when calculating metrics for each knot.
- `radius.mult`: The value to multiply each point's radius by when removing other points from the list of knot candidates.

**Return:**
- A dataframe containing the generated knot metrics.

**Example:**
```
knot.metrics = gkm_parallel(df.points=test.data, n.neighbors=5, radius.mult=3)
```

## Examples

For more examples on how to use the entropy-based knot selection algorithm, you can refer to the `cksm_v6.R` file in the root directory. This file contains sample code demonstrating how to apply the `entropy_max` and `gkm_parallel` functions on real data.
