# Entropy-Based Knot Selection Algorithm

This repository contains R code for an entropy-based knot selection algorithm. The purpose of this algorithm is to select knots from a spatial dataframe containing an additional continuous value column, which allows for faster and more accurate results than other methods that only factor in the spatial part. The algorithm's code is located in the `ka_v5_functions.R` file and includes two main functions: the base algorithm and the Grid Search method.

## Data

["8hour_44201_2021.csv"](https://aqs.epa.gov/aqsweb/airdata/download_files.html#eighthour)

## Usage

### Importing the Algorithm

```
source("ka_v5_functions.R")
```

**Note:**

A spatial dataframe, in this context, is a dataframe that contains the columns `x` and `y`.

### Function: Entropy Maximization

**Parameters:**
- `df.points`: A spatial dataframe containing the column `signal`.
- `n.neighbors`: The number of neighbors to consider when calculating metrics for each knot.
- `radius.mult`: The value to multiply each point's radius by when removing other points from the list of knot candidates.
- `max.knots`: The maximum number of iterations to pick knots (not guaranteed to always reach the max).

**Return:**
- The final list of knots: A dataframe containing the `x`, `y`, `signal` (continuous value), and any other metrics used to determine the knots.

**Example:**
```
knots.entropy = entropy_max(df.points=test.data, n.neighbors=5, radius.mult=3, max.knots=25)
```

### Function: Entropy Maximization Grid Search

**Parameters:**
- `df.points`: A spatial dataframe containing the column signal.
- `seq.nn`: A sequence of positive integers of the number of neighbors to consider when calculating metrics for each knot.
- `seq.rm`: A sequence of positive numeric values to multiply each point's radius by when removing other points from the list of knot candidates.
- `max.knots`: The maximum number of iterations to pick knots (not guaranteed to always reach the max).

**Return:**
- A list containing the list of best hyperparameters found, the grid search results, and the final list of knots.

**Example:**
```
emax.gs = entropy_max_gs(df.points=test.data, seq.nn=3:7, seq.rm=seq(from=1, to=2, by=.25), max.knots=n.knots)
knots.entropy = emax.gs$knots
```

**Accessing the Results:**

- `emax.gs$args$n_neighbors`: Best number of neighbors found.
- `emax.gs$args$radius_mult`: Best radius multiplier found.
- `emax.gs$grid`: The grid search results, containing all combinations of hyperparameters and their performance.
- `emax.gs$knots`: The final list of knots.

## Examples

For more examples on how to use the entropy-based knot selection algorithm, you can refer to the `cksm_v5.R` file in the root directory. This file contains sample code demonstrating how to apply the `entropy_max` and `entropy_max_gs` functions on real data.
