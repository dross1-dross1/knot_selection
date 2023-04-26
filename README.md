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

**Return:**
- A dataframe containing the generated knot metrics.

**Example:**
```
knot.metrics = gkm_parallel(df.points=test.data, n.neighbors=5)
```


### Function: Fit ellipses
Fits anisotropic ellipses to the data using a Gaussian Process model. Fit the model by fitting the `.stan` model supplied. Extract the ellipse data using the `extract_ellipses()` function, documented below.

`extract_ellipses()`
**Parameters:**
- `df`: A dataframe containing the locations on which to fit the data. Must contain columns `x` and `y.
- `df.all` (Optional): A dataframe containing all locations, including locations on which the model _was not_ fit. Useful only for graphing purposes when fitting the model on a subset of the data.
- `fit`: A `cmdstanr` model fit object which is the fitted model described above.
- `plot`: A bool of whether or not to plot the fitted ellipses.
- `scale_ellipses`: An integer indicating how much to scale the major and minor axes of each ellipse. So a value of 2 (for example) would make the ellipse half as long on both the major and minor axes. Useful because the model fits large values of the ellipse (for statistical reasons). Note that this _is only used for graphing and does not alter the actual major/minor axes of the ellipse_.
- `num_ellipses`: An integer specifying how many ellipses to plot. If `knots` are supplied, this argument is ignored. Otherwise, a random select of `num_ellipses` ellipses are selected for plotting.
- `knots`: (Optional) A dataframe giving the `x` and `y` coordinates of each knot.
- `psi_suffix`: Just set this to the string `"_all"`.
- `return_df`: A bool indicating whether or not to return the dataframe of ellipse info. Only set to `FALSE` if you just want to plot the data and aren't interested in manipulating it.

## Examples

For more examples on how to use the entropy-based knot selection algorithm, you can refer to the `cksm_v6.R` file in the root directory. This file contains sample code demonstrating how to apply the `entropy_max` and `gkm_parallel` functions on real data.

### Example of fitting ellipses
```
# Functions to load EPA data for testing/visualization purposes
source('load_epa_data.R')

# Load and plot EPA data
epa.df = load_epa_data()
ggplot(epa.df) +
  geom_point(aes(x=x, y=y, fill=pm2_5), size=3, pch=21, color='black') +
  scale_fill_gradient2(low='blue', high='red', midpoint=0)

# Compile stan model
model.knots = cmdstan_model('aniso_process_axes_latent_knots_spectral_profile.stan')

# Pick a sample of the EPA data for testing purposes
data.knots = prepare_data(epa.df, num_knots=15, sample_pct=0.25, plot=TRUE)

# Fit the stan model
fit.knots = model.knots$sample(data=data.knots,  # Data for the stan model
                                parallel_chains=4,
                                iter_warmup=2000, # Change as desired
                                max_treedepth=10, # Change to higher number only if sampling complains about iterations hitting max treedepth
                                init=function() list(sigma_z=0.1, # Initial values for parameters for more efficient sampling
                                                     ell_z=0.1,
                                                     sigma_interp=0.1,
                                                     ell_interp=0.1,
                                                     lambda_y=0.1))

# Extract locations and data
data.knots.spatial.df = cbind(data.knots$spatial_locs, data.knots$y_spatial)
names(data.knots.spatial.df) = c('x', 'y', 'pm2_5')
data.knots.knots.df = data.knots$knot_locs

# Plot the ellipses and save ellipse data
ellipse_df = extract_ellipses(data.knots.spatial.df, # Data showing locations and knots
                               fit=fit.knots, # Stan model fit
                               scale_ellipses=7, # Pick a value that makes it so that you can actually see the ellipses
                               # knots=data.knots.knots.df, # Uncomment if you only want to show ellipses at the knots
                               psi_suffix="_all") # Necessary when showing _all_ ellipses
                               
# Get major/minor axes (radii in x/y directions) and rotations from ellipses
ellipse_locs = ellipse_df[, c('x', 'y')]
ellipse_rotations = ellipse_df$alpha
ellipse_x_y_radii = ellipse_df[, c('a', 'b')]
```
