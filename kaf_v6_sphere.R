# Knot Algorithm Functions Version 6 for Spheres Only

# Setup -----------------------------------------------------------------------------------------------------------

# Imports
library(dplyr)    # piping and QoL functions; %>%
library(tidyr)    # for dropping NA valuse; drop_na()
library(scales)   # normalizing lists between 0 and 1; rescale()

# Functions (Formulas and Math) -----------------------------------------------------------------------------------

#' @title Euclidean Distance 3D
#' @description Finds the distance between 2 3D points using euclidean geometry.
#' @param x1 A numeric scalar for the x coordinate of the first point.
#' @param y1 A numeric scalar for the y coordinate of the first point.
#' @param z1 A numeric scalar for the z coordinate of the first point.
#' @param x2 A numeric scalar for the x coordinate of the second point.
#' @param y2 A numeric scalar for the y coordinate of the second point.
#' @param z2 A numeric scalar for the z coordinate of the second point.
#' @return A numeric scalar.
#' @example
#' euclid_dist_3d(3, 4, 5, 0, 0, 0)
euclid_dist_3d = function(x1, y1, z1, x2, y2, z2) { sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2) }

#' @title Pad Edges
#' @description Repeats the first and last elements of a vector `pad.length` times.
#' @param vector The target vector to pad.
#' @param pad.length The amount of repetitions to pad with
#' @return The padded vector.
#' @example
#' pad_edges(c(1, 2, 3, 4, 5), 3)
pad_edges = function(vector, pad.length) {
  first = vector[1]
  last = vector[length(vector)]
  c(rep(first, each=pad.length), vector, rep(last, each=pad.length))
}

#' @title Ebrahimi Entropy
#' @description Calculates the entropy of a vector using Ebrahimi's method.
#' @param signal The vector to find the entropy of.
#' @return A numeric scalar.
#' @example
#' ebrahimi_entropy(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5))
#' ebrahimi_entropy(c(0, 1000, 32, -400, 87, 23, 76, 10000, -1000, 1))
#' ebrahimi_entropy(c(1000, 1000, 1000, 1000, 999))
#' ebrahimi_entropy(c(1000, 1000, 0, 1000, 1000))
#' ebrahimi_entropy(c(1000, 1000, 1000, 1000, 1000))
ebrahimi_entropy = function(signal) {
  sorted_signal = sort(signal)
  n = length(signal)
  m = floor(sqrt(n) + 0.5)
  X = pad_edges(sorted_signal, m)
  differences = tail(X, length(X) - (2 * m)) - head(X, length(X) - (2 * m))
  i = seq(from = 1, to = n, by = 1)
  ci = rep(2, length(i))
  ci[i <= m] = 1 + (i[i <= m] - 1) / m
  ci[i >= n - m + 1] = 1 + (n - i[i >= n - m + 1]) / m
  logs = log(n * differences / (ci * m))
  return(mean(logs))
}

# Functions (Base Algorithm) --------------------------------------------------------------------------------------

#' @title Generate Knot Metrics in Parallel
#' @description Creates metrics for specific locations from a dataframe with spatial information. (In parallel)
#' @param df.points A spatial dataframe containing the column `signal`, which contains numeric data.
#' @param n.neighbors The number of neighbors to consider when calculating metrics for each knot.
#' @return A spatial dataframe, including only the columns `x`, `y`, `z`, `signal`, along with the metric columns.
gkm_parallel = function(df.points, n.neighbors) {
  # Load necessary libraries
  library(doParallel)
  library(foreach)

  # Initialize parallel backend
  n.cores = parallel::detectCores() - 1
  my.cluster = parallel::makeCluster(n.cores, type="PSOCK")
  registerDoParallel(my.cluster)

  # Export helper functions to parallel workers
  parallel::clusterExport(my.cluster, c("euclid_dist_3d", "ebrahimi_entropy", "pad_edges"))

  # Load required packages within parallel workers
  parallel::clusterEvalQ(my.cluster, {
    library(dplyr)
    library(tidyr)
    library(scales)
  })

  # Extract only the relevant columns
  metric.df = data.frame(df.points) %>% subset(select=c(x, y, z, signal))

  # Use foreach and %dopar% to parallelize the loop
  results = foreach(i=1:nrow(df.points), .combine="rbind", .multicombine=T) %dopar% {
    df.points$dist_from_point = euclid_dist_3d(df.points[i, "x"], df.points[i, "y"], df.points[i, "z"], df.points$x, df.points$y, df.points$z)
    df.points.subset = df.points[order(df.points$dist_from_point), ] %>% head(n.neighbors + 1)  # subset with only points within radius
    radius = df.points.subset$dist_from_point %>% max  # find the radius of the point
    entropy = df.points.subset$signal %>% ebrahimi_entropy  # find the entropy of the point
    c(radius, entropy)
  }

  # Stop cluster
  parallel::stopCluster(my.cluster)

  # Add the results to metric.df
  metric.df$radius = results[, 1]
  metric.df$entropy = results[, 2]

  metric.df
}

#' @title Generate Knot Metrics
#' @description Creates metrics for specific locations from a dataframe with spatial information.
#' @param df.points A spatial dataframe containing the column `signal`, which contains numeric data.
#' @param n.neighbors The number of neighbors to consider when calculating metrics for each knot.
#' @return A spatial dataframe, including only the columns `x`, `y`, `z`, `signal`, along with the metric columns.
generate_knot_metrics = function(df.points, n.neighbors) {
  metric.df = data.frame(df.points) %>% subset(select=c(x, y, z, signal))  # extract only the relavent columns
  for (i in 1:nrow(df.points)) {
    df.points$dist_from_point = euclid_dist_3d(df.points[i, "x"], df.points[i, "y"], df.points[i, "z"], df.points$x, df.points$y, df.points$z)
    df.points.subset = df.points[order(df.points$dist_from_point), ] %>% head(n.neighbors + 1)  # subset with only points within radius
    metric.df[i, "radius"] = df.points.subset$dist_from_point %>% max  # find the radius of the point
    metric.df[i, "entropy"] = df.points.subset$signal %>% ebrahimi_entropy  # find the entropy of the point
  }
  metric.df
}

#' @title Anti-Sphere Subset
#' @description Removes all points in a spatial dataframe within the radius of a specified sphere.
#' @param df.points A spatial dataframe.
#' @param x0 The x coordinate of the sphere's center.
#' @param y0 The y coordinate of the sphere's center.
#' @param z0 The z coordinate of the sphere's center.
#' @param r The radius of the sphere
#' @return A spatial dataframe, including all of the columns from `df.points`.
anti_sphere_subset = function(df.points, x0, y0, z0, r) { subset(df.points, euclid_dist_3d(x0, y0, z0, x, y, z) > r) }

#' @title Entropy Maximization
#' @description Given a spatial dataframe with a `signal` column that represents some continuous value, this function will return knots based on the other parametrs.
#' An iterative way of determining knot candidates in relation to their distance from each other.
#' @param df.points A spatial dataframe containing the column `signal`.
#' @param n.neighbors The number of neighbors to consider when calculating metrics for each knot.
#' @param radius.mult The value to multiply each point's radius by when removing other points from the list of knot candidates.
#' @param max.knots The maximum number of iterations to pick knots (not guaranteed to always reach the max).
#' @return The final list of knots.
entropy_max = function(df.points, n.neighbors, radius.mult, max.knots) {
  df.knots = data.frame()
  for(i in 1:max.knots) {
    df.points = generate_knot_metrics(df.points, n.neighbors)  # recompute entropy every iteration
    df.points = df.points[order(df.points$entropy, decreasing=T), ]  # sort df by entropy
    df.knots = df.knots %>% rbind(df.points[1, ])  # append knot with highest entropy to list of knots
    df.points = df.points %>% anti_sphere_subset(df.knots[i, "x"], df.knots[i, "y"], df.knots[i, "z"], df.knots[i, "radius"] * radius.mult)  # remove the points "too close" to the appended knot
  }
  return(df.knots %>% drop_na)
}
