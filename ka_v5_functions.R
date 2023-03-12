# Setup -----------------------------------------------------------------------------------------------------------

# Imports
library(dplyr)    # piping and QoL functions; %>%
library(tidyr)    # for dropping NA valuse; drop_na()
library(scales)   # normalizing lists between 0 and 1; rescale()
library(mgcv)     # grid search; gam()

# Functions (Formulas and Math) -----------------------------------------------------------------------------------

#' @title Euclidean Distance 2D
#' @description Finds the distance between 2 2D points using euclidean geometry.
#' @param x1 A numeric scalar for the x coordinate of the first point.
#' @param y1 A numeric scalar for the y coordinate of the first point.
#' @param x2 A numeric scalar for the x coordinate of the second point.
#' @param y2 A numeric scalar for the y coordinate of the second point.
#' @return A numeric scalar.
#' @example
#' euclid_dist_2d(3, 4, 0, 0)
euclid_dist_2d = function(x1, y1, x2, y2) { sqrt((x1-x2)^2 + (y1-y2)^2) }

#' @title Pad Edges
#' @description Repeats the first and last elements of a vector `pad.length` times.
#' @param vector The target vector to pad.
#' @param pad.length The amount of repetitions to pad with
#' @return The padded vector.
pad_edges = function(vector, pad.length) {
  first = vector[1]
  last = vector[length(vector)]
  c(rep(first, each = pad.length), vector, rep(last, each=pad.length))
}

#' @title Ebrahimi Entropy
#' @description Calculates the entropy of a vector using Ebrahimi's method.
#' @param signal The vector to find the entropy of.
#' @return A numeric scalar.
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
  mean(logs)
}

# Functions (Base Algorithm) --------------------------------------------------------------------------------------

# spatial dataframe: A dataframe that contains the columns `x` and `y`.

#' @title Generate Knot Metrics
#' @description Creates metrics for specific locations from a dataframe with spatial information.
#' @param df.points A spatial dataframe containing the column `signal`, which contains numeric data.
#' @param n.neighbors The number of neighbors to consider when calculating metrics for each knot.
#' @return A spatial dataframe, including only the columns `x`, `y`, `signal`, along with the metric columns.
generate_knot_metrics = function(df.points, n.neighbors) {
  metric.df = data.frame(df.points) %>% subset(select=c(x, y, signal))
  for (i in 1:nrow(df.points)) {
    df.points$dist_from_point = euclid_dist_2d(df.points[i, "x"], df.points[i, "y"], df.points$x, df.points$y)
    df.points.subset = df.points[order(df.points$dist_from_point), ] %>% head(n.neighbors + 1)  # subset with only points within radius
    metric.df[i, "radius"] = df.points.subset$dist_from_point %>% max  # find the radius of the point
    metric.df[i, "entropy"] = df.points.subset$signal %>% ebrahimi_entropy  # find the entropy of the point
  }
  metric.df
}

#' @title Anti-Circle Subset
#' @description Removes all points in a spatial dataframe within the radius of a specified circle.
#' @param df.points A spatial dataframe.
#' @param x0 The x coordinate of the circle's center.
#' @param y0 The y coordinate of the circle's center.
#' @param r The radius of the circle.
#' @return A spatial dataframe, including all of the columns from `df.points`.
anti_circle_subset = function(df.points, x0, y0, r) { subset(df.points, euclid_dist_2d(x0, y0, x, y) > r) }

#' @title Spaced Pruning
#' @description An iterative way of determining knot candidates in relation to their distance from each other.
#' @param df.points A spatial dataframe containing the column `entropy`.
#' @param radius.mult The value to multiply each point's radius by when removing other points from the list of knot candidates.
#' @param max.knots The maximum number of iterations to pick knots (not guaranteed to always reach the max).
#' @return The final list of knots.
spaced_prune = function(df.points, radius.mult, max.knots) {
  df.knots = data.frame()
  df.points = df.points[order(df.points$entropy, decreasing=T), ]  # sort df by entropy
  for(i in 1:max.knots) {
    df.knots = df.knots %>% rbind(df.points[1, ])  # append knot with highest entropy to list of knots
    df.points = df.points %>% anti_circle_subset(df.knots[i, "x"], df.knots[i, "y"], df.knots[i, "radius"] * radius.mult)  # remove the points "too close" to the appended knot
  }
  df.knots %>% drop_na()
}

#' @title Entropy Maximization
#' @description Given a spatial dataframe with a `signal` column that represents some continuous value, this function will return knots (representative points) based on the other parametrs.
#' @param df.points A spatial dataframe containing the column `signal`.
#' @param n.neighbors The number of neighbors to consider when calculating metrics for each knot.
#' @param radius.mult The value to multiply each point's radius by when removing other points from the list of knot candidates.
#' @param max.knots The maximum number of iterations to pick knots (not guaranteed to always reach the max).
#' @return The final list of knots.
entropy_max = function(df.points, n.neighbors, radius.mult, max.knots) {
  df.metrics = generate_knot_metrics(df.points, n.neighbors)
  df.metrics %>% spaced_prune(radius.mult, max.knots)
}

# Functions (Grid Search) -----------------------------------------------------------------------------------------

#' @title Evaluate Knots via Mean Squared Error
#' @description Given a spatial dataframe with a `signal` column and a list of knots, this function will fit a GAM to determine how well the knots represent the data.
#' @param df.points A spatial dataframe containing the column `signal`, which contains numeric data.
#' @param df.knots A subset of the spatial dataframe (columns may not match) containing the knots for `df.points`.
#' @return A scalar; the Mean Squared Error.
eval_knots_mse = function(df.points, df.knots) {
  gam.eval = gam(signal ~ te(x, y, k=3, bs="gp"),data=df.knots,method="REML", family=gaussian)
  df.points[, "pkn"] = predict(gam.eval, newdata=df.points)
  mean(df.points$signal - df.points[, "pkn"])^2  # mse
}

#' @title Entropy Maximization Grid Search
#' @description Uses grid search cross-validation to select the optimal hyperparameters for the Entropy Maximization algorithm from a list of input ranges.
#' @param df.points A spatial dataframe containing the column `signal`.
#' @param seq.nn A sequence of positive integers of the number of neighbors to consider when calculating metrics for each knot.
#' @param seq.rm A sequence of postive numeric values to multiply each point's radius by when removing other points from the list of knot candidates.
#' @param max.knots The maximum number of iterations to pick knots (not guaranteed to always reach the max).
#' @return A list containing the list of best hyperparameters found, the grid search results, and the final list of knots.
entropy_max_gs = function(df.points, seq.nn, seq.rm, max.knots) {
  results = data.frame(k_actual=integer(), nn=integer(), rm=double(), mse=double())
  for (nn in seq.nn) {
    for (rm in seq.rm) {
      knots = entropy_max(df.points, n.neighbors=nn, radius.mult=rm, max.knots=max.knots)
      mse = eval_knots_mse(df.points, knots)
      results = results %>% rbind(data.frame(k_actual=nrow(knots), nn=nn, rm=rm, mse=mse))
    }
  }
  best.mse = filter(results, mse == min(mse))[1, ]
  best.args = list(n_neighbors=best.mse$nn, radius_mult=best.mse$rm)
  best.knots = entropy_max(df.points, best.args$n_neighbors, best.args$radius_mult, max.knots)
  list(args=best.args, grid=results, knots=best.knots)
}