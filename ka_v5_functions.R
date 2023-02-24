# Setup -----------------------------------------------------------------------------------------------------------

# Imports
library(dplyr)    # piping and QoL functions; %>%
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
#' euclid_dist(3, 4, 0, 0)
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
  # logs = logs %>% replace(logs == -Inf, 0)
  mean(logs)
}

# Functions (Base Algorithm) --------------------------------------------------------------------------------------

# spatial dataframe: A dataframe that contains the columns `x` and `y`.

#' @title Generate Knot Metrics
#' @description Creates metrics for specific locations from a dataframe with spatial information.
#' @param df.points A spatial dataframe containing the column `signal`, which is normalized between 0 and 1.
#' @param n.neighbors The number of neighbors to consider when calculating metrics for each knot.
#' @return A dataframe, matching the row count of `df.locations` that gives metrics on `df.points`.
generate_knot_metrics = function(df.points, n.neighbors) {
  metric.df = data.frame(df.points) %>% subset(select=c(x, y, signal))
  for (i in 1:nrow(df.points)) {
    df.points$dist_from_point = euclid_dist_2d(df.points[i, "x"], df.points[i, "y"], df.points$x, df.points$y)
    df.points.subset = df.points[order(df.points$dist_from_point), ] %>% head(n.neighbors + 1)
    metric.df[i, "radius"] = df.points.subset$dist_from_point %>% max
    metric.df[i, "avg_signal"] = df.points.subset$signal %>% mean
    metric.df[i, "entropy"] = df.points.subset$signal %>% ebrahimi_entropy
  }
  metric.df
}

# takes center and applies function to the subset
modify_metrics = function(df.km, x0, y0, r) {

}

ss_new = function(df.km, radius.mult, max.knots, cols.to.sort) {
  df.knots = data.frame()
  for (i.knot in 1:max.knots) {
    for (i.col in length(cols.to.sort):1) { df.km = df.km[order(-df.km[, c(cols.to.sort[i.col])]), ] }  # sort descending
    # print(df.km %>% head)
    # print(view_world(df.km))
    # update df.km with new signal
    is.neighbor = euclid_dist_2d(df.km$x, df.km$y, df.km[1, "x"], df.km[1, "y"]) <= df.km[1, "radius"] * radius.mult
    df.km$signal[is.neighbor] = df.km$signal[is.neighbor] - df.km[1, "avg_signal"] #.5
    df.km = generate_knot_metrics(df.km, 10)
    df.knots = df.knots %>% rbind(df.km[1, ])  # append first row to df.knots
  }
  df.knots
}

# test.data %>% generate_knot_metrics(15) %>% ss_new(2, 15, c("entropy")) %>% view_knots_r(test.data, "Testing")

# TODO
anti_circle_subset = function(df.points, x0, y0, r) { subset(df.points, euclid_dist_2d(x0, y0, x, y) > r) }

# TODO
subset_spaced = function(df.points, radius.mult) {
  df.points = df.points[order(df.points$entropy), ]
  new.df = data.frame() %>% rbind(df.points[1, ])  # create new df with starting row of index 1
  while(nrow(df.points) > 0) {
    i = nrow(new.df)
    df.points = anti_circle_subset(df.points, new.df[i, "x"], new.df[i, "y"], new.df[i, "radius"] * radius.mult)
    # df.points$signal = df.points$signal - mean(df.points$signal)
    new.df = rbind(new.df, df.points[1, ])  # append subset
  }
  new.df[-c(nrow(new.df)), ]  # last row is always NA, this removes it
}

# TODO
ss_sorted = function(df.points, radius.mult, max.knots, cols.to.sort) {
  df.points = df.points %>% subset_spaced(radius.mult)
  for(i in length(cols.to.sort):1) {
    df.points = df.points[order(-df.points[, c(cols.to.sort[i])]), ]  # select highest entropy
  }
  df.points %>% head(max.knots)
}

# TODO
#' @title Variable Knot Radius Base with no Argument List
vkr_base_no_list = function(df.points, n.neighbors, radius.mult, max.knots, cols.to.sort) {
  df.metrics = generate_knot_metrics(df.points, n.neighbors)
  df.metrics %>% ss_sorted(radius.mult, max.knots, cols.to.sort)
}

# TODO
#' @title Variable Knot Radius Base
vkr_base = function(df.points, args.list) {
  vkr_base_no_list(
    df.points = df.points,
    n.neighbors = args.list$n_neighbors,
    radius.mult = args.list$radius_mult,
    max.knots = args.list$max_knots,
    cols.to.sort = args.list$cols_to_sort
  )
}

# VKR Grid Search -------------------------------------------------------------------------------------------------

eval_knots_mse = function(df.points, df.knots, gam.k=3) {
  gam.eval = gam(signal ~ te(x, y, k=gam.k, bs="gp"),data=df.knots,method="REML", family=gaussian)
  df.points[, "pkn"] = predict(gam.eval, newdata=df.points)
  (df.points$signal - df.points[, "pkn"])^2 %>% mean  # mse
}

vkr_gs = function(df.points, seq.nn, seq.rm, n.knots, cols.to.sort=c("entropy"), gam.k=3) {
  results = data.frame(nn=integer(), rm=double(), mse=double())
  for (nn in seq.nn) {
    for (rm in seq.rm) {
      knots = vkr_base(df.points, list(n_neighbors=nn, radius_mult=rm, max_knots=n.knots, cols_to_sort=cols.to.sort))
      # paste0("Currently Doing: knots=", knots %>% nrow, ", nn=", nn, ", rm=", rm) %>% print
      mse = eval_knots_mse(df.points, knots, gam.k)
      results = results %>% rbind(data.frame(nn=nn, rm=rm, mse=mse))
    }
  }
  best.mse = results %>% filter(mse == min(mse))
  best.list.args = list(n_neighbors=best.mse$nn, radius_mult=best.mse$rm, max_knots=n.knots, cols_to_sort=cols.to.sort)
  list(ls.args=best.list.args, grid=results)
}