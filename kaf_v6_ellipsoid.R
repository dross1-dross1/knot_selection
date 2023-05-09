# Knot Algorithm Functions Version 6 for Ellipsoids

# Setup -----------------------------------------------------------------------------------------------------------

# Imports
library(dplyr)    # piping and QoL functions; %>%
library(tidyr)    # for dropping NA valuse; drop_na()
library(scales)   # normalizing lists between 0 and 1; rescale()

# Functions (Formulas and Math) -----------------------------------------------------------------------------------

#' @title Check if a Point is Inside a Ellipsoid
#' @description Determines if a given point is inside, on the surface, or outside a specified ellipsoid.
#' @param a A numeric scalar for the semi-major axis length of the ellipsoid.
#' @param b A numeric scalar for the semi-intermediate axis length of the ellipsoid.
#' @param c A numeric scalar for the semi-minor axis length of the ellipsoid.
#' @param x0 A numeric scalar for the x coordinate of the center of the ellipsoid.
#' @param y0 A numeric scalar for the y coordinate of the center of the ellipsoid.
#' @param z0 A numeric scalar for the z coordinate of the center of the ellipsoid.
#' @param alpha A numeric scalar for the rotation angle around the x-axis in degrees.
#' @param beta A numeric scalar for the rotation angle around the y-axis in degrees.
#' @param gamma A numeric scalar for the rotation angle around the z-axis in degrees.
#' @param px A numeric scalar for the x coordinate of the test point.
#' @param py A numeric scalar for the y coordinate of the test point.
#' @param pz A numeric scalar for the z coordinate of the test point.
#' @return A boolean value indicating whether the point is inside the ellipsoid (TRUE) or not (FALSE).
#' @example
#' is_point_inside_ellipsoid(3, 2, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0)  # Inside (TRUE)
#' is_point_inside_ellipsoid(3, 2, 1, 0, 0, 0, 0, 0, 0, 3, 0, 0)  # On surface (TRUE)
#' is_point_inside_ellipsoid(3, 2, 1, 0, 0, 0, 0, 0, 0, 4, 0, 0)  # Outside (FALSE)
#' is_point_inside_ellipsoid(3, 2, 1, 0, 0, 0, 30, 45, 60, 0.5, 0.5, 0.5)  # Inside with rotation (TRUE)
#' is_point_inside_ellipsoid(5, 3, 2, 2, -1, 1, 0, 0, 0, 6, -1, 1)  # On surface with translation (TRUE)
#' is_point_inside_ellipsoid(3, 2, 1, 0, 0, 0, 30, 45, 60, 4, 1, 1)  # Outside with rotation (FALSE)
#' is_point_inside_ellipsoid(3, 2, 1, 0, 0, 0, 30, 45, 60, -1, -3, 1)  # Outside with rotation (FALSE)
is_point_inside_ellipsoid = function(a, b, c, x0, y0, z0, alpha, beta, gamma, px, py, pz) {
  # Convert angles from degrees to radians
  alpha = alpha * (pi / 180)
  beta = beta * (pi / 180)
  gamma = gamma * (pi / 180)

  # Calculate the rotation matrix
  R = matrix(
    c(
      cos(alpha) * cos(beta),
      sin(alpha) * cos(beta),
      -sin(beta),
      cos(alpha) * sin(beta) * sin(gamma) - sin(alpha) * cos(gamma),
      sin(alpha) * sin(beta) * sin(gamma) + cos(alpha) * cos(gamma),
      cos(beta) * sin(gamma),
      cos(alpha) * sin(beta) * cos(gamma) + sin(alpha) * sin(gamma),
      sin(alpha) * sin(beta) * cos(gamma) - cos(alpha) * sin(gamma),
      cos(beta) * cos(gamma)
    ),
    nrow = 3,
    ncol = 3,
    byrow = TRUE
  )

  # Translate the point
  pt = c(px - x0, py - y0, pz - z0)

  # Rotate the point
  pt_rot = solve(R) %*% pt

  # Check if the rotated and translated point is inside the ellipsoid
  inside = (pt_rot[1]^2 / a^2) + (pt_rot[2]^2 / b^2) + (pt_rot[3]^2 / c^2) <= 1

  return(inside)
}

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
  if (length(signal) == 0) { return (NA) }
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

#' @title Generate Knot Metrics for Ellipsoids in Parallel
gkm_ellipsoid_parallel = function(df.points) {
  # Load necessary libraries
  library(doParallel)
  library(foreach)

  # Initialize parallel backend
  n.cores = parallel::detectCores() - 1
  my.cluster = parallel::makeCluster(n.cores, type="PSOCK")
  registerDoParallel(my.cluster)

  # Export helper functions to parallel workers
  parallel::clusterExport(my.cluster, c("is_point_inside_ellipsoid", "ebrahimi_entropy", "pad_edges"))

  # Load required packages within parallel workers
  parallel::clusterEvalQ(my.cluster, {
    library(dplyr)
    library(tidyr)
    library(scales)
  })

  # Extract only the relevant columns
  metric.df = data.frame(df.points) %>% subset(select=c(a, b, c, x0, y0, z0, alpha, beta, gamma, signal))

  # Use foreach and %dopar% to parallelize the loop
  results = foreach(i=1:nrow(df.points), .combine="rbind", .multicombine=T) %dopar% {
    point_is_inside = apply(df.points, 1, function(row) {
      is_point_inside_ellipsoid(
        a = df.points[i, "a"],
        b = df.points[i, "b"],
        c = df.points[i, "c"],
        x0 = df.points[i, "x0"],
        y0 = df.points[i, "y0"],
        z0 = df.points[i, "z0"],
        alpha = df.points[i, "alpha"],
        beta = df.points[i, "beta"],
        gamma = df.points[i, "gamma"],
        px = row["x0"],
        py = row["y0"],
        pz = row["z0"]
      )
    })
    df.points.subset = df.points[point_is_inside, ]  # These are the subset of points inside the ellipsoid
    entropy = df.points.subset$signal %>% ebrahimi_entropy  # find the entropy of the point
  }

  # Stop cluster
  parallel::stopCluster(my.cluster)

  # Add the results to metric.df
  metric.df$entropy = results[, 1]

  metric.df
}

#' @title Generate Knot Metrics for Ellipsoids
gkm_ellipsoid = function(df.points) {
  metric.df = data.frame(df.points) %>% subset(select=c(a, b, c, x0, y0, z0, alpha, beta, gamma, signal))  # extract only the relavent columns
  for (i in 1:nrow(df.points)) {
    point_is_inside = apply(df.points, 1, function(row) {
      is_point_inside_ellipsoid(
        a = df.points[i, "a"],
        b = df.points[i, "b"],
        c = df.points[i, "c"],
        x0 = df.points[i, "x0"],
        y0 = df.points[i, "y0"],
        z0 = df.points[i, "z0"],
        alpha = df.points[i, "alpha"],
        beta = df.points[i, "beta"],
        gamma = df.points[i, "gamma"],
        px = row["x0"],
        py = row["y0"],
        pz = row["z0"]
      )
    })
    df.points.subset = df.points[point_is_inside, ]  # These are the subset of points inside the ellipsoid
    metric.df[i, "entropy"] = df.points.subset$signal %>% ebrahimi_entropy  # find the entropy of the point
  }
  metric.df
}

#' @title Anti Ellipsoid Subset
anti_ellipsoid_subset = function(df, a, b, c, x0, y0, z0, alpha, beta, gamma) {
  indices_to_keep = c()

  for (i in 1:nrow(df)) {
    row = df[i,]
    px = row$x0
    py = row$y0
    pz = row$z0

    if (!is_point_inside_ellipsoid(a, b, c, x0, y0, z0, alpha, beta, gamma, px, py, pz)) {
      indices_to_keep = c(indices_to_keep, i)
    }
  }

  return(df[indices_to_keep, ])
}

#' @title Entropy Maximization
#' @description Given a spatial dataframe with a `signal` column that represents some continuous value, this function will return knots based on the other parametrs.
#' An iterative way of determining knot candidates in relation to their distance from each other.
#' @param df.points A spatial dataframe containing the column `signal`.
#' @param scale.factor The value to multiply each point's ? by when removing other points from the list of knot candidates.
#' @param max.knots The maximum number of iterations to pick knots (not guaranteed to always reach the max).
#' @return The final list of knots.
entropy_max = function(df.points, scale.factor, max.knots) {
  df.knots = data.frame()
  for(i in 1:max.knots) {
    if (nrow(df.points) == 0) { break }
    df.points = gkm_ellipsoid(df.points)  # recompute entropy every iteration
    df.points = df.points[order(df.points$entropy, decreasing=T), ]  # sort df by entropy
    df.knots = df.knots %>% rbind(df.points[1, ])  # append knot with highest entropy to list of knots

    # remove the points "too close" to the appended knot
    df.points = df.points %>% anti_ellipsoid_subset(
      a=df.knots[i, "a"] * scale.factor,
      b=df.knots[i, "b"] * scale.factor,
      c=df.knots[i, "c"] * scale.factor,
      x0=df.knots[i, "x0"],
      y0=df.knots[i, "y0"],
      z0=df.knots[i, "z0"],
      alpha=df.knots[i, "alpha"],
      beta=df.knots[i, "beta"],
      gamma=df.knots[i, "gamma"]
    )
  }
  return(df.knots %>% drop_na)
}

# entropy_max(df.points = test_df, scale.factor = 1, max.knots = 50)
