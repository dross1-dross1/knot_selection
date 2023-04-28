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

#' @title Check if a Point is Inside a 3D Ellipsoid
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

#' @title Scale a 3D Ellipsoid
#' @description Scales a specified ellipsoid by a scalar value.
#' @param a A numeric scalar for the semi-major axis length of the ellipsoid.
#' @param b A numeric scalar for the semi-intermediate axis length of the ellipsoid.
#' @param c A numeric scalar for the semi-minor axis length of the ellipsoid.
#' @param scale_factor A numeric scalar for the scaling factor.
#' @return A named list containing the scaled ellipsoid semi-axes.
#' @example
#' scaled_axes = scale_ellipsoid_axes(3, 2, 1, 2)
#' print(scaled_axes$a)  # 6
#' print(scaled_axes$b)  # 4
#' print(scaled_axes$c)  # 2
scale_ellipsoid_axes = function(a, b, c, scale_factor) {
  # Scale the ellipsoid axes
  a_scaled = a * scale_factor
  b_scaled = b * scale_factor
  c_scaled = c * scale_factor

  # Return the scaled ellipsoid axes as a named list
  return(list(a = a_scaled, b = b_scaled, c = c_scaled))
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

# spatial dataframe: A dataframe that contains the columns `x`, `y`, and `z`.

# for every row in the world:
# find n nearest neighbors
# each iteration should generate the modified row of the old world

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
    df.points = df.points %>% anti_sphere_subset(df.knots[i, "x"], df.knots[i, "y"], df.knots[i, "z"], df.knots[i, "radius"] * radius.mult)  # remove the points "too close" to the appended knot
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
