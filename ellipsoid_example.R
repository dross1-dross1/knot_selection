# Setup -----------------------------------------------------------------------------------------------------------

# Set WD to directory of current script location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# UI and misc
if (!is.null(dev.list()["RStudioGD"])) { dev.off(dev.list()["RStudioGD"]) } # clear plots
rm(list = ls())                                                             # remove all variables
set.seed(100)                                                               # reproducable randomness
cat("\014")                                                                 # ctrl+L (clear console output)

# Imports (Graphing)
library(ggplot2)  # general; ggplot() + ...
library(ggforce)  # circles; geom_circle()
library(rgl)      # ellipsoids; plot3d()

# Imports
library(dplyr)    # piping and QoL functions; %>%
library(scales)   # normalizing vectors between 0 and 1; rescale()

# Code ------------------------------------------------------------------------------------------------------------

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

# Testing Code ----------------------------------------------------------------------------------------------------

# Ellipsoid parameters
a = 3
b = 2
c = 1
x0 = 0
y0 = 0
z0 = 0
alpha = 13
beta = 525
gamma = 16

# Test points
point_inside = c(1, 0, 0) # A point inside the ellipsoid
point_on_surface = c(3, 0, 0) # A point on the surface of the ellipsoid
point_outside = c(4, 0, 0) # A point outside the ellipsoid

# More test points
p_inside_1 = c(0, 0, 0)
p_inside_2 = c(2, 1, 0)
p_inside_3 = c(-2, -1, 0)
p_outside_1 = c(3, 2, 1)
p_outside_2 = c(-3, -2, -1)
p_outside_3 = c(0, 3, 0)

# Function calls
inside = is_point_inside_ellipsoid(a, b, c, x0, y0, z0, alpha, beta, gamma, point_inside[1], point_inside[2], point_inside[3])
on_surface = is_point_inside_ellipsoid(a, b, c, x0, y0, z0, alpha, beta, gamma, point_on_surface[1], point_on_surface[2], point_on_surface[3])
outside = is_point_inside_ellipsoid(a, b, c, x0, y0, z0, alpha, beta, gamma, point_outside[1], point_outside[2], point_outside[3])

# Print results
cat("Inside:", inside, "\n") # Should return TRUE
cat("On surface:", on_surface, "\n") # Should return TRUE
cat("Outside:", outside, "\n") # Should return FALSE

# More results
print("Inside:", is_point_inside_ellipsoid(a, b, c, x0, y0, z0, alpha, beta, gamma, p_inside_1[1], p_inside_1[2], p_inside_1[3]))
print("Inside:", is_point_inside_ellipsoid(a, b, c, x0, y0, z0, alpha, beta, gamma, p_inside_2[1], p_inside_2[2], p_inside_2[3]))
print("Inside:", is_point_inside_ellipsoid(a, b, c, x0, y0, z0, alpha, beta, gamma, p_inside_3[1], p_inside_3[2], p_inside_3[3]))
print("Outside:", is_point_inside_ellipsoid(a, b, c, x0, y0, z0, alpha, beta, gamma, p_outside_1[1], p_outside_1[2], p_outside_1[3]))
print("Outside:", is_point_inside_ellipsoid(a, b, c, x0, y0, z0, alpha, beta, gamma, p_outside_2[1], p_outside_2[2], p_outside_2[3]))
print("Outside:", is_point_inside_ellipsoid(a, b, c, x0, y0, z0, alpha, beta, gamma, p_outside_3[1], p_outside_3[2], p_outside_3[3]))

# Randomly Generated Points ---------------------------------------------------------------------------------------

random_points <- matrix(runif(1000*3, min = -3, max = 3), ncol = 3)

# Function to get color for each point based on whether it's inside or outside the ellipsoid
get_point_color <- function(point) {
  if (is_point_inside_ellipsoid(a, b, c, x0, y0, z0, alpha, beta, gamma, point[1], point[2], point[3])) {
    return("green")
  } else {
    return("red")
  }
}

# Get the colors for the random points
point_colors <- apply(random_points, 1, get_point_color)

# Plot the random points with corresponding colors
points3d(random_points, col = point_colors, size = 5)

# Axes labels
axes3d(col = "black")
title3d(xlab = "X", ylab = "Y", zlab = "Z")
