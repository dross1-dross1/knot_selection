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

# Imports
library(dplyr)    # piping and QoL functions; %>%
library(scales)   # normalizing vectors between 0 and 1; rescale()

# World Generation ------------------------------------------------------------------------------------------------

# EPA data
get_epa = function(month) {
  epa.data = readRDS("data/df_data_12list.RDS")[[month]]
  epa.data = epa.data %>% subset(select = c(x, y, health))
  names(epa.data) = c("x", "y", "signal")
  epa.data$signal = epa.data$signal %>% rescale
  epa.data
}
test.data = get_epa(month=5)

# a 3d ellipse will need:
# x0, y0, z0, a, b, c, angle

# view world
ggplot(test.data) + geom_density(aes(x=signal)) + theme_dark()
view_world = function(df.points) {
  ellipse_center = c(0.5, 0.5)
  ellipse_axes = c(0.3, 0.2)
  ellipse_angle = pi / 6
  ggplot(df.points) + geom_point(aes(x=x, y=y, color=signal)) + scale_color_gradient2(midpoint=.5) + theme_dark() +
    xlim(c(0, 1)) + ylim(c(0, 1)) + labs(title="World", subtitle="subtitle", caption="caption") +
    geom_ellipse(aes(x0=ellipse_center[1], y0=ellipse_center[2], a=ellipse_axes[1], b=ellipse_axes[2], angle=ellipse_angle))
}
view_world(test.data)
test.data %>% summary
test.data %>% nrow

# Knot Selection Methods ------------------------------------------------------------------------------------------

view_knots = function(df.points, df.knots, title="Default Title") {
  ggplot(df.points) + geom_point(aes(x=x, y=y, color=signal)) + scale_color_gradient2(midpoint=.5) +
    geom_point(data=df.knots, aes(x=x, y=y), size=5, alpha=.5) + theme_dark() +
    labs(title=title, subtitle="subtitle", caption="caption")
}

view_knots_r = function(df.knots, df.points, title="Default Title") {
  ggplot(df.points) + geom_point(aes(x=x, y=y, color=signal)) + scale_color_gradient2(midpoint=.5) + theme_dark() +
  geom_circle(data=df.knots, aes(x0=x, y0=y, r=radius)) +
  labs(title=title, subtitle="subtitle", caption="caption")
}

n.knots = 25

#### Knot Selection: Entropy Maximization ####
test.data$z = 1
source("kaf_v6_sphere.R")
knots.entropy = entropy_max(df.points=test.data, n.neighbors=5, radius.mult=2, max.knots=25)
view_knots_r(knots.entropy, test.data, "Entropy Maximization")
view_knots(test.data, knots.entropy, "Entropy Maximization")

# Testing Ellipsoid Code ------------------------------------------------------------------------------------------

# uncomment this source statement when running the ellipsoid code below
# source("kaf_v6_ellipsoid.R")

generate_test_dataframe = function(N) {
  # Randomly generate values for each column
  a = runif(N, min = 1, max = 5)
  b = runif(N, min = 1, max = 5)
  c = runif(N, min = 1, max = 5)
  x0 = runif(N, min = -5, max = 5)
  y0 = runif(N, min = -5, max = 5)
  z0 = runif(N, min = -5, max = 5)
  alpha = runif(N, min = 0, max = 360)
  beta = runif(N, min = 0, max = 360)
  gamma = runif(N, min = 0, max = 360)
  signal = runif(N, min = 0, max = 1)

  # Create a dataframe with the generated values
  df = data.frame(a, b, c, x0, y0, z0, alpha, beta, gamma, signal)

  return(df)
}

# Example usage
N = 100
test_df = generate_test_dataframe(N)

thingy1 = gkm_ellipsoid(test_df)
thingy1 %>% head(15)

thingy2 = gkm_ellipsoid_parallel(test_df)
thingy2 %>% head(15)
