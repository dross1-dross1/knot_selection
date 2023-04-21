# Setup -----------------------------------------------------------------------------------------------------------

# Set WD to directory of current script location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# UI and misc
if (!is.null(dev.list()["RStudioGD"])) { dev.off(dev.list()["RStudioGD"]) } # clear plots
rm(list = ls())                                                             # remove all variables
# set.seed(100)                                                               # reproducable randomness
cat("\014")                                                                 # ctrl+L (clear console output)

# Imports (Graphing)
library(ggplot2)  # graphing; ggplot() + ...
library(ggforce)  # graphing circles; geom_circle()
library(plotly)   # 3d scatter plots; plot_ly()

# Imports (Knot Selection)
library(tbart)    # for Teitz and Bart's p-Median algorithm; allocate()
library(nebcs)    # for cover.design knot selection; make.knots()

# Imports (Modeling)
library(mgcv)     # evaluating knots via GAMs; gam(), vis.gam()

# Imports
library(dplyr)    # piping and QoL functions; %>%
library(scales)   # normalizing lists between 0 and 1; rescale()
library(stringr)  # manipulating strings

# World Generation ------------------------------------------------------------------------------------------------

# EPA data
get_epa = function(month) {
  epa.data = readRDS("df_data_12list.RDS")[[month]]
  epa.data = epa.data %>% subset(select = c(x, y, health))
  names(epa.data) = c("x", "y", "signal")
  epa.data$signal = epa.data$signal %>% rescale
  epa.data
}
test.data = get_epa(month=5)

# view world
ggplot(test.data) + geom_density(aes(x=signal)) + theme_dark()
view_world = function(df.points) {
  ggplot(df.points) + geom_point(aes(x=x, y=y, color=signal)) + scale_color_gradient2(midpoint=.5) + theme_dark() +
    xlim(c(0, 1)) + ylim(c(0, 1)) + labs(title="World", subtitle="subtitle", caption="caption")
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
source("ka_v5_functions_parallel.R")
knots.entropy = entropy_max(df.points=test.data, n.neighbors=5, radius.mult=2, max.knots=25)
view_knots_r(knots.entropy, test.data, "Entropy Maximization")
view_knots(test.data, knots.entropy, "Entropy Maximization")

# view grid
ggplot(emax.gs$grid, aes(x=nn, y=rm)) + geom_tile(aes(fill=mse), colour="black") + scale_fill_gradient(low="blue", high="white")

# Functions (Knots Evaluation) ------------------------------------------------------------------------------------

eval_knots = function(df.points, df.knots, pred.knots.name="default_pred_knot_name", k=3, vis=F) {
  # gam fitting
  gam.eval = gam(
    signal ~ te(x, y, k=k, bs="gp"),  # , k=12, bs="so"
    data=df.knots,
    method="REML",
    family=gaussian
  )
  if (vis == T) { vis.gam(gam.eval, type="response", plot.type="contour", color="gray") }  # gam visualization
  df.points[, pred.knots.name] = predict(gam.eval, newdata=df.points)  # predictions
  df.points
}

eval_knots_mse = function(df.points, df.knots, gam.k=3) {
  gam.eval = gam(signal ~ te(x, y, k=gam.k, bs="gp"),data=df.knots,method="REML", family=gaussian)
  df.points[, "pkn"] = predict(gam.eval, newdata=df.points)
  (df.points$signal - df.points[, "pkn"])^2 %>% mean  # mse
}

eval_knots_metrics = function(df.points.evaled) {
  result.metrics = list()
  n.knot.types = (names(df.points.evaled) %>% length) - 3
  for(i.knot.type in 1:n.knot.types) {
    i.name = names(df.points.evaled)[3+i.knot.type]
    result.metrics[paste0(i.name, "_mse")] = (df.points.evaled$signal - df.points.evaled[, i.name])^2 %>% mean
    result.metrics[paste0(i.name, "_cor2")] = cor(df.points.evaled$signal, df.points.evaled[, i.name])^2
  }
  result.metrics
}

# Testing ---------------------------------------------------------------------------------------------------------

gam_k = 3
test.data = eval_knots(test.data, test.data, "all_data_preds", gam_k, vis=T)
test.data = eval_knots(test.data, knots.entropy, "entropy_knot_preds", gam_k, vis=T)
test.data = eval_knots(test.data, knots.cover.design, "cd_knot_preds", gam_k, vis=T)
test.knot.metrics = eval_knots_metrics(test.data)
tkm_results = function(kmetrics) {
  print(paste0("all:          ", kmetrics$all_data_preds_mse))
  print(paste0("entropy:      ", kmetrics$entropy_knot_preds_mse))
  print(paste0("cover.design: ", kmetrics$cd_knot_preds_mse))
}
tkm_results(test.knot.metrics)

alpha = .3
size = 3
ggplot(test.data) + geom_abline(slope=1, intercept=0) + theme_dark() +
  geom_point(aes(x=signal, y=all_data_preds, color="all"), alpha=alpha, size=size) +
  geom_point(aes(x=signal, y=entropy_knot_preds, color="entropy"), alpha=alpha, size=size)
ggplot(test.data) + geom_abline(slope=1, intercept=0) + theme_dark() +
  geom_point(aes(x=signal, y=all_data_preds, color="all"), alpha=alpha, size=size) +
  geom_point(aes(x=signal, y=cd_knot_preds, color="cd"), alpha=alpha, size=size)

# Automation ------------------------------------------------------------------------------------------------------

