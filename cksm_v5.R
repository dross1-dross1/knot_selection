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
# library(ggdark)   # dark theme; dark_theme_gray()

# Imports (Knot Selection)
library(tbart)    # for Teitz and Bart's p-Median algorithm; allocate()
library(nebcs)    # for cover.design knot selection; make.knots()

# Imports (Modeling)
library(mgcv)     # evaluating knots via GAMs; gam(), vis.gam()

# Imports
library(dplyr)    # piping and QoL functions; %>%
library(scales)   # normalizing lists between 0 and 1; rescale()
library(stringr)  # manipulating strings

# Knot Evaluation
source("fuentes_GP_model.R")

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

n.knots = 25

view_knots = function(df.points, df.knots, title="Default Title") {
  ggplot(df.points) + geom_point(aes(x=x, y=y, color=signal)) + scale_color_gradient2(midpoint=.5) +
    xlim(c(0, 1)) + ylim(c(0, 1)) +
    geom_point(data=df.knots, aes(x=x, y=y), size=5, alpha=.5) + theme_dark() +
    labs(title=title, subtitle="subtitle", caption="caption")
}

view_knots_r = function(df.knots, df.points, title="Default Title") {
  ggplot(df.points) + geom_point(aes(x=x, y=y, color=signal)) + scale_color_gradient2(midpoint=.5) + theme_dark() +
  geom_circle(data=df.knots, aes(x0=x, y0=y, r=radius)) +
  labs(title=title, subtitle="subtitle", caption="caption")
}

#### Knot Selection: Entropy Maximization ####
source("ka_v5_functions.R")
# knots.entropy = vkr_base(test.data, list(n_neighbors=10, radius_mult=4, max_knots=n.knots, cols_to_sort=c("entropy")))
vkr.gs = vkr_gs(df.points=test.data, seq.nn=1:10, seq.rm=seq(from=2, to=3, by=.25), n.knots=n.knots, cols.to.sort=c("entropy"), gam.k=3)
knots.entropy = vkr_base(test.data, vkr.gs$ls.args)
view_knots_r(knots.entropy, test.data, "Entropy Maximization")
view_knots(test.data, knots.entropy, "Entropy Maximization")

#### Knot Selection: Cover.Design ####
generate_cd_knots = function(df.points, N_k) {
  df.points.cd = data.frame(df.points) %>% subset(select=c(x, y, signal))
  names(df.points.cd) = c("X_COORD", "Y_COORD", "signal")
  cd.knots.prep = make_knots(df.points.cd, N_k)
  test.data[cd.knots.prep %>% rownames %>% as.numeric, ] %>% subset(select=c(x, y, signal))
}
knots.cover.design = generate_cd_knots(test.data, n.knots)
view_knots(test.data, knots.cover.design, "Cover.Design")

# Functions (Knots Evaluation) ------------------------------------------------------------------------------------

eval_knots = function(df.points, df.knots, pred.knots.name="default_pred_knot_name", k=3, vis=F) {
  # gam fitting
  gam.eval = gam(
    signal ~ te(x, y, k=k, bs="gp"),  # , k=12, bs="so"
    data=df.knots,
    method="REML",
    family=gaussian
  )
  # gam visualization
  if (vis == T) { vis.gam(gam.eval, type="response", plot.type="contour", color="gray") }
  # predictions
  df.points[, pred.knots.name] = predict(gam.eval, newdata=df.points)
  # mean squared error; lower is better
  # mse = (df.points$signal - df.points[, pred.knots.name])^2 %>% mean
  # print(paste0("Mean Squared Error (lower better): ", mse))
  # correlation squared; higher is better
  # cor2 = cor(df.points$signal, df.points[, pred.knots.name])^2
  # print(paste0("Correlation Squared (higher better): ", cor2))
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

# Fuentes Testing -------------------------------------------------------------------------------------------------
source("fuentes_GP_model.R")

# Training data with X columns x and y, and y column pm2_5
td = test.data[c("x", "y", "signal")]
names(td) = c("x", "y", "pm2_5")
td.train = td %>% sample_frac(.7)
td.test = td %>% anti_join(td.train)

td.X.train = td.train[, c("x", "y")]
td.y.train = td.train[, c("pm2_5")]

td.X.test = td.test[, c("x", "y")]
td.y.test = td.test[, c("pm2_5")]

# entropy fitting
params.entropy = fit.fuentes.gp(knots.entropy, td.X.train, td.y.train)
preds.df.entropy = predict.fuentes.gp(test.knots, td.X.test, params.entropy)

# Automation ------------------------------------------------------------------------------------------------------

cd_mse = function(n.knots, n.iter) {
  mses = c()
  for (i in 1:n.iter) {
    knots = generate_cd_knots(test.data, n.knots)
    # paste0("Currently Doing: i=", i) %>% print
    mses = c(mses, eval_knots_mse(test.data, knots))
  }
  mses
}

# TODO: create train-test split for more accurate MSEs
# TODO: add random knot selection and Fuentes knot seleciton
# TODO: use Fuentes model and compare to GAMs

auto_mse = function(seq.n.knots) {
  results = data.frame(nk=integer(), entropy=double(), cd=double())
  for (i.n.knots in seq.n.knots) {
    paste0("Currently Doing: n_knots=", i.n.knots) %>% print
    # gs.ent = gs_entropy(seq.nn=1:5, seq.rm=seq(from=1, to=3, by=.25), i.n.knots)
    gs.ent = vkr_gs(df.points=test.data, seq.nn=1:10, seq.rm=seq(from=2, to=3, by=.25), n.knots=n.knots, cols.to.sort=c("entropy"), gam.k=3)
    mse.entropy = gs.ent$grid$mse %>% min
    mse.cd = cd_mse(n.knots=i.n.knots, n.iter=10) %>% mean
    results = results %>% rbind(data.frame(n_k=i.n.knots, mse_entropy=mse.entropy, mse_cd=mse.cd))
    print(results)
  }
  results
}
test.mse = auto_mse(seq(from=10, to=20, by=10))
test.mse

ggplot(test.mse) + theme_dark() +
  labs(title="MSE for each knot selection method", subtitle="subtitle", caption="caption") +
  geom_line(aes(x=n_k, y=mse_entropy, color="mse_entropy"), size=2) + geom_line(aes(x=n_k, y=mse_cd, color="mse_cd"), size=2)

ggplot() + theme_dark() +
  labs(title="MSE for each knot selection method", subtitle="subtitle", caption="caption") +
  geom_line(data=gs.entropy, aes(x=nn, y=mse, color=as.factor(rm)), size=2) + geom_line(data=test.mse, aes(x=n_k, y=mse_cd, color="mse_cd"), size=2)

