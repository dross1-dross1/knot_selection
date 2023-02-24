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
vkr.gs = vkr_gs(df.points=test.data, seq.nn=1:20, seq.rm=seq(from=3, to=4.5, by=.25), n.knots=n.knots, cols.to.sort=c("entropy"), gam.k=3)
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

# Grid Search -----------------------------------------------------------------------------------------------------

gs_entropy = function(seq.nn, seq.rm, n.knots, cols.to.sort=c("entropy")) {
  results = data.frame(nn=integer(), rm=double(), mse=double())
  for (nn in seq.nn) {
    for (rm in seq.rm) {
      knots = vkr_base(test.data, list(n_neighbors=nn, radius_mult=rm, max_knots=n.knots, cols_to_sort=cols.to.sort))
      # paste0("Currently Doing: knots=", knots %>% nrow, ", nn=", nn, ", rm=", rm) %>% print
      mse = eval_knots_mse(test.data, knots)
      results = results %>% rbind(data.frame(nn=nn, rm=rm, mse=mse))
    }
  }
  results
}
gs.entropy = gs_entropy(seq.nn=1:5, seq.rm=seq(from=1, to=3, by=.25), 10)
gs.entropy[order(gs.entropy$mse), ][1, "mse"]
plot_ly(gs.entropy, x=~nn, y=~rm, z=~mse, color=~mse) %>% add_markers

cd_mse = function(n.knots, n.iter) {
  mses = c()
  for (i in 1:n.iter) {
    knots = generate_cd_knots(test.data, n.knots)
    # paste0("Currently Doing: i=", i) %>% print
    mses = c(mses, eval_knots_mse(test.data, knots))
  }
  mses
}
cd.mse = cd_mse(n.knots=25, n.iter=25)
cd.mse
cd.mse %>% mean

# Automation ------------------------------------------------------------------------------------------------------

auto_mse = function(seq.n.knots) {
  results = data.frame(nk=integer(), entropy=double(), cd=double())
  for (i.n.knots in seq.n.knots) {
    paste0("Currently Doing: n_knots=", i.n.knots) %>% print
    gs.ent = gs_entropy(seq.nn=1:5, seq.rm=seq(from=1, to=3, by=.25), i.n.knots)
    mse.entropy = gs.ent[order(gs.ent$mse), ][1, "mse"]
    mse.cd = cd_mse(n.knots=i.n.knots, n.iter=10)
    results = results %>% rbind(data.frame(nk=i.n.knots, entropy=mse.entropy, cd=mse.cd))
    print(results)
  }
  results
}
test.mse = auto_mse(seq(from=10, to=100, by=10))
test.mse

# month = 1:12

# gam_k = 3
# n.knots = seq(from=5, to=20, by=5)

# gam_k = 4
# n.knots = seq(from=20, to=40, by=10)

eval_knots_automated = function(seq.month, seq.n.knots, gam.k) {
  # results = data.frame(month=integer(), n_knots=integer(), all=double(), entropy=double(), random=double(), c_design=double(), tbart=double())
  results = data.frame()
  for (i.month in seq.month) {
    for (i.n.knots in seq.n.knots) {
      print(paste0("Currently Doing: month=", i.month, ", n_knots=", i.n.knots))
      # get data
      df.epa = get_epa(i.month)
      # k.entropy
      args.list = list(n_neighbors=10, lambda=0, radius_mult=NA, max_knots=1000, cols_to_sort=c("dist_from_median"))
      k.entropy = vkr_full(df.epa, args.list, i.n.knots)
      # k.cover.design
      k.cover.design = generate_cd_knots(df.epa, i.n.knots)
      # eval knots
      df.epa = eval_knots(df.epa, df.epa, "all_data_preds", gam.k)
      df.epa = eval_knots(df.epa, k.entropy, "entropy_knot_preds", gam.k)
      df.epa = eval_knots(df.epa, k.cover.design, "cd_knot_preds", gam.k)
      k.metrics = eval_knots_metrics(df.epa)
      tkm_results(k.metrics)
      results = results %>% rbind(data.frame(month=i.month, n_knots=i.n.knots, gam_k=gam.k) %>% cbind(k.metrics))
    }
  }
  results
}

plot_knot_eval = function(df.k.results, target.month, type="mse") {
  # digest results: mse (lower better)
  if (type == "mse") {
    ggplot(df.k.results %>% filter(month == target.month)) + theme_dark() +
      labs(title="GAM Model Accuracy (MSE) by Knot Selection Type", subtitle=paste0("Month = ", target.month), caption="caption") +
      geom_line(aes(x=n_knots, y=all_data_preds_mse, color="all_data_preds_mse"), size=2) +
      geom_line(aes(x=n_knots, y=entropy_knot_preds_mse, color="entropy_knot_preds_mse"), size=2) +
      geom_line(aes(x=n_knots, y=cd_knot_preds_mse, color="cd_knot_preds_mse"), size=2)
  } else if (type == "cor2") {
  # digest results: cor2 (higher better)
  ggplot(df.k.results %>% filter(month == target.month)) + theme_dark() +
    labs(title="GAM Model Accuracy (cor^2) by Knot Selection Type", subtitle=paste0("Month = ", target.month), caption="caption") +
    geom_line(aes(x=n_knots, y=all_data_preds_cor2, color="all_data_preds_cor2"), size=2) +
    geom_line(aes(x=n_knots, y=entropy_knot_preds_cor2, color="entropy_knot_preds_cor2"), size=2) +
    geom_line(aes(x=n_knots, y=cd_knot_preds_cor2, color="cd_knot_preds_cor2"), size=2)
  } else { print("Invalid knot eval type!") }
}

# knot.results = eval_knots_automated(seq.month=1:12, seq.n.knots=c(10, 15, 20, 25), gam.k=3)
# knot.results = eval_knots_automated(seq.month=1:12, seq.n.knots=10:20, gam.k=3)
knot.results = eval_knots_automated(seq.month=c(2, 5, 6, 10), seq.n.knots=seq(from=10, to=100, by=10), gam.k=3)
# knot.results = eval_knots_automated(seq.month=c(1, 6, 7, 12), seq.n.knots=c(20, 25, 30), gam.k=4)
write.csv(knot.results, paste0("kresults_", Sys.time() %>% str_replace_all("[:-]+", ".") %>% str_replace_all("\\s", "_"), "_", nrow(knot.results), "rows", ".csv"), row.names=F)
knot.results %>% filter(n_knots == 10)

plot_knot_eval(knot.results, target.month=10, type="mse")

# load previous runs
kr.prev = read.csv("kresults_2023.02.13_08.54.50_132rows.csv")
for (i in 1:12) { plot_knot_eval(kr.prev, target.month=i, type="mse") %>% print }
