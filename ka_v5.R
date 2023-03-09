# Setup -------------------------------------------------------------------

#if (!is.null(dev.list()["RStudioGD"])) { dev.off(dev.list()["RStudioGD"]) }
#rm(list=ls())                     # remove all variables
cat("\014")                       # ctrl+L
set.seed(100)

# Imports
library(ggplot2)
library(ggforce)
library(readr)      # don't use tidyverse
library(dplyr)
# library(bayesplot)

# Set WD to directory of current script location
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# TODO: Avoid groupings of knots where the area is very homogeneous (meaning the neighborhood of each knot essentially looks identical)
# TODO: Adaptive knot selection based on density of population so that sparsely populated areas aren't ignored
# TODO: Try using data tables as it is more efficient than data frames

# Functions (Base Algorithm) ----------------------------------------------

# in: X (numeric vector) & m (amount of padding to add)
# out: padded vector for computing the rolling window difference
pad_along_last_axis = function(X, m) {
  first = X[1]
  last = X[length(X)]
  c(rep(first, each=m), X, rep(last, each=m))
}

# signal: (numeric vector) list of values used as the signal to compute the Ebrahimi differential entropy
ebrahimi_entropy = function(signal) {
  sorted_signal = sort(signal)
  n = length(signal)
  m = floor(sqrt(n) + 0.5)
  X = pad_along_last_axis(sorted_signal, m)
  differences = tail(X, length(X) - (2*m)) - head(X, length(X) - (2*m))
  i = seq(from=1, to=n, by=1)
  ci = rep(2, length(i))
  ci[i <= m] = 1 + (i[i <= m] - 1)/m
  ci[i >= n - m + 1] = 1 + (n - i[i >= n-m+1])/m
  logs = log(n * differences / (ci * m))
  logs = logs %>% replace(logs == -Inf, 0)
  mean(logs)
}

# x: (numeric vector) list of values to be normalized between 0 and 1
normalize = function(x) {(x-min(x))/(max(x)-min(x))}

# in: 2 coordinate points
# out: real number
euclid_dist_2d = function(x1, y1, x2, y2) {sqrt((x1-x2)^2 + (y1-y2)^2)}

# in: df (contains numeric columns "x", "y") & coordinate circle
# out: df with points inside circle
circle_subset = function(in.df, in.x, in.y, in.r) {subset(in.df, euclid_dist_2d(in.x, in.y, x, y) <= in.r)}

# in: real number between 0 and 1
# out: "similarity" (inputs near 0 or 1 go to 0, inputs near .5 go to 1)
info_gain = function(P) {log10(1 / (P^P * (1-P)^(1-P))) / (log10(2))}
info_gain_flipped = function(P) {1 - info_gain(P)}
igf_lambda = function(P, lambda) {(info_gain_flipped(P) + lambda * P) / (lambda + 1)}

# in: df (contains numeric columns "x", "y", and a continuous column "signal") & df (contains numeric columns "x", "y") & radius
# out: new df with statistics about each location in the list of locations inputted (empty locations are dropped)
generate_metrics = function(in.df, in.loc.df, in.r) {
  metric.df = data.frame(in.loc.df)
  metric.df$n_prox = NA
  metric.df$entropy = NA
  metric.df$mean_signal = NA
  for(i in 1:nrow(in.loc.df)) {
    i.cs = circle_subset(in.df, in.loc.df[i,]$x, in.loc.df[i,]$y, in.r)
    if(nrow(i.cs) > 0) {
      metric.df[i, ]$n_prox = nrow(i.cs)
      # instead of doing IG on pct_case, now we're doing it on a normalized mean_health value, the rest should be the same
      metric.df[i, ]$entropy = ebrahimi_entropy(i.cs$signal)
      metric.df[i, ]$mean_signal = i.cs$signal %>% mean
    }
  }
  # EXTRA METRICS (ADD OR REMOVE AS NEEDED)
  metric.df$entropy = normalize(metric.df$entropy)
  metric.df$loc_ig = info_gain(metric.df$entropy)
  metric.df$igf_lambda = igf_lambda(metric.df$entropy, 1)
  # addition edit for continuous metrics
  metric.df$mhn_ig = info_gain(metric.df$entropy)
  metric.df$igfl_mhn = igf_lambda(metric.df$entropy, 1)
  # final metric for scoring "value", lowest is best
  metric.df$score = normalize(
    metric.df$igfl_mhn * -1 +
    # metric.df$mhn_ig * 1 +
    normalize(metric.df$n_prox) * -1.15
  )
  metric.df
}

# in: df (contains integer column "n_prox") & real number
# out: df containing only rows where n_prox is greater than mean * a multiplier for offsetting mean
proximity_subset = function(in.df, in.mean.mult=1) {
  n_prox.mean = mean(in.df$n_prox)
  subset(in.df, in.df$n_prox > n_prox.mean * in.mean.mult)
}

# in: df (contains numeric columns "x", "y") & coordinate circle
# out: df with points outside circle
anti_circle_subset = function(in.df, in.x, in.y, in.r) {subset(in.df, euclid_dist_2d(in.x, in.y, x, y) > in.r)}

# in: df (contains numeric column "x", "y") & real number (min distance between points) & integer (first row to use when starting loop)
# out: df containing rows (including in.start.row) that are spaced apart according to in.r
subset_spaced = function(in.loc.df, in.r, in.start.row=1) {
  new.df = data.frame()
  new.df = rbind(new.df, in.loc.df[in.start.row,])
  while(nrow(in.loc.df) > 0) {
    in.loc.df = anti_circle_subset(in.loc.df, new.df[nrow(new.df),]$x, new.df[nrow(new.df),]$y, in.r)
    new.df = rbind(new.df, in.loc.df[1,])
  }
  new.df[-c(nrow(new.df)),]  # INSPECT THIS IN MORE DETAIL (MAYBE SOMETIMES YOU DON'T NEED TO DO THIS)
}

# in: df (contains numeric column "x", "y", "loc_ig"(for sorting)) &
#     real number (min distance between points) &
#     integer (top n to consider when optimizing) &
#     character vector (name of columns to optimize (last to sort goes first))
# out: df containing rows that are spaced apart according to in.r that are also globally optimized to a numeric column
subset_spaced_sorted = function(in.loc.df, in.r, in.max.knots, in.optimized.cols) {
  for(i in length(in.optimized.cols):1) {in.loc.df = in.loc.df[order(in.loc.df[, c(in.optimized.cols[i])]),]}
  head(subset_spaced(in.loc.df, in.r), in.max.knots)
}

# in: df (contains numeric columns "x", "y", and boolean column "case") "df" &
#     df (contains numeric columns "x", "y") "loc_df" &
#     positive real number "knot_radius" &
#     real number "min_circle_proximity_radius_multiplier" &
#     positive integer "max_knots" &
#     character (name of column to optimize) "col_names_to_optimize" &
#     real number "proximity_subset_mean_multiplier" &
# out: df of locations that are optimized to a specific numeric column
generate_optimized_knots = function(df, loc_df, knot_radius, min_circle_proximity_radius_multiplier, max_knots, col_names_to_optimize, proximity_subset_mean_multiplier=1) {
  df.knots = generate_metrics(df, loc_df[, c("x", "y")], knot_radius)
  df.knots = proximity_subset(df.knots, proximity_subset_mean_multiplier)
  df.knots = subset_spaced_sorted(df.knots, knot_radius * min_circle_proximity_radius_multiplier, max_knots, col_names_to_optimize)
  head(df.knots, max_knots)
}

gok_list_args = function(arg_list) {
  generate_optimized_knots(
    # df; people's location and their cases (cols "x", "y", and "case")
    arg_list$df,

    # df; all potential knot locations (cols "x" and "y")
    arg_list$loc_df,

    # positive numeric; radius of all knots.
    arg_list$knot_radius,

    # positive numeric; how close you want the knots to be to each other.
    # (1=knots overlap no more than radius, 2=knots can't overlap)
    # WARNING: setting too close to 0 will result in VERY long run times! (0=infinite run time)
    arg_list$min_circle_proximity_radius_multiplier,

    # positive integer; algorithm will try to optimize only considering the top N that you set.
    arg_list$max_knots,

    # character vector; the names of numeric cols that are minimized during optimization.
    # the cols inputted are ordered from most important to least important.
    arg_list$col_names_to_optimize,

    # positive numeric; for "proximity_subset" function.
    # to make sure that the algorithm doesn't select knots that are too far from populated areas,
    # this selects a subset of points where n_prox is greater than the mean * some multiplier.
    # (0=0*mean, .5=.5*mean, 1=mean, 2=2*mean)
    # WARNING: setting too close to 0 will result in longer run times!
    arg_list$proximity_subset_mean_multiplier
  )
}

# Functions (Automatic Knot Radius) ---------------------------------------

# in: list of arguments (complies with code declared in gok_list_args()) &
#     positive real number
# out: df with single row containing summary statistics of the knot algorithm
algo_iter_single = function(in.arg.list, in.knot.radius) {
  in.arg.list$knot_radius = in.knot.radius
  df.knots.iter = gok_list_args(in.arg.list)
  data.frame(
    knot_radius=in.knot.radius,
    mean_n_prox=df.knots.iter$n_prox %>% mean,
    mean_mhn_ig=mean(df.knots.iter$mhn_ig),
    mean_igfl_mhn_norm=mean(df.knots.iter$igfl_mhn),
    pct_from_max_knots=nrow(df.knots.iter) / in.arg.list$max_knots,
    mean_igfl=sum(df.knots.iter$igf_lambda %>% normalize) / nrow(df.knots.iter)
  )
}

# in: list of arguments (complies with code declared in gok_list_args()) &
#     sequence object (contains all knot radii to consider)
# out: df containing summary statistics about the results of each knot algo run (1 row = 1 knot radius)
generate_akr_df = function(in.arg.list, in.radius.seq) {
  df.akr.stats = data.frame()
  print("Generating AKR Data Frame...")
  for(i in in.radius.seq) {
    df.akr.stats = rbind(df.akr.stats, algo_iter_single(in.arg.list, i))
    #print(paste("Completion percentage:", round(i/max(in.radius.seq), 4)))
  }
  # line below is important for automatic radius selection
  df.akr.stats$dist_mli_mc =
    euclid_dist_2d(df.akr.stats$mean_mhn_ig, 0, 0, 0) +
    euclid_dist_2d(df.akr.stats$mean_igfl, 0, .5, 0)*2
  print("Done!")
  df.akr.stats
}

# in: df (complies with output of generate_akr_df())
# out: finds knot_radius in the input df with smallest dist_mli_mc
find_best_knot_radius = function(in.akr.df) {in.akr.df[order(in.akr.df[, c("dist_mli_mc")]),][1,]$knot_radius}

# generate.optimum.knots_list.args_automatic.knot.radius
# in: list of arguments (complies with code declared in gok_list_args()) &
#     sequence object (contains all knot radii to consider)
# out: df of locations that are optimized to a specific numeric column with automatically optimized knot radius
gok_la_akr = function(in.arg.list, in.radius.seq) {
  df.akr = generate_akr_df(in.arg.list, in.radius.seq)
  in.arg.list$knot_radius = find_best_knot_radius(df.akr)
  gok_list_args(in.arg.list)
}

# Functions (Report) ------------------------------------------------------

# show all knot locations on the map with their respective case
graph_algo_results = function(in.df.algo, in.df.knots, in.arg.list) {
  ggplot(data=in.df.algo, label="test") +
    geom_point(aes(x=x, y=y, color=signal)) +
    scale_color_gradient2() +
    geom_circle(data=in.df.knots, aes(x0=x, y0=y, r=in.arg.list$knot_radius, color=entropy)) +
    geom_text(data=in.df.knots, aes(x=x, y=(y + .015), label=round(mhn_ig, 2)), color="black", size=5) +
    geom_text(data=in.df.knots, aes(x=x, y=(y - .015), label=round(igfl_mhn, 2)), color="#550077", size=5) +
    labs(title="BLACK(top) is 'mhn_ig' | PURPLE(bottom) is 'igfl_mhn'")
}

generate_knot_report = function(in.df.algo, in.df.knots, in.arg.list) {
  # show all final knot locations
  print(in.df.knots)

  # sometimes, the algorithm can't make as many knots as you specified.
  # this shows the percentage from max_knots it was able to make. the closer to 1, the better
  print(paste("'pct_from_max_knots':", nrow(in.df.knots) / in.arg.list$max_knots))

  # shows the average mhn_ig, which was minimized in the final step of the algorithm
  print(paste("Mean 'mhn_ig':", mean(in.df.knots$mhn_ig)))

  # show all knot locations on the map with their respective case
  graph_algo_results(in.df.algo, in.df.knots, in.arg.list)
}

graph_akr_results = function(in.akr.df) {
  best_knot_radius = find_best_knot_radius(in.akr.df)
  ggplot(data=in.akr.df) +
    geom_line(aes(x=knot_radius, y=mean_mhn_ig, color="mean_mhn_ig")) +
    # geom_line(aes(x=knot_radius, y=mean_igfl_mhn_norm, color="mean_igfl_mhn_norm")) +
    geom_line(aes(x=knot_radius, y=pct_from_max_knots, color="pct_from_max_knots")) +
    geom_line(aes(x=knot_radius, y=mean_igfl, color="mean_igfl")) +
    geom_vline(xintercept=best_knot_radius) + geom_text(aes(x=best_knot_radius+.005, y=.75, label=round(best_knot_radius, 4)))
}

generate_akr_report = function(in.akr.df) {
  print(in.akr.df)
  print(paste("Best knot radius:", find_best_knot_radius(in.akr.df)))
  graph_akr_results(in.akr.df)
}

# # Algorithm Testing -------------------------------------------------------
#
# # df.algo = data.frame(df.real)
# # df.algo = read.csv("med.csv")
#
# if(F) {
#   # Create a copy of the original data and work with that copy
#   df = data.frame(df.orig)
#   # Custom elevation of odds ratios to make this work better with a small number of points
#   df[df$ELEVATED == TRUE, 'POL'] = 2.5 * df[df$ELEVATED == TRUE, 'POL']
#   df$CASE_CNTL = rbinom(n=nrow(df), size=1, prob=c(df$POL, 1-df$POL))
#
#   df.algo = data.frame(df)
#   df.algo$x = df.algo$X_COORD
#   df.algo$y = df.algo$Y_COORD
#   df.algo$case = df.algo$CASE_CNTL
# }
#
# # World Generation --------------------------------------------------------
#
# generate_world_old = function(N) {
#   # generate random points with random health value
#   df = data.frame(
#     x=runif(N),
#     y=runif(N),
#     health=rnorm(N, 4, .5)) # sample(0:1, size=N, replace=T)
#   # include a "hot-spot" of increased change in health values
#   df$health = df$health - euclid_dist_2d(.25, .25, df$x, df$y)*10 + .005  # x, y; + r
#   # add "case" column
#   df$case = ifelse(df$health > 0, 1, 0)
#   df
# }
#
# N = 2500
#
# df = generate_world_old(N)
#
# # Experimental Testing ----------------------------------------------------
#
# df.algo = data.frame(df)
# ggplot(data=df) + geom_point(aes(x=x, y=y, color=health)) + scale_color_gradient2()
#
knot_algo_args = list(
  df=df.algo,
  loc_df=df.algo,
  knot_radius=.05,
  min_circle_proximity_radius_multiplier=1.5,
  max_knots=25,
  col_names_to_optimize=c("score"),
  proximity_subset_mean_multiplier=0
)
#
# df.knots = gok_list_args(knot_algo_args)
# generate_knot_report(df.algo, df.knots, knot_algo_args)
#
akr_seq = seq(from=.051, to=.15, length.out=100)
#
# # start
df.akr = generate_akr_df(knot_algo_args, akr_seq)
generate_akr_report(df.akr)
#
knot_algo_args$knot_radius = find_best_knot_radius(df.akr)
df.knots = gok_list_args(knot_algo_args)
generate_knot_report(df.algo, df.knots, knot_algo_args)
#
# # idea: have each point weighted by the

# print(paste0("Now saving month: ", ccnt, "/12"))
# df.kav4.12list[[ccnt]] = df.knots

# Extra -----------------------------------------------------------------------------------------------------------

library(mgcv)
library(ggeffects)

ggplot() + geom_point(aes(x=seq(from=-1, to=2, length.out=100), y=seq(from=-1, to=2, length.out=100) %>% igf_lambda(1)))

df.algo %>% head



url = "https://raw.githubusercontent.com/m-clark/generalized-additive-models/master/data/pisasci2006.csv"
pisa = read.csv(url)

ggplot(pisa) + geom_point(aes(x=Income, y=Health)) + geom_smooth(aes(x=Income, y=Health), method="gam")

mod_gam0 = gam(Income ~ Health, data=pisa)
mod_gam0 %>% summary()

mod_gam1 = gam(Income ~ s(Health), data=pisa)
mod_gam1 %>% summary()

# ggplot(pisa) + geom_point(aes(x=Income, y=Health)) + geom_ribbon(data=ggpredict(mod_gam0), aes(ymin=conf.low, ymax=conf.high), alpha=.1)
# plot(mod_gam1)


# full data (no knots)
full_gam = gam(
  signal ~ te(x, y),
  data=df.algo,
  method="REML",
  family=gaussian
)
full_gam %>% summary()

vis.gam(full_gam, type="response", plot.type="contour")

# knot data
knot_gam = gam(
  mean_signal ~ te(x, y),
  data=df.knots,
  method="REML",
  family=gaussian
)
knot_gam %>% summary()

vis.gam(knot_gam, type="response", plot.type="contour")