# Setup -----------------------------------------------------------------------------------------------------------

# set.seed(1)

# Set WD to directory of current script location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# UI and misc
# if (!is.null(dev.list()["RStudioGD"])) { dev.off(dev.list()["RStudioGD"]) } # clear plots
# rm(list = ls())                                                             # remove all variables
# set.seed(100)                                                               # reproducable randomness
# cat("\014")                                                                 # ctrl+L (clear console output)

# Imports (Graphing)
library(ggplot2)  # graphing; ggplot() + ...
library(ggforce)  # graphing circles; geom_circle()

# Imports (Modeling)
library(cmdstanr)
library(loo)
library(posterior)

# Imports
library(dplyr)    # piping and QoL functions; %>%
library(scales)   # normalizing lists between 0 and 1; rescale()
library(stringr)  # manipulating strings
library(mvtnorm)

# Imports (Knot Selection)
library(fields)   # cover.design; cover.design()

# Knot Evaluation
source("scripts/load_epa_data.R")
source("kaf_v6_sphere.R")
source("kaf_v6_ellipsoid.R")
source("scripts/stan_data_prep.R")

# TODO: [x] save summary.stats to file
# TODO: [x] add graphing
# TODO: [x] also save the files to the directory
# TODO: [ ] get hyperparameter testing
# TODO: [ ] maybe parallelize hyperparameter testing
# TODO: [ ]

# Data Loading ----------------------------------------------------------------------------------------------------

epa.df = load_epa_data()  # entire US

# Take a sample of the data for runtime purposes (California-ish)
epa.df = epa.df[(epa.df$x < 0.15) & (epa.df$y < 0.75), ]

# show data
ggplot(epa.df) +
  geom_point(aes(x=x, y=y, fill=pm2_5), size=3, pch=21, color='black') +
  scale_fill_gradient2(low='blue', high='red', midpoint=0)

epa.sample = data.frame(epa.df)
# Get a sample of the rows as data
# epa.sample.idx = sample(1:nrow(epa.df), 250)
# epa.sample = epa.df[epa.sample.idx, ]
# rownames(epa.sample) = NULL


# Functions -------------------------------------------------------------------------------------------------------

# draws is the fit$draws() aka posterior samples
# returns a vector with a size matching the nrow of the original dataset
get_chain_preds_from_draws = function(draws, chain.number) {
  # get the draws and turn it into a data frame
  draws.df = draws %>% as_draws_df

  # find columns that match the pattern "R_space_knot[*,n]"
  cols = grep(paste0("R_space_knot\\[\\d+,", chain.number, "\\]"), colnames(draws.df))
  df_filtered = draws.df[, cols] %>% as.data.frame

  # get the mean from each column
  chain_means = colMeans(df_filtered) %>% as.vector

  return(chain_means)
}

# df.points contains 3 columns: `x`, `y`, `pm2_5`
# df.knots contains 3 columns: `x`, `y`, `pm2_5`
# gives you the fit object
fit_model = function(df.points, df.knots) {
  # Fitting Predictive Process Model
  model = cmdstan_model('stan/predictive_process.stan')

  data = list(
    N_knots = nrow(df.knots),
    knot_locs = df.knots[, c('x', 'y')],
    N_spatial = nrow(df.points),
    spatial_locs = df.points[, c('x', 'y')],
    y_spatial = df.points$pm2_5 %>% as.vector,
    return_predictions = 1,
    return_log_likelihoods = 1
  )

  fit = model$sample(
    data = data,
    parallel_chains = 4,
    iter_warmup = 2000,
    max_treedepth = 10,
    init = function() list(
      sigma_z=0.1,
      ell_z=0.1,
      sigma_interp=0.1,
      ell_interp=0.1,
      lambda_y=0.1
    )
  )

  return(fit)
}

# also try closeAllConnections() if console output is not working
save_fit = function(fit.knots, preds.knots, summary.stats, timestamp, name) {
  file.path = paste0("results/", timestamp, "/")
  if (!dir.exists(file.path)) { dir.create(file.path, recursive=T) }
  sink(paste0(file.path, name, "_results.txt"))
  print(fit.knots$loo())
  print(summary.stats)
  sink()
  saveRDS(fit.knots, paste0(file.path, name, "_fit.RDS"))
  write.csv(preds.knots, paste0(file.path, name, "_preds.csv"), row.names=F)
}

get_preds = function(df.points, fit) {
  # Extract the predictions (these are slow to extract, so just be patient)
  fit.preds = fit$summary(variables = "y_spatial_sim")

  # Store the locations, true values, and predictions
  fit.preds.df = df.points[, c('x', 'y', 'pm2_5')]
  fit.preds.df$median_pred = fit.preds$median
  fit.preds.df$mean_pred = fit.preds$mean
  fit.preds.df$sd_pred = fit.preds$sd

  # Store the difference between the actual and predicted values
  fit.preds.df$diff_median = fit.preds.df$pm2_5 - fit.preds.df$median_pred
  fit.preds.df$diff_mean = fit.preds.df$pm2_5 - fit.preds.df$mean_pred

  # Also store the percentage difference
  fit.preds.df$pct_diff_median = (fit.preds.df$pm2_5 - fit.preds.df$median_pred) / fit.preds.df$pm2_5
  fit.preds.df$pct_diff_mean = (fit.preds.df$pm2_5 - fit.preds.df$mean_pred) / fit.preds.df$pm2_5

  return(fit.preds.df)
}

get_summary_statistics = function(preds.knots) {
  # Mean squared error
  mse = sum((preds.knots$pm2_5 - preds.knots$median_pred)**2) / nrow(preds.knots)
  # Mean absolute error
  mae = sum(abs(preds.knots$pm2_5 - preds.knots$median_pred)) / nrow(preds.knots)
  # R^2
  r2 = cor(preds.knots$pm2_5, preds.knots$median_pred)

  return(list(mse=mse, mae=mae, r2=r2[1, 1]))
}

# Create Different Kinds of Knots ---------------------------------------------------------------------------------

n.knots = 5

# Random Knot Selection
epa.sample.knots.rand.idx = sample(1:nrow(epa.sample), n.knots)
epa.sample.knots.rand = epa.sample[epa.sample.knots.rand.idx, ]

# Uniform Grid Knot Selection
generate_grid = function(num.x, num.y) { expand.grid(x=seq(from=0, to=1, length.out=num.x), y=seq(from=0, to=1, length.out=num.y)) }
epa.sample.knots.grid = generate_grid(num.x=2, num.y=2)

# Cover Design Knot Selection
epa.sample.knots.cd.locs = epa.sample[, c("x", "y")]
epa.sample.knots.cd = cover.design(epa.sample.knots.cd.locs, n.knots)$design %>% as.data.frame
epa.sample.knots.cd = merge(epa.sample.knots.cd, epa.sample, by=c("x", "y"), sort=F)

# Isotropic Entropy Maximization Knot Selection
epa.sample.prepped = data.frame(epa.sample)
names(epa.sample.prepped) = c("x", "y", "signal")
epa.sample.prepped$z = 0
epa.sample.knots.sphere_entropy = entropy_max.sphere(epa.sample.prepped, 5, 1, n.knots)
epa.sample.knots.sphere_entropy = epa.sample.knots.sphere_entropy %>% select(-z, -radius, -entropy) %>% rename(pm2_5=signal)

# Anisotropic Entropy Maximization Knot Selection
get_ellipsoid_knots = function(df.points, max_knots.sphere, radius_mult.sphere, num_neighbors.sphere, spectral_components, max_knots.ellipsoid, radius_mult.ellipsoid) {
  # Circle algorithm to select knot candidates ------------------------------
  entropy.sphere.df = entropy_max.sphere(df.points, num_neighbors.sphere, radius_mult.sphere, max_knots.sphere)

  # Prepare data for Stan model ---------------------------------------------
  # Rename columns as needed
  knots.sphere.df = entropy.sphere.df[, c('x', 'y', 'signal')]
  knots.sphere.df$pm2_5 = knots.sphere.df$signal
  df.points$pm2_5 = df.points$signal

  # Prepare data for stan fitting
  data = prepare_data(df.points, M=spectral_components, plot=F, df.knots=knots.sphere.df)

  # Stan model fitting (fit ellipses) -----------------------------------------
  model.ellipse = cmdstan_model('stan/aniso_process_axes_latent_knots_spectral_profile.stan')

  fit.ellipse = model.ellipse$sample(
  data=data,
  parallel_chains=4,
  iter_warmup=2000,
  max_treedepth=10,
  init=function() list(
      sigma_z=0.1,
      ell_z=0.1,
      sigma_interp=0.1,
      ell_interp=0.1,
      lambda_y=0.1
    )
  )

  # Extract the locations and PM2.5 data
  spatial.df = cbind(data$spatial_locs, data$y_spatial)
  names(spatial.df) = c('x', 'y', 'pm2_5')
  knots.df = data$knot_locs

  # Extract ellipses --------------------------------------------------------
  # Get the length and rotation of the ellipses
  ellipse.df = extract_ellipses(
    spatial.df,
    fit=fit.ellipse,
    plot=T,
    scale_ellipses=10, # Turn on if plotting
    psi_suffix="_all",
    return_df=TRUE
  )

  # Rename columns for ellipse algo
  ellipse.df$signal = ellipse.df$pm2_5
  ellipse.df$x0 = ellipse.df$x
  ellipse.df$y0 = ellipse.df$y
  ellipse.df$z0 = 0
  ellipse.df$c = 0

  ellipse.df$x = NULL
  ellipse.df$y = NULL
  ellipse.df$z = NULL
  ellipse.df$pm2_5 = NULL

  ellipse.df$alpha = atan(ellipse.df$b / ellipse.df$a)
  ellipse.df$beta = 0
  ellipse.df$gamma = 0

  ellipse.df = ellipse.df[, c('a', 'b', 'c', 'x0', 'y0', 'z0', 'alpha', 'beta', 'gamma', 'signal')]

  # Ellipse algo ------------------------------------------------------------
  knots.ellipsoid.df = entropy_max(ellipse.df, radius_mult.ellipsoid, max_knots.ellipsoid)
  knots.ellipsoid.df
}

epa.sample.knots.ellipsoid_entropy = get_ellipsoid_knots(
  df.points=epa.sample.prepped,
  max_knots.sphere=25,
  radius_mult.sphere=1,
  num_neighbors.sphere=5,
  spectral_components=5,
  max_knots.ellipsoid=n.knots,
  radius_mult.ellipsoid=1/25
)

epa.sample.knots.ellipsoid_entropy = epa.sample.knots.ellipsoid_entropy %>% select(x0, y0) %>% rename(x = x0, y = y0)

# Trying It Out ---------------------------------------------------------------------------------------------------

run.timestamp = Sys.time() %>% format("%Y-%m-%d_%H_%M_%S")

# Random Knot Selection
fit.knots.rand = fit_model(epa.sample, epa.sample.knots.rand)
preds.knots.rand = get_preds(epa.sample, fit.knots.rand)
preds.knots.rand$type = "random"
summary.stats.rand = get_summary_statistics(preds.knots.rand)
save_fit(fit.knots=fit.knots.rand, preds.knots=preds.knots.rand, summary.stats=summary.stats.rand, timestamp=run.timestamp, name="random")

# Uniform Grid Knot Selection
fit.knots.grid = fit_model(epa.sample, epa.sample.knots.grid)
preds.knots.grid = get_preds(epa.sample, fit.knots.grid)
preds.knots.grid$type = "grid"
summary.stats.grid = get_summary_statistics(preds.knots.grid)
save_fit(fit.knots=fit.knots.grid, preds.knots=preds.knots.grid, summary.stats=summary.stats.grid, timestamp=run.timestamp, name="grid")

# Cover Design Knot Selection
fit.knots.cd = fit_model(epa.sample, epa.sample.knots.cd)
preds.knots.cd = get_preds(epa.sample, fit.knots.cd)
preds.knots.cd$type = "cover_design"
summary.stats.cd = get_summary_statistics(preds.knots.cd)
save_fit(fit.knots=fit.knots.cd, preds.knots=preds.knots.cd, summary.stats=summary.stats.cd, timestamp=run.timestamp, name="cover_design")

# Isotropic Entropy Maximization Knot Selection
fit.knots.sphere_entropy = fit_model(epa.sample, epa.sample.knots.sphere_entropy)
preds.knots.sphere_entropy = get_preds(epa.sample, fit.knots.sphere_entropy)
preds.knots.sphere_entropy$type = "sphere_entropy"
summary.stats.sphere_entropy = get_summary_statistics(preds.knots.sphere_entropy)
save_fit(fit.knots=fit.knots.sphere_entropy, preds.knots=preds.knots.sphere_entropy, summary.stats=summary.stats.sphere_entropy, timestamp=run.timestamp, name="sphere_entropy")

# Anisotropic Entropy Maximization Knot Selection
fit.knots.ellipsoid_entropy = fit_model(epa.sample, epa.sample.knots.ellipsoid_entropy)
preds.knots.ellipsoid_entropy = get_preds(epa.sample, fit.knots.ellipsoid_entropy)
preds.knots.ellipsoid_entropy$type = "ellipsoid_entropy"
summary.stats.ellipsoid_entropy = get_summary_statistics(preds.knots.ellipsoid_entropy)
save_fit(fit.knots=fit.knots.ellipsoid_entropy, preds.knots=preds.knots.ellipsoid_entropy, summary.stats=summary.stats.ellipsoid_entropy, timestamp=run.timestamp, name="ellipsoid_entropy")

# Viewing Results -------------------------------------------------------------------------------------------------

summary.stats.rand
summary.stats.grid
summary.stats.cd
summary.stats.sphere_entropy
summary.stats.ellipsoid_entropy

# compare predictions between all sets of knots
preds.knots.all.df = rbind(preds.knots.rand, preds.knots.grid, preds.knots.cd, preds.knots.sphere_entropy, preds.knots.ellipsoid_entropy)

ggplot(preds.knots.all.df, aes(x=pm2_5, y=median_pred, color=type)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  geom_abline(slope=1, intercept=0, color='blue')
