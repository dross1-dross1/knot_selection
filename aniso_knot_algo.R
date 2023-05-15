# Setup -----------------------------------------------------------------------------------------------------------
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

# Knot Evaluation
source("scripts/load_epa_data.R")
source("kaf_v6_sphere.R")
source("kaf_v6_ellipsoid.R")
source("scripts/stan_data_prep.R")


# Parameters --------------------------------------------------------------
# Circle algo
max_knots.sphere = 25
radius_mult.sphere = 1
num_neighbors.sphere = 5

# Stan ellipse model
spectral_components = 5

# Ellipse algo
max_knots.ellipsoid = 15
radius_mult.ellipsoid = 1/25


# ------ === Select knots using elliptical neighborhood entropy === ----------
# Data loading (testing, replace with your selected data) -----------------
epa.df = load_epa_data()

# Take a sample of the data for runtime purposes (California-ish)
epa.df = epa.df[(epa.df$x < 0.15) & (epa.df$y < 0.75), ]

ggplot(epa.df) +
  geom_point(aes(x=x, y=y, fill=pm2_5), size=3, pch=21, color='black') +
  scale_fill_gradient2(low='blue', high='red', midpoint=0)

# Prepare data for algorithm
epa.df$signal = epa.df$pm2_5
epa.df$z = 0

# Circle algorithm to select knot candidates ------------------------------
entropy.sphere.df = entropy_max.sphere(epa.df, num_neighbors.sphere, radius_mult.sphere, max_knots.sphere)

# Show knots selected based on entropy in circular neighborhoods
ggplot(epa.df) +
  geom_point(aes(x=x, y=y, fill=pm2_5), size=3, pch=21, color='black') +
  geom_point(data=entropy.sphere.df, aes(x=x, y=y), shape=23, size=5, fill='black') +
  scale_fill_gradient2(low='blue', high='red', midpoint=0)


# Prepare data for Stan model ---------------------------------------------
# Rename columns as needed
knots.sphere.df = entropy.sphere.df[, c('x', 'y', 'signal')]
knots.sphere.df$pm2_5 = knots.sphere.df$signal
epa.df$pm2_5 = epa.df$signal

# Prepare data for stan fitting
data = prepare_data(epa.df, M=spectral_components, plot=F, df.knots=knots.sphere.df)


# Stan model fitting (fit ellipses) -----------------------------------------
model.ellipse = cmdstan_model('stan/aniso_process_axes_latent_knots_spectral_profile.stan')

fit.ellipse = model.ellipse$sample(data=data,
                                   parallel_chains=4,
                                   iter_warmup=2000,
                                   max_treedepth=10,
                                   init=function() list(sigma_z=0.1,
                                                        ell_z=0.1,
                                                        sigma_interp=0.1,
                                                        ell_interp=0.1,
                                                        lambda_y=0.1))

# Extract the locations and PM2.5 data
spatial.df = cbind(data$spatial_locs, data$y_spatial)
names(spatial.df) = c('x', 'y', 'pm2_5')
knots.df = data$knot_locs


# Extract ellipses --------------------------------------------------------
# Get the length and rotation of the ellipses
ellipse.df = extract_ellipses(spatial.df,
                              fit=fit.ellipse,
                              plot=T,
                              scale_ellipses=10, # Turn on if plotting
                              psi_suffix="_all",
                              return_df=TRUE)

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


# ------ === Model analysis === ----------
# Now that you've fit the knots, here is some starter code to run analysis.


# Plotting ellipses --------------------------------------------------
# If you want to plot the ellipses, you can do so easily. Here I plot _all_ ellipses,
# scaled by the scaling factor set in the parameters.
ggplot(ellipse.df) +
  geom_point(aes(x=x0, y=y0, fill=signal), size=3, pch=21, color='black') +
  scale_fill_gradient2(low='blue', high='red', midpoint=0) +
  geom_ellipse(aes(x0=x0, y0=y0, a=radius_mult.ellipsoid*a, b=radius_mult.ellipsoid*b, angle=alpha)) +
  geom_point(data=knots.ellipsoid.df, aes(x=x0, y=y0), shape=23, size=5, fill='black') +
  coord_fixed() +
  ggtitle("Ellipses")

# This can be a little messy, so if you want, you can also just plot the ellipses
# at the knots selected.
ggplot(ellipse.df) +
  geom_point(aes(x=x0, y=y0, fill=signal), size=3, pch=21, color='black') +
  scale_fill_gradient2(low='blue', high='red', midpoint=0) +
  geom_ellipse(data=knots.ellipsoid.df, aes(x0=x0, y0=y0, a=radius_mult.ellipsoid*a, b=radius_mult.ellipsoid*b, angle=alpha)) +
  geom_point(data=knots.ellipsoid.df, aes(x=x0, y=y0), shape=23, size=5, fill='black') +
  coord_fixed() +
  ggtitle("Ellipses")


# Making predictions ----------------------------------------------------
# The ellipses are nice to look at, but what we really care about are making predictions.
# Here I'll fit that "predictive process" model I talked about in Slack.
model.pp = cmdstan_model('stan/predictive_process.stan')

# Extract the knots we just fit in the ellipse algo
data.pp = list(N_knots=nrow(knots.ellipsoid.df),
              knot_locs=knots.ellipsoid.df[, c('x0', 'y0')],
              N_spatial=nrow(epa.df),
              spatial_locs=epa.df[, c('x', 'y')],
              y_spatial=epa.df$pm2_5 %>% as.vector)

# Fit the predictive process model
fit.pp = model.pp$sample(data=data.pp,
                         parallel_chains=4,
                         iter_warmup=2000,
                         max_treedepth=10,
                         init=function() list(sigma_z=0.1,
                                              ell_z=0.1,
                                              sigma_interp=0.1,
                                              ell_interp=0.1,
                                              lambda_y=0.1))

# Extract the predictions (these are slow to extract, so just be patient)
pp.preds = fit.pp$summary(variables = 'y_spatial_sim')

# Store the locations, true values, and predictions
preds.df = epa.df[, c('x', 'y', 'pm2_5')]
preds.df$median_pred = pp.preds$median
preds.df$mean_pred = pp.preds$mean
preds.df$sd_pred = pp.preds$sd

# Store the difference between the actual and predicted values
preds.df$diff_median = preds.df$pm2_5 - preds.df$median_pred
preds.df$diff_mean = preds.df$pm2_5 - preds.df$mean_pred

# Also store the percentage difference
preds.df$pct_diff_median = (preds.df$pm2_5 - preds.df$median_pred) / preds.df$pm2_5
preds.df$pct_diff_mean = (preds.df$pm2_5 - preds.df$mean_pred) / preds.df$pm2_5
  

# Plot the predictions (percent difference using the median, but you can change this)
ggplot(preds.df) +
  geom_point(data=knots.ellipsoid.df, aes(x=x0, y=y0), shape=23, size=5, fill='black') + # Also plot the knots
  geom_point(aes(x=x, y=y, fill=diff_median), size=3, pch=21, color='black') +
  scale_fill_gradient2(low='white', high='red')
  
  

# Error metrics ------------------------------------------------
# Mean squared error
mse = sum((preds.df$pm2_5 - preds.df$median_pred)**2) / nrow(preds.df)
# Mean absolute error
mae = sum(abs(preds.df$pm2_5 - preds.df$median_pred)) / nrow(preds.df)
# R^2
r2 = cor(preds.df$pm2_5, preds.df$median_pred)

# Print output (or sink it to a file). "cat" respects "\n", whereas print does not.
cat(paste0("MSE = ", round(mse, 4), "\nMAE = ", round(mae, 4), "\nR^2 = ", round(r2, 4)))

# R2 graph
ggplot(preds.df) +
  geom_point(aes(x=pm2_5, y=median_pred)) +
  geom_abline(slope=1, intercept=0, color='blue')


# LOO ---------------------------------------------------------------------
# Here's how you do that "loo" thing I was talking about, it's just a one-liner.
# It prints out a bunch of stuff. It probably makes sense to just save the entire output to a file.
# You can do that using "sink"

sink("my_file.txt")
fit.pp$loo()
sink()


# Predictions using "circular" knots --------------------------------------
# In addition, we need to check if this whole ellipse approach produces better
# results than just working with your original algorithm (with circles). It would
# be more efficient to turn the code above into functions and re-run them with
# different knots. But for the sake of simplicity, I'm just going to copy-paste
# a few things.

# Re-format the data to use the original (non-elliptical) knots
data.pp.sphere = list(N_knots=nrow(knots.sphere.df),
               knot_locs=knots.sphere.df[, c('x', 'y')],
               N_spatial=nrow(epa.df),
               spatial_locs=epa.df[, c('x', 'y')],
               y_spatial=epa.df$pm2_5 %>% as.vector)

# Fit the predictive process model
fit.pp.sphere = model.pp$sample(data=data.pp.sphere,
                                 parallel_chains=4,
                                 iter_warmup=2000,
                                 max_treedepth=10,
                                 init=function() list(sigma_z=0.1,
                                                      ell_z=0.1,
                                                      sigma_interp=0.1,
                                                      ell_interp=0.1,
                                                      lambda_y=0.1))

# Extract the predictions (these are slow to extract, so just be patient)
pp.sphere.preds = fit.pp.sphere$summary(variables = 'y_spatial_sim')

# Store the locations, true values, and predictions
preds.sphere.df = epa.df[, c('x', 'y', 'pm2_5')]
preds.sphere.df$median_pred = pp.sphere.preds$median
preds.sphere.df$sd_pred = pp.sphere.preds$sd

# Mean squared error
mse.sphere = sum((preds.sphere.df$pm2_5 - preds.sphere.df$median_pred)**2) / nrow(preds.sphere.df)
# Mean absolute error
mae.sphere = sum(abs(preds.sphere.df$pm2_5 - preds.sphere.df$median_pred)) / nrow(preds.sphere.df)
# R^2
r2.sphere = cor(preds.sphere.df$pm2_5, preds.sphere.df$median_pred)

# Print output (or sink it to a file). "cat" respects "\n", whereas print does not.
cat(paste0("MSE = ", round(mse.sphere, 4), "\nMAE = ", round(mae.sphere, 4), "\nR^2 = ", round(r2.sphere, 4)))

# Finally, let's compare predictions between both sets of knots
preds.df$type = 'elliptical'
preds.sphere.df$type = 'sphere'
preds.all.df = rbind(preds.df[, names(preds.sphere.df)], preds.sphere.df)

ggplot(preds.all.df, aes(x=pm2_5, y=median_pred, color=type)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  geom_abline(slope=1, intercept=0, color='blue')
