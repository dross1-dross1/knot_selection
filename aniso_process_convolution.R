# Setup -----------------------------------------------------------------------------------------------------------
# Imports (Graphing)
library(ggplot2)  # graphing; ggplot() + ...
library(ggforce)  # graphing circles; geom_circle()

# Imports (Modeling)
library(cmdstanr)

# Imports
library(dplyr)    # piping and QoL functions; %>%
library(scales)   # normalizing lists between 0 and 1; rescale()
library(stringr)  # manipulating strings
library(mvtnorm)

# Knot Evaluation
source("scripts/load_epa_data.R")


# Data loading ------------------------------------------------------------
epa.df = load_epa_data()

ggplot(epa.df) +
  geom_point(aes(x=x, y=y, fill=pm2_5), size=3, pch=21, color='black') +
  scale_fill_gradient2(low='blue', high='red', midpoint=0)


# Stan model fitting ------------------------------------------------------
# Prepare data for the stan model
prepare_data = function(df, num_knots=NULL, sample_pct=1, plot=FALSE, fixed_rank_dim=NULL) {
  # Reduce sample size if needed
  if (sample_pct < 1) {
    df.sample.idx = sample(1:nrow(df), as.integer(sample_pct*nrow(df)))
    df.sample = df[df.sample.idx, ]
    rownames(df.sample) = 1:nrow(df.sample)
  } else {
    df.sample = df
  }
  
  # Spatial process sites and data
  spatial_locs = df.sample[, c('x', 'y')]
  y_spatial = df.sample$pm2_5 %>% as.vector
  N_spatial = nrow(spatial_locs)
  
  # Package data
  data = list(N_spatial=N_spatial,
              spatial_locs=spatial_locs,
              y_spatial=y_spatial)
  
  # Include knots if needed
  if (!is.null(num_knots)) {
    # Randomly select num_knots
    df.knots.idx = sample(1:nrow(df.sample), size=num_knots)
    df.knots = df.sample[df.knots.idx, ]
    
    data$N_knots = nrow(df.knots)
    data$knot_locs = df.knots[, c('x', 'y')]
  }
  
  # Fixed rank dimension (if using)
  if (!is.null(fixed_rank_dim)) {
    data$fixed_rank_dim = fixed_rank_dim
  }
  
  if (plot) {
    g = ggplot(df) +
      geom_point(aes(x=x, y=y, fill=pm2_5), size=3, pch=21, color='black', alpha=0.1) +
      geom_point(data=df.sample, aes(x=x, y=y, fill=pm2_5), size=3, pch=21, color='black') +
      scale_fill_gradient2(low='blue', high='red', midpoint=0) +
      ggtitle('US (sample points highlighted)')
    
    if (!is.null(num_knots)) {
      g = g + geom_point(data=df.knots, aes(x=x, y=y), shape=23, size=5, fill='black') +
        coord_fixed()
    }
    show(g)
  }
  
  return(data)
}


# Extract ellipses for plotting ------------------------------------------
# After the model is fit, extract (and optionally plot) the ellipses
extract_ellipses = function(df, df.all=NULL, fit, plot=TRUE, scale_ellipses=1, num_ellipses=NULL, knots=NULL, psi_suffix='', return_df=FALSE) {
  psi_x_col = paste0('psi_x', psi_suffix)
  psi_y_col = paste0('psi_y', psi_suffix)
  
  a = fit$summary(variables = psi_x_col)$median
  b = fit$summary(variables = psi_y_col)$median
  rotation = atan(b / a)
  
  # Append rotations and ellipse major/minor axes to original data
  df$a = a
  df$b = b
  df$alpha = rotation
  
  # If knots are supplied, extract their axes/rotation for plotting
  if (!is.null(knots)) {
    knots.df = merge(knots, df, by=c('x', 'y'), all.x=TRUE)
  }
  
  # Plot ellipses on map
  if (plot) {
    g = ggplot(df) +
          geom_point(aes(x=x, y=y, fill=pm2_5), size=3, pch=21, color='black') +
          scale_fill_gradient2(low='blue', high='red', midpoint=0) +
          ggtitle('US (sample points highlighted)') +
          coord_fixed()
    
    if (!is.null(df.all)) {
      g = g + geom_point(aes(x=x, y=y, fill=pm2_5), size=3, pch=21, color='black', alpha=0.1)
    }
    
    if (!is.null(num_ellipses)) {
      ellipse.idx = sample(1:nrow(df), size=num_ellipses)
      g = g + geom_ellipse(data=df[ellipse.idx,], aes(x0=x, y0=y, a=a/scale_ellipses, b=b/scale_ellipses, angle=alpha)) +
              geom_point(data=df[ellipse.idx,], aes(x=x, y=y), shape=13)
    } else if (!is.null(knots)) {
      g = g + geom_ellipse(data=knots.df, aes(x0=x, y0=y, a=a/scale_ellipses, b=b/scale_ellipses, angle=alpha)) +
        geom_point(data=knots.df, aes(x=x, y=y), shape=13)
    } else {
      g = g + geom_ellipse(data=df, aes(x0=x, y0=y, a=a/scale_ellipses, b=b/scale_ellipses, angle=alpha))
    }
    show(g)
  }
  
  if (return_df) {
    return(df)
  }
}

# Fitting knot model ------------------------------------------------------
model.knots = cmdstan_model('stan/aniso_process_axes_latent_knots_spectral_profile.stan')

# Data
data.knots = prepare_data(epa.df, num_knots=15, sample_pct=0.32, plot=T, fixed_rank_dim=as.integer(sqrt(nrow(epa.df))))

fit.knots = model.knots$sample(data=data.knots,
                                parallel_chains=4,
                                iter_warmup=2000,
                                max_treedepth=10,
                                init=function() list(sigma_z=0.1,
                                                     ell_z=0.1,
                                                     sigma_interp=0.1,
                                                     ell_interp=0.1,
                                                     lambda_y=0.1))

# Extract the locations and data
data.knots.spatial.df = cbind(data.knots$spatial_locs, data.knots$y_spatial)
names(data.knots.spatial.df) = c('x', 'y', 'pm2_5')
data.knots.knots.df = data.knots$knot_locs

# Plot the ellipses
extract_ellipses(data.knots.spatial.df, 
                 fit=fit.knots, 
                 scale_ellipses=7, 
                 # knots=data.knots.knots.df,
                 psi_suffix="_all")

# Extract the data frame with ellipse data
ellipse_df = extract_ellipses(data.knots.spatial.df, 
                 fit=fit.knots, 
                 scale_ellipses=5, 
                 psi_suffix="_all",
                 return_df=TRUE)

