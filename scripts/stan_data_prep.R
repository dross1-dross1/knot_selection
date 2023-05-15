library(cmdstanr)
library(dplyr)
library(ggplot2)


# Prepare data for the stan model
prepare_data = function(df, M, num_knots=NULL, sample_pct=1, plot=FALSE, return_ellipses=TRUE, return_predictions=FALSE, return_log_likelihoods=TRUE, df.sample=NULL, df.knots=NULL) {
  # Reduce sample size if needed
  if (is.null(df.sample)) {
    if (sample_pct < 1) {
      df.sample.idx = sample(1:nrow(df), as.integer(sample_pct*nrow(df)))
      df.sample = df[df.sample.idx, ]
      rownames(df.sample) = 1:nrow(df.sample)
    } else {
      df.sample = df
    }
  }
  
  # Spatial process sites and data
  spatial_locs = df.sample[, c('x', 'y')]
  y_spatial = df.sample$pm2_5 %>% as.vector
  N_spatial = nrow(spatial_locs)
  
  # Package data
  data = list(N_spatial=N_spatial,
              spatial_locs=spatial_locs,
              y_spatial=y_spatial,
              M=M,
              return_ellipses=return_ellipses,
              return_predictions=return_predictions,
              return_log_likelihoods=return_log_likelihoods)
  
  # Include knots if needed
  if (!is.null(num_knots)) {
    # Randomly select num_knots
    df.knots.idx = sample(1:nrow(df.sample), size=num_knots)
    df.knots = df.sample[df.knots.idx, ]
  }
  
  data$N_knots = nrow(df.knots)
  data$knot_locs = df.knots[, c('x', 'y')]
  
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
