# Z(x) = Z_0(x) + sqrt(alpha) * sum[K(x-r_i)*Z_i(x)]
# Z_i = stationary GP with Matern covariance
# alpha = real, positive
# r_i = knot
# K = Epanechnikov kernel

# Parameters
# Matern:
# - tau2 = diagonal variance (nugget)
# - sigma = scale
# - rho = lengthscale
# - eta = variance (smoothness)
#
# Kernel:
# - h = bandwidth parameter

library(rSPDE)
library(fields)


# Helper functions --------------------------------------------------------
# Combining GPs -----------------------------------------------------------
# Epanechnikov kernel (distance function)
epan.kernel = function(u, h) {
  (2 / pi) * (1 / h**2) * (1-(u/h)**2)
}

# Create background and local spatial processes
generate.gp = function(X.knots, X, y) {
  N.KNOTS = nrow(X.knots)
  RADIUS = 0.2

  # Empty list to hold all local spatial processes
  gp.list = list()

  for (knot.id in 1:N.KNOTS) {
    # Get knot coordinates
    r = X.knots[knot.id, ] %>% matrix(ncol=ncol(X.knots))

    # Find distance from all points to current knot
    D = distance(r, X)

    # Get points within circle of radius RADIUS
    r.idx = which(D < RADIUS)
    r.nbhd = X[r.idx, ]
    r.nbhd.y = y[r.idx]

    # Fit spatial processes
    gp.fit = spatialProcess(x=r.nbhd, y=r.nbhd.y)

    # Store fitted spatial process
    gp.list[[knot.id]] = gp.fit
  }

  # Fit background GP
  gp.background = spatialProcess(x=X, y=y)

  return(list("gp.background"=gp.background,
              "gp.list"=gp.list))
}

# Generate simulations from all spatial processes (local + background)
simulate.gp = function(gp.list,
                       gp.background,
                       X) {

  # Number of samples (used for estimating variability/calculating entropy)
  num_samples = 100

  # Simulate from background spatial process
  sim.background = sim.spatialProcess(gp.background, X, M=num_samples)

  # Simulate from local spatial processes
  sim.local = list()
  for (i in 1:length(gp.list)) {
    sim.gp = sim.spatialProcess(gp.list[[i]], X, M=num_samples)
    sim.local[[i]] = sim.gp
  }

  return(list("sim.background"=sim.background,
              "sim.local"=sim.local))
}

# Generate predictions using total spatial process (background + local)
gp.preds = function(sim.background, sim.local.list, X.knots, X.obs, sqrt_alpha, ell) {
  # Dimensions:
  # sim.background: (num locations) x (num samples generated)
  # sim.local: same
  # epan.kernel(X.knots, X.obs): (num knots) x (num locations)

  # Holds predictions for each sample
  sim.total = matrix(NA, nrow=nrow(sim.background), ncol=ncol(sim.background))

  # Go sample-by-sample
  for (j in 1:ncol(sim.background)) {
    # Just the j'th sample from background spatial process
    sim.background.j = sim.background[, j]
    # Empty vector to hold the sum of all local spatial processes
    sim.local.sum.j = 0 * sim.background.j
    # Loop over spatial processes
    for (i in 1:length(sim.local.list)) {
      # Just the j'th sample of the i'th local spatial process
      sim.local.i.j = sim.local.list[[i]][, j]
      sim.local.sum.j = sim.local.sum.j +
        epan.kernel(u=distance(matrix(X.knots[i, ], ncol=2), X.obs), h=ell) * sim.local.i.j
    }
    # Combine local and background predictions into one row (all locations simultaneously)
    sim.total[, j] = sim.background.j + sqrt_alpha * sim.local.sum.j
  }

  return(sim.total)
}

# Function passed to optim to fit sqrt_alpha and ell
fit.params = function(sim.background, sim.local.list, X.knots, X.train, y.train, par) {
  # Extract parameters of interest
  sqrt_alpha = par[1]
  ell = par[2]

  sim.total = gp.preds(sim.background=sim.background,
                       sim.local.list=sim.local.list,
                       X.knots=X.knots,
                       X.obs=X.train,
                       sqrt_alpha=sqrt_alpha,
                       ell=ell)

  # Replicate actual values along each column to compute MSE
  y.train.matrix = matrix(y.train, nrow=length(y.train), ncol=ncol(sim.background))

  # Compute MSE
  mse = sum((sim.total - y.train.matrix)**2 / length(y.train.matrix))
  return(mse)
}

# Fitting local GPs using training data ------------------------------------------------------------
fit.fuentes.gp = function(X.knots, X.train, y.train) {
  N.KNOTS = nrow(X.knots)

  gp.all = generate.gp(X.knots, X.train, y.train)
  gp.background = gp.all$gp.background
  gp.list = gp.all$gp.list

  # Generate predictions from all processes
  sim.train.values = simulate.gp(gp.list, gp.background, X.train)
  # Save background and local processes separately for easy reference
  sim.train.background = sim.train.values$sim.background
  sim.train.local.list = sim.train.values$sim.local

  params = optim(par=c(1, 1), fn=fit.params,
                 method="L-BFGS",
                 lower=c(0, 0),
                 upper=c(10, 1),
                 sim.background=sim.train.background,
                 sim.local.list=sim.train.local.list,
                 X.knots=X.knots,
                 X.train=X.train,
                 y.train=y.train)

  sqrt_alpha = params$par[1]
  ell = params$par[2]

  return(list("sqrt_alpha"=sqrt_alpha,
              "ell"=ell,
              "gp.background"=gp.background,
              "gp.list"=gp.list))
}

predict.fuentes.gp = function(X.knots, X, params) {
  gp.background = params$gp.background
  gp.list = params$gp.list

  # Generate predictions from all processes
  sim.values = simulate.gp(gp.list, gp.background, X)
  # Save background and local processes separately for easy reference
  sim.background = sim.values$sim.background
  sim.local.list = sim.values$sim.local

  best_fit_preds = gp.preds(sim.background=sim.background,
                           sim.local.list=sim.local.list,
                           X.knots=X.knots,
                           X.obs=X,
                           sqrt_alpha=params$sqrt_alpha,
                           ell=params$ell)

  best_fit_preds.mean = apply(best_fit_preds, 1, mean)
  best_fit_preds.median = apply(best_fit_preds, 1, median)
  # Bootstrap standard error
  best_fit_preds.se = apply(best_fit_preds, 1, sd) / sqrt(ncol(sim.background))

  best_fit_preds.df = data.frame(mean_pred=best_fit_preds.mean,
                                 median_pred=best_fit_preds.median,
                                 se_pred=best_fit_preds.se)
  return(best_fit_preds.df)
}
