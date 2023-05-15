functions {
  real generalized_inverse_gaussian_lpdf(real x, int p,
                                        real a, real b) {
    return p * 0.5 * log(a / b)
      - log(2 * modified_bessel_second_kind(p, sqrt(a * b)))
      + (p - 1) * log(x)
      - (a * x + b / x) * 0.5;
  }
}

data {
  int<lower=0> N_knots;
  array[N_knots] vector[2] knot_locs;
  
  int<lower=0> N_spatial; // Number of sites at which spatial process is measured
  array[N_spatial] vector[2] spatial_locs; // x-y coordinates of spatial process
  array[N_spatial] real y_spatial; // Measured value of spatial process at each site
}

transformed data {
  // Compute the maximum distance between any two points. This is used
  // for spatial dependence priors.
  real max_dist = 1;
}

parameters {
  // Intercept
  real mu;
  
  // Spatial process
  real<lower=0> nugget_psi; // Precision of spatial error term
  real<lower=0> sigma_psi; // Variance of spatial process
  real<lower=0> ell_psi; // Range of spatial process
  
  // Overall spatial process
  vector[N_knots] eta;
  real<lower=0> lambda_y; // Precision of model error term
}

transformed parameters {
  // cov.space_knot = cov_exp_quad(x, sigma, ell, knots)
  // cov.knot_knot = cov_exp_quad(knots, sigma, ell)
  // cov.knot_knot.inv = solve(cov.knot_knot)
  
  // # Low-rank covariance matrix
  // cov.pp = cov.space_knot %*% cov.knot_knot.inv %*% t(cov.space_knot)
  
  // Covariance matrix for foci
  matrix[N_spatial, N_knots] R_space_knot;
  matrix[N_knots, N_knots] R_knot_knot;
  matrix[N_spatial, N_spatial] R_low_rank;
  
  // // -----------
  // // latent gp at knots
  // for (i in 1:(m-1)) {
  //   for (j in (i + 1):m) {
  //     Cstar[i, j] = eta_sq * exp(-D_star[i, j] * phi);
  //     Cstar[j, i] = Cstar[i, j];
  //   }
  // }
  // 
  // for (k in 1:m) Cstar[k, k] = eta_sq + sig_sq;
  // inv_Cstar = inverse(Cstar);
  // w_star = cholesky_decompose(Cstar) * w_z;
  // 
  // // latent gp at sample locations
  // C_site_star = eta_sq * exp(-D_site_star * phi);
  // C_ss_inv_Cstar = C_site_star * inv_Cstar;
  // w = C_site_star * inv_Cstar * w_star;
  // // -----------
  
  vector[N_spatial] f;
  vector[N_knots] w_star;
  profile("R_psi_covariance_mat") {
    R_knot_knot = gp_exp_quad_cov(knot_locs, sigma_psi, ell_psi);
    R_knot_knot = add_diag(R_knot_knot, nugget_psi);
    matrix[N_knots, N_knots] R_knot_knot_chol = cholesky_decompose(R_knot_knot);
    w_star = R_knot_knot_chol * eta;
    R_space_knot = gp_exp_quad_cov(spatial_locs, knot_locs, sigma_psi, ell_psi);
    f = R_space_knot / R_knot_knot * w_star;
  }

}

model {
  // Priors
  profile("simple_priors") {
    target += std_normal_lpdf(mu);
  
    // Reference for priors: https://mc-stan.org/docs/stan-users-guide/fit-gp.html
    // nugget_psi ~ lognormal(0, 1);
    target += lognormal_lpdf(nugget_psi | 0, 1);
    // sigma_psi ~ std_normal();
    target += std_normal_lpdf(sigma_psi);
    // ell_psi ~ normal(0, max_dist/3);
    target += normal_lpdf(ell_psi | 0, max_dist/3);
    // tau_psi ~ std_normal();
    // target += std_normal_lpdf(tau_psi);
    
    // eta ~ std_normal();
    target += std_normal_lpdf(eta);
    // lambda_y ~ std_normal();
    target += std_normal_lpdf(lambda_y);
  }
  
  // === Sampling from linear model ===
  profile("model_sampling") {
    // y_spatial ~ normal(mu + f, lambda_y);
    target += normal_lpdf(y_spatial | mu + f, lambda_y);
  }
}

generated quantities {
  vector[N_spatial] y_spatial_sim;
  for (i in 1:N_spatial) {
    y_spatial_sim[i] = normal_rng(mu + f[i], lambda_y);
  }

  vector[N_spatial] log_lik;
  for (i in 1:N_spatial) {
    log_lik[i] = normal_lpdf(y_spatial[i] | mu + f[i], lambda_y);
  }
}
