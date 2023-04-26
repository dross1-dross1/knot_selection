functions {
  real generalized_inverse_gaussian_lpdf(real x, int p,
                                        real a, real b) {
    return p * 0.5 * log(a / b)
      - log(2 * modified_bessel_second_kind(p, sqrt(a * b)))
      + (p - 1) * log(x)
      - (a * x + b / x) * 0.5;
 }
 
 matrix approx_L(int M, real scale, real[] xt, real sigma, real l) {
    int N = size(xt);
  
    real epsilon = sqrt(1 / (2 * l^2));
    real alpha = 1 / scale;
    real beta = (1.0 + (2.0 * epsilon / alpha)^2)^0.25;
    real delta = sqrt(alpha^2 * (beta^2 - 1.0) / 2.0);
    
    vector[N] Ht[M];
    vector[N] x = to_vector(xt);
    matrix[N, M] Htt;
    vector[N] xp = alpha * beta * x;
    real f = sqrt(epsilon^2 / (alpha^2 + delta^2 + epsilon^2));
    
    Ht[1] = sqrt(sqrt(alpha^2 / (alpha^2 + delta^2 + epsilon^2))) * sqrt(beta) * exp(-delta^2 * x .* x);
    Ht[2] = f * sqrt(2.0) * xp .* Ht[1];
    for(n in 3:M) {
      Ht[n] = f * sqrt(2.0 / (n - 1.0)) * xp .* Ht[n - 1] - f^2 * sqrt((n - 2.0) / (n - 1.0)) * Ht[n - 2];
    }
    
    for(n in 1:M) {
      Htt[:, n] = Ht[n];
    }
    
    return sigma * Htt;
  }
}

data {
  int<lower=0> N_knots;
  array[N_knots] vector[2] knot_locs;
  
  int<lower=0> fixed_rank_dim; // Dimension to project covariance matrix into
  
  int<lower=0> N_spatial; // Number of sites at which spatial process is measured
  array[N_spatial] vector[2] spatial_locs; // x-y coordinates of spatial process
  array[N_spatial] real y_spatial; // Measured value of spatial process at each site
}

transformed data {
  // Compute the maximum distance between any two points. This is used
  // for spatial dependence priors.
  real max_dist = 0;
  real temp_dist;
  for (i in 1:N_spatial) {
    for (j in (i+1):N_spatial) {
      temp_dist = distance(spatial_locs[i], spatial_locs[j]);
      max_dist = max([max_dist, temp_dist]);
    }
  }
}

parameters {
  // Intercept
  real mu;
  
  // Kernel (ellipse) process
  // real<lower=0> nugget_psi;
  real<lower=0> sigma_psi; // Variance of foci process
  real<lower=0> tau_psi; // Range of ellipse dependence
  real<lower=0> ell_psi; // Range of ellipse dependence
  real<lower=0> nugget_psi;
  vector<lower=0>[N_knots] psi_x; // x component of foci
  vector<lower=0>[N_knots] psi_y; // y component of foci
  
  // Kernel interpolation process
  real<lower=0> sigma_interp;
  real<lower=0> ell_interp;
  
  // (Latent) spatial process
  // real<lower=0> nugget_z; // Precision of spatial error term
  real<lower=0> sigma_z; // Variance of spatial process
  real<lower=0> ell_z; // Range of spatial process
  
  real<lower=0> global_scale; // Low-rank expansion scale factor
  
  // Overall spatial process
  vector[10] eta;
  real<lower=0> lambda_y; // Precision of model error term
}

transformed parameters {
  // Covariance matrix for foci
  matrix[N_knots, N_knots] R_psi;
  
  profile("R_psi_covariance_mat") {
    R_psi = gp_exp_quad_cov(knot_locs, sigma_psi, ell_psi);
    // Used in multi_normal priors for foci
    for (i in 1:N_knots) {
      R_psi[i, i] = R_psi[i, i] + nugget_psi;
    }
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
    target += std_normal_lpdf(tau_psi);
    
    // sigma_interp ~ std_normal();
    target += std_normal_lpdf(sigma_interp);
    // ell_interp ~ inv_gamma(3, 0.1);
    target += inv_gamma_lpdf(ell_interp | 1, 0.1);
    
    // sigma_z ~ std_normal();
    target += normal_lpdf(sigma_z | 0, 0.1);
    // ell_z ~ inv_gamma(3, 0.1);
    target += inv_gamma_lpdf(ell_z | 1, 0.1);
    
    target += std_normal_lpdf(global_scale);
    
    // eta ~ std_normal();
    target += std_normal_lpdf(eta);
    // lambda_y ~ std_normal();
    target += std_normal_lpdf(lambda_y);
  }
  
  matrix[N_knots, N_knots] R_psi_chol;

  profile("psi_multi_normal_cholesky") {
    R_psi_chol = cholesky_decompose(R_psi);
    // psi_x ~ multi_normal_cholesky(rep_vector(0., N_knots), R_psi_chol);
    target += multi_normal_cholesky_lpdf(psi_x | rep_vector(0., N_knots), R_psi_chol);
    // psi_y ~ multi_normal_cholesky(rep_vector(0., N_knots), R_psi_chol);
    target += multi_normal_cholesky_lpdf(psi_y | rep_vector(0., N_knots), R_psi_chol);
  }
  
  matrix[N_spatial, N_knots] W_interp;
  vector[N_spatial] psi_x_all;
  vector[N_spatial] psi_y_all;
  
  profile("W_interp") {
    W_interp = gp_exp_quad_cov(spatial_locs, knot_locs, sigma_interp, ell_interp);
    psi_x_all = W_interp * psi_x;
    psi_y_all = W_interp * psi_y;
  }
  
  // Compute rotations (sampling both psi_x and psi_y near zero is problematic)
  vector[N_spatial] alpha = atan(psi_y_all ./ psi_x_all);
  for (i in 1:N_spatial) {
    if (is_nan(alpha[i])) {
      alpha[i] = 0;
    }
  }
  
  // Compute rotation matrices
  array[N_spatial] matrix[2, 2] Sigma_sqrts;
  profile("Sigma_sqrts") {
    for (i in 1:N_spatial) {
      matrix[2, 2] rotation = [[cos(alpha[i]), sin(alpha[i])], [-sin(alpha[i]), cos(alpha[i])]];
      row_vector[2] ellipse_scale = [psi_x_all[i], psi_y_all[i]];
      // // Component of covariance matrix that handles ellipse rotation
      // matrix[2, 2] rotation = [[cos(alpha[i]), sin(alpha[i])], [-sin(alpha[i]), cos(alpha[i])]];
      // Full kernel covariance matrix (square root)
      Sigma_sqrts[i] = tau_psi * diag_post_multiply(rotation, ellipse_scale);
    }
  }

  // Construct kernel covariance matrices at each site
  array[N_spatial] real spatial_locs_transformed_norm;
  profile("spatial_transform") {
    for (i in 1:N_spatial) {
      // // Component of covariance matrix that handles ellipse scaling
      // row_vector[2] ellipse_scale = [psi_x_all[i], psi_y_all[i]];
      // // Component of covariance matrix that handles ellipse rotation
      // matrix[2, 2] rotation = [[cos(alpha[i]), sin(alpha[i])], [-sin(alpha[i]), cos(alpha[i])]];
      // Full kernel covariance matrix (square root)
      // matrix[2, 2] Sigma_sqrt = tau_psi * diag_post_multiply(rotations[i], ellipse_scale);
      // matrix[2, 2] Sigma_sqrt = tau_psi * [[sqrt(psi_x_all[i])*cos(alpha[i]), sqrt(psi_y_all[i])*sin(alpha[i])], [-sqrt(psi_x_all[i])*sin(alpha[i]), sqrt(psi_y_all[i])*cos(alpha[i])]];
      // Full kernel covariance matrix (matrix times its transpose)
      // Equivalent to Sigma_sqrt * Sigma_sqrt'
      // matrix[2, 2] Sigma = tcrossprod(Sigma_sqrt);
      // matrix[2, 2] Sigma = tcrossprod(tau_psi * [[sqrt(psi_x_all[i])*cos(alpha[i]), sqrt(psi_y_all[i])*sin(alpha[i])], [-sqrt(psi_x_all[i])*sin(alpha[i]), sqrt(psi_y_all[i])*cos(alpha[i])]]);
      // Transform each location according to its elliptical covariance matrix
      // spatial_locs_transformed[i] = tau_psi * [[sqrt(psi_x_all[i])*cos(alpha[i]), sqrt(psi_y_all[i])*sin(alpha[i])], [-sqrt(psi_x_all[i])*sin(alpha[i]), sqrt(psi_y_all[i])*cos(alpha[i])]] * spatial_locs[i];
      // spatial_locs_transformed[i] = tau_psi * [sqrt(psi_x_all[i])*cos(alpha[i])*spatial_locs[i][1] + sqrt(psi_y_all[i])*sin(alpha[i])*spatial_locs[i][2], -sqrt(psi_x_all[i])*sin(alpha[i])*spatial_locs[i][1] + sqrt(psi_y_all[i])*cos(alpha[i])*spatial_locs[i][2]]'
      spatial_locs_transformed_norm[i] = norm2(Sigma_sqrts[i] * spatial_locs[i]);
    }
  }

  // Reference: // Reference: https://discourse.mc-stan.org/t/approximate-gps-with-spectral-stuff/1041
  vector[N_spatial] f = approx_L(10, global_scale, spatial_locs_transformed_norm, sigma_z, ell_z) * eta;

  // // Latent variable formulation
  // // Reference: https://mc-stan.org/docs/stan-users-guide/fit-gp.html
  // matrix[N_spatial, N_spatial] R_z;
  // matrix[N_spatial, N_spatial] R_z_chol;
  // matrix[N_spatial, N_spatial] C_lp;
  // matrix[N_spatial, N_spatial] C_lp_chol;
  // vector[N_spatial] f;
  // {
  //   // Compute covariance matrix between transformed coordinates (simplified version)
  //   profile("R_z") {
  //     R_z = gp_exp_quad_cov(spatial_locs_transformed, 0.1, ell_z);
  //     R_z = add_diag(R_z, 1e-10);
  //   }
  // 
  //   // // === Banerjee linear projection method ===
  //   // // Project down to lower dimension
  //   // profile("C_lp") {
  //   //   // C_lp = (Omega_Phi * R_z)' / quad_form(R_z, Omega_Phi') * (Omega_Phi * R_z);
  //   //   C_lp = quad_form(inverse(quad_form(R_z, Omega_Phi')), Omega_Phi * R_z);
  //   // }
  //   // // matrix[N_spatial, N_spatial] C_lp = quad_form(Omega_Phi * R_z, solve(Omega_Phi * R_z * Omega_Phi'));
  //   // // matrix[N_spatial, fixed_rank_dim] C_lp = (Omega_Phi * R_z)' * solve(Omega_Phi * R_z * Omega_Phi') * Omega_Phi * R_z;
  //   // 
  //   // // Add a nugget term (correct for under-estimating variance)
  //   // profile("C_lp_nugget") {
  //   //   C_lp = add_diag(C_lp, diagonal(R_z) - diagonal(C_lp));
  //   //   for (n in 1:N_spatial) {
  //   //     C_lp[n, n] = C_lp[n, n] + 0.01;
  //   //   }
  //   // }
  //   // 
  //   // // Cholesky decomposition for later efficiency
  //   // profile("C_lp_chol") {
  //   //   C_lp_chol = cholesky_decompose(C_lp);
  //   // }
  //   // === End linear projection method ===
  // 
  //   // Latent variable formulation uses Cholesky decomposition (use this _or_ Banerjee linear projection)
  //   profile("R_z_chol") {
  //     R_z_chol = cholesky_decompose(R_z);
  //   }
  //   f = R_z_chol * eta;
  // }

  // // === Sampling from linear model ===
  profile("model_sampling") {
    // y_spatial ~ normal(mu + f, lambda_y);
    target += normal_lpdf(y_spatial | mu + f, lambda_y);
  }
  // // 
  // // // int grainsize = 1;
  // // // profile("reduce_sum_model") {
  // // //   target += reduce_sum(partial_sum_lpdf, y_spatial, 
  // // //                       grainsize,
  // // //                       mu, f, lambda_y);
  // // // }
}

generated quantities {
  matrix[N_spatial, N_knots] W_interp = gp_exp_quad_cov(spatial_locs, knot_locs, sigma_interp, ell_interp);

  vector[N_spatial] psi_x_all = W_interp * psi_x;
  vector[N_spatial] psi_y_all = W_interp * psi_y;
}
