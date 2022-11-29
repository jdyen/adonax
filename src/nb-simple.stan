data {
  
  // main indices
  int<lower=0> N;
  int<lower=0> Q;
  int<lower=0> K;

  // response data  
  int<lower=0> nflat;
  int<lower=0> yflat[nflat];

  // scaling factor for response variable (pseudo-effort)
  real scale_factor;

  // predictor variables  
  matrix[K, N] X;
  row_vector[N] arundo;

  // random effects
  int<lower=0> nbasin;
  int<lower=0> nblock;
  int<lower=0> nsite;
  int<lower=1,upper=nbasin> basin[N];
  int<lower=1,upper=nblock> block_id[N];
  int<lower=1,upper=nsite> site[N];
  
  // scale of main effects and overall variance
  real<lower=0> sigma_fixed;
  real<lower=0> sigma_resid;
  real<lower=0> sigma_random;

}

transformed data {
  matrix[Q, N] log_scale_factor;

  // expand and log-transform scale_factor
  log_scale_factor = rep_matrix(log(scale_factor), Q, N);
  
}

parameters {
  
  // linear predictor
  real zalpha_main;
  vector[Q] zalpha;
  row_vector[K] zbeta_main;
  matrix[Q, K] zbeta;
  real ztheta_main;
  vector[Q] ztheta;
  
  // pooling variances for coefficients
  real<lower=0> sigma_alpha;
  // row_vector<lower=0>[K] sigma_beta;
  real<lower=0> sigma_theta;

  // random effects
  matrix[Q, nbasin] gamma_basin;
  matrix[Q, nblock] gamma_block;
  matrix[Q, nsite] gamma_site;

  // random effects variances
  real<lower=0> sigma_main_basin;
  real<lower=0> sigma_main_block;
  real<lower=0> sigma_main_site;
  vector<lower=0>[Q] sigma_basin;
  vector<lower=0>[Q] sigma_block;
  vector<lower=0>[Q] sigma_site;

  // over-dispersion parameter
  vector<lower=0>[Q] phi;

}

transformed parameters {
  vector[Q] alpha;
  matrix[Q, K] beta;
  vector[Q] theta;
  matrix[Q, N] mu;
  vector[nflat] mu_flat;
  vector[nflat] phi_flat;

  // rescale z-transformed covariates
  alpha = sigma_fixed * zalpha_main + sigma_random * sigma_alpha * zalpha;
  beta = sigma_fixed * rep_matrix(zbeta_main, Q) +
    sigma_random * rep_matrix(sigma_beta, Q) .* zbeta;
  theta = sigma_fixed * ztheta_main + sigma_random * sigma_theta * ztheta;

  // calculate linear predictor
  mu = rep_matrix(alpha, N) +
    // beta * X +
    theta * arundo +
    sigma_random *
    (sigma_main_basin * rep_matrix(sigma_basin, N) .* gamma_basin[, basin] +
     sigma_main_block * rep_matrix(sigma_block, N) .* gamma_block[, block_id] +
     sigma_main_site * rep_matrix(sigma_site, N) .* gamma_site[, site]
    ) +
    log_scale_factor;

  // flatten mu and theta_zero for likelihood calcs
  mu_flat = to_vector(mu);
  phi_flat = to_vector(rep_matrix(phi, N));

}

model {

  // z-scaled priors for regression terms
  zalpha_main ~ std_normal();
  zalpha ~ std_normal();
  zbeta_main ~ std_normal();
  to_vector(zbeta) ~ std_normal();
  ztheta_main ~ std_normal();
  ztheta ~ std_normal();

  // regression term pooling variances
  sigma_alpha ~ std_normal();
  sigma_beta ~ std_normal();
  sigma_theta ~ std_normal();
  
  // priors for random effects
  to_vector(gamma_basin) ~ std_normal();
  to_vector(gamma_block) ~ std_normal();
  to_vector(gamma_site) ~ std_normal();

  // and their variances
  sigma_basin ~ std_normal();
  sigma_block ~ std_normal();
  sigma_site ~ std_normal();

  // prior for phi
  real log_half = -0.693147180559945286;  
  target += student_t_lpdf(phi | 3.0, 0, 1) - log_half;

  // calculate likelihood of response (NB model)
  yflat ~ neg_binomial_2_log(mu_flat, phi_flat);

}

generated quantities {
  int<lower=0> ypred[nflat];
  
  // random draws from the posterior of mu and phi for posterior checks
  for (i in 1:nflat) {
    real mu_pred = mu_flat[i]; 
    if (mu_pred > 15.) mu_pred = 15.;
    ypred[i] = neg_binomial_2_log_rng(mu_pred, phi_flat[i]);
  }
  
}
