data {
  
  // main indices
  int<lower=0> N;
  int<lower=0> Q;
  int<lower=0> K;

  // response data  
  int<lower=0> nflat;
  int<lower=0> yflat[nflat];
  int<lower=0> nzero;
  int<lower=0> notzero;
  int<lower=0> zero_idx[nzero];
  int<lower=0> nonzero_idx[notzero];
  
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
  int ones_vec_obs[nflat];
  int zeros_vec_obs[nflat];
  int ones_vec_zero_obs[nzero];
  int zeros_vec_zero_obs[nzero];
  int zeros_vec_nonzero_obs[notzero];
  int<lower=0> yflat_zero[nzero];
  int<lower=0> yflat_nonzero[notzero];
  matrix[Q, N] log_scale_factor;

  // define expanded ones and zeros vectors for all obs
  for (i in 1:nflat) {
    ones_vec_obs[i] = 1;
    zeros_vec_obs[i] = 0;
  }
  ones_vec_zero_obs = ones_vec_obs[zero_idx];
  zeros_vec_zero_obs = zeros_vec_obs[zero_idx];
  zeros_vec_nonzero_obs = zeros_vec_obs[nonzero_idx];
  
  // pull out zero and nonzero obs from yflat
  yflat_zero = yflat[zero_idx];
  yflat_nonzero = yflat[nonzero_idx];
  
  // expand and log-transform scale_factor
  log_scale_factor = rep_matrix(log(scale_factor), Q, N);
  
}

parameters {
  
  // linear predictor
  vector[Q] zalpha;
  matrix[Q, K] zbeta;
  vector[Q] ztheta;

  // residual term
  matrix[Q, N] eps;
  vector<lower=0>[Q] zsigma;
  
  // random effects
  matrix[Q, nbasin] gamma_basin;
  matrix[Q, nblock] gamma_block;
  matrix[Q, nsite] gamma_site;

  // random effects variances
  vector<lower=0>[Q] sigma_basin;
  vector<lower=0>[Q] sigma_block;
  vector<lower=0>[Q] sigma_site;

  // zero inflation parameter
  vector<lower=0,upper=1>[Q] theta_zero_phi;
  real<lower=0> theta_zero_lambda;
  matrix<lower=0,upper=1>[Q, nbasin] theta_zero;
  
}

transformed parameters {
  vector[Q] alpha;
  matrix[Q, K] beta;  
  vector[Q] theta;
  vector<lower=0>[Q] sigma;
  matrix[Q, N] mu;
  vector[nflat] mu_flat;
  matrix[Q, N] theta_zero_mat;
  vector[nflat] theta_zero_flat;
  vector[nzero] theta_zero_flat_zero;
  vector[notzero] theta_zero_flat_nonzero;
  vector[nzero] mu_flat_zero;
  vector[notzero] mu_flat_nonzero;

  // rescale z-transformed covariates
  alpha = sigma_fixed * zalpha;
  beta = sigma_fixed * zbeta;
  theta = sigma_fixed * ztheta;
  sigma = sigma_resid * zsigma;

  // calculate linear predictor
  mu = rep_matrix(alpha, N) +
    beta * X +
    theta * arundo +
    sigma_random *
    (rep_matrix(sigma_basin, N) .* gamma_basin[, basin] +
     rep_matrix(sigma_block, N) .* gamma_block[, block_id] +
     rep_matrix(sigma_site, N) .* gamma_site[, site]
    ) +
    log_scale_factor +
    rep_matrix(sigma, N) .* eps;

  // pull out zero-inflation param for each obs
  theta_zero_mat = theta_zero[, basin];

  // flatten mu and theta_zero for likelihood calcs
  mu_flat = to_vector(mu);
  theta_zero_flat = to_vector(theta_zero_mat);
  theta_zero_flat_zero = theta_zero_flat[zero_idx];
  theta_zero_flat_nonzero = theta_zero_flat[nonzero_idx];
  mu_flat_zero = mu_flat[zero_idx];
  mu_flat_nonzero = mu_flat[nonzero_idx];

}

model {

  // z-scaled priors for regression terms
  zalpha ~ std_normal();
  to_vector(zbeta) ~ std_normal();
  ztheta ~ std_normal();

  // priors for random effects
  to_vector(gamma_basin) ~ std_normal();
  to_vector(gamma_block) ~ std_normal();
  to_vector(gamma_site) ~ std_normal();

  // and their variances
  sigma_basin ~ std_normal();
  sigma_block ~ std_normal();
  sigma_site ~ std_normal();

  // priors for covariance matrix
  zsigma ~ std_normal();
  to_vector(eps) ~ std_normal();

  // define zero inflation parameter from river-level mean
  theta_zero_lambda ~ pareto(0.1, 1.5);
  theta_zero_phi ~ beta(1., 1.);
  for (i in 1:nbasin)
    theta_zero[, i] ~ beta(theta_zero_lambda * theta_zero_phi, theta_zero_lambda * (1 - theta_zero_phi));

  // calculate likelihood of response (zero-inflataed Poisson)
  target += log_sum_exp(
    bernoulli_lpmf(ones_vec_zero_obs | theta_zero_flat_zero),
    bernoulli_lpmf(zeros_vec_zero_obs | theta_zero_flat_zero) +
    poisson_log_lpmf(yflat_zero | mu_flat_zero)
    );
    target += bernoulli_lpmf(zeros_vec_nonzero_obs | theta_zero_flat_nonzero) +
    poisson_log_lpmf(yflat_nonzero | mu_flat_nonzero);
   
}
