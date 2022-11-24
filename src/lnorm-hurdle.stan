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
  vector[Q] zalpha;
  matrix[Q, K] zbeta;
  vector[Q] ztheta;

  // random effects
  matrix[Q, nbasin] gamma_basin;
  matrix[Q, nblock] gamma_block;
  matrix[Q, nsite] gamma_site;

  // random effects variances
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
  alpha = sigma_fixed * zalpha;
  beta = sigma_fixed * zbeta;
  theta = sigma_fixed * ztheta;

  // calculate linear predictor
  mu = rep_matrix(alpha, N) +
    beta * X +
    theta * arundo +
    sigma_random *
    (rep_matrix(sigma_basin, N) .* gamma_basin[, basin] +
     rep_matrix(sigma_block, N) .* gamma_block[, block_id] +
     rep_matrix(sigma_site, N) .* gamma_site[, site]
    ) +
    log_scale_factor;

  // flatten mu and theta_zero for likelihood calcs
  mu_flat = to_vector(mu);
  phi_flat = to_vector(rep_matrix(phi, N));

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

  // and observation variance
  sigma_main ~ std_normal();

  // calculate likelihood of response (zero-inflataed Poisson)
  target += lognormal_lpdf(yp1 | mu, sigma_main);

}
