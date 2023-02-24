data {
  
  // main indices
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> Q;
  
  // response data  
  vector<lower=0>[N] y;

  // predictor variables  
  matrix[N, K] X;
  vector[N] arundo;

  // random effects
  int<lower=0> norigin;
  int<lower=0> nbasin;
  int<lower=0> nsite;
  int<lower=0> nblock;
  int<lower=1,upper=norigin> origin[N];
  int<lower=1,upper=nbasin> basin[N];
  int<lower=1,upper=nsite> site[N];
  int<lower=1,upper=nblock> block_id[N];

  // taxonomic info
  int<lower=0> ngenus;
  int<lower=0> nfamily;
  int<lower=0> norder;
  int<lower=1,upper=Q> species[N];
  int<lower=1,upper=ngenus> genus[N];
  int<lower=1,upper=nfamily> family[N];
  int<lower=1,upper=norder> order[N];

  // scale of main effects and overall variance
  real<lower=0> sigma_fixed;
  real<lower=0> sigma_random;
  
}

parameters {
  
  // linear predictor
  real zalpha_main;
  vector[Q] zalpha;
  vector[K] zbeta_main;
  matrix[Q, K] zbeta;  
  vector[K] zdelta_main;
  matrix[Q, K] zdelta;
  real ztheta_main;
  vector[Q] ztheta;

  // random effects
  vector[norigin] gamma_origin;
  vector[nbasin] gamma_basin;
  vector[nsite] gamma_site;
  vector[nblock] gamma_block;
  vector[norder] gamma_order;
  vector[nfamily] gamma_family;
  vector[ngenus] gamma_genus;

  // random effects variances
  real<lower=0> sigma_origin;
  real<lower=0> sigma_basin;
  real<lower=0> sigma_site;
  real<lower=0> sigma_block;
  real<lower=0> sigma_order;
  real<lower=0> sigma_family;
  real<lower=0> sigma_genus;
  
  // and for regression coefs
  real<lower=0> sigma_alpha;
  vector<lower=0>[K] sigma_beta;
  vector<lower=0>[K] sigma_delta;
  real<lower=0> sigma_theta;
  
  // observation variance
  real<lower=0> sigma_main;

}

transformed parameters {
  real alpha_main;
  vector[Q] alpha;
  vector[N] beta_term;
  vector[N] delta_term;
  real theta_main;
  vector[Q] theta;
  vector[N] mu;

  // rescale z-transformed covariates
  alpha_main = sigma_fixed * zalpha_main;
  alpha = sigma_random * sigma_alpha * zalpha;
  theta_main = sigma_fixed * ztheta_main;
  theta = sigma_random * sigma_theta * ztheta;

  // rescale and calculate beta (covariate) effects
  beta_term = rep_vector(0.0, N);
  delta_term = rep_vector(0.0, N);
  for (k in 1:K) {
    beta_term += (sigma_fixed * zbeta_main[k] +
      sigma_random * sigma_beta[k] * zbeta[species, k]) .* X[, k];
    delta_term += (sigma_fixed * zdelta_main[k] +
      sigma_random * sigma_delta[k] * zdelta[species, k]) .* arundo .* X[, k];
  }

  // calculate linear predictor
  mu = alpha_main + alpha[species] +
    beta_term +
    delta_term +
    (theta_main + theta[species]) .* arundo +
    sigma_random * 
    (sigma_origin * gamma_origin[origin] +
     sigma_basin * gamma_basin[basin] +
     sigma_site * gamma_site[site] +
     sigma_block * gamma_block[block_id] +
     sigma_order * gamma_order[order] + 
     sigma_family * gamma_family[family] + 
     sigma_genus * gamma_genus[genus]);

}

model {

  // z-scaled priors for regression terms
  zalpha_main ~ std_normal();
  zalpha ~ std_normal();
  zbeta_main ~ std_normal();
  to_vector(zbeta) ~ std_normal();
  zdelta_main ~ std_normal();
  to_vector(zdelta) ~ std_normal();
  ztheta_main ~ std_normal();
  ztheta ~ std_normal();

  // priors for random effects
  gamma_origin ~ std_normal();
  gamma_basin ~ std_normal();
  gamma_site ~ std_normal();
  gamma_block ~ std_normal();
  gamma_order ~ std_normal();
  gamma_family ~ std_normal();
  gamma_genus ~ std_normal();

  // and their variances
  sigma_origin ~ std_normal();
  sigma_basin ~ std_normal();
  sigma_site ~ std_normal();
  sigma_block ~ std_normal();
  sigma_order ~ std_normal();
  sigma_family ~ std_normal();
  sigma_genus ~ std_normal();
  
  // and vars for regression coefs
  sigma_alpha ~ std_normal();
  sigma_beta ~ std_normal();
  sigma_delta ~ std_normal();
  sigma_theta ~ std_normal();
  
  // and observation variance
  sigma_main ~ std_normal();

  // calculate likelihood of response (lognormal)
  target += lognormal_lpdf(y | mu, sigma_main);
   
}

generated quantities {
  vector<lower=0>[N] ypred;
  
  // random draws from the posterior of mu and phi for posterior checks
  for (i in 1:N) {
    ypred[i] = lognormal_rng(mu[i], sigma_main);
  }
  
}
