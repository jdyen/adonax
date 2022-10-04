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
  vector[N] ftle;
  int<lower=0> norigin;
  int<lower=0,upper=norigin> origin[N];
  
  // random effects
  int<lower=0> nbasin;
  int<lower=0> nsite;
  int<lower=0> nblock;
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

transformed data {
  vector[N] yp1;
  
  // add 1 to y
  yp1 = y + 1;
  
}

parameters {
  
  // linear predictor
  real zalpha_main;
  vector[Q] zalpha;
  vector[K] zbeta_main;
  matrix[Q, K] zbeta;  
  vector[K] zbeta_main_arundo;
  matrix[Q, K] zbeta_arundo;
  real ztheta_main;
  vector[Q] ztheta;
  real zgamma;
  real zgamma_arundo;
  vector[norigin] zdelta_main;
  matrix[Q, norigin] zdelta;

  // random effects
  vector[nbasin] gamma_basin;
  vector[nsite] gamma_site;
  vector[nblock] gamma_block;
  vector[norder] gamma_order;
  vector[nfamily] gamma_family;
  vector[ngenus] gamma_genus;

  // random effects variances
  real<lower=0> sigma_basin;
  real<lower=0> sigma_site;
  real<lower=0> sigma_block;
  real<lower=0> sigma_order;
  real<lower=0> sigma_family;
  real<lower=0> sigma_genus;
  
  // observation variance
  real<lower=0> sigma_main;

}

transformed parameters {
  real alpha_main;
  vector[Q] alpha;
  vector[K] beta_main;
  matrix[Q, K] beta;  
  vector[K] beta_main_arundo;
  matrix[Q, K] beta_arundo;
  vector[N] beta_term;
  real theta_main;
  vector[Q] theta;
  real gamma;
  real gamma_arundo;
  vector[norigin] delta_main;
  matrix[Q, norigin] delta;
  vector[N] delta_term;
  vector[N] mu;

  // rescale z-transformed covariates
  alpha_main = sigma_fixed * zalpha_main;
  alpha = sigma_fixed * zalpha;
  beta_main = sigma_fixed * zbeta_main;
  beta = sigma_fixed * zbeta;
  beta_main_arundo = sigma_fixed * zbeta_main_arundo;
  beta_arundo = sigma_fixed * zbeta_arundo;
  theta_main = sigma_fixed * ztheta_main;
  theta = sigma_fixed * ztheta;
  gamma = sigma_fixed * zgamma;
  gamma_arundo = sigma_fixed * zgamma_arundo;
  delta_main = sigma_fixed * zdelta_main;
  delta = sigma_fixed * zdelta;

  // calculate species-specific covariate effects
  beta_term = rep_vector(0.0, N);
  for (k in 1:K) {
    beta_term += (
      beta_main[k] + beta[species, k] +
      arundo .* (beta_main_arundo[k] + beta_arundo[species, k])
    ) .* X[, k];
  }

  // calculate species-specific origin term
  for (i in 1:N)
    delta_term[i] = delta_main[origin[i]] + delta[species[i], origin[i]];

  // calculate linear predictor
  mu = alpha_main + alpha[species] +
    beta_term +
    (theta_main + theta[species]) .* arundo +
    (gamma + arundo * gamma_arundo) .* ftle +
    delta_term +
    sigma_random * 
    (sigma_basin * gamma_basin[basin] +
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
  zbeta_main_arundo ~ std_normal();
  to_vector(zbeta_arundo) ~ std_normal();
  ztheta_main ~ std_normal();
  ztheta ~ std_normal();
  zgamma ~ std_normal();
  zgamma_arundo ~ std_normal();
  zdelta_main ~ std_normal();
  to_vector(zdelta) ~ std_normal();
  
  // priors for random effects
  to_vector(gamma_basin) ~ std_normal();
  to_vector(gamma_site) ~ std_normal();
  to_vector(gamma_block) ~ std_normal();
  to_vector(gamma_order) ~ std_normal();
  to_vector(gamma_family) ~ std_normal();
  to_vector(gamma_genus) ~ std_normal();

  // and their variances
  sigma_basin ~ std_normal();
  sigma_site ~ std_normal();
  sigma_block ~ std_normal();
  sigma_order ~ std_normal();
  sigma_family ~ std_normal();
  sigma_genus ~ std_normal();
  
  // and observation variance
  sigma_main ~ std_normal();

  // calculate likelihood of response (zero-inflataed Poisson)
  target += lognormal_lpdf(yp1 | mu, sigma_main);
   
}
