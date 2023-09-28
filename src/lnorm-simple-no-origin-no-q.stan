data {
  
  // main indices
  int<lower=0> N;
  int<lower=0> K;

  // response data  
  vector<lower=0>[N] y;

  // predictor variables  
  matrix[N, K] X;
  vector[N] arundo;

  // random effects
  int<lower=0> nbasin;
  int<lower=0> nblock;
  int<lower=0> nsite;
  int<lower=1,upper=nbasin> basin[N];
  int<lower=1,upper=nblock> block_id[N];
  int<lower=1,upper=nsite> site[N];
  
  // scale of main effects and overall variance
  real<lower=0> sigma_fixed;
  real<lower=0> sigma_random;

}

parameters {
  
  // linear predictor
  real zalpha;
  vector[K] zbeta;
  vector[K] zdelta;
  real ztheta;

  // random effects
  vector[nbasin] gamma_basin;
  vector[nblock] gamma_block;
  vector[nsite] gamma_site;

  // random effects variances
  real<lower=0> sigma_basin;
  real<lower=0> sigma_block;
  real<lower=0> sigma_site;

  // residual variance parameter
  real<lower=0> sigma_main;
  
}

transformed parameters {
  real alpha;
  vector[K] beta;
  vector[K] delta;
  real theta;
  vector[N] mu;

  // rescale z-transformed covariates
  alpha = sigma_fixed * zalpha;
  beta = sigma_fixed * zbeta;
  delta = sigma_fixed * zdelta;
  theta = sigma_fixed * ztheta;

  // calculate linear predictor
  mu = alpha +
    X * beta +
    (rep_matrix(arundo, K) .* X) * delta +
    arundo * theta +
    sigma_random * 
    (sigma_basin * gamma_basin[basin] +
     sigma_block * gamma_block[block_id] +
     sigma_site * gamma_site[site]);

}

model {

  // z-scaled priors for regression terms
  zalpha ~ std_normal();
  zbeta ~ std_normal();
  zdelta ~ std_normal();
  ztheta ~ std_normal();

  // priors for random effects
  gamma_basin ~ std_normal();
  gamma_block ~ std_normal();
  gamma_site ~ std_normal();

  // and their variances
  sigma_basin ~ std_normal();
  sigma_block ~ std_normal();
  sigma_site ~ std_normal();
  
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
