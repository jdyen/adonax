data {
  
  // main indices
  int<lower=0> N;
  int<lower=0> K;
  
  // response data  
  int<lower=0,upper=1> y[N];

  // predictor variables  
  matrix[N, K] X;
  
  // random effects
  int<lower=0> nbasin;
  int<lower=1,upper=nbasin> basin[N];

  // scale of main effects and overall variance
  real<lower=0> sigma_fixed;
  real<lower=0> sigma_random;

}

parameters {
  
  // linear predictor
  real zalpha;
  vector[K] zbeta;

  // random effects
  vector[nbasin] gamma_basin;

  // random effects variances
  real<lower=0> sigma_basin;

}

transformed parameters {
  real alpha;
  vector[K] beta;  
  vector[N] mu;

  // rescale z-transformed covariates
  alpha = sigma_fixed * zalpha;
  beta = sigma_fixed * zbeta;

  // calculate linear predictor
  mu = rep_vector(alpha, N) +
    X * beta +
    sigma_random * sigma_basin * gamma_basin[basin];

}

model {

  // z-scaled priors for regression terms
  zalpha ~ std_normal();
  to_vector(zbeta) ~ std_normal();

  // priors for random effects
  to_vector(gamma_basin) ~ std_normal();

  // and their variances
  sigma_basin ~ std_normal();

  // calculate likelihood of response (zero-inflataed Poisson)
  target += bernoulli_logit_lpmf(y | mu);
   
}

generated quantities {
  int<lower=0,upper=1> ypred[N];
  
  // random draws from the posterior of mu and phi for posterior checks
  for (i in 1:N) {
    ypred[i] = bernoulli_logit_rng(mu[i]);
  }
  
}
