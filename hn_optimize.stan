// hn_optimize.stan
// Plain half-normal detection function, distance only, no covariate.
// Fit via rstan::optimizing() (penalized MLE / MAP), not full MCMC sampling.
// Closed-form P_c via erf -- no numerical integration needed.

functions {
  real hn_Pc(real sigma, real w) {
    return (sigma * sqrt(pi()/2) * erf(w / (sigma * sqrt(2)))) / w;
  }
}

data {
  int<lower=1> N;
  vector<lower=0>[N] y;
  real<lower=0> w;
}

parameters {
  real<lower=0> sigma;
}

model {
  // weakly informative prior, same across every sim/stratum -- keeps optimizing()
  // regularized in weak-signal cases rather than running off to a boundary
  sigma ~ normal(0, w/2);

  target += -0.5 * sum(square(y / sigma)) - N * log(hn_Pc(sigma, w));
}

generated quantities {
  real P_c = hn_Pc(sigma, w);
}
