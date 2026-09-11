// hn_y_posterior_sampling.stan
// Plain half-normal detection function, distance only, no covariate.

data {
  real<lower=10, upper=500> w;
  int<lower=1> N;
  vector<lower=0>[N] y;

}

parameters {
  real<lower=0> sigma;
}

model {
  // weakly informative prior, same across every sim/stratum -- keeps optimizing()
  // regularized in weak-signal cases rather than running off to a boundary
  sigma ~ normal(0, w/2);

  y ~normal(0, sigma) T[0, w];
}
