// hn_y_posterior_sampling.stan
// Plain half-normal detection function, distance only, no covariate.

data {
  int<lower=1> N;
  vector[N] theta;

}

parameters {
  real mu_theta;
  real<lower=0> sigma;
}

model {
  mu_theta ~ normal(0.5, 1);
  sigma_theta ~ normal(0, 1);

  y ~normal(mu_theta, sigma_theta) T[0, w];
}
