// hn_y_posterior_sampling.stan
// Plain half-normal detection function, distance only, no covariate.

data {
  int<lower=1> N;
  vector[N] theta;

}

parameters {
  real mu_p;
  real<lower=0> sigma_p;
}

model {
  mu_p ~ normal(0.5, 1);
  sigma_p ~ normal(0, 1);

  theta ~normal(mu_p, sigma_p);
}
