// beta_y_MCDS_posterior_sampling.stan  (STABILITY-FIXED)
//
// Same as before, but kappa (now per-individual, driven by theta) gets the
// same floor treatment. mu is shared across individuals, so kappa_floor is
// a single scalar; kappa itself stays a per-individual vector, floored
// elementwise. This mirrors the R fix directly -- and fixes a real bug
// that existed in the R port of this model, where the floor accidentally
// OVERWROTE the per-individual kappa instead of being added to it,
// silently disabling theta's effect entirely. Here the floor is correctly
// added on top of the real per-individual term.

data {
  real<lower=10, upper=500> w;
  int<lower=1> N;
  vector<lower=0, upper=500>[N] y;
  vector[N] theta;
  real<lower=0> min_shape;   // pass 1 from R
}

transformed data {
  vector<lower=0, upper=1>[N] z;
  z = y / w;
}

parameters {
  real<lower=0, upper=1> mu;
  real gamma_0;
  real gamma_1;
}

transformed parameters {
  real kappa_floor = min_shape / fmin(mu, 1 - mu);
  vector<lower=0>[N] kappa = kappa_floor + exp(gamma_0 + gamma_1 * theta);
  vector<lower=0>[N] alpha_p = mu * kappa;
  vector<lower=0>[N] beta_p  = (1 - mu) * kappa;
}

model {
  mu ~ beta(2, 2);
  gamma_0 ~ normal(log(10), 1);
  gamma_1 ~ normal(0, 1);
  z ~ beta(alpha_p, beta_p);
}
