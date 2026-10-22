// beta_y_posterior_sampling.stan  (STABILITY-FIXED)
//
// Same beta density on distance as before, but kappa is now reparameterized
// with a floor so alpha_p, beta_p can never drop below min_shape (pass 1
// from R). This is the exact fix worked out in R: kappa = kappa_floor +
// kappa_excess, where kappa_floor = min_shape / min(mu, 1-mu) is the
// smallest kappa that keeps BOTH alpha_p and beta_p >= min_shape.
//
// Without this, alpha_p or beta_p could drop below 1, sending the density's
// true maximum to +Inf near a boundary (or, via the naive mode formula,
// producing a nonsensical out-of-[0,1] mode) -- this was the root cause of
// P_c_beta coming back as Inf.

data {
  real<lower=0.5, upper=500> w;
  int<lower=1> N;
  vector<lower=0>[N] y;
  real<lower=0> min_shape;   // pass 1 from R
}

transformed data {
  vector<lower=0, upper=1>[N] z;
  z = y / w;
}

parameters {
  real<lower=0, upper=1> mu;
  real<lower=0> kappa_excess;   // kappa = kappa_floor + kappa_excess >= kappa_floor
}

transformed parameters {
  real kappa_floor = min_shape / fmin(mu, 1 - mu);
  real kappa = kappa_floor + kappa_excess;
  real alpha_p = mu * kappa;
  real beta_p  = (1 - mu) * kappa;
}

model {
  mu ~ beta(2, 2);
  kappa_excess ~ lognormal(log(10), 1);
  z ~ beta(alpha_p, beta_p);
}
