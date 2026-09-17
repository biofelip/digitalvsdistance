// beta_theta_posterior_sampling.stan  (STABILITY-FIXED)
//
// Two fixes here, both needed:
//
// 1. Same kappa floor as the other beta models (alpha_p, beta_p >= min_shape).
//
// 2. z = (theta-min)/(max-min) ALWAYS produces an exact 0 and an exact 1
//    (whichever animals have the min/max theta). Once alpha_p/beta_p are
//    floored >= 1, the beta density is EXACTLY zero at those boundary
//    points -- so those two observations contribute log(0) = -Inf to the
//    likelihood for nearly every candidate the optimizer tries. This is
//    what made this specific model return NULL (BFGS/L-BFGS can't take a
//    sensible step when the objective is -Inf almost everywhere). Fixed by
//    squeezing z slightly inside (0,1) before ever touching it.

data {
  int<lower=1> N;
  vector[N] theta;
  real<lower=0> min_shape;            // pass 1 from R
  real<lower=0, upper=0.01> eps;      // pass e.g. 1e-6 from R
}

transformed data {
  real theta_min = min(theta);
  real theta_max = max(theta);
  vector[N] z_raw = (theta - theta_min) / (theta_max - theta_min);
  vector<lower=0, upper=1>[N] z = z_raw * (1 - 2 * eps) + eps;  // never exactly 0 or 1
}

parameters {
  real<lower=0, upper=1> mu;
  real<lower=0> kappa_excess;
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
