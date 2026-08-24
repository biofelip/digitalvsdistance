// hn_flexible_optimize.stan
//
// Flexible detection-function shape for distance, using a beta DENSITY
// (not CDF) as the raw shape. Unlike the half-normal, this can be
// monotonic decreasing, monotonic increasing, or peak anywhere in the
// interior of [0,w] -- matching what actually broke the half-normal model.
//
// Anchored by the GENERALIZED assumption "wherever detection was highest,
// call that 1" -- not "detection is highest at y=0". The anchor's value M
// cancels out of the conditional-on-detection likelihood exactly (it's a
// constant multiplying both g(y) and P_c identically), so it only re-enters
// when reporting P_c itself, in generated quantities.

functions {
  // unnormalized beta-density shape on [0,w], via u = y/w
  real beta_shape(real y, real w, real alfa, real beta_) {
    real u = y / w;
    return exp((alfa - 1) * log(u) + (beta_ - 1) * log1m(u) - lbeta(alfa, beta_));
  }
}

data {
  int<lower=1> N;
  vector<lower=0>[N] y;
  real<lower=0> w;
  int<lower=2> n_grid;   // resolution for numerically finding the max and P_c
}

parameters {
  // constrained >=1 so the shape stays well-behaved (finite density
  // everywhere, no blow-up at the boundaries)
  real<lower=1> alfa;
  real<lower=1> beta_;
}

transformed parameters {
  vector[N] g_raw;
  for (i in 1:N) g_raw[i] = beta_shape(y[i], w, alfa, beta_);

  // (1/w) * integral_0^w g_raw(y) dy, via a Riemann sum -- this is what
  // survives in the likelihood after the anchor M cancels
  real Pc_numerator;
  {
    real total = 0;
    for (k in 1:n_grid) {
      real y_k = w * (k - 0.5) / n_grid;
      total += beta_shape(y_k, w, alfa, beta_);
    }
    Pc_numerator = total / n_grid;
  }
}

model {
  alfa  ~ normal(1, 3);
  beta_ ~ normal(1, 3);

  // f(y|detected) = (1/w) * g_raw(y) / Pc_numerator -- anchor M cancelled
  // out algebraically; -N*log(w) is a constant, dropped
  target += sum(log(g_raw)) - N * log(Pc_numerator);
}

generated quantities {
  // find M = max shape value over [0,w] -- this is where the anchor
  // re-enters, needed for P_c but not for the fit itself
  real M;
  {
    real best = negative_infinity();
    for (k in 1:n_grid) {
      real y_k = w * (k - 0.5) / n_grid;
      real val = beta_shape(y_k, w, alfa, beta_);
      if (val > best) best = val;
    }
    M = best;
  }
  real P_c = Pc_numerator / M;
}
