// hn_cov_optimize.stan
// Half-normal with theta entering through the scale parameter (log link),
// matching the standard MCDS structure. Uses the "conditioning" approach --
// each individual's own observed theta_i is used directly; f(theta) is never
// needed, matching how mrds/MCDS handles covariates.
// Fit via rstan::optimizing().

functions {
  real hn_Pc(real sigma, real w) {
    return (sigma * sqrt(pi()/2) * erf(w / (sigma * sqrt(2)))) / w;
  }
}

data {
  int<lower=1> N;
  vector<lower=0>[N] y;
  vector[N] theta;
  real<lower=0> w;
}

parameters {
  real beta0;
  real beta1;
}

transformed parameters {
  vector<lower=0>[N] sigma_i = exp(beta0 + beta1 * theta);
}

model {
  beta0 ~ normal(0, 5);
  beta1 ~ normal(0, 5);

  for (i in 1:N)
    target += -0.5 * square(y[i] / sigma_i[i]) - log(hn_Pc(sigma_i[i], w));
}

generated quantities {
  // per-individual detection probability -- averaged/summed in R afterward
  // for the sample-conditional P_c and Horvitz-Thompson N_hat
  vector[N] P_c_i;
  for (i in 1:N) P_c_i[i] = hn_Pc(sigma_i[i], w);
}
