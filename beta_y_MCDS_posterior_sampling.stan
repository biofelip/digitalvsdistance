// This model fits a  beta distribution to the distance data but allows the theta variaiton to enter trouhg a linear predictor that controls the parameter kappa (precision) of the the beta
data {
  real<lower=10, upper=500> w;
  int<lower=1> N;
  vector<lower=0, upper=500>[N] y;
  vector[N] theta;
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
  vector<lower=0>[N] kappa;
  vector<lower=0>[N] alpha_p;
  vector<lower=0>[N] beta_p;

  kappa = exp(gamma_0 + gamma_1 * theta);
  alpha_p = mu * kappa;
  beta_p = (1 - mu) * kappa;
}

model {
  mu ~ beta(2, 2);

  gamma_0 ~ normal(log(10), 1);
  gamma_1 ~ normal(0, 1);

  z ~ beta(alpha_p, beta_p);
}
