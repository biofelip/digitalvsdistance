
data {
  real<lower=10, upper=500> w;
  int<lower=1> N;
  vector<lower=0>[N] y;

}

transformed data {
   vector<lower=0, upper=1>[N] z;

   z = y/w;
}

parameters {
  real<lower=0, upper=1> mu;
  real<lower=0> kappa;
}

transformed parameters {
   real alpha_p = mu * kappa;
   real beta_p = (1-mu)* kappa;
}
model {
  mu ~ beta(2,2);
  kappa ~ lognormal(log(10),1);

  z ~beta(alpha_p, beta_p);
}