//code for a model that fits a beta distribution to the theta component and ignores the distance

data{
    int<lower=1> N;
    vector[N] theta;
}

transformed data {
  vector<lower=0, upper=1>[N] z;
  real theta_min = min(theta);
  real theta_max = max(theta);
  z = (theta - theta_min) / (theta_max - theta_min);
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