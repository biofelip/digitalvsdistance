
data {
  real<lower=10, upper=500> w;
  int<lower=1> N;
  vector<lower=0>[N] y;
  vector[N] theta;

}

transformed data {
  real<lower=0> var_y = variance(y);
  real nm1_over2 = 0.5 * (N - 1);
  real sqrt_N = sqrt(N);
}
parameters {
  real beta_0;
  real beta_1;  
}

transformed parameters {
    vector<lower=0>[N] sigma;
    sigma = exp(beta_0 + beta_1*theta);
}
model {
  var_y ~ gamma(nm1_over2, nm1_over2 / sigma^2);
}