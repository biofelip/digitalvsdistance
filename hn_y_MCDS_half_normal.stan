// hn_y_MCDS_half_normal.stan
//  half-normal with scale parameter contolled by a linear predictor of 
// theta 

data {
  real<lower=0.5, upper=500> w;
  int<lower=1> N;
  vector<lower=0>[N] y;
  vector[N] theta;

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
  beta_0 ~ normal(log(w/2),1);
  beta_1 ~ normal(0,1);

  y ~normal(0, sigma) T[0, w];
}

// generated quantities {
//   array[N] real y_rep;
//   real lower_cdf;
//   real upper_cdf;
//   real u;

//   lower_cdf = normal_cdf(0 | 0, sigma);
//   upper_cdf = normal_cdf(w | 0, sigma);

//   for (n in 1:N) {
//     u = uniform_rng(lower_cdf, upper_cdf);
//     y_rep[n] = sigma * inv_Phi(u);
//   }
// }

