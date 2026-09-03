// full_detection_model.stan
//
// Fits the FULL joint detection function g(y, theta) -- distance AND covariate
// together -- to distance-sampling data, conditional on detection.
//
// Design choice: uses an UNCENTERED basis (b0, by, bth, bint instead of
// g0, a_y, a_th, a_int). This mirrors the standard MCDS "conditioning" trick
// (Marques & Buckland 2004): theta enters the likelihood only through each
// individual's OWN observed value, never through the population distribution
// f(theta). That means we never need to integrate over f(theta) -- which is
// good, because doing so would require mean_g_theta = E_theta[pgamma(theta,k,1)],
// and that mean depends on k (a parameter being estimated), so it CANNOT be
// precomputed as fixed data. The uncentered basis sidesteps this entirely.
//
// Key derivation (see comments in the model block): the shoulder-condition
// rescale -- which enforces g(0,theta)=1 for every theta -- cancels exactly
// out of the conditional-on-detection likelihood. It re-enters only in
// generated quantities, where it's needed to get p_hat (individual detection
// probability) for Horvitz-Thompson abundance estimation.

functions {
  // Half-normal shape, evaluated at a vector of distances y (unnormalized
  // over [0,w] -- normalization is handled later via the shoulder rescale)
  vector hn_density(vector y, real sigma) {
    return sqrt(2 / (pi() * square(sigma))) * exp(-square(y) / (2 * square(sigma)));
  }

  // Same half-normal shape, evaluated at the single point y = 0 (the "shoulder")
  real hn_density0(real sigma) {
    return sqrt(2 / (pi() * square(sigma)));
  }

  // E_y[hn_density(y,sigma)] under y ~ Uniform(0,w).
  // Closed form via erf -- this is the ONLY integral this model ever needs,
  // because f(y) is uniform and known. There is no equivalent integral over
  // f(theta) anywhere in this model.
  real hn_mean(real sigma, real w) {
    return erf(w / (sigma * sqrt(2))) / w;
  }
}

data {
  int<lower=1> N;                  // number of DETECTED animals in this sim
  vector<lower=0>[N] y;            // detected distances
  vector<lower=0>[N] theta;        // covariate value for each detected animal
  real<lower=0> w;                 // strip half-width (truncation distance)
}

parameters {
  real b0;                         // intercept, uncentered basis (not literally P_c)
  real by;                         // coefficient on the half-normal shape in y
  real bth;                        // coefficient on the gamma-CDF shape in theta
  real bint;                       // coefficient on the y*theta interaction
  real<lower=0> sigma;             // half-normal scale parameter
  real<lower=0> k;                 // gamma shape controlling the theta functional form
}

transformed parameters {
  real m_y  = hn_mean(sigma, w);       // E_y[h_y(y)] -- closed form, no theta involved
  real h_y0 = hn_density0(sigma);      // h_y(0), i.e. the shoulder value of the y-shape

  vector[N] h_y = hn_density(y, sigma);

  vector[N] h_theta;
  for (i in 1:N) {
    h_theta[i] = gamma_cdf(theta[i], k, 1);   // conditions on each animal's OWN theta_i
  }

  // Unrescaled joint detection surface
  vector[N] g_raw = b0 + by * h_y + bth * h_theta + bint * (h_y .* h_theta);

  // g_raw evaluated AT y=0, per animal -- this is the shoulder-condition denominator
  vector[N] shoulder = b0 + by * h_y0 + bth * h_theta + bint * h_y0 * h_theta;

  // What survives integrating g_raw over y (before dividing by shoulder).
  // Only m_y appears -- no integral over f(theta) anywhere, by construction.
  vector[N] Pc_numerator = b0 + by * m_y + h_theta * (bth + bint * m_y);
}

model {
  // Same weakly-informative priors for every simulated dataset, as agreed
  b0   ~ normal(0.5, 0.5);
  by   ~ normal(0, 1);
  bth  ~ normal(0, 1);
  bint ~ normal(0, 1);
  sigma ~ normal(0, w / 2);
  k     ~ normal(0, 5);

  // Conditional-on-detection density: f(y|theta,detected) = f_y(y)*g(y,theta)/P_c(theta)
  // where g(y,theta) = g_raw/shoulder and P_c(theta) = Pc_numerator/shoulder.
  // The shoulder(theta) term is common to both and cancels exactly:
  //   f(y|theta,detected) = g_raw(y,theta) / Pc_numerator(theta) / w
  // The constant -N*log(w) is dropped since it doesn't affect the posterior.
  target += sum(log(g_raw) - log(Pc_numerator));
}

generated quantities {
  // Per-individual detection probability. This is where the shoulder rescale
  // re-enters -- needed for abundance estimation even though it cancelled
  // out of the fitting likelihood above.
  vector[N] p_hat = Pc_numerator ./ shoulder;
}
