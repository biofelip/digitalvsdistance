# fit_full_model.R
#
# Compiles full_detection_model.stan ONCE, then fits it against each
# simulated dataset's detected-animal data. Returns a one-row summary per
# sim, matching the naming pattern of your existing results_pc columns
# (P_c, P_c_cov, P_c_unif / N_hat, N_hat_cov, N_hat_unif).
#
# Adds divergences and rhat_max -- these are the Bayesian analogue of
# mrds's convergence failures, but instead of crashing (like ddf() does in
# the hard regime), they degrade gracefully into a quantifiable signal of
# "this simulation's variance-ratio combination is hard to estimate."

library(rstan)
rstan_options(auto_write = TRUE)

# If your machine restricts spawning multiple R processes (same restriction
# that blocked plan(multisession, ...) earlier), parallel chains via
# mc.cores may also be blocked. Set mc.cores = 1 as the safe default;
# raise it if/when parallel execution is actually available to you.
options(mc.cores = 1)

# ---- Compile once, reuse across every simulated dataset ----
full_model <- stan_model("full_detection_model.stan")

# ---- Fit one simulation's detected-animal data ----
# fit_data: data.frame with columns `distance` and `theta` (detected animals only)
# w: strip half-width for that sim
fit_full_model <- function(fit_data, w, iter = 2000, chains = 4, warmup = 1000) {

  stan_data <- list(
    N     = nrow(fit_data),
    y     = fit_data$distance,
    theta = fit_data$theta,
    w     = w
  )

  fit <- tryCatch(
    sampling(full_model, data = stan_data,
             iter = iter, warmup = warmup, chains = chains,
             refresh = 0, show_messages = FALSE, open_progress = FALSE),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(data.frame(P_c_full = NA_real_, N_hat_full = NA_real_,
                       divergences = NA_integer_, rhat_max = NA_real_))
  }

  draws <- tryCatch(rstan::extract(fit), error = function(e) NULL)
  if (is.null(draws)) {
    return(data.frame(P_c_full = NA_real_, N_hat_full = NA_real_,
                       divergences = NA_integer_, rhat_max = NA_real_))
  }

  # posterior mean of each detected animal's p_hat, then summarize the same
  # way you already do for the mrds covariate model (mean / Horvitz-Thompson sum)
  p_hat_mean <- colMeans(draws$p_hat)

  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  n_divergent <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
  rhat_max <- suppressWarnings(max(summary(fit)$summary[, "Rhat"], na.rm = TRUE))

  data.frame(
    P_c_full    = mean(p_hat_mean),
    N_hat_full  = sum(1 / p_hat_mean),
    divergences = n_divergent,   # >0 flags a hard-to-estimate regime
    rhat_max    = rhat_max       # >1.01 flags poor mixing / weak identifiability
  )
}

# ---- Example: fold into your existing per-sim loop ----
# fit_data <- detected_df %>%
#   filter(detected == 1, sim == s_pick) %>%
#   transmute(distance = y, theta = theta)
#
# full_result <- fit_full_model(fit_data, w = w)
#
# results_pc <- rbind(results_pc, cbind(sim = s_pick, full_result))
