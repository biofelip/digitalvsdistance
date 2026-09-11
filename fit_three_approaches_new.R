# fit_three_approaches.R
#
# Three approaches to estimating P_c from one sim's detected-animal data,
# all via rstan::optimizing() (penalized MLE, not full MCMC) -- fast enough
# to run across all nsim sims. Reserve full sampling() with divergence/Rhat
# diagnostics for a stratified subsample of hard (low r_y) sims separately.

library(rstan)
library(glue)
library(data.table)


# ---- compile all models ONCE, reused across every sim ----
hn_y_model <- stan_model("hn_y_posterior_sampling.stan")
hn_y_MCDS_model <- stan_model("hn_y_MCDS_half_normal.stan")
beta_y_model <- stan_model("beta_y_posterior_sampling.stan")
beta_y_MCDS_model <- stan_model("beta_y_MCDS_posterior_sampling.stan")
beta_theta_model <- stan_model("beta_theta_posterior_sampling.stan")
normal_theta_model <- stan_model("normal_thetha_posterior_sampling.stan")

# ---- Approach 1: plain distance-only half-normal ----
fit_approach1 <- function(y, w) {
  fit <- tryCatch(
    sampling(hn_y_model, data = list(N = length(y), y = y, w = w)),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(P_c = NA_real_, N_hat = NA_real_, converged = FALSE))
  }

  summary_sigma <- summary(
    fit,
    pars = "sigma",
    probs = c(0.025, 0.975)
  )$summary
  sigma_hat <- summary_sigma[, "mean"]

  # Distance grid
  y_grid <- seq(0, w, length.out = length(y))

  # Fitted observed-distance density:
  # truncated half-normal on [0, w]
  f_obs_hat <- dnorm(y_grid, mean = 0, sd = sigma_hat) /
    (pnorm(w / sigma_hat) - pnorm(0))

  # Value of the fitted observed-distance density at y = 0
  f_obs_hat_0 <- dnorm(0, mean = 0, sd = sigma_hat) /
    (pnorm(w / sigma_hat) - pnorm(0))

  # Rescale so that g_hat(0) = 1
  g_hat <- f_obs_hat / f_obs_hat_0
  # Estimated average detection probability
  Pc_hat <- mean(g_hat)
  list(P_c = P_c, N_hat = length(y) / P_c, converged = TRUE)
}

# ---- Approach 2: half-normal with theta on the scale (MCDS-style) ----
fit_approach2 <- function(y, theta, w) {
  fit <- tryCatch(
    optimizing(hn_cov_model, data = list(N = length(y), y = y, theta = theta, w = w)),
    error = function(e) NULL
  )
  if (is.null(fit) || fit$return_code != 0) {
    return(list(P_c_cov = NA_real_, N_hat_cov = NA_real_, converged = FALSE))
  }
  p_hat_i <- fit$par[grepl("^P_c_i\\[", names(fit$par))]
  list(P_c_cov = mean(p_hat_i), N_hat_cov = sum(1 / p_hat_i), converged = TRUE)
}

# ---- Approach 3: stratify theta into homogeneous cuts, fit approach-1 model per stratum ----
# n_strata: number of cuts. Default uses QUANTILE-based cuts (equal sample size
# per stratum) rather than equal-width bins, since theta can be quite skewed
# (beta with varying alfa/beta) -- equal-width bins risk near-empty, unstable
# strata. Switch to equal-width (seq(min,max,...)) if "homogeneous" should mean
# equal theta-range instead of equal sample count for your purposes.
fit_approach3 <- function(y, theta, w, n_strata = 4) {
  breaks <- quantile(theta, probs = seq(0, 1, length.out = n_strata + 1))
  breaks[1] <- breaks[1] - 1e-8 # ensure the lowest value falls inside the first interval
  strat <- cut(theta, breaks = breaks, include.lowest = TRUE, labels = FALSE)

  strata_results <- lapply(seq_len(n_strata), function(k) {
    idx <- which(strat == k)
    if (length(idx) < 5) {
      return(list(n = length(idx), P_c = NA_real_))
    } # too few to fit
    res <- fit_approach1(y[idx], w)
    list(n = length(idx), P_c = res$P_c)
  })

  n_k <- sapply(strata_results, `[[`, "n")
  P_c_k <- sapply(strata_results, `[[`, "P_c")

  valid <- !is.na(P_c_k) & P_c_k > 0
  if (!any(valid)) {
    return(list(P_c_strat = NA_real_, N_hat_strat = NA_real_, n_strata_failed = n_strata))
  }

  N_hat_strat <- sum(n_k[valid] / P_c_k[valid])
  P_c_strat <- sum(n_k[valid]) / N_hat_strat # pooled P_c implied by the stratified N_hat

  list(
    P_c_strat = P_c_strat, N_hat_strat = N_hat_strat,
    n_strata_failed = sum(!valid)
  )
}

# ---- Approach 4: true shape recovered exactly from the full simulated population ----
# Not a fitted model -- uses the actual known g(y,theta) surface (averaged over
# theta) from detected_df directly, then applies the classical g(0)=1 anchor.
# Isolates the bias from the anchoring CONVENTION itself, with zero estimation
# error, since nothing here is fit from noisy detected-only data.
library(mgcv)

fit_true_shape <- function(full_sim_df, w, n_detected) {
  gam_fit <- tryCatch(
    mgcv::bam(g ~ s(y), data = full_sim_df),
    error = function(e) NULL
  )
  if (is.null(gam_fit)) {
    return(list(P_c_true_shape = NA_real_, N_hat_true_shape = NA_real_))
  }

  y_grid <- seq(min(full_sim_df$y), w, length.out = 200)
  g_y_marginal <- predict(gam_fit, newdata = data.frame(y = y_grid))

  anchor <- max(g_y_marginal, na.rm = TRUE)
  P_c_true_shape <- mean(g_y_marginal / anchor)

  list(
    P_c_true_shape = P_c_true_shape,
    N_hat_true_shape = n_detected / P_c_true_shape
  )
}

# ---- Approach 5: flexible beta-density shape, anchored at its own max ----
# Generalizes "g(0)=1" to "wherever detection peaked, call that 1" -- does
# NOT assume the peak is at y=0, unlike approaches 1-4. Use when the shape
# is known/suspected not to be half-normal-like.
fit_approach5 <- function(y, w, n_grid = 100) {
  fit <- tryCatch(
    optimizing(hn_flex_model, data = list(N = length(y), y = y, w = w, n_grid = n_grid)),
    error = function(e) NULL
  )
  if (is.null(fit) || fit$return_code != 0) {
    return(list(P_c = NA_real_, N_hat = NA_real_, converged = FALSE))
  }
  P_c <- fit$par[["P_c"]]
  list(P_c = P_c, N_hat = length(y) / P_c, converged = TRUE)
}

# ---- combined wrapper: run all five, return one summary row ----
fit_all_approaches <- function(fit_data, full_sim_df, w, true_P_c, n_strata = 4) {
  n_detected <- nrow(fit_data)

  r1 <- fit_approach1(fit_data$distance, w)
  r2 <- fit_approach2(fit_data$distance, fit_data$theta, w)
  r3 <- fit_approach3(fit_data$distance, fit_data$theta, w, n_strata = n_strata)
  r4 <- fit_true_shape(full_sim_df, w, n_detected)
  r5 <- fit_approach5(fit_data$distance, w)

  P_c_uniform <- 1
  N_hat_uniform <- n_detected / P_c_uniform

  data.frame(
    n_detected = n_detected,
    true_P_c = true_P_c, true_N_hat = n_detected / true_P_c,
    P_c_plain = r1$P_c, N_hat_plain = r1$N_hat,
    P_c_cov = r2$P_c_cov, N_hat_cov = r2$N_hat_cov,
    P_c_strat = r3$P_c_strat, N_hat_strat = r3$N_hat_strat,
    P_c_uniform = P_c_uniform, N_hat_uniform = N_hat_uniform,
    P_c_true_shape = r4$P_c_true_shape, N_hat_true_shape = r4$N_hat_true_shape,
    P_c_flexible = r5$P_c, N_hat_flexible = r5$N_hat,
    bias_plain = r1$P_c - true_P_c,
    bias_cov = r2$P_c_cov - true_P_c,
    bias_strat = r3$P_c_strat - true_P_c,
    bias_uniform = P_c_uniform - true_P_c,
    bias_flexible = r5$P_c - true_P_c,
    converged_plain = r1$converged, converged_cov = r2$converged,
    converged_flex = r5$converged,
    n_strata_failed = r3$n_strata_failed
  )
}

# ---- example: fold into your per-sim loop ----

# ---- example: fold into your per-sim loop ----
nsim <- 10000

list_detected_df <- c(
  "E:/digitalvsdistance_results/new_results/1_1_detected_df_list_high_r_y_hn_g_y_smooth_scale.parquet",
  "E:/digitalvsdistance_results/new_results/1_detected_df_list_low_r_y_hn_g_y_smooth_scale.parquet",
  "E:/digitalvsdistance_results/new_results/2_2_detected_df_list_high_r_y_beta_g_y_smooth_scale.parquet",
  "E:/digitalvsdistance_results/new_results/2_detected_df_list_low_r_y_beta_g_y_smooth_scale.parquet",
  "E:/digitalvsdistance_results/new_results/3_detected_df_list_smooth_params_uniform_g_y.parquet"
)

list_sim_summary <- c(
  "E:/digitalvsdistance_results/new_results/1_1_sim_summary_list_high_r_y_hn_g_y_smooth_scale.rds",
  "E:/digitalvsdistance_results/new_results/1_sim_summary_list_low_r_y_hn_g_y_smooth_scale.rds",
  "E:/digitalvsdistance_results/new_results/2_2_sim_summary_list_high_r_y_beta_g_y_smooth_scale.rds",
  "E:/digitalvsdistance_results/new_results/2_sim_summary_list_low_r_y_beta_g_y_smooth_scale.rds",
  "E:/digitalvsdistance_results/new_results/3_sim_summary_list_smooth_params_uniform_g_y.rds"
)
# i <- sample(1:5, 1)
# nsim <- 10


library(profvis)
# f <- function(nsim) {
for (i in 1:5) {
  detected_df <- arrow::read_parquet(list_detected_df[i])
  sim_summary <- readRDS(list_sim_summary[i])
  sim_summary <- data.table::rbindlist(sim_summary)
  detected_df <- data.table::data.table(detected_df)
  data.table::setkey(detected_df, sim)
  data.table::setkey(sim_summary, sim)
  results_list <- vector("list", nsim)
  for (spick in 1:nsim) {
    sim_data <- detected_df[sim == spick]
    fit_data <- sim_data[detected == 1, .(distance = y, theta = theta)]
    real_Pc <- sim_data[, mean(g)]
    w_this_sim <- sim_summary[sim == spick]$w
    result_row <- fit_all_approaches(fit_data, sim_data, w = w_this_sim, n_strata = 4, true_P_c = real_Pc)
    results_list[[spick]] <- cbind(sim = spick, result_row)
  }
  results_pc <- data.table::rbindlist(results_list)
  saveRDS(results_pc, file = glue("E:/digitalvsdistance_results/new_results/results_pc_{tools::file_path_sans_ext(basename(list_sim_summary[i]))}.rds"))
  gc()
}
# }

p <- profvis(f(10))
# plot the difference between the observed and the real Pc  agains the proportion of
# r_Y that they explain
# fist join with the sim_summary
results_pc <- results_pc[sim_summary, on = "sim"]

results_pc |>
  pivot_longer(c("P_c_plain", "P_c_cov", "P_c_strat"), names_to = "approach", values_to = "Pc") |>
  ggplot(aes(x = realized_r_y, y = abs(Pc - P_c_real), color = approach)) +
  geom_smooth()


fit_data <- detected_df[sim == spick & detected == 1, .(distance = y, theta = theta)]
w_this_sim <- sim_summary[sim == spick]$w
result_row <- fit_all_approaches(fit_data, w = w_this_sim, n_strata = 4, real_Pc = real_Pc)
