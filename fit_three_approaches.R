# fit_three_approaches.R
#
# Three approaches to estimating P_c from one sim's detected-animal data,
# all via rstan::optimizing() (penalized MLE, not full MCMC) -- fast enough
# to run across all nsim sims. Reserve full sampling() with divergence/Rhat
# diagnostics for a stratified subsample of hard (low r_y) sims separately.

library(rstan)
rstan_options(auto_write = TRUE)

# ---- compile both models ONCE, reused across every sim ----
hn_model     <- stan_model("hn_optimize.stan")
hn_cov_model <- stan_model("hn_cov_optimize.stan")

# ---- Approach 1: plain distance-only half-normal ----
fit_approach1 <- function(y, w) {
  fit <- tryCatch(
    optimizing(hn_model, data = list(N = length(y), y = y, w = w)),
    error = function(e) NULL
  )
  if (is.null(fit) || fit$return_code != 0) {
    return(list(P_c = NA_real_, N_hat = NA_real_, converged = FALSE))
  }
  P_c <- fit$par[["P_c"]]
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
  breaks[1] <- breaks[1] - 1e-8   # ensure the lowest value falls inside the first interval
  strat <- cut(theta, breaks = breaks, include.lowest = TRUE, labels = FALSE)

  strata_results <- lapply(seq_len(n_strata), function(k) {
    idx <- which(strat == k)
    if (length(idx) < 5) return(list(n = length(idx), P_c = NA_real_))  # too few to fit
    res <- fit_approach1(y[idx], w)
    list(n = length(idx), P_c = res$P_c)
  })

  n_k   <- sapply(strata_results, `[[`, "n")
  P_c_k <- sapply(strata_results, `[[`, "P_c")

  valid <- !is.na(P_c_k) & P_c_k > 0
  if (!any(valid)) return(list(P_c_strat = NA_real_, N_hat_strat = NA_real_, n_strata_failed = n_strata))

  N_hat_strat <- sum(n_k[valid] / P_c_k[valid])
  P_c_strat   <- sum(n_k[valid]) / N_hat_strat   # pooled P_c implied by the stratified N_hat

  list(P_c_strat = P_c_strat, N_hat_strat = N_hat_strat,
       n_strata_failed = sum(!valid))
}

# ---- combined wrapper: run all three, return one summary row ----
fit_all_approaches <- function(fit_data, w, n_strata = 4, real_Pc = NULL) {
  n_detected <- nrow(fit_data)

  r1 <- fit_approach1(fit_data$distance, w)
  r2 <- fit_approach2(fit_data$distance, fit_data$theta, w)
  r3 <- fit_approach3(fit_data$distance, fit_data$theta, w, n_strata = n_strata)

  # Approach 4: uniform / strip-transect assumption -- P_c=1 by fiat, not fit
  P_c_uniform   <- 1
  N_hat_uniform <- n_detected / P_c_uniform   # = n_detected exactly

  data.frame(
    n_detected      = n_detected,
    P_c_real        = real_Pc,
    P_c_plain       = r1$P_c,        N_hat_plain     = r1$N_hat,
    P_c_cov         = r2$P_c_cov,    N_hat_cov       = r2$N_hat_cov,
    P_c_strat       = r3$P_c_strat,  N_hat_strat     = r3$N_hat_strat,
    converged_plain = r1$converged,  converged_cov   = r2$converged,
    n_strata_failed = r3$n_strata_failed
  )
}

# ---- example: fold into your per-sim loop ----
spick=1
nsim=1000
detected_df <- arrow::read_parquet(r"(L:\Documents\digitalvsdistance_results\1_detected_df_list_hn_g_y.parquet)")
sim_summary <- readRDS(r"(L:\Documents\digitalvsdistance_results\1_sim_summary_list_hn_g_y.rds)")
sim_summary <- data.table::rbindlist(sim_summary)
deteced_df <- data.table::data.table(detected_df)
results_pc <- data.table::data.table()
for (spick in 1:nsim) {
  fit_data <- detected_df[sim == spick & detected == 1, .(distance = y, theta = theta)]
  real_Pc = detected_df[sim == spick, mean(g)]
  w_this_sim <- sim_summary[sim == spick]$w
  result_row <- fit_all_approaches(fit_data, w = w_this_sim, n_strata = 4, real_Pc = real_Pc)
  results_pc <- rbind(results_pc, cbind(sim = spick, result_row))
}
fit_data <- detected_df[sim == spick & detected == 1, .(distance = y, theta = theta)]
w_this_sim <- sim_summary[sim == spick]$w
result_row <- fit_all_approaches(fit_data, w = w_this_sim, n_strata = 4, real_Pc = real_Pc)

