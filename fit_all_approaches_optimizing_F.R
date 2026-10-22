# fit_all_approaches_optimizing.R
# Fits every approach to each simulated dataset with rstan::optimizing()
# and saves P_c and bias per simulation.

library(rstan)
library(glue)
library(data.table)

hn_y_model <- stan_model("hn_y_posterior_sampling.stan") # approach 1 (and 3)
hn_y_MCDS_model <- stan_model("hn_y_MCDS_half_normal.stan") # approach 2
beta_y_model <- stan_model("beta_y_posterior_sampling_C.stan") # approach 4
beta_y_MCDS_model <- stan_model("beta_y_MCDS_posterior_sampling_C.stan") # approach 5
beta_theta_model <- stan_model("beta_theta_posterior_sampling_C.stan") # approach 6
normal_theta_model <- stan_model("normal_thetha_posterior_sampling.stan") # approach 7

MIN_SHAPE <- 1
EPS <- 1e-6
UNIFORM_TOL <- 1e-3

# ---- Approach 1: half-normal on distance ----
fit_approach1 <- function(y, w) {
  fit <- tryCatch(
    optimizing(hn_y_model, data = list(N = length(y), y = y, w = w), init = list(sigma = w / 4)),
    error = function(e) NULL
  )
  if (is.null(fit) || fit$return_code != 0) {
    return(list(P_c = NA_real_, N_hat = NA_real_, sigma_hat = NA_real_))
  }
  sigma_hat <- fit$par[["sigma"]]
  p_hat_i <- exp(-y^2 / (2 * sigma_hat^2))
  N_hat <- sum(1 / p_hat_i)
  P_c <- length(y) / N_hat
  list(P_c = P_c, N_hat = N_hat, sigma_hat = sigma_hat)
}

# ---- Approach 2: half-normal, theta on the scale ----
fit_approach2 <- function(y, theta, w) {
  fit <- tryCatch(
    optimizing(hn_y_MCDS_model, data = list(N = length(y), y = y, theta = theta, w = w), init = list(beta_0 = log(w / 4), beta_1 = 0)),
    error = function(e) NULL
  )
  if (is.null(fit) || fit$return_code != 0) {
    return(list(Pc_hat = NA_real_, N_hat_cov = NA_real_, beta_0 = NA_real_, beta_1 = NA_real_))
  }
  beta_0_hat <- fit$par[["beta_0"]]
  beta_1_hat <- fit$par[["beta_1"]]
  sigma_hat <- exp(beta_0_hat + beta_1_hat * theta)
  p_hat_i <- exp(-y^2 / (2 * sigma_hat^2))
  N_hat <- sum(1 / p_hat_i)
  Pc_hat <- length(y) / N_hat
  list(Pc_hat = Pc_hat, N_hat_cov = N_hat, beta_0 = beta_0_hat, beta_1 = beta_1_hat)
}

# ---- Approach 3: half-normal stratified on theta ----
fit_approach3 <- function(y, theta, w, n_strata = 4) {
  breaks <- quantile(theta, probs = seq(0, 1, length.out = n_strata + 1))
  breaks[1] <- breaks[1] - 1e-8
  strat <- cut(theta, breaks = breaks, include.lowest = TRUE, labels = FALSE)

  strata_results <- lapply(seq_len(n_strata), function(k) {
    idx <- which(strat == k)
    if (length(idx) < 5) {
      return(list(n = length(idx), P_c = NA_real_))
    }
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
  P_c_strat <- sum(n_k[valid]) / N_hat_strat
  list(P_c_strat = P_c_strat, N_hat_strat = N_hat_strat, n_strata_failed = sum(!valid))
}

# ---- Approach 4: Beta on distance ----
fit_approach4 <- function(y, w, min_shape = MIN_SHAPE, uniform_tol = UNIFORM_TOL) {
  fit <- tryCatch(
    optimizing(beta_y_model, data = list(N = length(y), y = y, w = w, min_shape = min_shape), init = list(mu = 0.5, kappa_excess = 10)),
    error = function(e) NULL
  )
  if (is.null(fit) || fit$return_code != 0) {
    return(list(Pc_hat = NA_real_, N_hat = NA_real_, alpha = NA_real_, beta = NA_real_))
  }
  alpha_hat <- fit$par[["alpha_p"]]
  beta_hat <- fit$par[["beta_p"]]

  if (alpha_hat + beta_hat <= 2 + uniform_tol) {
    return(list(Pc_hat = 1, N_hat = length(y), alpha = alpha_hat, beta = beta_hat))
  }

  z <- y / w
  mode <- (alpha_hat - 1) / (alpha_hat + beta_hat - 2)
  anchor <- dbeta(mode, alpha_hat, beta_hat)
  p_hat_i <- dbeta(z, alpha_hat, beta_hat) / anchor
  N_hat <- sum(1 / p_hat_i)
  list(Pc_hat = length(y) / N_hat, N_hat = N_hat, alpha = alpha_hat, beta = beta_hat)
}


# ---- Approach 5: Beta on distance, theta on the concentration ----
fit_approach5 <- function(y, w, theta, min_shape = MIN_SHAPE, uniform_tol = UNIFORM_TOL) {
  fit <- tryCatch(
    optimizing(beta_y_MCDS_model, data = list(N = length(y), y = y, w = w, theta = theta, min_shape = min_shape), init = list(mu = 0.5, kappa_excess = 10)),
    error = function(e) NULL
  )
  if (is.null(fit) || fit$return_code != 0) {
    return(list(Pc_hat = NA_real_, N_hat = NA_real_, mu = NA_real_, gamma_0 = NA_real_, gamma_1 = NA_real_))
  }
  mu_hat <- fit$par[["mu"]]
  gamma_0_hat <- fit$par[["gamma_0"]]
  gamma_1_hat <- fit$par[["gamma_1"]]

  alpha_hat <- fit$par[grepl("^alpha_p\\[", names(fit$par))]
  beta_hat <- fit$par[grepl("^beta_p\\[", names(fit$par))]

  z <- y / w
  near_uniform <- (alpha_hat + beta_hat) <= (2 + uniform_tol)
  mode_i <- ifelse(near_uniform, NA, (alpha_hat - 1) / (alpha_hat + beta_hat - 2))
  anchor_i <- ifelse(near_uniform, 1, dbeta(mode_i, alpha_hat, beta_hat))
  p_hat_i <- ifelse(near_uniform, 1, dbeta(z, alpha_hat, beta_hat) / anchor_i)

  N_hat <- sum(1 / p_hat_i)
  list(Pc_hat = length(y) / N_hat, N_hat = N_hat, mu = mu_hat, gamma_0 = gamma_0_hat, gamma_1 = gamma_1_hat)
}

# ---- Approach 6: Beta on theta ----
fit_approach6 <- function(theta, min_shape = MIN_SHAPE, eps = EPS, uniform_tol = UNIFORM_TOL) {
  fit <- tryCatch(
    optimizing(beta_theta_model,
      data = list(
        N = length(theta), theta = theta,
        min_shape = min_shape, eps = eps
      ),
      init = list(mu = 0.5, kappa_excess = 10)
    ),
    error = function(e) NULL
  )
  if (is.null(fit) || fit$return_code != 0) {
    return(list(Pc_hat = NA_real_, N_hat = NA_real_, alpha = NA_real_, beta = NA_real_))
  }
  alpha_hat <- fit$par[["alpha_p"]]
  beta_hat <- fit$par[["beta_p"]]

  if (alpha_hat + beta_hat <= 2 + uniform_tol) {
    return(list(Pc_hat = 1, N_hat = length(theta), alpha = alpha_hat, beta = beta_hat))
  }

  z_raw <- (theta - min(theta)) / (max(theta) - min(theta))
  z <- z_raw * (1 - 2 * eps) + eps
  mode <- (alpha_hat - 1) / (alpha_hat + beta_hat - 2)
  anchor <- dbeta(mode, alpha_hat, beta_hat)
  p_hat_i <- dbeta(z, alpha_hat, beta_hat) / anchor
  N_hat <- sum(1 / p_hat_i)
  list(Pc_hat = length(theta) / N_hat, N_hat = N_hat, alpha = alpha_hat, beta = beta_hat)
}

# ---- Approach 7: Normal on theta ----
fit_approach7 <- function(theta) {
  fit <- tryCatch(
    optimizing(normal_theta_model,
      data = list(N = length(theta), theta = theta),
      init = list(mu_p = mean(theta), sigma_p = sd(theta))
    ),
    error = function(e) NULL
  )
  if (is.null(fit) || fit$return_code != 0) {
    return(list(Pc_hat = NA_real_, N_hat = NA_real_, mu = NA_real_, sigma = NA_real_))
  }
  mu_hat <- fit$par[["mu_p"]]
  sigma_hat <- fit$par[["sigma_p"]]
  anchor <- dnorm(mu_hat, mu_hat, sigma_hat)
  p_hat_i <- dnorm(theta, mu_hat, sigma_hat) / anchor
  N_hat <- sum(1 / p_hat_i)
  Pc_hat <- length(theta) / N_hat
  list(Pc_hat = Pc_hat, N_hat = N_hat, mu = mu_hat, sigma = sigma_hat)
}

# ---- Approach 8: true shape from the full population (binned mean of g over y) ----
fit_true_shape <- function(full_sim_df, n_detected, n_bins = 100) {
  bins <- cut(full_sim_df$y, breaks = n_bins)
  g_y_marginal <- as.vector(tapply(full_sim_df$g, bins, mean))
  g0 <- g_y_marginal[1]
  gmax <- max(g_y_marginal, na.rm = TRUE)
  Pc_anchor0 <- mean(g_y_marginal / g0)
  Pc_anchormax <- mean(g_y_marginal / gmax)
  list(
    P_c_true_shape_anchor0 = Pc_anchor0, P_c_true_shape_anchormax = Pc_anchormax,
    N_hat_true_shape_anchor0 = n_detected / Pc_anchor0,
    N_hat_true_shape_anchormax = n_detected / Pc_anchormax
  )
}

param_string <- function(lst) {
  if (is.null(lst)) {
    return(NA_character_)
  }
  drop <- c("Pc_hat", "P_c", "N_hat", "N_hat_cov", "N_hat_strat", "P_c_strat", "n_strata_failed")
  keep <- setdiff(names(lst), drop)
  if (length(keep) == 0) {
    return(NA_character_)
  }
  vals <- vapply(keep, function(nm) {
    v <- lst[[nm]]
    if (is.null(v) || length(v) == 0 || all(is.na(v))) {
      return(NA_real_)
    }
    if (length(v) > 1) v <- mean(v, na.rm = TRUE)
    round(v, 4)
  }, numeric(1))
  paste(paste0(keep, "=", vals), collapse = "; ")
}

safe_fit <- function(expr) tryCatch(expr, error = function(e) NULL)
g <- function(lst, nm) if (is.null(lst)) NA else lst[[nm]]

fit_all_approaches <- function(fit_data, full_sim_df, w, true_P_c, sim, scenario,
                               n_strata = 4, csv_path = NULL) {
  n_detected <- nrow(fit_data)
  y <- fit_data$distance
  theta <- fit_data$theta

  r1 <- safe_fit(fit_approach1(y, w))
  r2 <- safe_fit(fit_approach2(y, theta, w))
  r3 <- safe_fit(fit_approach3(y, theta, w, n_strata = n_strata))
  r4 <- safe_fit(fit_approach4(y, w))
  r5 <- safe_fit(fit_approach5(y, w, theta))
  r6 <- safe_fit(fit_approach6(theta))
  r7 <- safe_fit(fit_approach7(theta))
  r8 <- safe_fit(fit_true_shape(full_sim_df, n_detected))

  P_c_uniform <- 1
  N_hat_uniform <- n_detected / P_c_uniform

  row <- data.table::data.table(
    scenario = scenario, sim = sim,
    n_detected = n_detected, true_P_c = true_P_c, true_N_hat = n_detected / true_P_c,
    P_c_plain = g(r1, "P_c"), N_hat_plain = g(r1, "N_hat"), bias_plain = g(r1, "P_c") - true_P_c,
    params_plain = param_string(r1),
    P_c_cov = g(r2, "Pc_hat"), N_hat_cov = g(r2, "N_hat_cov"), bias_cov = g(r2, "Pc_hat") - true_P_c,
    params_cov = param_string(r2),
    P_c_strat = g(r3, "P_c_strat"), N_hat_strat = g(r3, "N_hat_strat"),
    bias_strat = g(r3, "P_c_strat") - true_P_c, n_strata_failed = g(r3, "n_strata_failed"),
    P_c_beta = g(r4, "Pc_hat"), N_hat_beta = g(r4, "N_hat"), bias_beta = g(r4, "Pc_hat") - true_P_c,
    params_beta = param_string(r4),
    P_c_beta_cov = g(r5, "Pc_hat"), N_hat_beta_cov = g(r5, "N_hat"),
    bias_beta_cov = g(r5, "Pc_hat") - true_P_c, params_beta_cov = param_string(r5),
    P_c_beta_theta = g(r6, "Pc_hat"), N_hat_beta_theta = g(r6, "N_hat"),
    bias_beta_theta = g(r6, "Pc_hat") - true_P_c, params_beta_theta = param_string(r6),
    P_c_normal_theta = g(r7, "Pc_hat"), N_hat_normal_theta = g(r7, "N_hat"),
    bias_normal_theta = g(r7, "Pc_hat") - true_P_c, params_normal_theta = param_string(r7),
    P_c_uniform = P_c_uniform, N_hat_uniform = N_hat_uniform, bias_uniform = P_c_uniform - true_P_c,
    P_c_true_anchor0 = g(r8, "P_c_true_shape_anchor0"), N_hat_true_anchor0 = g(r8, "N_hat_true_shape_anchor0"),
    P_c_true_anchormax = g(r8, "P_c_true_shape_anchormax"), N_hat_true_anchormax = g(r8, "N_hat_true_shape_anchormax")
  )

  return(row)
}

nsim <- 10000

list_detected_df <- c(
  "C:/Users/247404/Documents/2024/digital_vs_distance/simulation_results/1_2_detected_df_list_high_r_y_halfnormal_g_y_smooth_scale.parquet",
  "C:/Users/247404/Documents/2024/digital_vs_distance/simulation_results/1_detected_df_list_low_r_y_halfnormal_g_y_smooth_scale.parquet",
  "C:/Users/247404/Documents/2024/digital_vs_distance/simulation_results/2_2_detected_df_list_high_r_y_beta_g_y_smooth_scale.parquet",
  "C:/Users/247404/Documents/2024/digital_vs_distance/simulation_results/2_detected_df_list_low_r_y_beta_g_y_smooth_scale.parquet",
  "C:/Users/247404/Documents/2024/digital_vs_distance/simulation_results/3_detected_df_list_uniform_g_y_smooth_scale.parquet"
)


list_sim_summary <- c(
  "C:/Users/247404/Documents/2024/digital_vs_distance/simulation_results/1_2_sim_summary_list_high_r_y_halfnormal_g_y_smooth_scale.rds",
  "C:/Users/247404/Documents/2024/digital_vs_distance/simulation_results/1_sim_summary_list_low_r_y_halfnormal_g_y_smooth_scale.rds",
  "C:/Users/247404/Documents/2024/digital_vs_distance/simulation_results/2_2_sim_summary_list_high_r_y_beta_g_y_smooth_scale.rds",
  "C:/Users/247404/Documents/2024/digital_vs_distance/simulation_results/2_sim_summary_list_low_r_y_beta_g_y_smooth_scale.rds",
  "C:/Users/247404/Documents/2024/digital_vs_distance/simulation_results/3_sim_summary_list_uniform_g_y_smooth_scale.rds"
)


for (i in 1:5) {
  detected_df <- arrow::read_parquet(list_detected_df[i])
  sim_summary <- data.table::rbindlist(readRDS(list_sim_summary[i]))
  detected_df <- data.table::data.table(detected_df)
  data.table::setkey(detected_df, sim)
  data.table::setkey(sim_summary, sim)

  scenario_name <- tools::file_path_sans_ext(basename(list_sim_summary[i]))
  csv_path <- glue("C:/Users/247404/Documents/2024/digital_vs_distance/simulation_results/results_pc_{scenario_name}5.csv")

  sim_data_1 <- detected_df[sim == 1]
  template <- fit_all_approaches(
    sim_data_1[detected == 1, .(distance = y, theta = theta)],
    sim_data_1,
    w = 1, n_strata = 4,
    true_P_c = sim_data_1[, mean(g)], sim = 1, scenario = scenario_name
  )

  results_table <- template[rep(1, nsim)]
  data.table::set(results_table, 1L, names(template), template)

  for (spick in 2:nsim) {
    sim_data <- detected_df[sim == spick]
    fit_data <- sim_data[detected == 1, .(distance = y, theta = theta)]
    real_Pc <- sim_data[, mean(g)]
    w_this_sim <- 1

    row <- fit_all_approaches(fit_data, sim_data,
      w = w_this_sim, n_strata = 4,
      true_P_c = real_Pc, sim = spick, scenario = scenario_name
    )

    data.table::set(results_table, spick, names(row), row)
  }

  data.table::fwrite(results_table, csv_path)
  gc()
}
