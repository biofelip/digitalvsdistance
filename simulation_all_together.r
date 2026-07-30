library(future.apply)
library(data.table)
library(mrds)

plan(multisession, workers = parallel::detectCores() - 1)

nsim <- 10000
n <- 1e4

run_one_sim <- function(s, n) {
  V  <- runif(1, 0, 1)
  g0 <- runif(1, 0, 1)
  w  <- runif(1, 1, 50) * 10

  r <- rgamma(3, shape = 1); r <- r / sum(r); names(r) <- c("y", "theta", "int")

  f_theta_shape <- runif(1, 1, 5)
  f_theta_scale <- runif(1, 1, 5)
  f_theta <- rgamma(n, shape = f_theta_shape, scale = f_theta_scale)

  g_theta_shape <- runif(1, 1, 5)
  g_theta_raw <- pgamma(f_theta, shape = g_theta_shape, rate = 1)
  g_theta <- g_theta_raw - mean(g_theta_raw)

  f_y <- c(0.01, runif(n - 1, 0, w))
  g_y1_scale <- runif(1, w / 4, w / 2)
  g_y1_raw <- sqrt(2) / (g_y1_scale * sqrt(pi)) * exp(-f_y^2 / (2 * g_y1_scale^2))
  g_y1 <- g_y1_raw - mean(g_y1_raw)
  g_y1_theta <- g_y1 * g_theta

  v_y <- var(g_y1); v_th <- var(g_theta); v_int <- var(g_y1_theta)
  target <- V * r
  a_y   <- sqrt(target["y"]     / v_y)
  a_th  <- sqrt(target["theta"] / v_th)
  a_int <- sqrt(target["int"]   / v_int)

  g <- g0 + a_y*g_y1 + a_th*g_theta + a_int*g_y1_theta
  margin <- min(g0, 1 - g0)
  excess <- max(abs(g - g0)) / margin
  k <- if (excess > 1) 1 / excess else 1
  a_y <- a_y*k; a_th <- a_th*k; a_int <- a_int*k

  g_raw <- g0 + a_y*g_y1 + a_th*g_theta + a_int*g_y1_theta
  g_raw <- pmin(pmax(g_raw, 0), 1)             # BUG FIX: was clamping the un-shrunk `g`
  g <- g_raw / g_raw[match(0.01, f_y)]         # shoulder rescale, g(0,theta)=1
  g <- pmin(pmax(g, 0), 1)
  # hard safety net: if anything is still NA/NaN/Inf, this sim is unusable -- bail out cleanly
  if (anyNA(g) || any(!is.finite(g))) {
    return(data.frame(sim = s, error = "invalid_g", 
                      V = V, g0_true = g0, w = w))  # minimal row, rest NA by absence
  }
  v_y_real <- var(a_y*g_y1); v_th_real <- var(a_th*g_theta); v_int_real <- var(a_int*g_y1_theta)
  v_total <- v_y_real + v_th_real + v_int_real

  detected <- rbinom(n, 1, g)
  n_detected <- sum(detected)

  fit_data <- data.frame(object = seq_len(n_detected),
                          distance = f_y[detected == 1],
                          theta = f_theta[detected == 1],
                          size = 1)

  fit_hn <- tryCatch(ddf(dsmodel = ~mcds(key="hn", formula=~1), data=fit_data,
                         method="ds", meta.data=list(width=w)), error=function(e) NULL)
  P_c <- tryCatch(if(!is.null(fit_hn)) summary(fit_hn)$average.p else NA_real_, error=function(e) NA_real_)
  if (length(P_c)==0 || !is.numeric(P_c)) P_c <- NA_real_
  N_hat <- if (!is.na(P_c)) n_detected/P_c else NA_real_

  fit_hn_cov <- tryCatch(ddf(dsmodel = ~mcds(key="hn", formula=~theta), data=fit_data,
                             method="ds", meta.data=list(width=w)), error=function(e) NULL)
  if (!is.null(fit_hn_cov)) {
    p_hat_cov <- tryCatch(predict(fit_hn_cov, integrate=TRUE)$fitted, error=function(e) NULL)
    P_c_cov <- if(!is.null(p_hat_cov)) mean(p_hat_cov) else NA_real_
    N_hat_cov <- if(!is.null(p_hat_cov)) sum(1/p_hat_cov) else NA_real_
  } else { P_c_cov <- NA_real_; N_hat_cov <- NA_real_ }

  fit_unif <- tryCatch(ddf(dsmodel = ~mcds(key="unif", formula=~1), data=fit_data,
                          method="ds", meta.data=list(width=w)), error=function(e) NULL)
  P_c_unif <- tryCatch(if(!is.null(fit_unif)) summary(fit_unif)$average.p else NA_real_, error=function(e) NA_real_)
  if (length(P_c_unif)==0 || !is.numeric(P_c_unif)) P_c_unif <- NA_real_
  N_hat_unif <- if (!is.na(P_c_unif)) n_detected/P_c_unif else NA_real_

  data.frame(
    sim = s, n_population = n, n_detected = n_detected,
    V = V, g0_true = g0, w = w,
    f_theta_shape = f_theta_shape, f_theta_scale = f_theta_scale,
    g_y1_scale = g_y1_scale, g_theta_shape = g_theta_shape,
    target_r_y = r[["y"]], target_r_theta = r[["theta"]], target_r_int = r[["int"]],
    realized_r_y = v_y_real/v_total, realized_r_theta = v_th_real/v_total, realized_r_int = v_int_real/v_total,
    V_realized = v_total, k_shrink = k,
    P_c = P_c, P_c_cov = P_c_cov, P_c_unif = P_c_unif,
    N_hat = N_hat, N_hat_cov = N_hat_cov, N_hat_unif = N_hat_unif,
    N_hat_strip = n_detected
  )
}

results_list <- future_lapply(seq_len(nsim), run_one_sim, n = n, future.seed = TRUE)
results_df <- rbindlist(results_list)
fwrite(results_df, "results_summary_all_in_one.csv")
