n <- 1e5
V <- 0.2
g0 <- 0.5
nsim <- 1000
w <- 50
detected_df_list <- vector("list", nsim) # y, theta for detected animals only

for (s in seq_len(nsim)) {
  r <- rgamma(3, shape = 1)
  r <- r / sum(r)
  names(r) <- c("y", "theta", "int")

  f_theta_shape <- runif(1, 1, 5)
  f_theta_scale <- runif(1, 1, 5)
  f_theta <- rgamma(n, shape = f_theta_shape, scale = f_theta_scale)

  g_theta_shape <- runif(1, 1, 5)
  g_theta_raw <- pgamma(f_theta, shape = g_theta_shape, rate = 1)
  g_theta <- g_theta_raw - mean(g_theta_raw)

  f_y <- runif(n, 0, w)

  g_y1_scale <- runif(1, w / 4, w / 2)
  g_y1_raw <- exp(-f_y^2 / (2 * g_y1_scale^2))
  g_y1 <- g_y1_raw - mean(g_y1_raw)

  g_y1_theta <- g_y1 * g_theta

  v_y <- var(g_y1)
  v_th <- var(g_theta)
  v_int <- var(g_y1_theta)

  target <- V * r

  a_y <- sqrt(target["y"] / v_y)
  a_th <- sqrt(target["theta"] / v_th)
  a_int <- sqrt(target["int"] / v_int)

  g <- g0 + a_y * g_y1 + a_th * g_theta + a_int * g_y1_theta

  margin <- min(g0, 1 - g0)
  excess <- max(abs(g - g0)) / margin
  k <- if (excess > 1) 1 / excess else 1

  a_y <- a_y * k
  a_th <- a_th * k
  a_int <- a_int * k

  g <- g0 + a_y * g_y1 + a_th * g_theta + a_int * g_y1_theta
  g <- pmin(pmax(g, 0), 1)

  v_y_real <- var(a_y * g_y1)
  v_th_real <- var(a_th * g_theta)
  v_int_real <- var(a_int * g_y1_theta)
  v_total <- v_y_real + v_th_real + v_int_real

  detected <- rbinom(n, 1, g)

  detected_df <- data.frame(
    sim = s,
    y = f_y,
    theta = f_theta,
    g = g,
    target_r_y = r["y"],
    target_r_theta = r["theta"],
    target_r_int = r["int"],
    realized_r_y = v_y_real / v_total,
    realized_r_theta = v_th_real / v_total,
    realized_r_int = v_int_real / v_total,
    V_realized = v_total,
    g_min = min(g),
    g_max = max(g),
    detected = detected
  )
  if (nrow(detected_df) > 0) {
    detected_df <- cbind(object = 1:nrow(detected_df), detected_df)
  }
  detected_df_list[[s]] <- detected_df
}


detected_df <- do.call(rbind, detected_df_list)
dim(detected_df)
