library(patchwork)
library(ggplot2)
library(glue)
plot_sim_diagnostics <- function(detected_df, sim_summary, sim_id,
                                  y_shape = c("beta", "halfnormal", "uniform"),
                                  n_grid = 300, n_theta_sample = 20000) {

  y_shape <- match.arg(y_shape)

  d <- detected_df[sim == sim_id]
  detected_d <- d[detected == 1]
  s_row <- sim_summary[sim == sim_id]

  w             <- s_row$w
  g0            <- s_row$g_0
  f_theta_alfa  <- s_row$f_theta_alfa
  f_theta_beta  <- s_row$f_theta_beta
  g_theta_alfa  <- s_row$g_theta_alfa
  g_theta_beta  <- s_row$g_theta_beta
  rho           <- s_row$copula_rho
  realized_r_y  <- s_row$realized_r_y
  realized_r_th <- s_row$realized_r_theta
  realized_r_int<- s_row$realized_r_int
  V_realized    <- s_row$V_realized

  y_grid <- seq(0, w, length.out = n_grid)
  theta_sample <- rbeta(n_theta_sample, f_theta_alfa, f_theta_beta)
  u_y_grid <- punif(y_grid, 0, w)

  g_y1_raw <- switch(y_shape,
    beta       = pbeta(u_y_grid, s_row$g_y1_alfa, s_row$g_y1_beta),
    halfnormal = { sc <- s_row$g_y1_scale; sqrt(2)/(sc*sqrt(pi)) * exp(-u_y_grid^2/(2*sc^2)) },
    uniform    = rep(1, n_grid)
  )
  g_y1_c <- g_y1_raw - mean(g_y1_raw)

  u_theta_sample <- pbeta(theta_sample, f_theta_alfa, f_theta_beta)
  g_theta_raw <- pbeta(u_theta_sample, g_theta_alfa, g_theta_beta)
  g_theta_c   <- g_theta_raw - mean(g_theta_raw)

  z_y     <- qnorm(u_y_grid)
  z_theta <- qnorm(u_theta_sample)
  Z_Y  <- matrix(z_y, nrow = n_grid, ncol = n_theta_sample)
  Z_TH <- matrix(z_theta, nrow = n_grid, ncol = n_theta_sample, byrow = TRUE)
  Z_TH_corr <- rho * Z_Y + sqrt(1 - rho^2) * Z_TH
  V_Y  <- pnorm(Z_Y); V_TH <- pnorm(Z_TH_corr)
  G_int_raw <- V_Y^2 * V_TH
  G_int_c   <- G_int_raw - mean(G_int_raw)

  v_y   <- var(g_y1_c)
  v_th  <- var(g_theta_c)
  v_int <- var(as.vector(G_int_c))

  a_y   <- if (y_shape == "uniform") 0 else sqrt(realized_r_y * V_realized / v_y)
  a_th  <- sqrt(realized_r_th  * V_realized / v_th)
  a_int <- sqrt(realized_r_int * V_realized / v_int)

  # --- Panel 1: distance component -- theoretical g_y1(y), averaged over theta, vs observed g~y ---
  theo_y_curve <- data.frame(y = y_grid, g = g0 + a_y * g_y1_c)

  p1 <- ggplot() +
    geom_point(data = d, aes(x = y, y = g), alpha = 0.05, size = 0.3) +
    geom_smooth(data = detected_d, aes(x = y, y = g, color = "Observed"), se = FALSE) +
    geom_line(data = theo_y_curve, aes(x = y, y = g, color = "Theoretical"), linewidth = 1) +
    scale_color_manual(values = c(Theoretical = "black", Observed = "firebrick")) +
    labs(title = "Distance component: g vs y", color = NULL) + theme_minimal()

  # --- Panel 2: theta component -- theoretical g_theta(theta), averaged over y, vs observed g~theta ---
  theta_grid <- seq(0.001, 0.999, length.out = n_grid)
  g_theta_grid_raw <- pbeta(theta_grid, g_theta_alfa, g_theta_beta)
  g_theta_grid_c   <- g_theta_grid_raw - mean(g_theta_grid_raw)
  theo_theta_curve <- data.frame(theta = theta_grid, g = g0 + a_th * g_theta_grid_c)

  p2 <- ggplot() +
    geom_point(data = d, aes(x = theta, y = g), alpha = 0.05, size = 0.3) +
    geom_smooth(data = detected_d, aes(x = theta, y = g, color = "Observed"), se = FALSE) +
    geom_line(data = theo_theta_curve, aes(x = theta, y = g, color = "Theoretical"), linewidth = 1) +
    scale_color_manual(values = c(Theoretical = "black", Observed = "firebrick")) +
    labs(title = "Theta component: g vs theta", color = NULL) + theme_minimal()

  # --- Panel 3: interaction -- theoretical, stratified by theta level, vs observed residual ---
  theta_levels <- quantile(theta_sample, c(0.1, 0.5, 0.9))

  theo_int <- do.call(rbind, lapply(seq_along(theta_levels), function(i) {
    th <- theta_levels[i]
    u_th <- pbeta(th, f_theta_alfa, f_theta_beta)
    z_th <- qnorm(u_th)
    z_th_corr <- rho * z_y + sqrt(1-rho^2) * z_th
    v_th_i <- pnorm(z_th_corr)
    v_y_i <- pnorm(z_y)
    g_int_raw_i <- v_y_i^2 * v_th_i - mean(G_int_raw)   # center using same constant as full surface
    data.frame(y = y_grid, interaction = a_int * g_int_raw_i,
               theta_level = names(theta_levels)[i])
  }))

  # observed: residual after removing the additive y- and theta-effects, evaluated at each animal's own values
  g_y1_obs <- switch(y_shape,
    beta       = pbeta(punif(d$y, 0, w), s_row$g_y1_alfa, s_row$g_y1_beta) - mean(g_y1_raw),
    halfnormal = { sc <- s_row$g_y1_scale; sqrt(2)/(sc*sqrt(pi))*exp(-punif(d$y,0,w)^2/(2*sc^2)) - mean(g_y1_raw) },
    uniform    = rep(0, nrow(d))
  )
  g_theta_obs <- pbeta(pbeta(d$theta, f_theta_alfa, f_theta_beta), g_theta_alfa, g_theta_beta) - mean(g_theta_raw)
  d$interaction_resid <- d$g - g0 - a_y*g_y1_obs - a_th*g_theta_obs
  d$theta_level <- cut(d$theta, breaks = c(-Inf, theta_levels[1], theta_levels[3], Inf),
                        labels = c("10%", "50%", "90%"))

  p3 <- ggplot() +
    geom_point(data = d[detected==1], aes(x = y, y = interaction_resid, color = theta_level),
               alpha = 0.1, size = 0.3) +
    geom_line(data = theo_int, aes(x = y, y = interaction, color = theta_level), linewidth = 1) +
    labs(title = "Interaction: contribution vs y, by theta level", color = "theta level") +
    theme_minimal()

  (p1 | p2) / p3 + plot_annotation(title = glue("Simulation diagnostics - sim {sim_id}, realized r_y = {round(realized_r_y,2)}, realized r_theta = {round(realized_r_th,2)}, realized r_int = {round(realized_r_int,2)}"))
}
