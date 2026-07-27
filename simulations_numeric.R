library(tidyverse)
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
    detected = detected,
    g_y1_scale = g_y1_scale,
    g_theta_shape = g_theta_shape
  )
  if (nrow(detected_df) > 0) {
    detected_df <- cbind(object = 1:nrow(detected_df), detected_df)
  }
  detected_df_list[[s]] <- detected_df
}


detected_df <- do.call(rbind, detected_df_list)
detected_df <- data.table(detected_df)


# one row per sim with its realized_r_y, to pick spread-out examples
sim_summary <- detected_df %>%
  distinct(sim, realized_r_y) %>%
  arrange(realized_r_y)

sample_sims <- sim_summary$sim[round(seq(1, nrow(sim_summary), length.out = 6))]

# filter to just those sims AND detected==1 -- small subset, safe to work with
plot_data <- detected_df %>%
  filter(sim %in% sample_sims, detected == 1) %>%
  select(sim, y, realized_r_y)

ggplot(detected_df %>% filter(detected == 1), aes(x = y, color = , color = factor(round(realized_r_y, 2)))) +
  geom_density(linewidth = 0.8) +
  labs(
    x = "Detected distance (y)", y = "Density",
    color = "realized\nr_y"
  ) +
  theme_minimal()



# grab one scale value per chosen sim
scale_lookup <- detected_df %>%
  filter(sim %in% sample_sims) %>%
  distinct(sim, g_y1_scale, realized_r_y)

# build theoretical half-normal density curves, normalized to integrate to 1 over [0,w]
theory_curves <- scale_lookup %>%
  rowwise() %>%
  reframe({
    y_grid <- seq(0, w, length.out = 256)
    raw <- exp(-y_grid^2 / (2 * g_y1_scale^2))
    norm_const <- integrate(function(y) exp(-y^2 / (2 * g_y1_scale^2)), 0, w)$value
    data.frame(
      sim = sim, realized_r_y = realized_r_y,
      y_grid = y_grid, dens = raw / norm_const
    )
  })

# empirical densities, same sims as before
plot_data <- detected_df %>%
  filter(sim %in% sample_sims, detected == 1) %>%
  select(sim, y, realized_r_y)

ggplot() +
  geom_density(
    data = plot_data, aes(
      x = y, group = sim,
      color = factor(round(realized_r_y, 2))
    ),
    linewidth = 0.8
  ) +
  geom_line(
    data = theory_curves, aes(
      x = y_grid, y = dens, group = sim,
      color = factor(round(realized_r_y, 2))
    ),
    linetype = "dashed", linewidth = 0.6
  ) +
  labs(
    x = "Detected distance (y)", y = "Density",
    color = "realized\nr_y",
    title = "Solid = observed detected distances, Dashed = theoretical half-normal"
  ) +
  theme_minimal()
