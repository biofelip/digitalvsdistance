library(tidyverse)
library(data.table)
n <- 1e4
V <- 0.2
g0 <- 0.5
nsim <- 100
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

  f_y <- c(0.01,runif(n-1, 0, w))

  g_y1_scale <- runif(1, w / 4, w / 2)
  g_y1_raw <- sqrt(2)/(g_y1_scale * sqrt(pi)) * exp(-f_y^2 / (2 * g_y1_scale^2))
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

  g_raw<- g0 + a_y * g_y1 + a_th * g_theta + a_int * g_y1_theta
  g_raw<- pmin(pmax(g, 0), 1)
  # rescale based on the assumption that g(0,theta)=1 for all theta
  g <- g_raw/g_raw[match(0.01, f_y)]
  # rescale based on the assumption that g(y,theta)=1 at the max point of the detection curve
  # g <- g_raw/max(g_raw)

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

### EXTREME EXAMPLES PLOT
# one row per sim with its realized_r_y, to pick spread-out examples
sim_summary <- detected_df %>%
  distinct(sim, realized_r_y, realized_r_theta, realized_r_int) %>%
  arrange(realized_r_y)

sample_sims <- sim_summary$sim[round(seq(1, nrow(sim_summary), length.out = 10))]

# filter to just those sims AND detected==1 -- small subset, safe to work with
plot_data <- detected_df %>%
  filter(sim %in% sample_sims, detected == 1) %>%
  select(sim, y, realized_r_y, g, g_y1_scale, g_theta_shape)
ggplot(plot_data,
       aes(x = y, y=g, group = sim, color = realized_r_y)) +
  geom_smooth(linewidth = 0.8, alpha = 0.5, se = FALSE) +
  scale_color_viridis_c(name = "realized\nr_y") +
  labs(x = "Detected distance (y)", y = "Density") +
  theme_minimal()

### PLOT the curves with their expected variance prorportion colores

by_var_y <- ggplot(detected_df %>% filter(detected == 1),
       aes(x = y, y=g, group = sim, color = realized_r_y)) +
  geom_smooth(linewidth = 1.5, alpha = 0.3, se = FALSE) +
  scale_color_viridis_c(name = "realized\nr_y") +
  labs(x = "Detected distance (y)", y = "Prob") +
  theme_minimal()

by_var_theta <- ggplot(detected_df %>% filter(detected == 1),
       aes(x = y, y=g, group = sim, color = realized_r_theta)) +
  geom_smooth(linewidth = 1.5, alpha = 0.5, se = FALSE) +
  scale_color_viridis_c(name = "realized\nr_theta") +
  labs(x = "Detected distance (y)", y = "Prob") +
  theme_minimal()

by_var_int <- ggplot(detected_df %>% filter(detected == 1),
       aes(x = y, y=g, group = sim, color = realized_r_int)) +
  geom_smooth(linewidth = 1.5, alpha = 0.5, se = FALSE) +
  scale_color_viridis_c(name = "realized\nr_int") +
  labs(x = "Detected distance (y)", y = "Prob") +
  theme_minimal()


cowplot::plot_grid(by_var_y, by_var_theta, by_var_int, ncol = 1)

## PLOT a tern plot with the residuals from  a half normal with the correct parameter scaleparamete . 

detected_df <- detected_df %>%
  mutate(
    hn_true    = sqrt(2)/(g_y1_scale * sqrt(pi)) * exp(-f_y^2 / (2 * g_y1_scale^2)),
    residual_hn = g - hn_true
  )

sim_summary <- detected_df %>%
  filter(detected == 1) %>%
  group_by(sim) %>%
  summarise(
    realized_r_y     = first(realized_r_y),
    realized_r_theta = first(realized_r_theta),
    realized_r_int   = first(realized_r_int),
    mean_residual_hn  = mean(residual_hn),
    rmse_hn           = sqrt(mean(residual_hn^2)),
    .groups = "drop"
  )

library(ggtern)

ggtern(detected_df, aes(x = realized_r_y, y = realized_r_theta, z = realized_r_int,
                         color = residual_hn)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c(name = "RMSE vs\ntrue HN") +
  theme_bw()

# export to a csv file
write.csv(detected_df, file = "results_sim_df.csv", row.names = FALSE)
