
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
ggplot(
  plot_data,
  aes(x = y, y = g, group = sim, color = realized_r_y)
) +
  geom_smooth(linewidth = 0.8, alpha = 0.5, se = FALSE) +
  scale_color_viridis_c(name = "realized\nr_y") +
  labs(x = "Detected distance (y)", y = "Density") +
  theme_minimal()

### PLOT the curves with their expected variance prorportion colores

by_var_y <- ggplot(
  detected_df %>% filter(detected == 1),
  aes(x = y, y = g, group = sim, color = realized_r_y)
) +
  geom_smooth(linewidth = 1.5, alpha = 0.3, se = FALSE) +
  scale_color_viridis_c(name = "realized\nr_y") +
  labs(x = "Detected distance (y)", y = "Prob") +
  theme_minimal()

by_var_theta <- ggplot(
  detected_df %>% filter(detected == 1),
  aes(x = y, y = g, group = sim, color = realized_r_theta)
) +
  geom_smooth(linewidth = 1.5, alpha = 0.5, se = FALSE) +
  scale_color_viridis_c(name = "realized\nr_theta") +
  labs(x = "Detected distance (y)", y = "Prob") +
  theme_minimal()

by_var_int <- ggplot(
  detected_df %>% filter(detected == 1),
  aes(x = y, y = g, group = sim, color = realized_r_int)
) +
  geom_smooth(linewidth = 1.5, alpha = 0.5, se = FALSE) +
  scale_color_viridis_c(name = "realized\nr_int") +
  labs(x = "Detected distance (y)", y = "Prob") +
  theme_minimal()


cowplot::plot_grid(by_var_y, by_var_theta, by_var_int, ncol = 1)

## PLOT a tern plot with the residuals from  a half normal with the correct parameter scaleparamete .

detected_df <- detected_df %>%
  mutate(
    hn_true = sqrt(2) / (g_y1_scale * sqrt(pi)) * exp(-f_y^2 / (2 * g_y1_scale^2)),
    residual_hn = g - hn_true
  )

sim_summary <- detected_df %>%
  filter(detected == 1) %>%
  group_by(sim) %>%
  summarise(
    realized_r_y = first(realized_r_y),
    realized_r_theta = first(realized_r_theta),
    realized_r_int = first(realized_r_int),
    mean_residual_hn = mean(residual_hn),
    rmse_hn = sqrt(mean(residual_hn^2)),
    .groups = "drop"
  )
sim_summary2 <- left_join(sim_summary, results_pc, by = "sim")

# plot  SE against realized_r_y
sim_summary2 %>%
  pivot_longer(cols = c(P_c, P_c_cov, P_c_unif), names_to = "key", values_to = "P_c") %>%
  ggplot(aes(x = realized_r_y, y = P_c, color = realized_r_theta)) +
  geom_point() +
  geom_smooth() +
  geom_hline(yintercept = .5, color = "red") +
  facet_wrap(~key) +
  scale_color_viridis_c(name = "realized\nr_theta")

sim_summary2 %>%
  select(-P_c_unif) %>%
  pivot_longer(cols = c(P_c, P_c_cov), names_to = "key", values_to = "P_c") %>%
  ggtern(aes(x = realized_r_y, y = realized_r_theta, z = realized_r_int, value = abs(P_c - 0.5))) +
  stat_interpolate_tern(
    geom = "polygon", aes(fill = after_stat(level)),
    formula = value ~ x + y, method = lm, n = 100,
    breaks = seq(0, 1, by = 0.05), # fewer, wider bands
    color = "black", linewidth = 0.2
  ) +
  # scale_fill_viridis_c(name = "Absolute error in P_c" ) +
  scale_fill_distiller(
    palette = "Spectral", direction = -1, name = "P_c",
    limits = c(0, 0.5)
  ) +
  theme_bw() +
  facet_wrap(~key)

ggtern(sim_summary2, aes(x = realized_r_y, y = realized_r_theta, z = realized_r_int)) +
  geom_point(aes(color = abs(P_c_unif - 0.5)), size = 10, color = "Red") +
  theme_bw()



sim_summary2_long <- sim_summary2 %>%
  select(-P_c_unif) %>%
  pivot_longer(cols = c(P_c, P_c_cov), names_to = "key", values_to = "P_c")


Ternary_Plot_custom <- function(simulation_summary) {
  abc <- with(
    simulation_summary %>% filter(!is.na(P_c)),
    rbind(realized_r_y, realized_r_theta, realized_r_int)
  )
  dim(abc)
  response <- abs(simulation_summary %>% filter(!is.na(P_c)) %>% pull(P_c) - 0.5)
  length(response)
  TernaryPlot(alab = "V(y)", blab = "V(theta)", clab = "V(int)")

  xy <- TernaryToXY(abc)

  # Use an inverse distance weighting to interpolate between measured points
  Predict <- function(predXY) {
    Distance <- function(a, b) {
      apply(a, 2, function(pt) sqrt(colSums((pt - b)^2)))
    }
    dists <- Distance(xy, predXY)
    id <- 1 / dists
    idw <- id / rowSums(id)

    # Return:
    colSums(response * t(idw))
  }

  # Predict at triangle centres
  tri <- TriangleCentres(resolution = 48L)
  # Adjust the resolution to suit your own dataset

  # Now we interpolate between our known values to generate a colour for each
  # of our tiles
  predicted <- Predict(tri[1:2, ])
  map <- rbind(
    x = tri["x", ], y = tri["y", ], z = predicted,
    down = tri["triDown", ]
  )

  # Place a semitransparent colour fill over grid lines:
  ColourTernary(map)

  # Calculate contours
  PredictABC <- function(a, b, c) Predict(TernaryToXY(rbind(a, b, c)))
  TernaryContour(PredictABC,
    resolution = 36L, legend = 6,
    legend... = list(bty = "n")
  )

  # Mark the points at which we took measurements
  TernaryPoints(abc, pch = 3, col = "#cc3333")
}

Ternary_Plot_custom(sim_summary2_long %>% filter(key == "P_c_cov"))
Ternary_Plot_custom(sim_summary2_long %>% filter(key == "P_c"))


