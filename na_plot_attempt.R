library(data.table)
library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)

list_sim_summary <- c(
  "E:/digitalvsdistance_results/new_results/1_2_sim_summary_list_high_r_y_halfnormal_g_y_smooth_scale.rds",
  "E:/digitalvsdistance_results/new_results/1_sim_summary_list_low_r_y_halfnormal_g_y_smooth_scale.rds",
  "E:/digitalvsdistance_results/new_results/2_2_sim_summary_list_high_r_y_beta_g_y_smooth_scale.rds",
  "E:/digitalvsdistance_results/new_results/2_sim_summary_list_low_r_y_beta_g_y_smooth_scale.rds",
  "E:/digitalvsdistance_results/new_results/3_sim_summary_list_high_r_y_halfnormal_g_y_smooth_scale.rds"
)

list_results <- c(
  "E:/digitalvsdistance_results/new_results/results_pc_1_2_sim_summary_list_high_r_y_halfnormal_g_y_smooth_scale4.csv",
  "E:/digitalvsdistance_results/new_results/results_pc_1_sim_summary_list_low_r_y_halfnormal_g_y_smooth_scale4.csv",
  "E:/digitalvsdistance_results/new_results/results_pc_2_2_sim_summary_list_high_r_y_beta_g_y_smooth_scale4.csv",
  "E:/digitalvsdistance_results/new_results/results_pc_2_sim_summary_list_low_r_y_beta_g_y_smooth_scale4.csv",
  "E:/digitalvsdistance_results/new_results/results_pc_3_sim_summary_list_high_r_y_halfnormal_g_y_smooth_scale4.csv"
)

results_l <- list()

for (i in 1:5) {
  sim_summary <- data.table::rbindlist(readRDS(list_sim_summary[i]))
  sim_summary <- data.table::data.table(sim_summary)
  results <- read.csv(list_results[i])
  results <- data.table::data.table(results)
  results_l[[i]] <- results[sim_summary, on = "sim"]
}

names(results_l) <- c("Half normal g_y, high r_y", "half normal g_y", "Beta g_y, high r_y", "Beta g_y", "Uniform g_y")


all_res_joined <- data.table::rbindlist(results_l, fill = TRUE, idcol = "Scenario")
dim(all_res_joined)
colnames(all_res_joined)
table(all_res_joined$Scenario)
# Erase all columns that are just NAs
all_res_joined %>%
  select_if(~ !all(is.na(.))) -> all_res_joined

# select the columns that start with P_c and the realized r_y, r_theta, and r_int
all_res_joined %>%
  select(Scenario, sim, n_detected, true_P_c, true_N_hat, starts_with("P_c"), starts_with("bias"), starts_with("realized_r")) -> all_res_joined

# cut the r_y column into bins
all_res_joined %>%
  mutate(r_y_bin = cut(realized_r_y, breaks = seq(0, 1, 0.1), include.lowest = TRUE)) -> all_res_joined

# plot the number of NAs per P_c column, per r_y bin
windows()
all_res_joined %>%
  select(Scenario, starts_with("P_c"), r_y_bin) %>%
  group_by(r_y_bin, Scenario) %>%
  reframe(across(starts_with("P_c"), ~ sum(is.na(.)))) %>%
  pivot_longer(-c(Scenario, r_y_bin), names_to = "method", values_to = "value") %>%
  ggplot(aes(x = as.factor(r_y_bin), y = value, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Scenario, scales = "free_x", ncol = 2) +
  xlab("r_y bin") +
  ylab("Number of NAs") +
  theme_minimal()


# wichi sims are uniform g_y
results_wchis <- all_res_joined %>%
  filter(Scenario == "Uniform g_y") %>%
  pull(sim)
# the uniform scenario  with the beta method produce the most amount of NAs with 4000 to 6000 simg showing NA we will explore the origin of this
windows()
all_res_joined %>%
  filter(Scenario == "Uniform g_y" & (is.na(P_c_beta) | is.na(P_c_beta_cov))) %>%
  select(realized_r_theta, realized_r_int) %>%
  pivot_longer(everything(), names_to = "component", values_to = "variance") %>%
  ggplot(aes(x = variance, color = component)) +
  geom_density(linewidth = 4) +
  lims(y = c(0, 2))
# facet_wrap(~component, scales = "free")

windows()
all_res_joined %>%
  filter(Scenario == "Uniform g_y" & (!is.na(P_c_beta) | !is.na(P_c_beta_cov))) %>%
  select(realized_r_theta, realized_r_int, bias_beta, bias_beta_cov) %>%
  pivot_longer(c(realized_r_int, realized_r_theta), names_to = "component", values_to = "variance") %>%
  pivot_longer(c(bias_beta, bias_beta_cov), names_to = "method", values_to = "bias") %>%
  ggplot(aes(x = variance, y = abs(bias), grup = component)) +
  geom_point(linewidth = 4) +
  geom_smooth(linewidth = 4) +
  facet_wrap(~method, scales = "free")

## some observation in this group had very low bias which where they?
windows()
all_res_joined %>%
  filter(Scenario == "Uniform g_y" & (!is.na(P_c_beta) | !is.na(P_c_beta_cov))) %>%
  mutate(n_detected_bin = cut(n_detected, breaks = quantile(n_detected, probs = seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)) %>%
  select(n_detected_bin, realized_r_theta, realized_r_int, bias_beta, bias_beta_cov) %>%
  ggplot(aes(x = realized_r_int, y = bias_beta, )) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~n_detected_bin, scales = "free")
# geom_smooth()

all_res_joined %>%
  filter(Scenario == "Uniform g_y") %>%
  ggplot(aes(x = realized_r_int, y = n_detected)) +
  geom_point(alpha = 0.2) +
  geom_smooth()


detected_df %>%
  mutate(y_bin = cut(y, breaks = 30)) %>%
  group_by(y_bin) %>%
  summarise(mean_g = mean(g)) %>%
  clipr::write_clip()

all_res_joined %>%
  filter(Scenario == "Uniform g_y") %>%
  arrange(desc(realized_r_int)) %>%
  head()

detected_df %>%
  filter(sim == 8308, detected == 1) %>%
  mutate(y_bin = cut(y, breaks = 30)) %>%
  group_by(y_bin) %>%
  summarise(n = n()) %>%
  clipr::write_clip()

set.seed(1)
results <- data.frame(n = numeric(), P_c_hat = numeric())

for (i in 1:200) {
  n_i <- sample(200:5000, 1) # vary sample size, like real sims do
  z <- runif(n_i) # perfectly, exactly uniform -- ground truth known
  fit <- fit_approach4(z * 100, w = 100, uniform_tol = 1) # or whatever your actual function signature is
  results <- rbind(results, data.frame(n = n_i, P_c_hat = fit$Pc_hat))
}

ggplot(results, aes(x = n, y = P_c_hat)) +
  geom_point() +
  geom_smooth()


## lets plot all the detection curves of this sims  with respect to distance
chosen_sims <- sample(unique(detected_df$sim), 100) # pick a handful of full curves

detected_df %>%
  filter(detected == 1, sim %in% chosen_sims) %>%
  arrange(sim, y) %>%
  ggplot(aes(x = y, y = g, group = sim)) +
  # geom_point(size = 0.5, alpha = 0.2) +
  geom_smooth(linewidth = 0.8, alpha = 0.5, se = FALSE) +
  scale_color_viridis_c(name = "sim") +
  labs(x = "distance (y)", y = "g") +
  theme_minimal()

detected_df %>%
  filter(sim == 7756, detected == 1) %>%
  mutate(y_bin = cut(y, breaks = 30)) %>%
  group_by(y_bin) %>%
  summarise(var_g = var(g), mean_g = mean(g)) %>%
  pull(var_g) |>
  clipr::write_clip()
for (i in 1:5) {
  r <- read.csv(list_results[i])
  s <- data.table::rbindlist(readRDS(list_sim_summary[i]))
  cat(list_results[i], ": results has", nrow(r), "/ sim_summary has", nrow(s), "\n")
}

all_res_joined %>%
  filter(r_y_bin == "(0.9,1]") %>%
  summarise(
    n = n(),
    na_uniform = sum(is.na(P_c_uniform)),
    na_plain = sum(is.na(P_c_plain)),
    na_normal_theta = sum(is.na(P_c_normal_theta))
  )
