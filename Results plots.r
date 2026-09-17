# Results plots
library(ggplot2)
library(dplyr)
library(magrittr)
library(tidyr)
library(glue)
# load all the results files
results_files <- c(
  "E:/digitalvsdistance_results/new_results/results_pc_1_2_sim_summary_list_high_r_y_halfnormal_g_y_smooth_scale3.csv",
  "E:/digitalvsdistance_results/new_results/results_pc_1_sim_summary_list_low_r_y_halfnormal_g_y_smooth_scale3.csv",
  "E:/digitalvsdistance_results/new_results/results_pc_2_2_sim_summary_list_high_r_y_beta_g_y_smooth_scale3.csv",
  "E:/digitalvsdistance_results/new_results/results_pc_2_sim_summary_list_low_r_y_beta_g_y_smooth_scale3.csv",
  "E:/digitalvsdistance_results/new_results/results_pc_3_sim_summary_list_high_r_y_halfnormal_g_y_smooth_scale3.csv"
)

summary_files <- c(
  "E:/digitalvsdistance_results/new_results/1_2_sim_summary_list_high_r_y_halfnormal_g_y_smooth_scale.rds",
  "E:/digitalvsdistance_results/new_results/1_sim_summary_list_low_r_y_halfnormal_g_y_smooth_scale.rds",
  "E:/digitalvsdistance_results/new_results/2_2_sim_summary_list_high_r_y_beta_g_y_smooth_scale.rds",
  "E:/digitalvsdistance_results/new_results/2_sim_summary_list_low_r_y_beta_g_y_smooth_scale.rds",
  "E:/digitalvsdistance_results/new_results/3_sim_summary_list_high_r_y_halfnormal_g_y_smooth_scale.rds"
)
summary_list <- lapply(summary_files, readRDS)
summary_list <- lapply(summary_list, data.table::rbindlist)
results_list <- lapply(results_files, read.csv)
names(results_list) <- names(summary_list) <- c("Half normal g_y, high r_y", "half normal g_y", "Beta g_y, high r_y", "Beta g_y", "Uniform g_y")

results_pc <- data.table::rbindlist(results_list, idcol = "Scenario")
summary_sim_pc <- data.table::rbindlist(summary_list, idcol = "Scenario", fill = TRUE)
results_pc <- results_pc[summary_sim_pc, on = c("Scenario", "sim")]
results_pc$Scenario %>% table()
results_pc$Scenario <- factor(results_pc$Scenario, levels = c("Half normal g_y, high r_y", "half normal g_y", "Beta g_y, high r_y", "Beta g_y", "Uniform g_y"))
dim(results_pc)

plot_param_dist <- function(results_pc, scenarios) {
  # g_y
  p1 <- results_pc %>%
    filter(Scenario %in% scenarios) %>%
    ggplot(aes(x = g_y1_scale)) +
    geom_density() +
    facet_wrap(~Scenario)
  # f_theta
  p2 <- results_pc %>%
    filter(Scenario %in% scenarios) %>%
    ggplot(aes(x = f_theta_alfa, y = f_theta_beta)) +
    geom_bin_2d() +
    facet_wrap(~Scenario)
  # g_theta
  p3 <- results_pc %>%
    filter(Scenario %in% scenarios) %>%
    ggplot(aes(x = g_theta_alfa, y = g_theta_beta)) +
    geom_bin_2d() +
    facet_wrap(~Scenario)
  pf <- cowplot::plot_grid(plotlist = list(p1, p2, p3)) +
    ggtitle(scenarios)
  return(pf)
}

# check for missing values in the columns
apply(results_pc, 2, function(x) sum(is.na(x))) %>% sort()
# check the distribution of parameters randomly generated
### Scenarios of the half normal: look good quite homogenous around the  range of values
windows()
plot_param_dist(results_pc, c("Half normal g_y, high r_y", "half normal g_y"))
plot_param_dist(results_pc, c("Beta g_y, high r_y", "Beta g_y"))
plot_param_dist(results_pc, c("Uniform g_y"))



# 1 Bias with change in r_y
# 1.1 Just r_y
windows()
results_pc |>
  filter(Scenario != "Uniform g_y") %>%
  pivot_longer(c(
    "P_c_plain", "P_c_cov", "P_c_strat", "P_c_beta", "P_c_beta_cov",
    "P_c_beta_theta", "P_c_normal_theta", "P_c_uniform", "P_c_true_anchor0",
    "P_c_true_anchormax"
  ), names_to = "approach", values_to = "Pc") |>
  mutate(approach = case_when(
    approach == "P_c_plain" ~ "DS traditional",
    approach == "P_c_cov" ~ "MCDS",
    approach == "P_c_strat" ~ "Stratification",
    approach == "P_c_beta" ~ "Beta shape on distance",
    approach == "P_c_beta_cov" ~ "Beta shape on distance with MCDS",
    approach == "P_c_beta_theta" ~ "Beta shape on theta",
    approach == "P_c_normal_theta" ~ "Normal shape on theta",
    approach == "P_c_uniform" ~ "Strip sampling",
    approach == "P_c_true_anchor0" ~ "Ideal estimation, anchor at g(0)",
    approach == "P_c_true_anchormax" ~ "Ideal estimation, anchor at g(y*)"
  )) %>%
  filter(approach != "Ideal estimation, anchor at g(0)") %>%
  mutate(approach = forcats::fct_reorder(approach, abs(Pc - true_P_c), .fun = median, .desc = FALSE)) %>%
  ggplot(aes(x = realized_r_y, y = abs(Pc - true_P_c), color = approach)) +
  geom_smooth() +
  xlab("Realized r_y") +
  ylab("Bias") +
  theme_minimal() +
  facet_wrap(~Scenario, scales = "free")

results_pc |>
  filter(Scenario != "Uniform g_y") %>%
  pivot_longer(c(
    "P_c_plain", "P_c_cov", "P_c_strat", "P_c_beta", "P_c_beta_cov",
    "P_c_beta_theta", "P_c_normal_theta", "P_c_uniform", "P_c_true_anchor0",
    "P_c_true_anchormax"
  ), names_to = "approach", values_to = "Pc") |>
  mutate(approach = case_when(
    approach == "P_c_plain" ~ "DS traditional",
    approach == "P_c_cov" ~ "MCDS",
    approach == "P_c_strat" ~ "Stratification",
    approach == "P_c_beta" ~ "Beta shape on distance",
    approach == "P_c_beta_cov" ~ "Beta shape on distance with MCDS",
    approach == "P_c_beta_theta" ~ "Beta shape on theta",
    approach == "P_c_normal_theta" ~ "Normal shape on theta",
    approach == "P_c_uniform" ~ "Strip sampling",
    approach == "P_c_true_anchor0" ~ "Ideal estimation, anchor at g(0)",
    approach == "P_c_true_anchormax" ~ "Ideal estimation, anchor at g(y*)"
  )) %>%
  filter(approach != "Ideal estimation, anchor at g(0)") %>%
  ggplot(aes(x = approach, y = abs(Pc - true_P_c), fill = approach)) +
  geom_boxplot(na.rm = TRUE) +
  theme_minimal()

results_pc |>
  filter(Scenario != "Uniform g_y") %>%
  pivot_longer(c(
    "P_c_plain", "P_c_cov", "P_c_strat", "P_c_beta", "P_c_beta_cov",
    "P_c_beta_theta", "P_c_normal_theta", "P_c_uniform", "P_c_true_anchor0",
    "P_c_true_anchormax"
  ), names_to = "approach", values_to = "Pc") |>
  mutate(approach = case_when(
    approach == "P_c_plain" ~ "DS traditional",
    approach == "P_c_cov" ~ "MCDS",
    approach == "P_c_strat" ~ "Stratification",
    approach == "P_c_beta" ~ "Beta shape on distance",
    approach == "P_c_beta_cov" ~ "Beta shape on distance with MCDS",
    approach == "P_c_beta_theta" ~ "Beta shape on theta",
    approach == "P_c_normal_theta" ~ "Normal shape on theta",
    approach == "P_c_uniform" ~ "Strip sampling",
    approach == "P_c_true_anchor0" ~ "Ideal estimation, anchor at g(0)",
    approach == "P_c_true_anchormax" ~ "Ideal estimation, anchor at g(y*)"
  )) %>%
  filter(approach != "Ideal estimation, anchor at g(0)") %>%
  group_by(Scenario, approach) %>%
  summarise(pct_na = mean(is.na(Pc)), .groups = "drop") %>%
  mutate(approach = forcats::fct_reorder(approach, pct_na, .fun = mean, .desc = TRUE)) %>%
  filter(approach != "Strip sampling") %>%
  ggplot(aes(x = approach, y = pct_na, fill = approach)) +
  geom_col() +
  facet_wrap(~Scenario) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "non-convergence rate", fill = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Max point of the detection curve against r_y and scenario
# ridge plot of max_y against r_y and scenario
# install.packages("ggridges")
library(ggridges)

# This ridges show some tail to the right despite the fact that the
# distribution of parameters for the betas was homogenous and did not seem
# to favor compbinations of lower alpha than beta (which will result in a increeasing function)
results_pc %>%
  filter(!is.na(realized_r_y)) %>%
  mutate(
    u_max_y = max_y / w, # normalize to [0,1] so it's comparable across sims with different w
    r_y_bin = cut(realized_r_y, breaks = seq(0, 1, 0.2), include.lowest = TRUE)
  ) %>%
  ggplot(aes(x = u_max_y, y = r_y_bin, fill = r_y_bin)) +
  geom_density_ridges(alpha = 0.7, scale = 1.2, na.rm = TRUE) +
  facet_wrap(~Scenario, scales = "free") +
  scale_fill_viridis_d(guide = "none") +
  labs(
    x = "location of max detection (as fraction of strip width)",
    y = "realized r_y bin"
  )

# Explore the missignes patterns for functions with low r_y to try to
# explain why the maximun of the cuntion accumulates more towards the end
windows()
# the pattern of of missines does not seem to accumulate on specific places
# of the distribution of theta
results_pc |>
  filter(realized_r_y < 0.5) %>%
  ggplot(aes(x = f_theta_alfa, y = f_theta_beta, group = n_strata_failed)) +
  geom_bin_2d() +
  facet_wrap(~Scenario, scales = "free")
# neither the paramters of the of the g_y,theta function
results_pc |>
  filter(realized_r_y < 0.5) %>%
  ggplot(aes(x = g_theta_alfa, y = g_theta_beta, group = n_strata_failed)) +
  geom_bin_2d() +
  facet_wrap(~Scenario, scales = "free")

# but as the proportion of r_y increases the porportion of strata that fail
# increases
results_pc %>%
  ggplot(aes(x = realized_r_y, y = n_strata_failed)) +
  geom_smooth()
# an maybe in those that failes we ended up with more common betas with decreasein shape
results_pc %>%
  filter(n_strata_failed > 0) %>%
  ggplot(aes(x = g_theta_alfa, y = g_theta_beta)) +
  geom_bin_2d()
## ok no that is not it, but then probably is because as the r_y decreases there
## more chances for the interacion to be more important and probably the rhos at that shape were
## going to positive.
windows()
results_pc %>%
  filter(n_strata_failed > 0) %>%
  ggplot(aes(x = realized_r_y, y = realized_r_int, color = Scenario)) +
  geom_point()

results_pc %>%
  mutate(true_Pc_bin = cut(true_P_c, breaks = seq(0, 1, 0.25), include.lowest = TRUE)) %>%
  pivot_longer(c(N_hat_plain, N_hat_cov, N_hat_strat, N_hat_uniform, N_hat_true_shape, N_hat_flexible),
    names_to = "approach", values_to = "Nhat"
  ) %>%
  mutate(approach = case_when(
    approach == "N_hat_plain" ~ "DS traditional",
    approach == "N_hat_cov" ~ "MCDS",
    approach == "N_hat_strat" ~ "Stratification",
    approach == "N_hat_uniform" ~ "Strip sampling",
    approach == "N_hat_true_shape" ~ "Ideal estimation",
    approach == "N_hat_flexible" ~ "Flexible shape model"
  )) %>%
  filter(!is.na(Nhat)) %>%
  mutate(approach = forcats::fct_reorder(approach, abs(Nhat - 10000), .fun = median, na.rm = TRUE)) %>%
  ggplot(aes(x = approach, y = Nhat, color = approach)) +
  geom_boxplot(na.rm = TRUE) +
  geom_hline(yintercept = 10000, linetype = "dashed", color = "gray40") +
  scale_y_log10() +
  facet_grid(true_Pc_bin ~ Scenario, scales = "free_y") +
  labs(x = NULL, y = expression(hat(N)), color = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Reverse the effect in other ways.
results_pc %>%
  mutate(true_N = n_detected / true_P_c) %>%
  pivot_longer(c(N_hat_plain, N_hat_cov, N_hat_strat, N_hat_uniform, N_hat_true_shape, N_hat_flexible),
    names_to = "approach", values_to = "Nhat"
  ) %>%
  mutate(approach = case_when(
    approach == "N_hat_plain" ~ "DS traditional",
    approach == "N_hat_cov" ~ "MCDS",
    approach == "N_hat_strat" ~ "Stratification",
    approach == "N_hat_uniform" ~ "Strip sampling",
    approach == "N_hat_true_shape" ~ "Ideal estimation",
    approach == "N_hat_flexible" ~ "Flexible shape model"
  )) %>%
  filter(!is.na(Nhat)) %>%
  mutate(bias = abs(Nhat - true_N)) %>%
  ggplot(aes(x = true_P_c, y = bias, color = Scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_smooth(se = FALSE, na.rm = TRUE) +
  scale_x_continuous(labels = rev) +
  facet_wrap(~approach, scales = "free") +
  labs(x = "true P_c", y = expression(hat(N) - N), color = NULL)
