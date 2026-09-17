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

list_results = c("E:/digitalvsdistance_results/new_results/results_pc_1_2_sim_summary_list_high_r_y_halfnormal_g_y_smooth_scale3.csv",
"E:/digitalvsdistance_results/new_results/results_pc_1_sim_summary_list_low_r_y_halfnormal_g_y_smooth_scale3.csv",
"E:/digitalvsdistance_results/new_results/results_pc_2_2_sim_summary_list_high_r_y_beta_g_y_smooth_scale3.csv",
"E:/digitalvsdistance_results/new_results/results_pc_2_sim_summary_list_low_r_y_beta_g_y_smooth_scale3.csv",
"E:/digitalvsdistance_results/new_results/results_pc_3_sim_summary_list_high_r_y_halfnormal_g_y_smooth_scale3.csv")

results_l = list()

for (i in 1:5) {
  sim_summary <- data.table::rbindlist(readRDS(list_sim_summary[i]))
  sim_summary <- data.table::data.table(sim_summary)
  results <- read.csv(list_results[i])
  results <- data.table::data.table(results)
  results_l[[i]] <- results[sim_summary, on = "sim"]
  }

all_res_joined <- data.table::rbindlist(results_l, fill = TRUE)
dim(all_res_joined)
colnames(all_res_joined)
all_res_
# Erase all columns that are just NAs
all_res_joined %>% 
  select_if(~ !all(is.na(.))) -> all_res_joined

# select the columns that start with P_c and the realized r_y, r_theta, and r_int
all_res_joined %>% 
  select(starts_with("P_c"), starts_with("realized_r")) -> all_res_joined

# cut the r_y column into bins
all_res_joined %>% 
  mutate(r_y_bin = cut(realized_r_y, breaks = seq(0, 1, 0.1), include.lowest = TRUE)) -> all_res_joined

# plot the number of NAs per P_c column, per r_y bin
all_res_joined %>% 
  select(starts_with("P_c"), r_y_bin) %>% 
  group_by(r_y_bin, Scenario) %>%
  summarise(across(starts_with("P_c"), ~ sum(is.na(.)))) %>% 
  pivot_longer(-r_y_bin) %>% 
  ggplot(aes(x =as.factor(r_y_bin) , y = value, fill = name)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal()


for (i in 1:5) {
  r <- read.csv(list_results[i])
  s <- data.table::rbindlist(readRDS(list_sim_summary[i]))
  cat(list_results[i], ": results has", nrow(r), "/ sim_summary has", nrow(s), "\n")
}

all_res_joined %>%
  filter(r_y_bin == "(0.9,1]") %>%
  summarise(n = n(),
            na_uniform = sum(is.na(P_c_uniform)),
            na_plain   = sum(is.na(P_c_plain)),
            na_normal_theta = sum(is.na(P_c_normal_theta)))