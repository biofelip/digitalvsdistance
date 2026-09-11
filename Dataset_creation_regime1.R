# Simulation with free varying variance contributions but a super flexible beta shape  for the
# g_y
library(tidyverse)
library(data.table)
library(glue)
library(logger)
source("utils.R")

log_appender(appender_file(".txt"))
nsim <- 10000
n <- 1e4
set.seed(0.101992)
high_r_y <- TRUE
sim_summary_list <- vector("list", nsim)
detected_df_list <- vector("list", nsim)
type <- "low_r_y"
number_file <- "1"
if (high_r_y) {
  type <- "high_r_y"
  number_file <- "1_2"
}

log_appender(appender_file(glue("log_step_1_decreasing_g_y_{type}.txt")))

gc()

for (s in 1:nsim) {
  V <- runif(1, 0, 1)
  g0 <- runif(1, 0, 1)
  w <- runif(1, 1, 50) * 10
  if (!high_r_y) {
    r <- rgamma(3, shape = 1)
    r <- r / sum(r)
    names(r) <- c("y", "theta", "int")
  } else {
    r <- rgamma(3, shape = 1)
    r <- (r / sum(r)) * 0.1
    r <- c(0.9 + r[1], r[2], r[3])
    names(r) <- c("y", "theta", "int")
  }

  f_theta_alfa <- runif(1, 1, 4)
  f_theta_beta <- runif(1, 1, 4)
  f_theta <- rbeta(n, shape1 = f_theta_alfa, shape2 = f_theta_beta)

  # Pit of f_theta
  u_theta <- pbeta(f_theta, shape1 = f_theta_alfa, shape2 = f_theta_beta)

  g_theta_alfa <- runif(1, 1, 4)
  g_theta_beta <- runif(1, 1, 4)
  g_theta_raw <- dbeta(u_theta, shape1 = g_theta_alfa, shape2 = g_theta_beta)
  g_theta <- g_theta_raw - mean(g_theta_raw)

  # fy is already a uniform
  f_y <- c(0.01, runif(n - 1, 0, w))
  u_y <- punif(f_y, 0, w)


  g_y1_scale <- runif(1, 0.4, 0.6)
  g_y1_raw <- sqrt(2) / (g_y1_scale * sqrt(pi)) * exp(-u_y^2 / (2 * g_y1_scale^2))
  g_y1 <- g_y1_raw - mean(g_y1_raw)


  ## the interaction will be defined as the product of two additional betas

  g_y_interaction_alfa <- runif(1, 1, 4)
  g_y_interaction_beta <- runif(1, 1, 4)

  g_theta_interaction_alfa <- runif(1, 1, 4)
  g_theta_interaction_beta <- runif(1, 1, 4)

  g_y_interaction_raw <- dbeta(u_y, shape1 = g_y_interaction_alfa, shape2 = g_y_interaction_beta)
  g_theta_interaction_raw <- dbeta(u_theta, shape1 = g_theta_interaction_alfa, shape2 = g_theta_interaction_beta)

  g_y_interaction <- g_y_interaction_raw - mean(g_y_interaction_raw)
  g_theta_interaction <- g_theta_interaction_raw - mean(g_theta_interaction_raw)
  # random contant multiplier
  constant <- runif(1, 0, 1)


  g_interaction <- constant * g_y_interaction * g_theta_interaction


  var_y <- var(g_y1)
  var_th <- var(g_theta)
  var_int <- var(g_interaction)

  target <- V * r

  a_y <- sqrt(target["y"] / var_y)
  a_th <- sqrt(target["theta"] / var_th)
  a_int <- sqrt(target["int"] / var_int)


  g_raw <- g0 + a_y * g_y1 + a_th * g_theta + a_int * g_interaction

  components_final <- cbind(
    y = a_y * g_y1,
    theta = a_th * g_theta,
    interaction = a_int * g_interaction
  )

  # calculate the sum of the off diagonal elements of the covaraince matrix
  # they need to be very small. If they are not throw a warning and log it in the
  # log file

  off_diag <- Map(apply(cov(components_final), 2, sum), diag(cov(components_final)), f = "-")
  if (any(off_diag > 0.01)) {
    mesage <- glue("simulation {s} the sum of the off diagonal elements of the covaraince matrix are y{off_diag[1]}, theta{off_diag[2]}, interaction{off_diag[3]}, the components might not be orthogonal ")
    log_warn(skip_formatter(mesage))
  }
  # rescale based on the  new assumption that g(0,theta)=1 for all theta
  # g <- g_raw / g_raw[match(0.01, f_y)]
  # rescale based on the assumption that g(y,theta)=1 at the max point of the detection curve
  # no min max sacling will be done because we dont want to assume that there are zero detection zones
  g <- g_raw / max(g_raw)
  g <- pmax(g, sd(g) / 10) # the min value allowed is a tenth of the standard deviation of detection
  # save the points that corresponde to max detection
  max_point <- which.max(g_raw)
  max_y <- f_y[max_point]
  max_theta <- f_theta[max_point]


  v_y_real <- var(a_y * g_y1)
  v_th_real <- var(a_th * g_theta)
  v_int_real <- var(a_int * g_interaction)
  v_total <- v_y_real + v_th_real + v_int_real

  detected <- rbinom(n, 1, g)

  sim_summary_list[[s]] <- data.table(
    # parameters of the simulation
    sim = s,
    w = w,
    g_0 = g0,

    # location of maximum detection (used to anchor the max-rescale)
    max_y = max_y,
    max_theta = max_theta,

    # hyper parameters of the distribution of theta
    f_theta_alfa = f_theta_alfa,
    f_theta_beta = f_theta_beta,

    # hyper parameters of the detection function by distance
    g_y1_scale = g_y1_scale,

    # hyper parameters of the detection function by theta
    g_theta_alfa = g_theta_alfa,
    g_theta_beta = g_theta_beta,

    # hyper parameters of the interaction
    g_y_interaction_alfa = g_y_interaction_alfa,
    g_y_interaction_beta = g_y_interaction_beta,
    g_theta_interaction_alfa = g_theta_interaction_alfa,
    g_theta_interaction_beta = g_theta_interaction_beta,
    interaction_constant = constant,

    # range of realized detection probability
    g_min = min(g),
    g_max = max(g),

    # values of the target variance and realized variance
    target_r_y = r["y"],
    target_r_theta = r["theta"],
    target_r_int = r["int"],
    realized_r_y = v_y_real / v_total,
    realized_r_theta = v_th_real / v_total,
    realized_r_int = v_int_real / v_total,
    V_realized = v_total
  )

  detected_df_list[[s]] <- data.table(
    # linking key back to sim_summary
    sim = s,
    object = 1:n,

    # values of the distribution of distances and its PIT
    y = f_y,
    u_y = u_y,

    # values of the distribution of theta and its PIT
    theta = f_theta,
    u_theta = u_theta,

    # probability of detection
    g = g,

    # detection yes or no
    detected = detected
  )
  log_info("Simulation {s} finished")
}

length(detected_df_list)
length(sim_summary_list)

# save the list of dataframes
# saveRDS(detected_df_list, file = "L:/Documents/digitalvsdistance_results/2_2_detected_df_list_high_variance_free_shape_g_y.rds")

# # one random example
# index= sample(1:nsim, 1)
# example_df= detected_df_list[[index]]
# example_summary= sim_summary_list[[index]]
# # plot the 3 distributions and the points on top
# plot_sim_diagnostics(example_df, example_summary,y_shape = "beta", sim_id = index)
# # save the dataframe with arrow because is too heavy right now
# # this is soooo much faster than rbind
detected_df_FULL <- data.table::rbindlist(detected_df_list)
dim(detected_df_FULL)


# # save the dataframe with arrow because is too heavy right now
saveRDS(sim_summary_list, file = glue("E:/digitalvsdistance_results/new_results/{number_file}_sim_summary_list_{type}_halfnormal_g_y_smooth_scale.rds"))

arrow::write_parquet(detected_df_FULL, glue("E:/digitalvsdistance_results/new_results/{number_file}_detected_df_list_{type}_halfnormal_g_y_smooth_scale.parquet"))
