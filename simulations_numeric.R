library(tidyverse)
library(data.table)
library(Ternary)
nsim <- 10000
n <- 1e4
detected_df_list <- vector("list", nsim) # y, theta for detected animals only
set.seed(0.101992)
for (s in 1:100) {
  V <- runif(1, 0, 1) # Variance of the detection function
  g0 <- runif(1, 0, 1) # Mean detection across al y and theta
  w <- runif(1, 1, 50) * 10 # transect width

  # Proportions of variance of each component of the detection function
  r <- rgamma(3, shape = 1)
  r <- r / sum(r)
  names(r) <- c("y", "theta", "int")

  # Paramaeters of the distribution of the covariate theta (pi(theta))
  # a random shape and scale parameter
  f_theta_shape <- runif(1, 1, 5)
  f_theta_scale <- runif(1, 1, 5)
  # draw from a random gamma distribution for theta, i.e the values of theta for all the animals
  # inside the sampling area.

  f_theta <- rgamma(n, shape = f_theta_shape, scale = f_theta_scale)

  # Parameters of the effect of theta on the detection function g(y,theta)
  g_theta_shape <- runif(1, 1, 5)
  g_theta_raw <- pgamma(f_theta, shape = g_theta_shape, rate = 1)
  # center the g_theta values to have mean 0 so its a valid sobol  component of the detection function g(y,theta)
  g_theta <- g_theta_raw - mean(g_theta_raw)


  # the distribution of distances of detected and undetected animals is uniform across the transect width,
  # so we can draw from a uniform distribution for the distances of all the animals inside the sampling area.
  f_y <- c(0.01, runif(n - 1, 0, w))

  # Parameters of the effect of y on the detection function g(y,theta)
  g_y1_scale <- runif(1, w / 4, w / 2) # scale of the half normal.
  g_y1_raw <- sqrt(2) / (g_y1_scale * sqrt(pi)) * exp(-f_y^2 / (2 * g_y1_scale^2))
  # again center so is a valid sobol component
  g_y1 <- g_y1_raw - mean(g_y1_raw)

  # Due to the fact that we need the detection function to be =1 at y=0 for all theta
  # we need to impose that g_y1(0) = 1- g_0, and at the same time that
  #

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

  g_raw <- g0 + a_y * g_y1 + a_th * g_theta + a_int * g_y1_theta
  g_raw <- pmin(pmax(g_raw, 0), 1)
  g_y1_at_0 <- g_y1[1]
  shoulder <- g0 + a_y * g_y1_at_0 + a_th * g_theta + a_int * g_y1_at_0 * g_theta
  g <- g_raw / shoulder
  # rescale based on the assumption that g(0,theta)=1 for all theta
  # g <- g_raw / g_raw[match(0.01, f_y)]
  # rescale based on the assumption that g(y,theta)=1 at the max point of the detection curve
  # g <- g_raw/max(g_raw)

  ### plot 3d
  # data.frame(f_y = f_y, f_theta = f_theta, g = g) %>%
  #   ggplot(aes(x = f_y, y = f_theta, z = g)) +
  #   stat_summary_hex() +
  #   labs(x = "y", y = "theta", fill = "g(y,theta)") +
  #   theme_minimal()
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


write.csv(detected_df, file = "results_sim_df.csv", row.names = FALSE)
