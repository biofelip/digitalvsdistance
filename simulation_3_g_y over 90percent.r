# Simulation script for cases in which g_y is half normal and varies abocve 0.9
library(tidyverse)
library(data.table)
nsim <- 10000
n <- 1e4
set.seed(0.101992)
detected_df_list <- vector("list", nsim) # y, theta for detected animals only
gc()
for (s in 1:nsim) {
  V <- runif(1, 0, 1)
  g0 <- runif(1, 0, 1)
  w <- runif(1, 1, 50) * 10
    r <- rgamma(3, shape = 1)
    r <- (r / sum(r))*0.1
    r <- c(0.9+r[1],r[2],r[3])
    names(r) <- c("y", "theta", "int")

    f_theta_shape <- runif(1, 1, 5)
    f_theta_scale <- runif(1, 1, 5)
    f_theta <- rgamma(n, shape = f_theta_shape, scale = f_theta_scale)

    #Pit of f_theta
    u_theta = pgamma(f_theta, shape = f_theta_shape, scale = f_theta_scale)

    g_theta_shape <- runif(1, 1, 5)
    g_theta_raw <- pgamma(u_theta, shape = g_theta_shape, rate = 1)
    g_theta <- g_theta_raw - mean(g_theta_raw)

    # fy is already a uniform
    f_y <- c(0.01, runif(n - 1, 0, w))
    u_y = punif(f_y, 0, w)

    g_y1_scale <- runif(1, w / 4, w / 2)
    g_y1_raw <- sqrt(2) / (g_y1_scale * sqrt(pi)) * exp(-u_y^2 / (2 * g_y1_scale^2))
    g_y1 <- g_y1_raw - mean(g_y1_raw)

    # the interaction will be defined with a gaussian copula between the two variables with a random
    # correlation coefficient
    rho <- runif(1, -0.9, 0.9)  # Kendall's tau range
    z_y <- qnorm(u_y)      # Convert u_y to standard normal
    z_theta <- qnorm(u_theta)  # Convert u_theta to standard normal

    Sigma <- matrix(c(1, rho, rho, 1), nrow=2)
    L <- chol(Sigma)  # Upper triangular
    z_corr <- t(L) %*% rbind(z_y, z_theta)
    z_y_corr <- z_corr[1, ]
    z_theta_corr <- z_corr[2, ]

    v_y <- pnorm(z_y_corr)
    v_theta <- pnorm(z_theta_corr)

    g_y1_theta_raw <- v_y**2 * v_theta
    g_y1_theta <- g_y1_theta_raw - mean(g_y1_theta_raw)

    var_y <- var(g_y1)
    var_th <- var(g_theta)
    var_int <- var(g_y1_theta)

    target <- V * r

    a_y <- sqrt(target["y"] / var_y)
    a_th <- sqrt(target["theta"] / var_th)
    a_int <- sqrt(target["int"] / var_int)


    g_raw <- g0 + a_y * g_y1 + a_th * g_theta + a_int * g_y1_theta
    # rescale based on the  new assumption that g(0,theta)=1 for all theta
    # g <- g_raw / g_raw[match(0.01, f_y)]
    # rescale based on the assumption that g(y,theta)=1 at the max point of the detection curve
    # no min max sacling will be done because we dont want to assume that there are zero detection zones
    g <- g_raw/max(g_raw)
    # save the points that corresponde to max detection
    max_point= which.max(g_raw )
    max_y= f_y[max_point]
    max_theta= f_theta[max_point]


    v_y_real <- var(a_y * g_y1)
    v_th_real <- var(a_th * g_theta)
    v_int_real <- var(a_int * g_y1_theta)
    v_total <- v_y_real + v_th_real + v_int_real

    detected <- rbinom(n, 1, g)

    detected_df <- data.frame(
      sim = s,
      y = f_y,
      u_y = u_y,
      max_y = max_y,
      theta = f_theta,
      u_theta = u_theta,
      max_theta = max_theta,
      copula_v_y=v_y, 
      copula_v_theta=v_theta,
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
      g_theta_shape = g_theta_shape, 
      copula_rho = rho
    )
    if (nrow(detected_df) > 0) {
      detected_df <- cbind(object = 1:nrow(detected_df), detected_df)
    }
    detected_df_list[[s]] <- detected_df
  }

length(detected_df_list)
# save the list of dataframes
saveRDS(detected_df_list, file = "detected_df_list_r_y_over90.rds")
# detected_df_list= readRDS("detected_df_list.rds")
# this is soooo much faster than rbind
detected_df <- data.table::rbindlist(detected_df_list)
dim(detected_df)


# # save the dataframe with arrow because is too heavy right now
# arrow::write_parquet(detected_df, "detected_df.parquet")
# saveRDS(detected_df, file = "detected_df.rds")
