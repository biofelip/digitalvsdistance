#Simulation with free uniform in g_y, meaning that the variance of g_y is zero. 
library(tidyverse)
library(data.table)
library(ggplot2)
nsim <- 10000
n <- 1e4
set.seed(0.101992)
sim_summary_list <- vector("list", nsim)
detected_df_list <- vector("list", nsim)
gc()
for (s in 1:nsim) {
  V <- runif(1, 0, 1)
  g0 <- runif(1, 0, 1)
  w <- runif(1, 1, 50) * 10
    r <- rgamma(2, shape = 1)
    r <- r / sum(r)
    r <- c(0, r[1], r[2])
    names(r) <- c("y", "theta", "int")

    f_theta_alfa <- runif(1, 1, 10)
    f_theta_beta<- runif(1, 1, 10)
    f_theta <- rbeta(n, shape1 = f_theta_alfa, shape2 = f_theta_beta)

    #Pit of f_theta
    u_theta = pbeta(f_theta, shape1 = f_theta_alfa, shape2 = f_theta_beta)

    g_theta_alfa <- runif(1, 1, 10)
    g_theta_beta<- runif(1, 1, 10)
    g_theta_raw <- pbeta(u_theta, shape1 = g_theta_alfa, shape2 = g_theta_beta)
    g_theta <- g_theta_raw - mean(g_theta_raw)

    # fy is already a uniform
    f_y <- c(0.01, runif(n - 1, 0, w))
    u_y = punif(f_y, 0, w)

    g_y1_raw <- rep(1, n)
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


    a_th <- sqrt(target["theta"] / var_th)
    a_int <- sqrt(target["int"] / var_int)


    g_raw <- g0 +  g_y1 + a_th * g_theta + a_int * g_y1_theta
    # rescale based on the  new assumption that g(0,theta)=1 for all theta
    # g <- g_raw / g_raw[match(0.01, f_y)]
    # rescale based on the assumption that g(y,theta)=1 at the max point of the detection curve
    # no min max sacling will be done because we dont want to assume that there are zero detection zones
    g <- g_raw/max(g_raw)
    g <- pmax(g, sd(g)/10) # the min value allowed is a tenth of the standard deviation of detection
    # save the points that corresponde to max detection
    max_point= which.max(g_raw )
    max_y= f_y[max_point]
    max_theta= f_theta[max_point]


    v_y_real <- 0
    v_th_real <- var(a_th * g_theta)
    v_int_real <- var(a_int * g_y1_theta)
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
        g_y1_alfa = NA,
        g_y1_beta = NA,

        # hyper parameters of the detection function by theta
        g_theta_alfa = g_theta_alfa,
        g_theta_beta = g_theta_beta,

        # hyper parameters of the interaction
        copula_rho = rho,

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

        # values of the distribution of the interaction
        copula_v_y = v_y,
        copula_v_theta = v_theta,

        # probability of detection
        g = g,

        # detection yes or no
        detected = detected
      )
  }

length(detected_df_list)
length(sim_summary_list)

# save the list of dataframes
saveRDS(detected_df_list, file = "L:/Documents/digitalvsdistance_results/3_detected_df_list_uniform_g_y.rds")
saveRDS(sim_summary_list, file = "L:/Documents/digitalvsdistance_results/3_sim_summary_list_uniform_g_y.rds")
# detected_df_list= readRDS("detected_df_list.rds")
# this is soooo much faster than rbind
# detected_df <- data.table::rbindlist(detected_df_list)
# dim(detected_df)

# one random example
index= sample(1:nsim, 1)
example_df= detected_df_list[[index]]
example_summary= sim_summary_list[[index]]
# plot the 3 distributions and the points on top
plot_sim_diagnostics(example_df, example_summary, y_shape = "uniform", sim_id = index)
# # save the dataframe with arrow because is too heavy right now
# arrow::write_parquet(detected_df, "detected_df.parquet")
# saveRDS(detected_df, file = "detected_df.rds")
