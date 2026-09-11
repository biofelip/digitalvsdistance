library(rstan)

list_sim_summary <- c(
    "E:/digitalvsdistance_results/new_results/1_2_sim_summary_list_high_r_y_halfnormal_g_y_smooth_scale.rds",
    "E:/digitalvsdistance_results/new_results/1_sim_summary_list_low_r_y_halfnormal_g_y_smooth_scale.rds",
    "E:/digitalvsdistance_results/new_results/2_2_sim_summary_list_high_r_y_beta_g_y_smooth_scale.rds",
    "E:/digitalvsdistance_results/new_results/2_sim_summary_list_low_r_y_beta_g_y_smooth_scale.rds",
    "E:/digitalvsdistance_results/new_results/3_sim_summary_list_high_r_y_halfnormal_g_y_smooth_scale.rds"
)

list_detected_df <- c(
    "E:/digitalvsdistance_results/new_results/1_2_detected_df_list_high_r_y_halfnormal_g_y_smooth_scale.parquet",
    "E:/digitalvsdistance_results/new_results/1_detected_df_list_low_r_y_halfnormal_g_y_smooth_scale.parquet",
    "E:/digitalvsdistance_results/new_results/2_2_detected_df_list_high_r_y_beta_g_y_smooth_scale.parquet",
    "E:/digitalvsdistance_results/new_results/2_detected_df_list_low_r_y_beta_g_y_smooth_scale.parquet",
    "E:/digitalvsdistance_results/new_results/3_detected_df_list_high_r_y_halfnormal_g_y_smooth_scale.parquet"
)


# select a random regime and nsim
i <- sample(1:5, 1)

detected_df <- arrow::read_parquet(list_detected_df[i])
detected_df <- data.table::data.table(detected_df)
sim_summary <- readRDS(list_sim_summary[i])
sim_summary <- data.table::rbindlist(sim_summary, fill = TRUE)
nsim <- sample(1:10000, 1)
# select only the specific simulation
detected_df_nsim <- detected_df[sim == nsim & detected == 1]
sim_summary_nsim <- sim_summary[sim == nsim]
# compile model


hn_model <- stan_model("hn_posterior_sampling.stan")

stan_data <- list(
    w = sim_summary_nsim$w,
    N = nrow(detected_df_nsim),
    y = detected_df_nsim$y
)

fit <- sampling(hn_model,
    data = stan_data,
    chains = 4,
    iter = 2000,
    warmup = 1000,
    refresh = 0,
    open_progress = FALSE,
    show_messages = FALSE
)


sigma_draws <- extract(fit, pars = "sigma")$sigma

# lets plot the results
posterior <- as.matrix(fit)

y_rep <- rstan::extract(fit, pars = "y_rep")$y_rep


sigma_draw_plot <- bayesplot::mcmc_areas(posterior,
    pars = c("sigma"),
    prob = 0.8
)

# posterior predictive check.
initial_draw <- sample(1:dim(y_rep)[2], 1)
ppc_plot <- bayesplot::ppc_dens_overlay(detected_df_nsim$y, y_rep[initial_draw:(initial_draw + 50), , drop = FALSE])


# plot the real distribution og g_Y of the detection function against
# the estimated distribution of  g_Y the detection function
hist(
    detected_df_nsim$y,
    breaks = 100,
    probability = TRUE
)

lines(
    density(detected_df_nsim$y),
    col = "blue",
    lwd = 2
)

curve(
    dnorm(x, mean = 0, sd = sigma_hat) /
        (pnorm(w, mean = 0, sd = sigma_hat) - pnorm(0, mean = 0, sd = sigma_hat)),
    from = 0,
    to = w,
    add = TRUE,
    col = "red",
    lwd = 2
)
legend(
    "topright",
    legend = c("Real", "Estimated"),
    col = c("blue", "red"),
    lwd = 2
)


library(ggplot2)

# Estimated truncated half-normal density
model_df <- data.frame(
    y = seq(0, w, length.out = 1000)
)

model_df$density_estimated <-
    dnorm(model_df$y, mean = 0, sd = sigma_hat) /
        (
            pnorm(w, mean = 0, sd = sigma_hat) -
                pnorm(0, mean = 0, sd = sigma_hat)
        )

# Choose a scale factor to display g on the density plot
density_scale <- max(model_df$density_estimated) / 1

ggplot(detected_df_nsim, aes(x = y)) +

    # Histogram on density scale
    geom_histogram(
        aes(y = after_stat(density)),
        bins = 100,
        fill = "grey85",
        colour = "black"
    ) +

    # Empirical density
    geom_density(
        aes(colour = "Empirical density"),
        linewidth = 1.2,
        adjust = 1
    ) +

    # Real g-values, rescaled to the density axis
    geom_smooth(
        data = detected_df_nsim,
        aes(
            x = y,
            y = g * density_scale,
            colour = "Real g(y)"
        ),
        method = "loess",
        se = FALSE,
        span = 0.75,
        linewidth = 1.1
    ) +

    # Estimated truncated half-normal density
    geom_line(
        data = model_df,
        aes(
            x = y,
            y = density_estimated,
            colour = "Estimated density"
        ),
        linewidth = 1.2
    ) +
    scale_y_continuous(
        name = "Distance density",
        sec.axis = sec_axis(
            ~ . / density_scale,
            name = "Detection probability g(y)"
        )
    ) +
    scale_colour_manual(
        values = c(
            "Real g(y)" = "black",
            "Estimated density" = "red",
            "Empirical density" = "blue"
        )
    ) +
    labs(
        x = "Distance",
        colour = NULL
    ) +
    theme_minimal()

library(ggplot2)

# Distance grid
g_plot_df <- data.frame(
    y = seq(0, w, length.out = 1000)
)

# Estimated half-normal detection probability
g_plot_df$g_estimated <- exp(
    -g_plot_df$y^2 / (2 * sigma_hat^2)
)


ggplot(detected_df_nsim, aes(x = y, y = g)) +

    # Stored real detection probabilities
    geom_point(
        aes(colour = "Real g(y)"),
        alpha = 0.25,
        size = 1
    ) +

    # Smoothed representation of the stored real g-values
    geom_smooth(
        aes(
            colour = "Real g(y)"
        ),
        method = "loess",
        formula = y ~ x,
        se = FALSE,
        span = 0.75,
        linewidth = 1.1
    ) +

    # Estimated half-normal detection function
    geom_line(
        data = g_plot_df,
        aes(
            x = y,
            y = g_estimated,
            colour = "Estimated g(y)"
        ),
        linewidth = 1.2
    ) +
    scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2)
    ) +
    scale_colour_manual(
        values = c(
            "Real g(y)" = "blue",
            "Estimated g(y)" = "red"
        )
    ) +
    labs(
        x = "Distance",
        y = "Probability of detection",
        colour = NULL
    ) +
    theme_minimal()


# if if do it grid search of the sigma in the close form posterior might be a lot faster
# according to thehta wise the lotg posterior of a truncated half normal with parameter sigma
# and know mean with a normal prior on sigma is:

log_posterior_hn <- function(sigma, y, w) {
    N <- length(y)
    p <- -N * log(sigma) - ((sum(y^2)) / (2 * sigma^2)) - N * log(pnorm(w / sigma) - (1 / 2)) - ((2 * sigma^2) / (w^2))
    return(p)
}

sigma_grid <- seq(1, 500, 0.1)

log_posterior_hn_grid <- sapply(sigma_grid, function(x) log_posterior_hn(sigma = x, y = detected_df_nsim$y, w = stan_data$w))

posterior_unnormalized <- exp(
    log_posterior_hn_grid - max(log_posterior_hn_grid)
)


# cdf by trapezoidal rule
dx <- diff(sigma_grid)


cdf <- c(
    0,
    cumsum(
        0.5 *
            (posterior_unnormalized[-1] +
                posterior_unnormalized[-length(posterior_unnormalized)]) *
            dx
    )
)

# Normalize the CDF
cdf <- cdf / max(cdf)
approx(
    x = cdf,
    y = sigma_grid,
    xout = c(0.1, 0.9)
)$y

plot(sigma_grid, posterior_unnormalized, xlim = c(80, 120))
abline(v = sigma_grid[which.max(posterior_unnormalized)], col = "red")
abline(v = c(89.37700, 93.98595), col = "blue")


# actual speed test

library(profvis)
# profiling with compilation time (not fair)
profvis({
    stan_data <- list(
        w = sim_summary_nsim$w,
        N = nrow(detected_df_nsim),
        y = detected_df_nsim$y
    )

    fit <- sampling(hn_model,
        data = stan_data,
        chains = 4,
        iter = 2000,
        warmup = 1000,
        refresh = 0,
        open_progress = FALSE,
        show_messages = FALSE
    )


    sigma_draws <- extract(fit, pars = "sigma")$sigma
})

profvis({
    sigma_grid <- seq(1, 500, 0.01)

    log_posterior_hn_grid <- sapply(sigma_grid, function(x) log_posterior_hn(sigma = x, y = detected_df_nsim$y, w = stan_data$w))

    posterior_unnormalized <- exp(
        log_posterior_hn_grid - max(log_posterior_hn_grid)
    )


    # cdf by trapezoidal rule
    dx <- diff(sigma_grid)


    cdf <- c(
        0,
        cumsum(
            0.5 *
                (posterior_unnormalized[-1] +
                    posterior_unnormalized[-length(posterior_unnormalized)]) *
                dx
        )
    )

    # Normalize the CDF
    cdf <- cdf / max(cdf)
    approx(
        x = cdf,
        y = sigma_grid,
        xout = c(0.1, 0.9)
    )$y
})
