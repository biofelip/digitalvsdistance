# library(mrds)

# results_pc=data.frame(sim = numeric(),
#      P_c=numeric(),
#      P_c_cov=numeric())
# for(s_pick in 1:nsim) {
#     fit_data <- detected_df %>%
#   filter(detected == 1, sim == s_pick) %>%   # one sim at a time
#   transmute(
#     object   = object,
#     distance = y,
#     theta    = theta,     # your covariate, if you want it in the model
#     size     = 1          # single-animal detections; required column name
#   )
# # plain half-normal, no covariates -- recovers a single pooled scale
# fit_hn <- ddf(dsmodel = ~mcds(key = "hn", formula = ~1),
#               data = fit_data, method = "ds",
#               meta.data = list(width = w))



# # half-normal with theta on the scale parameter -- MCDS covariate model
# fit_hn_cov <- ddf(dsmodel = ~mcds(key = "hn", formula = ~theta),
#                    data = fit_data, method = "ds",
#                    meta.data = list(width = w))

# results_pc <- rbind(results_pc,
#                     data.frame(sim = s_pick,
#                                P_c = summary(fit_hn)$average.p,
#                                P_c_cov = mean(predict(fit_hn_cov, integrate=TRUE)$fitted)))
# }

results_pc <- data.frame(sim = numeric(),
                          P_c = numeric(),
                          P_c_cov = numeric(),
                          N_hat = numeric(),
                          N_hat_cov = numeric())

for (s_pick in 1:nsim) {

  fit_data <- detected_df %>%
    filter(detected == 1, sim == s_pick) %>%
    transmute(object = object, distance = y, theta = theta, size = 1)

  n_detected  <- nrow(fit_data)
  N_hat_strip <- n_detected                      # <-- add this, was missing

  # --- plain half-normal ---
  fit_hn <- tryCatch(
    ddf(dsmodel = ~mcds(key = "hn", formula = ~1),
        data = fit_data, method = "ds", meta.data = list(width = w)),
    error = function(e) NULL
  )
  P_c   <- tryCatch(if (!is.null(fit_hn)) summary(fit_hn)$average.p else NA,
                     error = function(e) NA)
  if (length(P_c) == 0) P_c <- NA
  N_hat <- if (!is.na(P_c)) n_detected / P_c else NA

  # --- half-normal with theta covariate ---
  fit_hn_cov <- tryCatch(
    ddf(dsmodel = ~mcds(key = "hn", formula = ~theta),
        data = fit_data, method = "ds", meta.data = list(width = w)),
    error = function(e) NULL
  )
  if (!is.null(fit_hn_cov)) {
    p_hat_cov <- tryCatch(predict(fit_hn_cov, integrate = TRUE)$fitted, error = function(e) NULL)
    P_c_cov   <- if (!is.null(p_hat_cov)) mean(p_hat_cov) else NA
    N_hat_cov <- if (!is.null(p_hat_cov)) sum(1 / p_hat_cov) else NA
  } else {
    P_c_cov   <- NA
    N_hat_cov <- NA
  }

  # --- uniform key ---
  fit_unif <- tryCatch(
    ddf(dsmodel = ~mcds(key = "unif", formula = ~1),
        data = fit_data, method = "ds", meta.data = list(width = w)),
    error = function(e) NULL
  )
  P_c_unif <- tryCatch(if (!is.null(fit_unif)) summary(fit_unif)$average.p else NA,
                        error = function(e) NA)
  if (length(P_c_unif) == 0) P_c_unif <- NA
  N_hat_unif <- if (!is.na(P_c_unif)) n_detected / P_c_unif else NA

  results_pc <- rbind(results_pc,
                       data.frame(sim = s_pick,
                                  P_c = P_c, P_c_cov = P_c_cov,        # <-- direct, no cov_result$
                                  N_hat = N_hat, N_hat_cov = N_hat_cov,
                                  N_hat_strip = N_hat_strip,
                                  P_c_unif = P_c_unif, N_hat_unif = N_hat_unif))
}
#save results
saveRDS(results_pc, file = "results_pc.rds")

colMeans(results_pc, na.rm = TRUE)
