# simulate two scenario one in which the detection function is decreasing with g(0) = 1
# and one in which the detection function is uniform with respect to distance and g(0) = 1
# distance in independent from the other variables
# for this example theta is gonna be 4 variables
# Size of the animal
# Visibility
# distance from center
# correct training.
# I will assume for simplicity that the three of them are independent
library(ggplot2)
library(magrittr)
library(dplyr)
library(tidyr)
library(scales)
N <- 15000
set.seed(1)
n <- 10000


safe_rescale <- function(x, to = c(0, 1)) {
  # Remove NAs and infinite values to compute range
  x_finite <- x[is.finite(x)]

  # If no finite values, return NAs
  if (length(x_finite) == 0) {
    return(rep(NA_real_, length(x)))
  }

  # Compute min and max safely
  x_min <- min(x_finite)
  x_max <- max(x_finite)

  # If all values are equal, return constant vector at upper bound
  if (x_max == x_min) {
    return(rep(to[2], length(x)))
  }

  # Rescale safely, keeping original NAs
  scaled <- (x - x_min) / (x_max - x_min) * (to[2] - to[1]) + to[1]
  scaled[!is.finite(x)] <- NA_real_ # optional: keep NA/Inf as NA

  return(scaled)
}


# Prbability functions oof the different variables
# f_D_hn <- function(d, sigma) rescale(exp(-(d^2) / (2 * sigma^2)), to = c(0, 1)) # Half-normal: p = 1 at d=0, decreases
f_D_hn <- function(d, alfa, beta) safe_rescale(dbeta(rescale(d, to = c(0, 1)), alfa, beta), to = c(0, 1)) # just another function that can vary in many ways
# 1. Size (S): Use Normal CDF, centered lower to increase and plateau
f_S <- function(s, alfa, beta) safe_rescale(dbeta(rescale(S, to = c(0, 1)), alfa, beta), to = c(0, 1))
# 2. Visibility (V): Reverse Logistic CDF for decreasing shape with early plateau
f_V <- function(v, alfa, beta) safe_rescale(dbeta(rescale(v, to = c(0, 1)), alfa, beta), to = c(0, 1))
# 3. Distance to Center (D_c): Increasing with log-saturation, or just scaled CDF
f_D_c <- function(c, alfa, beta) safe_rescale(dbeta(rescale(c, to = c(0, 1)), alfa, beta), to = c(0, 1)) # Smooth increase, plateaus near 1
# Training effect: convert factor to numeric multiplier
f_Tr <- function(tr, alfa, beta) safe_rescale(dbeta(rescale(tr, to = c(0, 1)), alfa, beta), to = c(0, 1))


# vectors to store information
sigma_v <- numeric(N)
alpha_v <- beta_v <- matrix(NA, nrow = N, 5)
Pa_of_g_d <- numeric(N)
Pa_of_g_d_theta <- numeric(N)

D <- runif(n, 0, 100) # distribution of distances of animals detected and undetected
S <- rnorm(n, 10, 4) # size in pixels
V <- runif(n, 0, 1) # visibility simpe 0 to 1
D_c <- runif(n, 0, 25) # distance to center of the picture
Tr <- runif(n, 0, 1) # like bvisibility

## Scenario 1: The distance has a half normal with varying parameter sigma and all the others
## are betas with many different shape paremeter (included 1 and 1 so they are uniform)
# Define nonlinear functions for each covariate

for (i in 1:N) {
  ## defininf the parameters to iterate along
  sigma_v[i] <- sigmas <- sample(seq(0.1, 17, 0.5), 1)
  alpha_v[i, ] <- alphas <- sample(seq(0.1, 10, by = 0.1), 5) # runif(5, 0, 10)
  beta_v[i, ] <- betas <- sample(seq(0.1, 10, by = 0.1), 5) # runif(5, 0, 10)
  # covariate values

  # calculate the probabilities and scale them to 1 (TODO add a way in which we can change the range of the probs)
  # fs_D_rs <- rescale(f_D_hn(D, sigma = sigmas), to = c(0, 1))
  # fs_rs <- rescale(f_S(S, alfa = alphas[1], beta = betas[1]), to = c(0, 1))
  # fv_rs <- rescale(f_V(V, alfa = alphas[2], beta = betas[2]), to = c(0, 1))
  # fdc_rs <- rescale(f_D_c(D_c, alfa = alphas[3], beta = betas[3]), to = c(0, 1))
  # ftr_rs <- rescale(f_Tr(Tr, alfa = alphas[4], beta = betas[4]), to = c(0, 1))
  # fs_D_rs <- f_D_hn(D, alfa = 1, 1)
  fs_D_rs <- f_D_hn(D, alfa = alphas[1], betas[1])
  fs_rs <- f_S(S, alfa = alphas[2], beta = betas[2])
  fv_rs <- f_V(V, alfa = alphas[3], beta = betas[3])
  fdc_rs <- f_D_c(D_c, alfa = alphas[4], beta = betas[4])
  ftr_rs <- f_Tr(Tr, alfa = alphas[5], beta = betas[5])
  # Combine: assume independence and multiply for the detection function
  # p <- f_D(D) * fs_rs * fv_rs * fdc_rs * ftr_rs
  # Combine: assume the true proabbility of detection is a unweighted linear combination of the covariates
  p <- rowMeans(matrix(c(fs_D_rs, fs_rs, fv_rs, fdc_rs, ftr_rs), nrow = n, byrow = FALSE), na.rm = TRUE)
  Pa_of_g_d_theta[i] <- mean(p, na.rm = TRUE)
  Pa_of_g_d[i] <- mean(fs_D_rs, na.rm = TRUE)
}
results <- cbind(sigma_v, alpha_v, beta_v, Pa_of_g_d, Pa_of_g_d_theta)
colnames(results) <- c(
  "sigma_v", "alpha_D", "alpha_S", "alpha_V", "alpha_CD", "alpha_Tr",
  "beta_D", "beta_S", "beta_V", "beta_CD", "beta_Tr", "Pa_of_g_d",
  "Pa_of_g_d_theta"
)

N <- 1000

# vectors to store information
sigma_v <- numeric(N)
alpha_v <- beta_v <- matrix(NA, nrow = N, 5)
Pa_of_g_d <- numeric(N)
Pa_of_g_d_theta <- numeric(N)
diff_Pa <- Pa_of_g_d - Pa_of_g_d_theta
for (i in 1:N) {
  ## defininf the parameters to iterate along
  sigma_v[i] <- sigmas <- sample(seq(0.1, 17, 0.5), 1)
  alpha_v[i, ] <- alphas <- sample(seq(0.1, 10, by = 0.1), 5) # runif(5, 0, 10)
  beta_v[i, ] <- betas <- sample(seq(0.1, 10, by = 0.1), 5) # runif(5, 0, 10)
  # covariate values

  # calculate the probabilities and scale them to 1 (TODO add a way in which we can change the range of the probs)
  # fs_D_rs <- rescale(f_D_hn(D, sigma = sigmas), to = c(0, 1))
  # fs_rs <- rescale(f_S(S, alfa = alphas[1], beta = betas[1]), to = c(0, 1))
  # fv_rs <- rescale(f_V(V, alfa = alphas[2], beta = betas[2]), to = c(0, 1))
  # fdc_rs <- rescale(f_D_c(D_c, alfa = alphas[3], beta = betas[3]), to = c(0, 1))
  # ftr_rs <- rescale(f_Tr(Tr, alfa = alphas[4], beta = betas[4]), to = c(0, 1))
  fs_D_rs <- f_D_hn(D, alfa = 1, 1)
  # fs_D_rs <- f_D_hn(D, alfa = alphas[1], betas[1])
  fs_rs <- f_S(S, alfa = alphas[2], beta = betas[2])
  fv_rs <- f_V(V, alfa = alphas[3], beta = betas[3])
  fdc_rs <- f_D_c(D_c, alfa = alphas[4], beta = betas[4])
  ftr_rs <- f_Tr(Tr, alfa = alphas[5], beta = betas[5])
  # Combine: assume independence and multiply                                                                                                          for the detection function
  # p <- f_D(D) * fs_rs * fv_rs * fdc_rs * ftr_rs
  # Combine: assume the true proabbility of detection is a unweighted linear combination of the covariates
  p <- rowMeans(matrix(c(fs_D_rs, fs_rs, fv_rs, fdc_rs, ftr_rs), nrow = n, byrow = FALSE), na.rm = TRUE)
  Pa_of_g_d_theta[i] <- mean(p, na.rm = TRUE)
  Pa_of_g_d[i] <- mean(fs_D_rs, na.rm = TRUE)
}

results2 <- cbind(sigma_v, alpha_v, beta_v, Pa_of_g_d, Pa_of_g_d_theta)

# names if beta
colnames(results2) <- c(
  "sigma_v", "alpha_D", "alpha_S", "alpha_V", "alpha_CD", "alpha_Tr",
  "beta_D", "beta_S", "beta_V", "beta_CD", "beta_Tr", "Pa_of_g_d",
  "Pa_of_g_d_theta"
)
results2[, "alpha_D"] <- 1
results2[, "beta_D"] <- 1

results_f <- rbind(results, results2)
write.csv(results_f, "results_distance_final.csv")



## load the results for both simulations
setwd(r"(E:\bird_label_data\ECS_DS)")
results_hn <- read.csv("results_distance_hn2.csv")
results_beta <- read.csv("results_distance_final.csv")
# results_beta = results2

## Usando solo las betas porque la combinacion de betas
# aunque sean distintat tiene menos varianza que la combinacion de
# betas y otra joda mira voz

results_beta <- results_beta %>%
  as.data.frame() %>%
  mutate(
    shape = case_when(
      alpha_D < 1 & beta_D > 1 ~ "Decreasing (CDS)",
      alpha_D == 1 & beta_D == 1 ~ "Uniform",
      alpha_D < 1 & beta_D < 1 ~ "U-shaped",
      alpha_D > 1 & beta_D < 1 ~ "Increasing",
      alpha_D > 1 & beta_D > 1 ~ "Bell-Shaped",
      alpha_D > 1 & beta_D == 1 ~ "Increasing",
      alpha_D == 1 & beta_D > 1 ~ "Decreasing (CDS)",
    )
  )
table(results_beta$shape)
results_beta$res <- results_beta$Pa_of_g_d_theta - results_beta$Pa_of_g_d
results_beta$MSE <- results_beta$res**2

plot_sims <- results_beta %>%
  filter(shape %in% c("Decreasing (CDS)", "Uniform")) %>%
  ggplot(aes(x = MSE, fill = shape)) +
  geom_density(bw = 0.1, alpha = 0.5) +
  xlim(0, 0.7) +
  geom_vline(xintercept = 0, color = "red") +
  scale_fill_brewer(palette = "Blues") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(
  filename = "plot_sims2.svg", plot_sims,
  units = "cm", width = 15, height = 13, device = svglite::svglite
)


results_beta %>%
  filter(shape != "NA") %>%
  ggplot(aes(x = MSE, y = shape, color = shape)) +
  geom_density_ridges(fill = "white") +
  theme_minimal() +
  theme(legend.position = "none")
results_beta %>%
  group_by(shape) %>%
  summarise(
    mean(Pa_of_g_d_theta),
    mean(Pa_of_g_d),
    mean(MSE)
  ) %>%
  View()


## el resultado parece indicar que lo unico que importa es que
## la distancia tenga una forma lo peor es no tener forma
## siempre ser una gonorrea sin princpios.

# with(results_hn_grouped,
# {plot(Pa_of_g_d_theta, Pa_of_g_d_theta - Pa_of_g_d,
#      xlab = "Actual", ylab = "Residual (Actual - Approx)",
#      main = "Residual Plot")
# abline(h = 0, col = "red", lty = 2)})
#
# with(results_hn_grouped_alpha_beta,
#      {plot(Pa_of_g_d_theta, Pa_of_g_d_theta - Pa_of_g_d,
#            xlab = "Actual", ylab = "Residual (Actual - Approx)",
#            main = "Residual Plot beta")
#        abline(h = 0, col = "red", lty = 2)})


# plotear los datos crudos


## group by sigma in results 1 and by beta and alfa of

# set.seed(1)
# n <- 100000
# library(scales)

# # Covariates
# D <- runif(n, 0, 100) # distribution of distances of animals detected and undetected
# S <- rnorm(n, 10, 4) # size in pixels
# V <- runif(n, 0, 1) # visibility simpe 0 to 1
# D_c <- runif(n, 0, 25) # distance to center of the picture
# Tr <- runif(n, 0,1) # like bvisibility


# # Define nonlinear functions for each covariate
# f_D <- function(d) exp(-(d^2) / (2 * 30^2)) # Half-normal: p = 1 at d=0, decreases

# # 1. Size (S): Use Normal CDF, centered lower to increase and plateau
# f_S <- function(s) pnorm(s, mean = 20, sd = 5)

# # 2. Visibility (V): Reverse Logistic CDF for decreasing shape with early plateau
# f_V <- function(v) 1 - plogis((v - 0.5) * 10)

# # 3. Distance to Center (D_c): Increasing with log-saturation, or just scaled CDF
# f_D_c <- function(c) 1-plogis((c - 25) / 5)  # Smooth increase, plateaus near 1

# # Training effect: convert factor to numeric multiplier
# f_Tr <-  function(v) plogis((v - 0.7) * 10)

# # plot all the variables together in a nice plot
# fs_rs <- rescale(f_S(S), to=c(0, 1))
# fv_rs <- rescale(f_V(V), to = c(0,1))
# fdc_rs <- rescale(f_D_c(D_c), to = c(0,1))
# ftr_rs <-  rescale(f_Tr(Tr), to=c(0,1))

# df_y_theta <- data.frame(
#     "Variable" = rep(c("Distance", "Size", "Visibility", "Dist_center", "Training"), each = n),
#     "Value" = c(D, S, V, D_c, Tr),
#     "Probs" = c(f_D(D),
#                 fs_rs,
#                 fv_rs,
#                 fdc_rs,
#                 ftr_rs))
# df_y_theta$Variable = factor(df_y_theta$Variable, levels = c("Distance", "Size", "Visibility", "Dist_center", "Training"), ordered = TRUE)
# levels(df_y_theta$Variable)=c("Distance", "Size", "Visibility", "Dist_center", "Training")
# #densities
# pi_D <- dunif(D, 0, 100)
# pi_S <- dnorm(S, mean = 10, sd = 4)
# pi_V <- dunif(V, 0, 1)
# pi_D_c <- dunif(D_c, 0, 25)
# pi_Tr <- dunif(Tr, 0, 1)
# library(ggplot2)
# library(magrittr)
# library(dplyr)
# library(tidyr)

# df_y_theta[sample(n*5, 1000),]%>% ggplot(aes(x = Value, y = Probs, fill = Variable)) +
#     geom_jitter() +
#     geom_smooth() +
#     facet_wrap(~Variable, scales = "free_x")+theme_minimal()+
#   theme(legend.position = "none")
# # Combine: assume independence and multiply for the detection function
# p <- f_D(D) * fs_rs * fv_rs * fdc_rs * ftr_rs

# df_total <- data.frame(
#   "Distance"=D,"Size"=S, "Visibility"=V, "Dist_Center"=D_c,"Training"=Tr,
#   "detection_p"=p, "Prob_D"= f_D(D), "Prob_S"=fs_rs,
#   "Prov_V"=fv_rs,"Prob_D_c"=fdc_rs,"Prob_Tr"=ftr_rs
# )
# p2 = rowMeans(select(df_total,c(Prob_D, Prob_S, Prov_V,Prob_D_c, Prob_Tr)))

# ## esto parece estar causado solo por la fuerza de la relacion que
# # modele porque Tr tambien muestra una relacion muy visible
# #mientras que los otros se ven mal
# df_total[sample(n,10000),] %>% ggplot(aes(x=Distance, y=detection_p)) + geom_point()
# df_total[sample(n,10000),] %>% ggplot(aes(x=Size, y=detection_p)) +  geom_point()
# df_total[sample(n,10000),] %>% ggplot(aes(x=Visibility, y=detection_p)) +  geom_point()
# df_total[sample(n,10000),] %>% ggplot(aes(x=Dist_Center, y=detection_p))+ geom_point()
# df_total[sample(n,10000),] %>% ggplot(aes(x=Training, y=detection_p))+ geom_point()

# numerator <-p2 * pi_D *  pi_S * pi_V * pi_D_c * pi_Tr
# f <-  numerator / sum(numerator)
# hist(f)
# u <- runif(n)
# detected <-  u < p2
# plot(D[detected], f[detected], xlab = "Distance", ylab = "Normalized f(D, ...)", main = "Conditional density for detected animals", pch = 16, col = rgb(0,0,1,0.2))
# lines(lowess(D[detected], f[detected]), col = "red", lwd = 2)


# # De nuevo  como en el caso anterior aca  el efecto de la covariables es visible
# # en todos los casos porque las puse con variacion tan fuerte
# # una presupuesto que puedo discutir, si disminyo el rango de vraicion de la
# # probailidad el efecto de la mayria de variables se diluye
# #ESTO PODRIA USARSRE COMO UN ARGUMENTO A FAVOR YA QUE LA UNICA VARIABLE
# # QUE SIEMPRE TIENE QUE SER DECRECIENTE Y TENER PICO EN 1 ES DISTANCIA
# hist(D[detected])
# hist(S[detected])
# hist(V[detected])
# hist(D_c[detected])
# hist(Tr[detected])




# mean(p)
# mean(f_D(D))
# plot(density(rescale(p, to=c(0,1))))
