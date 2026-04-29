library(ggplot2)
library(dplyr)
library(magrittr)
library(tidyr)

data_folder <- r"(C:\Users\247404\ownCloud - juan.felipe.escobar.calderon@tiho-hannover.de@sync.academiccloud.de\HABITATWal\DATEN\SURVEYS\DIGI\data with coordinates)" # nolint

all_surveys_paths <- list.files(file.path(data_folder, "bird_label_data"))

# read all the files

all_surveys <- lapply(all_surveys_paths, function(x) {
  read.csv(file.path(data_folder, "bird_label_data", x))
})

surveys_df <- bind_rows(all_surveys)


## distribution of detection of all species in all frames across the width of the pictures

surveys_df %>%
  filter(Species == "Harbour Porpoise") %>%
  ggplot(aes(x = Frame.X)) +
  geom_histogram(bins = 12, color = "black", fill = "white") +
  ggtitle("Detection across the width of the picture") +
  xlim(range(surveys_df$Frame.X, na.rm = TRUE)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

surveys_df$D_c <- abs(surveys_df$Frame.X - 1643.5)
model_D_c <- lm(density(surveys_df2$D_c)$y ~ density(surveys_df2$D_c)$x)
surveys_df2 <- surveys_df[complete.cases(surveys_df$D_c), ]
hist(surveys_df2$D_c, probability = TRUE)
lines(density(surveys_df2$D_c), col = "red", lwd = 2)

## distribution of detection of all species in all frames across the height of the pictures
surveys_df %>%
  filter(Species == "Harbour Porpoise") %>%
  ggplot(aes(y = Frame.Y)) +
  geom_histogram(bins = 10, color = "black", fill = "white") +
  ggtitle("Detection across the height of the picture") +
  ylim(range(surveys_df$Frame.Y, na.rm = TRUE)) +
  theme_minimal()

## heatmap plot of the frequency of detections across the width and height of the picture
surveys_df %>%
  filter(Species == "Harbour Porpoise") %>%
  ggplot(aes(x = Frame.X, y = Frame.Y)) +
  geom_bin2d() +
  ggtitle("2D density plot of the frequency of
    detections across the width and height of the picture")



## make a 2D histogram with uniform distribution on both x and y with 10000 points
set.seed(123) # for reproducibility
x_uniform <- runif(10000, min = , max = 3287)
y_uniform <- runif(10000, min = 0, max = 1907)

data_frame <- data.frame(x_uniform, y_uniform)

data_frame %>%
  ggplot(aes(x = x_uniform, y = y_uniform)) +
  geom_bin2d() +
  ggtitle("2D Histogram with Uniform Distribution on X and Y") +
  theme_minimal()


## relative density of the deteciton depending on the camera
surveys_df$Camera <- as.factor(surveys_df$Camera)
surveys_df$Camera <- factor(surveys_df$Camera, levels = c(1, 5, 2, 6, 3, 7, 4, 8), ordered = TRUE)
surveys_df %>%
  # filter(Species == "Harbour Porpoise") %>%
  group_by(Camera) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = Camera, y = n)) +
  geom_bar(stat = "identity") +
  ggtitle("Relative density of the deteciton depending on the camera") +
  theme_minimal()
