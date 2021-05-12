# Title     : Visualisation of KS Stat for Poisson Distribution
# Objective : Implement a graphical representation to show what the Kolmogorov-Smirnov stat is and how its calculated
# Created by: Luke
# Created on: 01/04/2021

# Originally implemented when I was exploring the KS test.

if (!exists("poisson_generator", mode = "function")) source("poisson_generator.R")
library("ggplot2")

# N is relatively small to exacerbate and difference between the 2 ECDFs
N <- 150
lambda <- 120

my_pois <- poisson_generator(N, lambda)
inbuilt_pois <- rpois(N, lambda)

f.a <- ecdf(inbuilt_pois)
f.b <- ecdf(my_pois)

# x <- create a sequence between min and max of the samples with length = N
x <- seq(min(inbuilt_pois, my_pois), max(inbuilt_pois, my_pois), length.out = length(inbuilt_pois))
x0 <- mean(x[which(abs(f.a(x) - f.b(x)) == max(abs(f.a(x) - f.b(x))))]) # find x coord for greatest distance between ECDFs
y0 <- mean(f.a(x0)) # find y coordinates for line between ECDFs at our x0 coordinate
y1 <- mean(f.b(x0))

# Im sure there is a more elegant way of solving this, but this works
inbuilt_data <- data.frame(
  # Data for inbuilt function and setting colour/label for each
  x = inbuilt_pois,
  g = gl(2, N, labels = c("Theoretical Distribution", "Generated Variates"))
)

my_data <- data.frame(
  # Data for my pois generator
  xs = my_pois
)

ks_data <- data.frame(
  # Data required to plot the line between the 2 ECDFs
  xs = c(x0, x0),
  ys = c(y0, y1),
  xpair = c(x0, x0),
  ypair = c(y0, y1)
)


ggplot(inbuilt_data, aes(x, colour = "R Poisson Generator")) +
  stat_ecdf(geom = "step", alpha = 0.75) +
  stat_ecdf(data = my_data, aes(xs, colour = "My Poisson Generator"), geom = "step", alpha = 0.75) +
  geom_point(data = ks_data, aes(x = xs, y = ys, colour = "KS Stat"), shape = 10) +
  geom_line(data = ks_data, aes(x = xpair, y = ypair, colour = "KS Stat"), size = 0.5, lineend = "round", linetype = "dashed") +
  # This sets the color of the KS line and the 2 ECDFs
  scale_colour_manual(name = NULL, values = c('black', '#4c72b0', '#c44e52')) +
  labs(title = "A Visualisation of the KS Statistic", x = "X", y = "Pr(X <= x)", subtitle = "Illustrating the maximum distance between 2 ECDFs") +
  theme(text = element_text(size = 18, family = "serif"))

###############
# File Output #
###############

# FILENAME <- "KS_test.png"
# PATH <- "C:/Users/Luke/Sussex Code/Sussex R/Poisson Tests/Plot Images/PNGs/"
# ggsave(
#   FILENAME,
#   plot = last_plot(),
#   device = "png",
#   path = PATH,
#   scale = 1,
#   width = NA,
#   height = NA,
#   units = c("in", "cm", "mm"),
#   dpi = 600
# )