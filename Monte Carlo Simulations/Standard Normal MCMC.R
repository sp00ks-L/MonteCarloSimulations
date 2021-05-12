# Title     : Standard Normal MCMC Sampler
# Objective : Given the sum of N Poisson r.v.s (Z_n), as N -> Inf, the r.v Z_n will approach a normal distribution
#             Here I will plot that along side the standard normal distribution
# Created by: Luke
# Created on: 02/05/2021

# You will have to change the path to correctly look at my MCMC sampler
if (!exists("poisson_MCMC_sampler", mode = "function")) source("Sussex R/Final Code/Submission Code/Poisson MH_MCMC Sampler.R")
library("ggplot2")

'THIS SCRIPT DOES TAKE A LONG TIME TO EXECUTE WHEN THE VALUE OF N IS HIGH'


N <- 30000
lambda <- 100
# For my final plot I set num_of_variables to 500 - this takes a considerable amount of time
num_of_variables <- 50
new_lambda <- lambda * num_of_variables

pois_pmf <- function(x, lambda) {
  return(exp(x * log(lambda) - lambda - lgamma(x + 1)))
}

X <- rep(0, N)
for (i in 1:num_of_variables) {
  X <- X + poisson_MCMC_sampler(N, lambda, seed = lambda, 0, pois_pmf)
}


# Use this to extract the fd breaks
h <- hist(X, breaks = "fd", plot = FALSE)
sim_data <- data.frame(my_pois = X)

# Standard Normal Calculations
get_std_norm <- function(x) {
  UN <- rep(0, N)
  EX <- mean(X)
  VARX <- var(X)
  for (i in 1:N) {
    top <- X[i] - EX
    bottom <- sqrt(VARX)
    UN[i] <- (top / bottom)
  }
  return(UN)
}

UN <- get_std_norm(X)
# Use this to extract the fd breaks
h2 <- hist(UN, "fd", plot = FALSE)
un_data <- data.frame(my_norm = UN)

title <- paste("Normalised Poisson Distribution : lambda = ", lambda * num_of_variables)
ggplot(un_data, aes(x = my_norm)) +
  geom_histogram(aes(y = ..density.., color = "Poisson MCMC"), breaks = h2$breaks, fill = "steelblue", size = 0.2, alpha = 0.5) +
  scale_colour_manual(name = NULL, values = c("black", "red")) +
  labs(title = title, x = "X Normalised", y = "Density") +
  scale_x_continuous(breaks = seq(-3, 3, 1)) +
  theme(text = element_text(size = 18, family = "serif"))

title <- paste("Poisson Distribution : lambda = ", lambda * num_of_variables)
ggplot(sim_data, aes(x = my_pois)) +
  geom_histogram(aes(y = ..density.., color = "Poisson MCMC"), breaks = h$breaks, fill = "steelblue", size = 0.2, alpha = 0.5) +
  scale_colour_manual(name = NULL, values = c("black", "red")) +
  labs(title = title, x = "X", y = "Density") +
  scale_x_continuous(breaks = seq(49000, 50800, 400)) +
  theme(text = element_text(size = 18, family = "serif"))

# write.csv(data.frame(data = X), "C:\\Users\\Luke\\Sussex Code\\Sussex R\\Final Code\\Poisson\\NSUMresults.csv")


###############
# File Output #
###############
# FILENAME <- "Std Norm Plot.png"
# PATH <- "C:/Users/Luke/Sussex Code/Sussex R/Poisson Tests/Plot Images/PNGs/"
# ggsave(
#   FILENAME,
#   plot = std_norm_plot,
#   device = "png",
#   path = PATH,
#   scale = 1,
#   width = NA,
#   height = NA,
#   units = c("in", "cm", "mm"),
#   dpi = 600
# )
#
# FILENAME <- "MCMC Plot Big N.png"
# PATH <- "C:/Users/Luke/Sussex Code/Sussex R/Poisson Tests/Plot Images/PNGs/"
# ggsave(
#   FILENAME,
#   plot = mcmc_plot,
#   device = "png",
#   path = PATH,
#   scale = 1,
#   width = NA,
#   height = NA,
#   units = c("in", "cm", "mm"),
#   dpi = 600
# )