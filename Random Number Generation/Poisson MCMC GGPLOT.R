# Title     : Poisson MCMC GGPlot
# Objective : Use ggplot to illustrate importance of burnin time for MCMC
# Created by: Luke
# Created on: 30/04/2021

if (!exists("poisson_MCMC_sampler", mode = "function")) source("Sussex R/Final Code/Poisson/Poisson MH_MCMC Sampler.R")
library("ggplot2")

# This is extra - just to change plot fonts for report
# library(extrafont)
# loadfonts(device = "win")


N <- 100000
lambda <- 100

target <- function(x, lambda) {
  # This is the poisson PMF
  return(exp(x * log(lambda) - lambda - lgamma(x + 1)))
}

X <- poisson_MCMC_sampler(N, lambda, seed = 0, burnin = 2500, target)
inbuilt_X <- rpois(N, lambda)

k <- min(X):max(X)
# Use this to extract the fd breaks for histogram
h <- hist(X, breaks = "fd", plot = FALSE)

# Create Data frames for use in ggplot
sim_data <- data.frame(inbuilt = inbuilt_X, my_pois = X)
line_data <- data.frame(xs = k, ys = exp(k * log(lambda) - lambda - lgamma(k + 1)))

title <- paste("Poisson MCMC: Metropolis - Hastings Random Walk\nBurnin Period = 2500")
ggplot(sim_data, aes(x = my_pois)) +
  geom_histogram(aes(y = ..density.., color = "Poisson MCMC"), breaks = h$breaks, fill = "steelblue", size = 0.3, alpha = 0.5) +
  geom_line(data = line_data, aes(x = xs, y = ys, color = "Inbuilt rpois")) +
  scale_colour_manual(name = NULL, values = c("red", "black")) +
  labs(title = title, x = "X", y = "Density") +
  scale_x_continuous(breaks = seq(0, 130, 25))
# add line below with extra font imports to change plot font
# theme(text = element_text(size = 13, family = "CM Sans CE"))


###############
# File Output #
###############

# FILENAME <- "Burnin.pdf"
# PATH <- "C:/Users/Luke/Sussex Code/Sussex R/Poisson Tests/Plot Images/PDFS/"
# ggsave(
#   FILENAME,
#   plot = last_plot(),
#   device = "pdf",
#   path = PATH,
#   scale = 1,
#   width = NA,
#   height = NA,
#   units = c("in", "cm", "mm"),
#   dpi = 300
# )