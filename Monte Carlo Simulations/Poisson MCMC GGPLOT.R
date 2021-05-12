# Title     : Poisson MCMC GGPlot
# Objective : Use ggplot to plot my MCMC sampler and illustrate the importance of burnin time
# Created by: Luke
# Created on: 30/04/2021

if (!exists("poisson_MCMC_sampler", mode = "function")) source("Poisson MH_MCMC Sampler.R")
library("ggplot2")


N <- 100000
lambda <- 100

pois_pmf <- function(x, lambda) {
  # This is the poisson PMF
  return(exp(x * log(lambda) - lambda - lgamma(x + 1)))
}

# seed and burnin are not both required, either seed the sampler or use burnin
X <- poisson_MCMC_sampler(N, lambda, seed = 0, burnin = 0, pois_pmf)
inbuilt_X <- rpois(N, lambda)

k <- min(X):max(X)
# Use this to extract the fd breaks for histogram
h <- hist(X, breaks = "fd", plot = FALSE)

# Create Data frames for use in ggplot
sim_data <- data.frame(inbuilt = inbuilt_X, my_pois = X)
line_data <- data.frame(xs = k, ys = exp(k * log(lambda) - lambda - lgamma(k + 1)))

title <- paste("Poisson MCMC: Metropolis - Hastings Random Walk\nBurnin Period = 0")
ggplot(sim_data, aes(x = my_pois)) +
  geom_histogram(aes(y = ..density.., color = "Poisson MCMC"), breaks = h$breaks, fill = "steelblue", size = 0.3, alpha = 0.5) +
  geom_line(data = line_data, aes(x = xs, y = ys, color = "Builtin rpois")) +
  scale_colour_manual(name = NULL, values = c("red", "black")) +
  labs(title = title, x = "X", y = "Density") +
  # Line below will need to be changed depending if burnin is on and lambda is changed
  scale_x_continuous(breaks = seq(0, 130, 25)) +
  # Comment out line below if font issues
  theme(text = element_text(size = 18, family = "serif"))


###############
# File Output #
###############

# FILENAME <- "No_Burnin.png"
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