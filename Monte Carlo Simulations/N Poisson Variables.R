# Title     : Addition of N Poisson Random Variables
# Objective : To implement a MCMC sampler that is capable of sampling from an abritrary 'Zn' distribution
#             where Zn is the sum of N i.i.d poisson random variables. With Poisson that means if
#             if lambda = 3 and N = 10, the poisson distribution with be Pois(3 * 10 = 30)
# Created by: Luke
# Created on: 01/05/2021

if (!exists("poisson_MCMC_sampler", mode = "function")) source("Sussex R/Final Code/Submission Code/Poisson MH_MCMC Sampler.R")
library("ggplot2")

N <- 1000
lambda <- 3
num_of_variables <- 10

theoretical_pmf <- function(x, lambda) {
  # The theoretical PMF for the addition of N Poisson random variables
  # N given by the variable 'num_of_variables'
  new_lambda <- 0
  for (i in 1:num_of_variables) {
    new_lambda <- new_lambda + lambda
  }
  return(exp(x * log(new_lambda) - new_lambda - lgamma(x + 1)))
}

pois_pmf <- function(x, lambda) {
  # Computationally stable variation of the Poisson PMF
  # https://en.wikipedia.org/wiki/Poisson_distribution#Computational_methods
  return(exp(x * log(lambda) - lambda - lgamma(x + 1)))
}

theoretical_cdf <- function(x, lambda) {
  # The theoretical CDF for the addition of N poisson random variables
  # In essence, this is just the sum of the PMFs
  result <- 0
  for (i in 0:x) {
    result <- result + theoretical_pmf(i, lambda)
  }
  return(result)
}

# Where Pois(lam) + Pois(lam) = Pois(2*lam)  -->> Pois(N*lam)
new_lambda <- num_of_variables * lambda
X <- rep(0, N)
# Runs my MCMC sampler of Pois(lambda) N times, summing each iteration
# This was more stable than running extremely high values of lambda
for (i in 1:num_of_variables) {
  X <- X + poisson_MCMC_sampler(N, lambda, 0, 0, pois_pmf)
}

# Finding range to calculate the theoretical CDF. I have done max(X) + lambda to extend the plot horizontally (purely visual)
k <- min(X):max(X) + lambda
theo_pdf <- theoretical_pmf(k, lambda)
my_cdf <- 0
for (i in seq_along(k)) {
  my_cdf[i] <- theoretical_cdf(i, lambda)
}

# Use this to extract the fd breaks
h <- hist(X, breaks = "fd", plot = FALSE)
# GGPLOT requires data in dataframe
sim_data <- data.frame(my_pois = X)
line_data <- data.frame(xs = k, ys = exp(k * log(new_lambda) - new_lambda - lgamma(k + 1)))
cdf_data <- data.frame(cdf = my_cdf)

# title <- paste("Poisson MCMC PMF compared with the Theoretical PMF\nN = 10 : Pois(3)")
# ggplot(sim_data, aes(x = my_pois)) +
#   geom_histogram(aes(y = ..density.., color = "Poisson MCMC"), breaks = h$breaks, fill = "steelblue", size = 0.1, alpha = 0.5) +
#   geom_line(data = line_data, aes(x = xs, y = ys, color = "Theoretical PMF")) +
#   scale_colour_manual(name = NULL, values = c("black", "red")) +
#   labs(title = title, x = "X", y = "Density") +
#   scale_x_continuous(breaks = seq(0, 65, 5)) +
#   theme(text = element_text(size = 18, family = "serif"))


# Highlighted code below is for plotting the ECDF of my sample vs the theoretical CDF for my 10 Poisson variables

title <- paste("Poisson MCMC ECDF compared with the Theoretical CDF\nN = 10 : Pois(3) : 1,000 Samples")
ggplot(sim_data, aes(my_pois, colour = "MCMC ECDF")) +
  stat_ecdf(geom = "step", alpha = 0.75) +
  geom_step(data = cdf_data, aes(x = seq_along(theo_pdf), y = cdf, colour = "Theoretical CDF"), alpha = 0.75) +
  scale_colour_manual(name = NULL, values = c("black", "red")) +
  labs(title = title, x = "X", y = "Pr(X <= x)") +
  scale_x_continuous(breaks = seq(0, max(X) + 1, 5)) +
  theme(text = element_text(size = 18, family = "serif"))

###############
# File Output #
###############

# FILENAME <- "Sum of N Pois Vars CDF 1k.png"
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