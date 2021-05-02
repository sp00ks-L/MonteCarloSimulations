# Title     : The Convolution of 2 Poisson PMFs
# Objective : Implement a simple script illustrating the sum of 2 i.i.d poisson variables via convoltion of their PMFs
# Created by: Luke
# Created on: 30/04/2021


pmf <- function(x) {
  # Computationally stable PMF for Poisson Distribution
  return(exp(x * log(lambda) - lambda - lgamma(x + 1)))
}

convolve_pmf <- function(lambda) {
  # Set abritrary support to 2 * lambda
  support <- 0:(2 * lambda)
  # Set z support to 2 * support because we are summing 2 * poisson variables
  z_support <- 0:(2 * length(support))

  # I wrote the loop in this way so that I could write the following line
  # total <- total + (pmf(z - y) * pmf(y))
  # Whilst the rest of the loop is a bit cluttered, this line nicely
  # illustrates that actual operation being employed, that is
  # sum(  PMFx(z - y) * PMFy(y)  ) where z = x + y therefore x = z - y
  results <- 0
  ind <- 1
  for (z in z_support) {
    total <- 0
    for (y in support) {
      total <- total + (pmf(z - y) * pmf(y))
    }
    results[ind] <- total
    ind <- ind + 1
  }
  return(results)
}

####################
# Additional Exploration of the Poisson probability generating function (PGF)
####################

pgf <- function(k, s, lambda) {
  # The normal PGF is just e^lambda(s - 1)
  # This function allows the calculation of the kth derivative
  # When k = 0, this is just the normal PGF
  # Allows use within the get_prob function to extract probabilities from PGF
  return(lambda^k * exp(lambda * (s - 1)))
}

get_prob <- function(k, lambda) {
  # Extracting probabilities from the Poisson PGF where k
  # is the kth derivative of the PGF
  return(pgf(k, 0, lambda) / factorial(k))
}


lambda <- 34
results <- convolve_pmf(lambda)
cdf <- cumsum(results)
plot(results)
title(paste("PMF Convolution Plot : lambda = ", lambda, " : N = 2"))
plot(cdf)
title(paste("CDF Plot : lambda = ", lambda, " : N = 2"))


### PGF - alternate way to plot the PMF using the PGF

pgf_results <- get_prob(0:(lambda * 2), lambda)
plot(pgf_results)
title(paste("PGF plot : lambda = ", lambda))