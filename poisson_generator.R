# Title     : Poisson Generator
# Objective : Implement a file to encapsulate the Poisson Random Number Generator
# Created by: Luke
# Created on: 01/04/2021

sequential_search <- function(lambda) {
  # "Non-Uniform Random Variate Generation"
  # Devroye, Chapter. 10, page 505

  X <- 0
  total <- prod <- exp(-lambda) # for large lambda, this can cause overflows
  U <- runif(1)
  while (U > total) {
    X <- X + 1
    prod <- ((lambda / X) * prod)
    total <- (total + prod)
  }
  return(X)
}


logK <- function(K) {
  return(log(sqrt(2 * pi)) + (K + 1 / 2) * log(K) - K + (1 / 12 - 1 / (360 * K^2)) / K)
}

transformed_rejection <- function(lambda) {
  # "The transformed rejection method for generating Poisson random variables"
  # W. Hormann, 1993

  slam <- sqrt(lambda)
  loglam <- log(lambda)
  b <- 0.931 + 2.53 * slam
  a <- -0.059 + 0.02483 * b
  invalpha <- 1.1239 + 1.1328 / (b - 3.4)
  vr <- 0.9277 - 3.6224 / (b - 2)

  while (TRUE) {
    U <- runif(1) - 0.5
    V <- runif(1)
    us <- 0.5 - abs(U)
    K <- floor((2 * a / us + b) * U + lambda + 0.43)
    if (us >= 0.07 && V <= vr) {
      return(K)
    }
    if ((K < 0) || (us < 0.013 && V > us)) {
      next
    }
    # The LogK(K) function can be replaced with lgamma(K + 1)
    if ((log(V) + log(invalpha) - log(a / us^2 + b)) <= (-lambda + K * loglam - logK(K))) {
      return(K)
    }

  }
}

poisson_generator <- function(N, lambda) {
  # Generate N variables ~ Pois(lambda)
  rvs <- rep(0, N)
  if (lambda >= 60) {
    for (i in 1:N) {
      rvs[i] <- transformed_rejection(lambda)
    }
  }
  else if (lambda == 0) {
    return(0)
  }
  else {
    for (i in 1:N) {
      rvs[i] <- sequential_search(lambda)
    }
  }
  return(rvs)
}