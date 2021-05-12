# Title     : Poisson MCMC Sampler: Metropolis-Hastings Algorithm
# Objective : Create a MCMC sampler that answers Q3.2
# Created by: Luke
# Created on: 29/04/2021


get_proposal <- function(x) {
  #' Get Proposal
  #' @description For use in my MCMC sampler, this returns a proposal for the Metropolis-Hastings algorithm
  #' @param x The current MCMC value
  #' @output The proposed value for the next step in the chain
  if (x == 0) {
    return(x + rbinom(1, 1, 0.5))
  }
  else {
    return(x + 2 * (rbinom(1, 1, 0.5) - 0.5))
  }
}

poisson_MCMC_sampler <- function(N, lambda, seed, burnin, target) {
  #' MCMC Poisson Sampler
  #' @description Using MCMC methods, will sample the desired Poisson distribution
  #' @param N The sample size
  #' @param lambda The Poisson parameter
  #' @param seed Allows seeding of the sampler with a chosen value
  #' @param burnin Burn-in time for MCMC sampler
  #' @param target The target distribution function - here that is the Poisson PMF
  #' @details Seed and Burnin generally not used together.
  #' @output A sample of size N from the Poisson(lambda) distribution
  X <- seed
  proposed <- get_proposal(seed)
  if (burnin == 0) { N <- N + 1 }
  for (k in 2:(N + burnin)) {
    u <- runif(1)
    alpha <- target(proposed, lambda) / target(X[k - 1], lambda)
    if (u <= alpha) {
      X <- c(X, proposed)
    }
    else { X <- c(X, X[k - 1]) }
    proposed <- get_proposal(X[k])
  }
  return(X[-(1:burnin)])
}

