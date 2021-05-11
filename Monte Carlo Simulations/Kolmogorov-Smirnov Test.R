# Title     : KS Test Functions
# Objective : Implement a Kolmogorov-Smirnov Test that can be used with a discrete distribution
# Created by: Luke
# Created on: 17/02/2021

library("KSgeneral")

get_d_val <- function(alpha, n) {
  c_alpha <- -log(alpha / 2) * 1 / 2
  d_val <- c_alpha * sqrt((n + n) / (n * n))
  return(d_val)
}

KStest <- function(N, sample, reference, alpha) {
  # result <- ks.test(sample, reference)
  # Use the KSgeneral implementation for discrete KS test
  # It is quite slow for very large sample sizes
  result <- KSgeneral::disc_ks_test(sample, ecdf(reference), exact=TRUE)
  critical_d_val <- get_d_val(alpha, N)
  cat("Critical D value:\t", round(critical_d_val, 4), "\n")
  cat("Calculated D value:\t", round(result$statistic, 4), "\n")
  cat("Confidence level:\t", alpha, "\n")
  cat("Calculated p-value:\t", round(result$p.value, 4), "\n")

  # if (result$statistic > critical_d_val) {
  #   print("Reject Null: Distributions are not identical")
  # } else {
  #   print("Cannot Reject Null: Distributions are identical")
  # }

  if (result$p.value < alpha) {
    print("Reject Null: Distributions are not identical")
  } else {
    print("Cannot Reject Null: Distributions are identical")
  }
}

