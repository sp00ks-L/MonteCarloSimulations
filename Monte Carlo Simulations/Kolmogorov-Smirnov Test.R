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
  #' Discrete Kolmogorov-Smirnov Test
  #' @description Performs the Kolmogorov-Smirnov test for a discrete distribution using
  #' the R package KSgeneral
  #' @param N The size of the sample
  #' @param sample Your observed / generated distribution
  #' @param reference The theoretical distribution
  #' @param alpha The confidence level for calulcation of the D-statistic
  #' @output None - prints test result to console
  #' @details This is quite slow when the sample size is larger ~100k and upwards

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

