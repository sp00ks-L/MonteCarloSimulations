# Title     : Chi Squared Goodness-of-fit test
# Objective : Implement the above
# Created by: Luke
# Created on: 08/05/2021

if (!exists("poisson_generator", mode = "function")) source("poisson_generator.R")
library("plyr")

chisq_test <- function(X, alpha = 0.05, theoretical_lambda) {
  #' Chi-Squared Test
  #' @description  Performs a Chi-squared goodness-of-fit-test
  #' @param X  The data sample
  #' @param alpha  The confidence level for the statistical test. Default = 0.05
  #' @param theoretical_lambda  Your theorised Lambda value for the distribution of interest
  #' @output  None - function prints test result
  #' @details This function was developed specifically for use with the discrete Poisson distribution

  counts <- as.data.frame(count(X))
  probs <- dpois(counts$x, theoretical_lambda)
  expected <- probs * N

  ind <- which(expected >= 5)
  counts <- counts[ind,]
  expected <- expected[ind]

  result <- chisq.test(counts$freq, p = expected, rescale.p = TRUE)
  cat("Confidence Level:   ", alpha, "\n")
  cat("Calculated P-value: ", result$p.value, "\n")
  if (result$p.value < alpha) {
    print("Reject Null: Distributions are not identical")
  } else {
    print("Cannot Reject Null: Distributions are identical")
  }
}
