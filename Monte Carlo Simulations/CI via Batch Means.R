# Title     : Confidence Intervals via Batch Means
# Objective : Calculate Confidence Intervals for MCMC simulation based on batch means
# Created by: Luke
# Created on: 04/05/2021


# You will have to change the path to run this script
# You may also have to actually generate the CSV
results <- read.csv("../Poisson/NSUMresults.csv")
my_data <- results$data

# This is a splitting function used in the get_ci function
chunk2 <- function(my_data, b) split(my_data, cut(seq_along(my_data), b, labels = FALSE))


get_ci <- function(data, batch_num, conf_level = 0.05) {
  N <- length(data)
  # Split data into batches of size N / batch_num
  batches <- chunk2(data, batch_num)

  batch_means <- rep(0, N / batch_num)
  for (i in seq_along(batches)) {
    total <- 0
    for (num in batches[[i]]) {
      total <- total + num
    }
    batch_means[i] <- total / (N / batch_num)
  }

  # Calculate Sample Mean
  sample_mean <- sum(batch_means) / batch_num

  # Calculate Sample Variance
  sample_variance <- 0
  for (k in batch_means) {
    sample_variance <- sample_variance + (k - sample_mean)^2
  }
  sample_variance <- sample_variance / (batch_num - 1)

  # two tailed t-test. DoF = batch_num - 1
  alpha <- qt(conf_level / 2, batch_num - 1, lower.tail = FALSE)
  coeff <- alpha * sqrt(sample_variance / batch_num)

  return(c(round(sample_mean - coeff, 4), round(sample_mean + coeff, 4)))
}


get_ci(my_data, 500)
MCMC_standard_error <- sd(my_data) / sqrt(length(my_data))
print(round(MCMC_standard_error, 3))
