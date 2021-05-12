# Title     : Markov Monte Carlo GGPLOT
# Objective : Use GGPLOT to visualise my markov chains
# Created by: Luke
# Created on: 14/03/2021

library("Matrix")
library("ggplot2")
# Additional package for visualising the markov chain - not required
# library("markovchain")

get_alpha <- function(proposed, current, pi, Q) {
              #' get_alpha
              #' @description Calculates the acceptance probability using the Metropolis-Hastings ratio
              #' @param proposed The propsed value within the Markov Chain
              #' @param current The current value within the Markov Chain
              #' @param pi The invariant (target) distribution
              #' @param Q The transition matrix
              #' @return The acceptance probability
  top <- pi[proposed] * Q[proposed, current]
  bottom <- pi[current] * Q[current, proposed]
  return(min(1, top / bottom))
}

get_invariant <- function(p) {
              #' @title Get Invariant
              #' @description Calculates the invariant distribution for the provided transition matrix
              #' @param p The transition matrix for your Markov Chain
              #' @return The invariant distribution for P
  n <- nrow(p)
  # I create a local p so as not to change the global p variable
  local_p <- p
  I <- diag(n)
  local_p <- t(local_p) - I
  local_p <- rbind(local_p, c(1, 1, 1))
  b <- c(0, 0, 0, 1)
  # For further explaination, check the appendix of my dissertation
  pi <- drop(solve(t(local_p) %*% local_p, t(local_p) %*% b))

  return(pi)
}

# Transition probabilities, using Matrix package
probs <- c(1 / 3, 1 / 3, 1 / 3, 0.5, 0, 0.5, 0, 1, 0)
P <- matrix(probs, ncol = 3, nrow = 3, byrow = TRUE)

# If you want to plot the Markov Chain, uncomment the following code
# mc <- new("markovchain",
#           states = c("A", "B", "C"),
#           transitionMatrix = P,
#           name = "Weather")
# plot(mc)

PI <- get_invariant(P)

# Create state vector
s <- rep(0, 3)
X <- sample(3, size = 1) # Random init of starting node (1:A, 2:B, 3:C)
epochs <- 50000
# K is my chain
K <- rep(0, epochs)
K[1] <- X # Set initial state to equal X at time 0


for (i in 2:epochs) {
  # This samples my transition matrix with the corresponding probabilites
  k <- sample(ncol(P), size = 1, prob = P[X,])
  alpha <- get_alpha(k, X, PI, P)
  U <- runif(1)
  if (U <= alpha) {
    X <- k
  }
  K[i] <- X
}

# Get frequency of each state
for (j in 1:3) {
  s[j] <- length(which(K == j)) / epochs
}

# GGPLOT requires data in datafram
line_data <- data.frame(state = 1:3, freq = s)

# Calculating Expected Value - I originally used a moving average but decided it was not necessary
Mean <- seq(0, epochs)
for (i in 1:epochs) { Mean[i] <- mean(K[1:i]) }
m <- mean(K) # expected value of k
TMean <- rep(m, epochs)

# Calculating Variance
Var <- seq(0, epochs)
for (i in 1:epochs) { Var[i] <- var(K[1:i]) }
v <- var(K) # variance of k
TVar <- rep(v, epochs)

mean_var <- data.frame(xs = seq_along(Mean), m = Mean, v = Var)
mean_var_lines <- data.frame(tmean = TMean, tvar = TVar)

MCMC_standard_error <- sd(K) / sqrt(epochs)
chart_title <- paste(epochs, " Iterations : MCMC Error : ", round(MCMC_standard_error, 4))
plot_marker <- 13

ggplot(line_data, aes(x = state)) +
  geom_linerange(aes(color = "Theoretical Distribution"), x = 1:3, ymin = c(0, 0, 0), ymax = PI, size = 1.5) +
  geom_point(aes(y = freq, color = "Monte Carlo Distribution"), shape = plot_marker, size = 3) +
  labs(title = "The Probability of Each State at the Invariant Distribution", y = "Pr(state)", x = "State", subtitle = chart_title) +
  scale_x_continuous(breaks = c(1.0, 2.0, 3.0), labels = c("A", "B", "C")) +
  scale_colour_manual(name = NULL, values = c("black", "#4c72b0")) +
  guides(color = guide_legend(override.aes = list(shape = c(plot_marker, NA), size = c(3, 1.5), linetype = c(0, 1)))) +
  ylim(0, max(s) + 0.05) +
  theme(text = element_text(size = 18, family = "serif"))

# Used this to write data to a csv
# I used this to calculate the confidence intervals
# write.csv(data.frame(data = K), "C:\\Users\\Luke\\Sussex Code\\Sussex R\\Final Code\\Poisson\\results.csv")

ggplot(data = mean_var, aes(xs)) +
  geom_point(aes(y = Mean, color = "Expected State"), shape = 1, size = 1) +
  geom_line(data = mean_var_lines, aes(x = seq_along(TMean), y = TMean, color = "Mean Convergence"), size = 1.5, alpha = 0.5) +
  labs(title = "Mean Convergence of Monte Carlo Simulation", subtitle = chart_title, x = "Iteration", y = "Expected State") +
  # ylim(1, 3) +
  scale_y_continuous(breaks = c(1, 2, 3), limits = c(1, 3), labels = c("A", "B", "C")) +
  scale_colour_manual(name = NULL, values = c('black', '#4c72b0')) +
  guides(color = guide_legend(override.aes = list(shape = c(1, NA), linetype = c(0, 1)))) +
  theme(text = element_text(size = 18, family = "serif"))

ggplot(data = mean_var, aes(xs)) +
  geom_point(aes(y = Var, color = "State Variance"), shape = 1, size = 1) +
  geom_line(data = mean_var_lines, aes(x = seq_along(TVar), y = TVar, color = "Variance Convergence"), size = 1.5, , alpha = 0.5) +
  labs(title = "Variance Convergence of Monte Carlo Simulation", subtitle = chart_title, x = "Iteration", y = "State Variance") +
  scale_colour_manual(name = NULL, values = c('black', '#4c72b0')) +
  ylim(min(TVar) - 0.05, max(TVar) + 0.05) +
  guides(color = guide_legend(override.aes = list(shape = c(1, NA), linetype = c(0, 1)))) +
  theme(text = element_text(size = 18, family = "serif"))


###############
# File Output #
###############

# FILENAME <- paste("MCMC_", epochs, ".png")
# PATH <- "C:/Users/Luke/Sussex Code/Sussex R/Poisson Tests/Plot Images/PNGs/"
# ggsave(
#   FILENAME,
#   plot = state_plot,
#   device = "png",
#   path = PATH,
#   scale = 1,
#   width = NA,
#   height = NA,
#   units = c("in", "cm", "mm"),
#   dpi = 600
# )
# FILENAME <- paste("MCMC_mean_convergence.png")
# ggsave(
#   FILENAME,
#   plot = mean_plot,
#   device = "png",
#   path = PATH,
#   scale = 1,
#   width = NA,
#   height = NA,
#   units = c("in", "cm", "mm"),
#   dpi = 600
# )
# FILENAME <- paste("MCMC_var_convergence.png")
# ggsave(
#   FILENAME,
#   plot = var_plot,
#   device = "png",
#   path = PATH,
#   scale = 1,
#   width = NA,
#   height = NA,
#   units = c("in", "cm", "mm"),
#   dpi = 600
# )
