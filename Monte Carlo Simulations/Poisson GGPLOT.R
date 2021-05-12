# Title     : Poisson GGPLOT
# Objective : Plot the poisson distribution in GGPLOT using my RNG poisson_generator.
#             Compare to builtin rpois using Kolmogorov-Smirnov Goodness of Fit test and Chi-Squared
# Created by: Luke
# Created on: 01/04/2021

# You will have to change the paths in order to correctly execute the script

# This is the optional KS test function that I created
# if (!exists("KStest", mode = "function")) source("Kolmogorov-Smirnov Test.R")

if (!exists("poisson_generator", mode = "function")) source("Sussex R/Final Code/Submission Code/poisson_generator.R")
if (!exists("chisq_test", mode = "function")) source("Sussex R/Final Code/Submission Code/Chi_Squared_Test.R")
library("ggplot2")
library("tictoc")


N <- 100000
lambda <- 500

tic() # this times the execution of my poisson generator
X <- poisson_generator(N, lambda)
toc(log = TRUE, quiet = TRUE)
log.txt <- tic.log(format = TRUE)

inbuilt_X <- rpois(N, lambda)

# This is the Freedman-Diaconis rule for binwidth. I actually use plot 'h' below to extract this for ggplot
bw <- floor(2 * IQR(X) / length(X)^(1 / 3))
k <- min(X):max(X)
# Use this to extract the Freedman-Diaconis breaks
h <- hist(X, breaks = "fd", plot = FALSE)


# GGPLOT requires the data within a data.frame
sim_data <- data.frame(inbuilt = inbuilt_X, my_pois = X)
line_data <- data.frame(xs = k, ys = exp(k * log(lambda) - lambda - lgamma(k + 1)))

title <- paste("Poisson Distribution : ", tail(log.txt, n = 1)) # tail() is to retrieve the latest runtime in log.txt
ggplot(sim_data, aes(x = my_pois)) +
  geom_histogram(aes(y = ..density.., color = "PTRS"), breaks = h$breaks, fill = "steelblue", size = 0.3, alpha = 0.5) +
  geom_line(data = line_data, aes(x = xs, y = ys, color = "Builtin rpois")) +
  scale_colour_manual(name = NULL, values = c("red", "black")) +
  labs(title = title, x = "X", y = "Density") +
  theme(text = element_text(size = 18, family = "serif")) +
  # Will need to adjust the below scale depending on the value of lambda
  scale_x_continuous(breaks = seq(0, 70, 10))
# scale_x_continuous(breaks = seq(400, 600, 25))

# This is a seperate function / file that tests my generator and its goodness-of-fit using KS + Chi-squared
chisq_test(X, theoretical_lambda = lambda)
# KStest(N, X, inbuilt_X, 0.05) <- this is the optional KS test


###############
# File Output #
###############

# FILENAME <- "Poiss_high_Lambda.png"
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