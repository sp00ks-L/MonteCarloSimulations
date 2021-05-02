# Title     : Poisson GGPLOT
# Objective : Plot the poisson distribution in GGPLOT using my RNG poisson_generator. Compare to inbuilt rpois using Kolmogorov-Smirnov Goodness of Fit test
# Created by: Luke
# Created on: 01/04/2021

if (!exists("KStest", mode = "function")) source("Sussex R/Final Code/Poisson/Kolmogorov-Smirnov Test.R")
if (!exists("poisson_generator", mode = "function")) source("Sussex R/Final Code/Poisson/poisson_generator.R")
library("ggplot2")
library("tictoc")


# Added serif fonts for final plots. Not required to run script
# library(extrafont)
# loadfonts(device = "win")
# windowsFonts()

N <- 100000
lambda <- 34

tic() # this times the execution of my poisson generator
X <- poisson_generator(N, lambda)
toc(log = TRUE, quiet = TRUE)
log.txt <- tic.log(format = TRUE)

inbuilt_X <- rpois(N, lambda)

# This is the Freedman-Diaconis rule for binwidth. I actually use plot 'h' below to extract this for ggplot
bw <- floor(2 * IQR(X) / length(X)^(1 / 3))
k <- min(X):max(X)
# Use this to extract the fd breaks
h <- hist(X, breaks = "fd", plot = FALSE)

sim_data <- data.frame(inbuilt = inbuilt_X, my_pois = X)
line_data <- data.frame(xs = k, ys = exp(k * log(lambda) - lambda - lgamma(k + 1)))

title <- paste("Poisson Distribution : ", tail(log.txt, n = 1)) # tail() is to retrieve the latest runtime in log.txt
ggplot(sim_data, aes(x = my_pois)) +
  geom_histogram(aes(y = ..density.., color = "PTRS"), breaks = h$breaks, fill = "steelblue", size = 0.3, alpha = 0.5) +
  geom_line(data = line_data, aes(x = xs, y = ys, color = "Inbuilt rpois")) +
  scale_colour_manual(name = NULL, values = c("red", "black")) +
  labs(title = title, x = "X", y = "Density")

# Will need to adjust the below scale depending on the value of lambda
# scale_x_continuous(breaks = seq(400, 600, 25))

# Additional line below to alter the font on the plot
# theme(text = element_text(size = 13, family = "CM Sans CE")) +

# This is a seperate function / file that tests my generator and its goodness-of-fit using KS
# KStest(N, X, inbuilt_X, 0.05)

###############
# File Output #
###############

# FILENAME <- "Big lambda.pdf"
# PATH <- "C:/Users/Luke/Sussex Code/Sussex R/Poisson Tests/Plot Images/PDFS/"
# ggsave(
#   FILENAME,
#   plot = last_plot(),
#   device = "pdf",
#   path = PATH,
#   scale = 1,
#   width = NA,
#   height = NA,
#   units = c("in", "cm", "mm"),
#   dpi = 300
# )