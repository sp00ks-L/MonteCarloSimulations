# Title     : Goodness-of-fit Plotting
# Objective : Create a few visual aids to enhance my goodness-of-fit test using rootograms
# Created by: Luke
# Created on: 10/05/2021

if (!exists("poisson_generator", mode = "function")) source("Sussex R/Final Code/Submission Code/poisson_generator.R")
library("ggplot2")
library("lattice")
library("latticeExtra")
library("plyr")

N <- 100000
lambda <- 34

x <- poisson_generator(N, lambda)
# plot p if you want a single rootogram
p <- rootogram(~x, dfun = function(x) dpois(x, lambda = lambda), type = 'l', col = '#c44e52', lw = 3, lwd = 2)

# These lambda values need adjusting depending on your chosen lambda values
lambdav <- c(30, 32, 34, 36)
# lambdav <- c(495, 500, 505, 510)
update(p[rep(1, length(lambdav))],
       aspect = "xy",
       prepanel = function(x, ...) {
         tmp <-
           lapply(lambdav,
                  function(lambda) { prepanel.rootogram(x, dfun = function(x) dpois(x, lambda = lambda)) })
         list(xlim = range(sapply(tmp, "[[", "xlim")),
              ylim = range(sapply(tmp, "[[", "ylim")),
              dx = do.call("c", lapply(tmp, "[[", "dx")),
              dy = do.call("c", lapply(tmp, "[[", "dy")))
       },
       panel = function(x, ...) {
         panel.rootogram(x, dfun = function(x) dpois(x, lambda = lambdav[panel.number()]), col = '#c44e52')
         grid::grid.text(bquote(Pois(lambda == .(foo)),
                                where = list(foo = lambdav[panel.number()])),
                         y = 0.15,
                         gp = grid::gpar(cex = 1.5))
       },
       xlab = "",
       sub = paste("Random sample from Poisson(", lambda, ")\nusing my Poisson PRNG"), cex = 5)


