# Title     : Uniform Generation Counter
# Objective : Compare the number of generated uniform R.Vs in 2 different Poisson algorithms
# Created by: Luke
# Created on: 29/04/2021

N <- 100000
lambda <- 15

sequential_search <- function(lambda) {
  # Modified to count the number of uniform R.Vs
  # "Non-Uniform Random Variate Generation"
  # Devroye, Chapter. 10, page 505

  u_count <- 0
  X <- 0
  total <- prod <- exp(-lambda) # for large lambda, this can cause overflows
  U <- runif(1)
  u_count <- u_count + 1
  while (U > total) {
    X <- X + 1
    prod <- ((lambda / X) * prod)
    total <- (total + prod)
  }
  return(u_count)
}

uniform_multiplication <- function(lambda) {
  # Modified to count the number of uniform R.Vs
  # "Non-Uniform Random Variate Generation"
  # Devroye, Chapter. 10, page 504

  u_count <- 0
  X <- 0
  prod <- 1
  constant <- exp(-lambda)
  while (TRUE) {
    U <- runif(1)
    u_count <- u_count + 1
    prod <- prod * U
    if (prod > constant) {
      X <- X + 1
    }
    else { break }
  }
  return(u_count)
}

mult <- rep(0, N)
seq <- rep(0, N)
for (i in 1:N) {
  mult[i] <- uniform_multiplication(lambda)
  seq[i] <- sequential_search(lambda)
}


print(paste("Inversion Method      : ", sum(seq), "   uniforms : ~", round(sum(seq) / N), "  per Poisson r.v"))
print(paste("Multiplication Method : ", sum(mult), " uniforms : ~", round(sum(mult) / N), " per Poisson r.v"))
