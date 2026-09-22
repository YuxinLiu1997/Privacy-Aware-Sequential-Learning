

#Generate Figure 4: The expectation
if (!require(gsl)) install.packages("gsl")
library(gsl)

# Constants
sigma <-1
C1 <- 1
N <- 10000  # truncation point for zeta approx

# Define C(epsilon)
C_eps <- function(eps) {
  eps*sigma^2 / 2 * (exp(eps + (eps^2 * sigma^2)/2) - exp(-eps + (eps^2 * sigma^2)/2))
}

# Approximate zeta(s)
zeta_approx <- function(s, N) {
  sum((1:N)^(-s))
}

# Expected tau
E_tau <- function(eps) {
  s <- 2 / (eps * sigma^2)
  Cval <- C_eps(eps)
  if (Cval <= 0) return(NA)  # avoid taking negative number to power
  C1 * (Cval^(-s)) * zeta_approx(s, N)
}

# Plot
eps_values <- seq(0.5, 1, by = 0.01)
E_tau_values <- sapply(eps_values, E_tau)

plot(eps_values, E_tau_values, type = "l", col = "blue", lwd = 2,
     xlab = expression(epsilon), ylab = expression(E[tau]))
grid()

# Find minimum
min_index <- which.min(E_tau_values)
min_eps <- eps_values[min_index]
min_tau <- E_tau_values[min_index]

# Print result

cat("Minimum E[tau] is", round(min_tau, 4), "at epsilon =", round(min_eps, 4), "\n")
