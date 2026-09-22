# You may need this package for zeta function
if (!require(gsl)) install.packages("gsl")
library(gsl)

# Constants
sigma <- sqrt(2)
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





# Load required package
library(stats)

# Define parameters
epsilon <- 0.1
sigma <- 1
a <- 1 / (1 + exp(epsilon))

# Define the conditional probability integrands (with ln in exponent)
integrand_theta_pos1_part1 <- function(sn, ln) {
  dnorm(sn, mean = 1, sd = sigma) * (1 - a * exp(epsilon * (sn + sigma^2 * ln / 2)))
}

integrand_theta_pos1_part2 <- function(sn, ln) {
  dnorm(sn, mean = 1, sd = sigma) * a * exp(-epsilon * (sn + sigma^2 * ln / 2))
}

integrand_theta_neg1_part1 <- function(sn, ln) {
  dnorm(sn, mean = -1, sd = sigma) * (1 - a * exp(epsilon * (sn + sigma^2 * ln / 2)))
}

integrand_theta_neg1_part2 <- function(sn, ln) {
  dnorm(sn, mean = -1, sd = sigma) * a * exp(-epsilon * (sn + sigma^2 * ln / 2))
}

# Define full probability expressions
P_theta_pos1 <- function(ln) {
  upper <- -sigma^2 * ln / 2
  integrate(integrand_theta_pos1_part1, -Inf, upper, ln = ln)$value +
    integrate(integrand_theta_pos1_part2, upper, Inf, ln = ln)$value
}

P_theta_neg1 <- function(ln) {
  upper <- -sigma^2 * ln / 2
  integrate(integrand_theta_neg1_part1, -Inf, upper, ln = ln)$value +
    integrate(integrand_theta_neg1_part2, upper, Inf, ln = ln)$value
}

# Evaluate over range of ln
ln_values <- seq(-30, 30, length.out = 100)
P_pos1 <- sapply(ln_values, P_theta_pos1)
P_neg1 <- sapply(ln_values, P_theta_neg1)

# Plot
plot(ln_values, P_pos1, type = "l", col = "blue", lwd = 2,
     ylim = range(c(P_pos1, P_neg1)),
     xlab = expression(l[n]), ylab = "Probability",
     main = expression(paste("Conditional Probabilities vs. ", l[n])))
lines(ln_values, P_neg1, col = "red", lwd = 2)
legend("topright", legend = c(expression(P(x[n]==-1~"|"~l[n]~","~theta==+1)),
                              expression(P(x[n]==-1~"|"~l[n]~","~theta==-1))),
       col = c("blue", "red"), lwd = 2)
