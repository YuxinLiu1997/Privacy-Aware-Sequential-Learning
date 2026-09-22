#Genearate Figure K1: Evolution of l_n

# Load required package
library(pracma) # For numerical integration

# Define parameters
sigma <- sqrt(2)

# Define function a(epsilon)
a <- function(epsilon) {
  1 / (1 + exp(epsilon))
}

# Define the integrand functions for theta = -1
integrand1_theta_minus1 <- function(sn, ln, sigma) {
  (1 / (sqrt(2 * pi) * sigma)) * exp(-((sn + 1)^2) / (2 * sigma^2))
}

integrand2_theta_minus1 <- function(sn, ln, sigma, epsilon, a) {
  (1 / (sqrt(2 * pi) * sigma)) * exp(-((sn + 1)^2) / (2 * sigma^2)) * a * exp(epsilon * (sn + sigma^2 * ln / 2))
}

integrand3_theta_minus1 <- function(sn, ln, sigma, epsilon, a) {
  (1 / (sqrt(2 * pi) * sigma)) * exp(-((sn + 1)^2) / (2 * sigma^2)) * a * exp(-epsilon * (sn + sigma^2 * ln / 2))
}

# Probability function for theta = -1
P_xn_given_ln_theta_minus1 <- function(ln, sigma, epsilon) {
  a_val <- a(epsilon)
  term1 <- integral(function(sn) integrand1_theta_minus1(sn, ln, sigma), -Inf, -sigma^2 * ln / 2)
  term2 <- integral(function(sn) integrand2_theta_minus1(sn, ln, sigma, epsilon, a_val), -Inf, -sigma^2 * ln / 2)
  term3 <- integral(function(sn) integrand3_theta_minus1(sn, ln, sigma, epsilon, a_val), -sigma^2 * ln / 2, Inf)
  return(term1 - term2 + term3)
}

# Define the integrand functions for theta = +1
integrand1_theta_plus1 <- function(sn, ln, sigma) {
  (1 / (sqrt(2 * pi) * sigma)) * exp(-((sn - 1)^2) / (2 * sigma^2))
}

integrand2_theta_plus1 <- function(sn, ln, sigma, epsilon, a) {
  (1 / (sqrt(2 * pi) * sigma)) * exp(-((sn - 1)^2) / (2 * sigma^2)) * a * exp(epsilon * (sn + sigma^2 * ln / 2))
}

integrand3_theta_plus1 <- function(sn, ln, sigma, epsilon, a) {
  (1 / (sqrt(2 * pi) * sigma)) * exp(-((sn - 1)^2) / (2 * sigma^2)) * a * exp(-epsilon * (sn + sigma^2 * ln / 2))
}

# Probability function for theta = +1
P_xn_given_ln_theta_plus1 <- function(ln, sigma, epsilon) {
  a_val <- a(epsilon)
  term1 <- integral(function(sn) integrand1_theta_plus1(sn, ln, sigma), -Inf, -sigma^2 * ln / 2)
  term2 <- integral(function(sn) integrand2_theta_plus1(sn, ln, sigma, epsilon, a_val), -Inf, -sigma^2 * ln / 2)
  term3 <- integral(function(sn) integrand3_theta_plus1(sn, ln, sigma, epsilon, a_val), -sigma^2 * ln / 2, Inf)
  return(term1 - term2 + term3)
}

# Define the ratio function
ratio_function <- function(ln, sigma, epsilon) {
  P_minus1 <- P_xn_given_ln_theta_minus1(ln, sigma, epsilon)
  P_plus1 <- P_xn_given_ln_theta_plus1(ln, sigma, epsilon)
  return((1 - P_plus1) / (1 - P_minus1))
}

# Simulate l_n based on the recurrence relationship (integration-based)
simulate_ln <- function(ln_start, n_steps, sigma, epsilon) {
  ln_values <- numeric(n_steps)
  ln_values[1] <- ln_start
  
  for (i in 2:n_steps) {
    ratio <- ratio_function(ln_values[i - 1], sigma, epsilon)
    ln_values[i] <- ln_values[i - 1] + log(ratio)
  }
  
  return(ln_values)
}

# Parameters for simulation
ln_start <- 6
n_steps <- 200

# Epsilon values to compare
epsilon_values <- c(0.1, 0.5, 1)

# Generate l_n values for each epsilon (integration-based)
ln_curves <- lapply(epsilon_values, function(epsilon) simulate_ln(ln_start, n_steps, sigma, epsilon))

# Compute the constant C for the recursive formula
epsilon <- 1
C <- 1 / (exp(epsilon) + 1) * (exp(epsilon + epsilon^2 * sigma^2 / 2) - exp(-epsilon + epsilon^2 * sigma^2 / 2))

# Simulate l_n based on the new recurrence formula
simulate_ln_new <- function(ln_start, n_steps) {
  ln_values <- numeric(n_steps)
  ln_values[1] <- ln_start
  
  for (n in 2:n_steps) {
    ln <- ln_values[n - 1]
    ln_values[n] <- ln + C * exp(-epsilon * sigma^2 * ln / 2)
  }
  
  return(ln_values)
}

# Generate l_n values for non-private case (ε = 100)
epsilon_non_private <- 100
ln_non_private <- simulate_ln(ln_start, n_steps, sigma, epsilon_non_private)

# Reset the plotting layout and margins to default
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2)) 

# Left plot (a): Evolution of l_n for different epsilon values (integration-based)
plot(1:n_steps, ln_curves[[1]], type = "l", col = "blue", lwd = 2, 
     xlab = "Step (n)", ylab = expression(l[n]))
lines(1:n_steps, ln_curves[[2]], col = "red", lwd = 2)
lines(1:n_steps, ln_curves[[3]], col = "green", lwd = 2)

# Add the non-private function curve (dashed black line)
lines(1:n_steps, ln_non_private, col = "black", lwd = 2, lty = 2)

# Update the legend to include the non-private function
legend("topleft", legend = c(expression(epsilon == 0.1), expression(epsilon == 0.5), expression(epsilon == 1), 
                             "Non-private"), 
       col = c("blue", "red", "green", "black"), 
       lty = c(1, 1, 1, 2), lwd = 2)

mtext("(a)", side = 3, line = 0.5, adj = 0, cex = 1.2, font = 2) # Label for left plot
grid()

# Right plot (b): Comparison of integration-based and recursive formula for epsilon = 1
plot(1:n_steps, ln_curves[[3]], type = "l", col = "blue", lwd = 2, 
     xlab = "Step (n)", ylab = expression(l[n]))
lines(1:n_steps, ln_values_new, col = "red", lwd = 2)
legend("topleft", legend = c("Exact", "Approximate"),
       col = c("blue", "red"), lty = 1, lwd = 2)
mtext("(b)", side = 3, line = 0.5, adj = 0, cex = 1.2, font = 2) # Label for right plot
grid()

