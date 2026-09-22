#Generate Figure 3:Smooth randomized response and the signal likelihoods

# Load necessary library
library(ggplot2)
# Arrange plots side by side
library(gridExtra)
library(cowplot)
# Set parameters
mu_positive <- 1
mu_negative <- -1
sigma <- 1
epsilon <- 1
t_ln <- -3

# Generate data
x <- seq(-5, 5, by = 0.01)
f_positive <- dnorm(x, mean = mu_positive, sd = sigma)
f_negative <- dnorm(x, mean = mu_negative, sd = sigma)
u_n <- (1 / 2) * exp(-epsilon * (abs(x - t_ln)))


u_n_complement <- 1 - u_n

# Apply the (1 - threshold) to the blue and orange distributions
f_combined_negative <- f_negative * u_n
f_combined_positive <- f_positive * u_n

# Create a data frame for plotting
data <- data.frame(x = x, 
                   f_negative = f_negative, 
                   f_positive = f_positive, 
                   u_n = u_n, 
                   f_combined_negative = f_combined_negative, 
                   f_combined_positive = f_combined_positive)


# Plot p1: Only the blue and orange normal distributions with the threshold line and shading
p1 <- ggplot(data, aes(x = x)) + 
  geom_line(aes(y = f_negative), color = "blue") +
  geom_line(aes(y = f_positive), color = "orange") +
  geom_vline(xintercept = t_ln, color = "black", linetype = "dashed") +
  geom_area(aes(y = f_negative), data = data[data$x >= t_ln, ], fill = "blue", alpha = 0.2) + 
  geom_area(aes(y = f_positive), data = data[data$x >= t_ln, ], fill = "orange", alpha = 0.2) +
  labs(title = "(a)",  # Title added
       x = expression(s[n]), 
       y = "Density") +
  ylim(0, 0.5) +
  theme_minimal()

# Plot p2: Original functions with smooth randomized response and threshold
p2 <- ggplot(data, aes(x = x)) + 
  geom_line(aes(y = f_negative), color = "blue") +
  geom_line(aes(y = f_positive), color = "orange") +
  geom_line(aes(y = u_n), color = "black") +
  geom_vline(xintercept = t_ln, color = "black", linetype = "dashed") +
  labs(title = "(b)",  # Title added
       x = expression(s[n]), 
       y = "Density") +
  ylim(0, 0.5) +
  theme_minimal()

# Plot p3: Combined response with (1 - threshold) only for x >= t_ln, with shading and y-axis limit
p3 <- ggplot(data[data$x >= t_ln, ], aes(x = x)) + 
  geom_line(aes(y = f_combined_negative), color = "blue") +
  geom_line(aes(y = f_combined_positive), color = "orange") +
  geom_vline(xintercept = t_ln, color = "black", linetype = "dashed") +
  geom_area(aes(y = f_combined_negative), fill = "blue", alpha = 0.2) + 
  geom_area(aes(y = f_combined_positive), fill = "orange", alpha = 0.2) +
  labs(title = "(c)",  # Title added
       x = expression(s[n]), 
       y = "Density") +
  ylim(0, 0.1) + # Set y-axis limit
  theme_minimal()

# Arrange the three plots in a single row
grid.arrange(p1, p2, p3, ncol = 3)






# Load required package
library(pracma) # For numerical integration

# Define parameters
sigma <- 1

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



# Simulate first time to make correct action
simulate_first_n_with_sn <- function(ln_start, sigma, n_steps, repetitions, epsilon) {
  results <- numeric(repetitions)
  
  for (rep in 1:repetitions) {
    ln <- ln_start
    for (n in 1:n_steps) {
      # Sample s_n from N(1, sigma)
      s_n <- rnorm(1, mean = 1, sd = sigma)
      
      # Compute condition value (used to determine intended action)
      condition_value <- 2 * s_n / sigma^2 + ln
      
      # Compute probability of flipping the action
      u_n <- a(epsilon) * exp(-epsilon * abs(s_n + sigma^2 * ln / 2))
      
      # Determine intended action (+1 or -1)
      intended_action <- ifelse(condition_value > 0, 1, -1)
      
      # Simulate action flipping based on probability u_n
      final_action <- ifelse(runif(1) < u_n, -intended_action, intended_action)
      
      # If the final action is +1, store the time step and stop
      if (final_action == 1) {
        results[rep] <- n
        break
      }
      
      # Compute probabilities for updating l_n
      P_plus1 <- P_xn_given_ln_theta_plus1(ln, sigma, epsilon)
      P_minus1 <- P_xn_given_ln_theta_minus1(ln, sigma, epsilon)
      
      # Update l_n based on the original intended action
      if (intended_action == 1) {
        ln <- ln + log(P_plus1 / P_minus1)
      } else {
        ln <- ln + log((1 - P_plus1) / (1 - P_minus1))
      }
    }
  }
  
  return(results)
}


ln_start <- 0    # Initial belief value
sigma <- 5       # Noise level
n_steps <- 50   # Maximum steps allowed
repetitions <- 10000  # Number of trials
epsilon <- 10  # Privacy parameter

# Run the simulation
results <- simulate_first_n_with_sn(ln_start, sigma, n_steps, repetitions, epsilon)

# Print results
cat("Max first n of correct action:", max(results), "\n")

# Compute and print average
average_n <- mean(results, na.rm = TRUE)
cat("Average first n of correct action:", average_n, "\n")






#binary case



# Load the viridis package for better color contrast
# Load necessary library
library(viridis)

# Set up plot parameters with extra margin space on the right
par(mfrow = c(1, 1), mar = c(5, 5, 2, 8))  # Extend right margin for legends

# Define function to calculate w_k and epsilon_k
calculate_w_k <- function(k, alpha) {
  (1 - alpha ^ ((k - 2) / (k - 1))) / 
    (1 - alpha ^ ((k - 2) / (k - 1)) + alpha ^ (-1 / (k - 1)) - alpha)
}

calculate_epsilon_k <- function(w_k) {
  log((1 - w_k) / w_k)
}

# Function for probability of right cascade based on rho and epsilon
probability_right_cascade <- function(epsilon, k, p) {
  u <- 1 / (1 + exp(epsilon))
  u_tilde <- u * (1 - p) + p * (1 - u)
  rho <- (1 - u_tilde) / u_tilde
  (rho ^ k - 1) / (rho ^ (2 * k) - 1)
}

# Function to compute non-private probability
non_private_probability <- function(p) {
  p^2 / (2 * p^2 - 2 * p + 1)
}





# Set parameters
p_values <- c(0.90, 0.70, 0.55)
# Different p values
epsilon_values <- seq(0, 5, length.out = 1000)  # Epsilon range
k_values <- 2:5  # Range of k values
colors <- c("red", "green", "blue", "purple")  # Colors for different k values
line_types <- c("solid", "dashed", "dotted") # Different line styles for different p values

# Set up the plot
plot(NULL, xlim = range(epsilon_values), ylim = c(0.5, 1), type = "n",
     xlab = expression(epsilon), ylab = "Probability of Correct Cascade",
     cex.lab = 1.5, cex.axis = 1.2)

# Loop over different p values with different line styles
for (j in seq_along(p_values)) {
  p <- p_values[j]
  alpha <- (1 - p) / p  # Compute alpha
  
  for (i in seq_along(k_values)) {
    k <- k_values[i]
    
    # Compute thresholds
    w_k <- calculate_w_k(k, alpha)
    epsilon_k <- calculate_epsilon_k(w_k)
    
    w_k_plus_1 <- calculate_w_k(k + 1, alpha)
    epsilon_k_plus_1 <- calculate_epsilon_k(w_k_plus_1)
    
    # Filter epsilon values in the range (epsilon_k+1, epsilon_k]
    epsilon_range <- epsilon_values[epsilon_values <= epsilon_k & epsilon_values > epsilon_k_plus_1]
    prob_right_cascade <- sapply(epsilon_range, function(eps) probability_right_cascade(eps, k, p))
    
    # Plot with different line styles for different p values
    lines(epsilon_range, prob_right_cascade, col = colors[i], lwd = 2, lty = line_types[j])
    
    # Add vertical dashed lines at segment boundaries
    abline(v = epsilon_k_plus_1, col = "gray", lty = "dashed", lwd = 0.8)
  }
}

# Add non-private case horizontal lines (matching line styles of corresponding p-values)
non_private_probs <- sapply(p_values, non_private_probability)

for (j in seq_along(p_values)) {
  abline(h = non_private_probs[j], col = "black", lty = line_types[j], lwd = 2)
}

# Add legends **outside the plot**
par(xpd = TRUE)  # Allow plotting outside the normal plot region

# Legend for Non-private and k values (Non-private first)
legend("topright", inset = c(0.005, 0.4), 
       legend = c("Non-private", paste("k =", k_values)), 
       col = c("black", colors), 
       lwd = 2, 
       lty = c("solid", rep("solid", length(k_values))), 
       bty = "n")


# Legend for p values (placed outside)
legend("bottomright", inset = c(0.045, 0), 
       legend = paste("p =", p_values),
       col = "black", lwd = 2, lty = line_types, bty = "n")


# Reset xpd to default
par(xpd = FALSE)
