#Figure 1 Probability of correct cascade vs. privacy budget

# Load the viridis package for better color contrast
# Load necessary library
library(viridis)
library(ggplot2)
library(gridExtra)
library(cowplot)

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
