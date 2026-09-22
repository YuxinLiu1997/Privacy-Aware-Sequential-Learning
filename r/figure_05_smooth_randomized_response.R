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

