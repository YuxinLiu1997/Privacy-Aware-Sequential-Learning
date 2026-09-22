import numpy as np
import matplotlib.pyplot as plt
from math import erf, sqrt, log

# ============================================================
# Settings
# ============================================================

np.random.seed(123)

sigma = 1.0
theta_true = 1

N_max = 30000       # maximum simulation horizon
n_runs = 100        # Monte Carlo runs

B = 10.0            # fixed confidence threshold
c = 0.001           # stopping cost

eps_grid = np.linspace(0.05, 2.0, 30)

use_censoring = True
# If True: if not hit by N_max, count T_B=N_max.
# If False: ignore non-hit runs when estimating E[T_B].

# Numerical integration settings
l_min = -10
l_max = B + 30
n_l_grid = 2000


# ============================================================
# Basic functions
# ============================================================

def normal_cdf(x):
    return 0.5 * (1.0 + erf(x / sqrt(2.0)))


def gaussian_cdf(x, mean, sigma):
    return normal_cdf((x - mean) / sigma)


def decision_threshold(l, sigma):
    """
    Agent chooses +1 if s > t(l).
    """
    return -0.5 * sigma**2 * l


def prob_x_plus_smooth_rr(l, theta, sigma, beta):
    """
    Compute P(x=+1 | l, theta) under smooth randomized response:
        P(flip | s,l) = 0.5 exp(- beta |s - t(l)|)
    """
    t = decision_threshold(l, sigma)

    lower = theta - 8 * sigma
    upper = theta + 8 * sigma
    s_grid = np.linspace(lower, upper, 1400)

    density = (
        1.0 / (sigma * np.sqrt(2 * np.pi))
        * np.exp(-0.5 * ((s_grid - theta) / sigma) ** 2)
    )

    a_plus = s_grid > t
    distance = np.abs(s_grid - t)

    flip_prob = 0.5 * np.exp(-beta * distance)

    # If a=+1, x=+1 with prob 1-flip.
    # If a=-1, x=+1 with prob flip.
    px_plus_given_s = np.where(a_plus, 1.0 - flip_prob, flip_prob)

    return np.trapezoid(px_plus_given_s * density, s_grid)


def precompute_llr_increments(beta, sigma, l_min, l_max, n_l_grid):
    """
    Precompute LLR increments on a grid of l values.
    """
    l_grid = np.linspace(l_min, l_max, n_l_grid)

    p_pos = []
    p_neg = []

    for l in l_grid:
        pp = prob_x_plus_smooth_rr(l, theta=1, sigma=sigma, beta=beta)
        pn = prob_x_plus_smooth_rr(l, theta=-1, sigma=sigma, beta=beta)

        p_pos.append(pp)
        p_neg.append(pn)

    p_pos = np.clip(np.array(p_pos), 1e-12, 1 - 1e-12)
    p_neg = np.clip(np.array(p_neg), 1e-12, 1 - 1e-12)

    delta_plus = np.log(p_pos / p_neg)
    delta_minus = np.log((1 - p_pos) / (1 - p_neg))

    return l_grid, delta_plus, delta_minus


def interp_increment(l, x, l_grid, delta_plus, delta_minus):
    """
    Interpolate LLR increment.
    """
    l_clipped = np.clip(l, l_grid[0], l_grid[-1])

    if x == 1:
        return np.interp(l_clipped, l_grid, delta_plus)
    else:
        return np.interp(l_clipped, l_grid, delta_minus)


def generate_report(l, theta_true, sigma, beta):
    """
    Generate one privatized report x under smooth randomized response.
    """
    s = np.random.normal(loc=theta_true, scale=sigma)
    t = decision_threshold(l, sigma)

    a = 1 if s > t else -1

    flip_prob = 0.5 * np.exp(-beta * abs(s - t))

    if np.random.rand() < flip_prob:
        return -a
    else:
        return a


def simulate_one_stopping_time(
    B,
    N_max,
    theta_true,
    sigma,
    beta,
    l_grid,
    delta_plus,
    delta_minus
):
    """
    Simulate one path and return T_B.
    """
    l = 0.0

    for n in range(1, N_max + 1):
        x = generate_report(
            l=l,
            theta_true=theta_true,
            sigma=sigma,
            beta=beta
        )

        inc = interp_increment(
            l=l,
            x=x,
            l_grid=l_grid,
            delta_plus=delta_plus,
            delta_minus=delta_minus
        )

        l += inc

        if l >= B:
            return n

    if use_censoring:
        return N_max
    else:
        return np.nan


def estimate_expected_stopping_time(beta):
    """
    Estimate E[T_B] for one privacy budget beta.
    """
    l_grid, delta_plus, delta_minus = precompute_llr_increments(
        beta=beta,
        sigma=sigma,
        l_min=l_min,
        l_max=l_max,
        n_l_grid=n_l_grid
    )

    T_values = []

    for _ in range(n_runs):
        T = simulate_one_stopping_time(
            B=B,
            N_max=N_max,
            theta_true=theta_true,
            sigma=sigma,
            beta=beta,
            l_grid=l_grid,
            delta_plus=delta_plus,
            delta_minus=delta_minus
        )
        T_values.append(T)

    T_values = np.array(T_values, dtype=float)

    mean_T = np.nanmean(T_values)
    hit_rate = np.mean(~np.isnan(T_values))

    return mean_T, hit_rate


def posterior_accuracy(B):
    """
    Approximate final decision accuracy after reaching threshold B.
    """
    return np.exp(B) / (1.0 + np.exp(B))


# ============================================================
# Compute objective values
# ============================================================

accuracy = posterior_accuracy(B)

objective_values = []
mean_T_values = []
hit_rates = []

for eps in eps_grid:
    print(f"Simulating epsilon = {eps:.3f}...")

    mean_T, hit_rate = estimate_expected_stopping_time(beta=eps)

    J = accuracy - c * mean_T

    objective_values.append(J)
    mean_T_values.append(mean_T)
    hit_rates.append(hit_rate)

objective_values = np.array(objective_values)
mean_T_values = np.array(mean_T_values)
hit_rates = np.array(hit_rates)


# ============================================================
# Plot objective value vs privacy budget
# ============================================================

plt.figure(figsize=(9, 5.8))

plt.plot(
    eps_grid,
    objective_values,
    color="black",
    linewidth=2.5,
    marker="o",
    markersize=4
)

plt.xlabel(r"Privacy budget $\varepsilon$", fontsize=16)
plt.ylabel(r"Objective value $J(\varepsilon)$", fontsize=16)
plt.title(
    fr"Platform Objective vs. Privacy Budget, $B={B}$, $c={c}$",
    fontsize=16
)

plt.grid(alpha=0.25)
plt.tight_layout()
plt.show()


# ============================================================
# Optional: also plot expected stopping time
# ============================================================

plt.figure(figsize=(9, 5.8))

plt.plot(
    eps_grid,
    mean_T_values,
    color="blue",
    linewidth=2.5,
    marker="o",
    markersize=4
)

plt.xlabel(r"Privacy budget $\varepsilon$", fontsize=16)
plt.ylabel(r"Estimated expected stopping time $\widehat{E}[T_B]$", fontsize=16)
plt.title(
    fr"Expected Stopping Time vs. Privacy Budget, $B={B}$",
    fontsize=16
)

plt.grid(alpha=0.25)
plt.tight_layout()
plt.show()


# ============================================================
# Print summary
# ============================================================

print("\nSummary")
print("=======")
print(f"Fixed B = {B}")
print(f"Accuracy approximation = {accuracy:.6f}")
print(f"Cost c = {c}")

print("\nepsilon\tE[T_B]\t\tObjective\tHit rate")
for eps, mean_T, J, hit_rate in zip(eps_grid, mean_T_values, objective_values, hit_rates):
    print(f"{eps:.3f}\t{mean_T:.2f}\t\t{J:.6f}\t{hit_rate:.2f}")