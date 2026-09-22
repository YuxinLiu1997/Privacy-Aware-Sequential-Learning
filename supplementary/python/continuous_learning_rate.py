import numpy as np
import matplotlib.pyplot as plt
from math import erf, sqrt, exp, log
from matplotlib.lines import Line2D

# ============================================================
# Settings
# ============================================================

np.random.seed(123)

N = 10000  # number of agents
n_runs = 20  # number of simulation runs per regime
sigma = 1.0
theta_true = 1  # condition on theta = +1

# Homogeneous privacy levels
eps_low = 0.1
eps_medium = 0.5
eps_high = 1.0

# Heterogeneous privacy distribution
kappa = 1.0  # epsilon_n ~ Uniform[0, kappa]

# How many sample paths to show lightly
n_thin_paths = 8


# ============================================================
# Helper functions
# ============================================================

def normal_cdf(x):
    """Standard normal CDF using erf."""
    return 0.5 * (1.0 + erf(x / sqrt(2.0)))


def gaussian_cdf(x, mean, sigma):
    return normal_cdf((x - mean) / sigma)


def decision_threshold(l, sigma):
    """
    Gaussian sequential-learning threshold.
    Agent chooses +1 if s > t(l).

    For signals s ~ N(theta, sigma^2), theta in {-1,+1},
    threshold is:
        t(l) = - sigma^2 * l / 2
    """
    return -0.5 * sigma ** 2 * l


def prob_x_plus_nonprivate(l, theta, sigma):
    """
    Non-private case:
        x = a = +1 iff s > t(l)
    """
    t = decision_threshold(l, sigma)
    return 1.0 - gaussian_cdf(t, theta, sigma)


def prob_x_plus_smooth_rr(l, theta, sigma, beta):
    """
    Smooth randomized response:
        a = +1 if s > t(l), a = -1 otherwise.
        P(x != a | s, l) = 0.5 * exp(-beta * |s - t(l)|)

    This function computes:
        P(x = +1 | l, theta)

    by numerical integration over the signal distribution.
    """
    t = decision_threshold(l, sigma)

    # numerical integration grid around the Gaussian mass
    lower = theta - 8 * sigma
    upper = theta + 8 * sigma
    grid = np.linspace(lower, upper, 1200)

    density = (
            1.0 / (sigma * np.sqrt(2 * np.pi))
            * np.exp(-0.5 * ((grid - theta) / sigma) ** 2)
    )

    a_plus = grid > t
    flip_prob = 0.5 * np.exp(-beta * np.abs(grid - t))

    # If a=+1, x=+1 with prob 1-flip_prob.
    # If a=-1, x=+1 with prob flip_prob.
    px_plus_given_s = np.where(a_plus, 1.0 - flip_prob, flip_prob)

    return np.trapezoid(px_plus_given_s * density, grid)


def prob_x_plus_heterogeneous(l, theta, sigma, kappa, n_beta_grid=80):
    """
    Personalized privacy:
        epsilon_n ~ Uniform[0, kappa].

    We compute the marginal probability:
        P(x = +1 | l, theta)
        = E_epsilon[P(x = +1 | l, theta, epsilon)]
    """
    betas = np.linspace(1e-4, kappa, n_beta_grid)
    probs = np.array([
        prob_x_plus_smooth_rr(l, theta, sigma, beta)
        for beta in betas
    ])
    return np.mean(probs)


def update_llr(l, x, regime, sigma, beta=None, kappa=None):
    """
    Public belief update:
        l_{n+1} = l_n + log P(x | l_n, theta=+1) / P(x | l_n, theta=-1)
    """
    if regime == "nonprivate":
        p_plus_pos = prob_x_plus_nonprivate(l, theta=1, sigma=sigma)
        p_plus_neg = prob_x_plus_nonprivate(l, theta=-1, sigma=sigma)

    elif regime == "homogeneous":
        p_plus_pos = prob_x_plus_smooth_rr(l, theta=1, sigma=sigma, beta=beta)
        p_plus_neg = prob_x_plus_smooth_rr(l, theta=-1, sigma=sigma, beta=beta)

    elif regime == "heterogeneous":
        p_plus_pos = prob_x_plus_heterogeneous(l, theta=1, sigma=sigma, kappa=kappa)
        p_plus_neg = prob_x_plus_heterogeneous(l, theta=-1, sigma=sigma, kappa=kappa)

    else:
        raise ValueError("Unknown regime")

    # numerical safety
    eps = 1e-12
    p_plus_pos = np.clip(p_plus_pos, eps, 1 - eps)
    p_plus_neg = np.clip(p_plus_neg, eps, 1 - eps)

    if x == 1:
        increment = log(p_plus_pos / p_plus_neg)
    else:
        increment = log((1 - p_plus_pos) / (1 - p_plus_neg))

    return l + increment


def generate_report(l, theta_true, regime, sigma, beta=None, kappa=None):
    """
    Generate one privatized report x_n under the true state.
    """
    s = np.random.normal(loc=theta_true, scale=sigma)
    t = decision_threshold(l, sigma)

    a = 1 if s > t else -1

    if regime == "nonprivate":
        return a

    if regime == "homogeneous":
        beta_used = beta

    elif regime == "heterogeneous":
        beta_used = np.random.uniform(0.0, kappa)

    else:
        raise ValueError("Unknown regime")

    flip_prob = 0.5 * np.exp(-beta_used * abs(s - t))

    if np.random.rand() < flip_prob:
        return -a
    else:
        return a


def simulate_llr_path(N, theta_true, regime, sigma, beta=None, kappa=None):
    """
    Simulate one trajectory of l_n.
    """
    l_values = np.zeros(N + 1)
    l = 0.0

    for n in range(1, N + 1):
        x = generate_report(
            l=l,
            theta_true=theta_true,
            regime=regime,
            sigma=sigma,
            beta=beta,
            kappa=kappa
        )

        l = update_llr(
            l=l,
            x=x,
            regime=regime,
            sigma=sigma,
            beta=beta,
            kappa=kappa
        )

        l_values[n] = l

    return l_values


# ============================================================
# Run simulations
# ============================================================

regimes = {
    "Uniform[0,1]": {
        "regime": "heterogeneous",
        "beta": None,
        "kappa": kappa,
        "color": "purple",
        "linestyle": "-",
    },
    "Low": {
        "regime": "homogeneous",
        "beta": eps_low,
        "kappa": None,
        "color": "red",
        "linestyle": "--",
    },
    "Medium": {
        "regime": "homogeneous",
        "beta": eps_medium,
        "kappa": None,
        "color": "blue",
        "linestyle": ":",
    },
    "High": {
        "regime": "homogeneous",
        "beta": eps_high,
        "kappa": None,
        "color": "green",
        "linestyle": "-.",
    },
    "Non-private": {
        "regime": "nonprivate",
        "beta": None,
        "kappa": None,
        "color": "black",
        "linestyle": "--",
    },
}

all_paths = {}

for name, cfg in regimes.items():
    paths = []

    print(f"Simulating {name}...")

    for r in range(n_runs):
        path = simulate_llr_path(
            N=N,
            theta_true=theta_true,
            regime=cfg["regime"],
            sigma=sigma,
            beta=cfg["beta"],
            kappa=cfg["kappa"]
        )
        paths.append(path)

    all_paths[name] = np.array(paths)

# ============================================================
# Plot
# ============================================================

n_grid = np.arange(N + 1)

plt.figure(figsize=(9, 5.8))

for name, cfg in regimes.items():
    paths = all_paths[name]
    mean_path = paths.mean(axis=0)

    # Thin sample paths
    for i in range(min(n_thin_paths, n_runs)):
        plt.plot(
            n_grid,
            paths[i],
            color=cfg["color"],
            linestyle=cfg["linestyle"],
            alpha=0.18,
            linewidth=1
        )

    # Mean path
    plt.plot(
        n_grid,
        mean_path,
        color=cfg["color"],
        linestyle=cfg["linestyle"],
        linewidth=2.5,
        label=name
    )

plt.xlabel(r"$n$", fontsize=16)
plt.ylabel(r"$l_n$", fontsize=16)
plt.title("Log-likelihood Ratio Dynamics under Different Privacy Regimes", fontsize=16)

plt.grid(alpha=0.25)
plt.legend(fontsize=12, loc="lower right")

plt.tight_layout()
plt.show()