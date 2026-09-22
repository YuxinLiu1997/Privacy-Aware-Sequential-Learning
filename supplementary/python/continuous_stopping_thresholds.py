import numpy as np
import matplotlib.pyplot as plt
from math import erf, sqrt, log
from matplotlib.lines import Line2D

# ============================================================
# Settings
# ============================================================

np.random.seed(123)

sigma = 1.0
theta_true = 1

N_max = 20000          # maximum simulation horizon
n_runs = 100           # number of Monte Carlo runs

B_values = np.arange(2, 30, 1)   # confidence thresholds

# Homogeneous privacy levels
eps_list = [0.1, 0.5, 1.0]

# Personalized privacy: epsilon_n ~ Uniform[0, kappa]
kappa = 1.0

# If a path does not hit B within N_max:
# False: ignore that run for that B
# True: count it as N_max
use_censoring = False


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


def smooth_flip_prob(distance, beta):
    """
    Homogeneous smooth randomized response:
        P(flip | s,l) = 0.5 exp(- beta |s - t(l)|)
    """
    return 0.5 * np.exp(-beta * distance)


def personalized_avg_flip_prob(distance, kappa):
    """
    Personalized privacy with epsilon ~ Uniform[0,kappa].
    We average 0.5 exp(-epsilon * d) over epsilon.

    E[0.5 exp(-epsilon d)]
    =
    0.5 * (1 - exp(-kappa d)) / (kappa d)

    with limit 0.5 when d=0.
    """
    d = np.asarray(distance)
    out = np.empty_like(d, dtype=float)

    small = d < 1e-10
    out[small] = 0.5
    out[~small] = 0.5 * (1.0 - np.exp(-kappa * d[~small])) / (kappa * d[~small])

    return out


# ============================================================
# Precompute update probabilities
# ============================================================

def prob_x_plus_nonprivate(l, theta, sigma):
    """
    Non-private:
        x=+1 iff s > t(l)
    """
    t = decision_threshold(l, sigma)
    return 1.0 - gaussian_cdf(t, theta, sigma)


def prob_x_plus_privacy(l, theta, sigma, regime, beta=None, kappa=None):
    """
    Compute P(x=+1 | l, theta) by numerical integration.
    """
    t = decision_threshold(l, sigma)

    lower = theta - 8 * sigma
    upper = theta + 8 * sigma
    s_grid = np.linspace(lower, upper, 1600)

    density = (
        1.0 / (sigma * np.sqrt(2 * np.pi))
        * np.exp(-0.5 * ((s_grid - theta) / sigma) ** 2)
    )

    a_plus = s_grid > t
    distance = np.abs(s_grid - t)

    if regime == "homogeneous":
        flip_prob = smooth_flip_prob(distance, beta)
    elif regime == "personalized":
        flip_prob = personalized_avg_flip_prob(distance, kappa)
    else:
        raise ValueError("Unknown privacy regime")

    px_plus_given_s = np.where(a_plus, 1.0 - flip_prob, flip_prob)

    return np.trapezoid(px_plus_given_s * density, s_grid)


def precompute_llr_increments(regime, sigma, beta=None, kappa=None, l_min=-10, l_max=60, n_l_grid=2500):
    """
    Precompute:
        p_plus_pos(l) = P(x=+1 | l, theta=+1)
        p_plus_neg(l) = P(x=+1 | l, theta=-1)

    Then compute increments for x=+1 and x=-1:
        Delta_plus(l)  = log p_plus_pos / p_plus_neg
        Delta_minus(l) = log (1-p_plus_pos) / (1-p_plus_neg)
    """
    l_grid = np.linspace(l_min, l_max, n_l_grid)

    p_pos = []
    p_neg = []

    for l in l_grid:
        if regime == "nonprivate":
            pp = prob_x_plus_nonprivate(l, theta=1, sigma=sigma)
            pn = prob_x_plus_nonprivate(l, theta=-1, sigma=sigma)
        elif regime == "homogeneous":
            pp = prob_x_plus_privacy(l, theta=1, sigma=sigma, regime="homogeneous", beta=beta)
            pn = prob_x_plus_privacy(l, theta=-1, sigma=sigma, regime="homogeneous", beta=beta)
        elif regime == "personalized":
            pp = prob_x_plus_privacy(l, theta=1, sigma=sigma, regime="personalized", kappa=kappa)
            pn = prob_x_plus_privacy(l, theta=-1, sigma=sigma, regime="personalized", kappa=kappa)
        else:
            raise ValueError("Unknown regime")

        p_pos.append(pp)
        p_neg.append(pn)

    p_pos = np.clip(np.array(p_pos), 1e-12, 1 - 1e-12)
    p_neg = np.clip(np.array(p_neg), 1e-12, 1 - 1e-12)

    delta_plus = np.log(p_pos / p_neg)
    delta_minus = np.log((1 - p_pos) / (1 - p_neg))

    return l_grid, delta_plus, delta_minus


def interp_increment(l, x, l_grid, delta_plus, delta_minus):
    """
    Interpolate update increment at current l.
    """
    l_clipped = np.clip(l, l_grid[0], l_grid[-1])

    if x == 1:
        return np.interp(l_clipped, l_grid, delta_plus)
    else:
        return np.interp(l_clipped, l_grid, delta_minus)


# ============================================================
# Generate reports
# ============================================================

def generate_report(l, theta_true, sigma, regime, beta=None, kappa=None):
    """
    Generate x_n under the true state theta=+1.
    """
    s = np.random.normal(loc=theta_true, scale=sigma)
    t = decision_threshold(l, sigma)

    a = 1 if s > t else -1

    if regime == "nonprivate":
        return a

    elif regime == "homogeneous":
        flip_prob = 0.5 * np.exp(-beta * abs(s - t))

    elif regime == "personalized":
        eps_n = np.random.uniform(0.0, kappa)
        flip_prob = 0.5 * np.exp(-eps_n * abs(s - t))

    else:
        raise ValueError("Unknown regime")

    if np.random.rand() < flip_prob:
        return -a
    else:
        return a


# ============================================================
# Simulate one path and get stopping times for all B
# ============================================================

def simulate_stopping_times(
    B_values,
    N_max,
    theta_true,
    sigma,
    regime,
    l_grid,
    delta_plus,
    delta_minus,
    beta=None,
    kappa=None
):
    """
    Simulate one LLR path and record T_B for each B.
    """
    B_values = np.array(B_values)
    stopping_times = np.full(len(B_values), np.nan)

    l = 0.0
    not_hit = np.ones(len(B_values), dtype=bool)

    for n in range(1, N_max + 1):
        x = generate_report(
            l=l,
            theta_true=theta_true,
            sigma=sigma,
            regime=regime,
            beta=beta,
            kappa=kappa
        )

        inc = interp_increment(
            l=l,
            x=x,
            l_grid=l_grid,
            delta_plus=delta_plus,
            delta_minus=delta_minus
        )

        l += inc

        newly_hit = not_hit & (l >= B_values)
        stopping_times[newly_hit] = n
        not_hit[newly_hit] = False

        if not np.any(not_hit):
            break

    if use_censoring:
        stopping_times[np.isnan(stopping_times)] = N_max

    return stopping_times


# ============================================================
# Regimes to simulate
# ============================================================

regimes = []

regimes.append({
    "name": "Non-private",
    "regime": "nonprivate",
    "beta": None,
    "kappa": None,
    "color": "black",
    "linestyle": "--"
})

for eps in eps_list:
    regimes.append({
        "name": rf"$\varepsilon={eps}$",
        "regime": "homogeneous",
        "beta": eps,
        "kappa": None,
        "color": None,
        "linestyle": "-"
    })

regimes.append({
    "name": rf"$\varepsilon_n\sim U[0,{kappa}]$",
    "regime": "personalized",
    "beta": None,
    "kappa": kappa,
    "color": "purple",
    "linestyle": "-"
})

colors = ["red", "blue", "green"]
color_idx = 0
for cfg in regimes:
    if cfg["color"] is None:
        cfg["color"] = colors[color_idx]
        color_idx += 1


# ============================================================
# Run simulations
# ============================================================

results = {}

l_min = -10
l_max = max(B_values) + 30

for cfg in regimes:
    print(f"Precomputing increments for {cfg['name']}...")

    l_grid, delta_plus, delta_minus = precompute_llr_increments(
        regime=cfg["regime"],
        sigma=sigma,
        beta=cfg["beta"],
        kappa=cfg["kappa"],
        l_min=l_min,
        l_max=l_max,
        n_l_grid=2500
    )

    all_T = []

    print(f"Simulating {cfg['name']}...")

    for r in range(n_runs):
        T_vals = simulate_stopping_times(
            B_values=B_values,
            N_max=N_max,
            theta_true=theta_true,
            sigma=sigma,
            regime=cfg["regime"],
            l_grid=l_grid,
            delta_plus=delta_plus,
            delta_minus=delta_minus,
            beta=cfg["beta"],
            kappa=cfg["kappa"]
        )
        all_T.append(T_vals)

    all_T = np.array(all_T)

    mean_T = np.nanmean(all_T, axis=0)
    hit_rate = np.mean(~np.isnan(all_T), axis=0)

    results[cfg["name"]] = {
        "mean_T": mean_T,
        "hit_rate": hit_rate,
        "all_T": all_T,
        "cfg": cfg
    }


# ============================================================
# Plot expected stopping time vs B
# ============================================================

plt.figure(figsize=(9, 5.8))

for name, res in results.items():
    cfg = res["cfg"]
    mean_T = res["mean_T"]
    hit_rate = res["hit_rate"]

    # only plot B values with reasonable hitting rate
    valid = hit_rate >= 0.8

    plt.plot(
        B_values[valid],
        mean_T[valid],
        color=cfg["color"],
        linestyle=cfg["linestyle"],
        linewidth=2.5,
        marker="o",
        markersize=4,
        label=name
    )

plt.xlabel(r"Confidence threshold $B$", fontsize=16)
plt.ylabel(r"Estimated expected stopping time $\widehat{E}[T_B]$", fontsize=15)
plt.title("Expected High-Confidence Stopping Time vs. Confidence Threshold", fontsize=16)

plt.grid(alpha=0.25)
plt.legend(fontsize=11)
plt.tight_layout()
plt.show()


# ============================================================
# Optional: log-scale plot
# ============================================================

plt.figure(figsize=(9, 5.8))

for name, res in results.items():
    cfg = res["cfg"]
    mean_T = res["mean_T"]
    hit_rate = res["hit_rate"]

    valid = hit_rate >= 0.8

    plt.plot(
        B_values[valid],
        mean_T[valid],
        color=cfg["color"],
        linestyle=cfg["linestyle"],
        linewidth=2.5,
        marker="o",
        markersize=4,
        label=name
    )

plt.yscale("log")
plt.xlabel(r"Confidence threshold $B$", fontsize=16)
plt.ylabel(r"Estimated expected stopping time $\widehat{E}[T_B]$ (log scale)", fontsize=15)
plt.title("Expected High-Confidence Stopping Time vs. Confidence Threshold", fontsize=16)

plt.grid(alpha=0.25, which="both")
plt.legend(fontsize=11)
plt.tight_layout()
plt.show()


# ============================================================
# Print summary table
# ============================================================

print("\nSummary: mean stopping time and hit rate")
print("========================================")
for name, res in results.items():
    print(f"\n{name}")
    print("B\tMean T_B\tHit rate")
    for B, mean_T, hit_rate in zip(B_values, res["mean_T"], res["hit_rate"]):
        print(f"{B:.1f}\t{mean_T:.2f}\t\t{hit_rate:.2f}")