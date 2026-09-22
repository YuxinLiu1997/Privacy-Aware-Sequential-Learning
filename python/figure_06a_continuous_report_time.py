import numpy as np
import matplotlib.pyplot as plt
import math
import os
import warnings

# ============================================================
# Settings
# ============================================================

sigma = 1.0
theta_true = 1

B = 10.0                         # confidence threshold
eps_grid = np.linspace(0.1, 2.0, 25)

# Monte Carlo settings
n_runs = 300
N_max = 50000  # report-count cap; incomplete runs are reported, not discarded

# Bellman settings
bellman_grid_size = 401
bellman_max_iter = 100000
bellman_tol = 1e-8

# Two-sided report stopping rule: tau_B = inf{n >= 1: |l_n| >= B}.
# l_0 = 0; each update incorporates one submitted report.
# This comparison implements the two-sided theorem only.
two_sided_stop = True


# ============================================================
# Normal CDF
# ============================================================

try:
    from scipy.special import ndtr
    normal_cdf = ndtr
except Exception:
    def normal_cdf(x):
        x = np.asarray(x)
        return 0.5 * (1.0 + np.vectorize(math.erf)(x / np.sqrt(2.0)))


# ============================================================
# Model functions
# ============================================================

def decision_threshold(l, sigma):
    """
    Agent chooses +1 if s > t(l).
    """
    return -0.5 * sigma**2 * l


def prob_x_plus_smooth_rr_closed_form(l, theta, sigma, eps):
    """
    Closed-form P(x=+1 | l, theta) under smooth randomized response.

    Signal:
        s ~ N(theta, sigma^2)

    Initial recommendation:
        a = +1 if s > t(l), otherwise -1

    Smooth randomized response:
        P(x != a | s,l) = 0.5 * exp(- eps * |s - t(l)|)
    """
    l = np.asarray(l)
    t = decision_threshold(l, sigma)

    # P(s > t)
    z = (t - theta) / sigma
    p_s_gt_t = 1.0 - normal_cdf(z)

    # E[exp(-eps(s-t)) 1{s>t}]
    term_right = (
        np.exp(eps * t - eps * theta + 0.5 * eps**2 * sigma**2)
        * normal_cdf((theta - eps * sigma**2 - t) / sigma)
    )

    # E[exp(-eps(t-s)) 1{s<t}]
    term_left = (
        np.exp(-eps * t + eps * theta + 0.5 * eps**2 * sigma**2)
        * normal_cdf((t - theta - eps * sigma**2) / sigma)
    )

    p_x_plus = p_s_gt_t - 0.5 * term_right + 0.5 * term_left

    return np.clip(p_x_plus, 1e-12, 1.0 - 1e-12)


def llr_increments(l_values, sigma, eps):
    """
    Compute Delta_+(l) and Delta_-(l):
        Delta_+(l) = log P(x=+1|theta=+1,l) / P(x=+1|theta=-1,l)
        Delta_-(l) = log P(x=-1|theta=+1,l) / P(x=-1|theta=-1,l)
    """
    p_pos = prob_x_plus_smooth_rr_closed_form(
        l_values, theta=1, sigma=sigma, eps=eps
    )
    p_neg = prob_x_plus_smooth_rr_closed_form(
        l_values, theta=-1, sigma=sigma, eps=eps
    )

    delta_plus = np.log(p_pos / p_neg)
    delta_minus = np.log((1.0 - p_pos) / (1.0 - p_neg))

    return p_pos, delta_plus, delta_minus


def C_eps(eps, sigma):
    """Auxiliary C(eps) appearing in the definition of K(eps, sigma)."""
    if not np.isfinite(eps) or not np.isfinite(sigma) or eps <= 0 or sigma <= 0:
        raise ValueError("eps and sigma must be finite and strictly positive.")
    return (
        eps * sigma**2 / 4.0
        * np.exp(-eps + 0.5 * eps**2 * sigma**2)
        * np.expm1(2.0 * eps)
    )


def K_eps(eps, sigma):
    """
    Expected two-sided stopping-time constant:
        K(eps,sigma) = C(eps)
          - 0.5 exp(-eps + eps^2 sigma^2/2) (1 - exp(-eps^2 sigma^2)).
    """
    C = C_eps(eps, sigma)
    correction = (
        0.5 * np.exp(-eps + 0.5 * eps**2 * sigma**2)
        * (-np.expm1(-eps**2 * sigma**2))
    )
    K = C - correction
    if not np.isfinite(K) or K <= 0:
        raise FloatingPointError("K must be finite and positive; check parameter precision/range.")
    return K


def asymptotic_stopping_time(eps, B, sigma):
    """
    Revised Theorem 4.8, for fixed eps > 0 and sigma > 0 as B -> infinity:
        E[tau_B | theta] ~ exp(eps sigma^2 B / 2) / K(eps,sigma),
        tau_B = inf{n >= 1: |l_n| >= B}.
    This is a high-threshold approximation, not an exact finite-B value.
    """
    if not np.isfinite(B) or B <= 0:
        raise ValueError("B must be finite and strictly positive.")
    K = K_eps(eps, sigma)
    return np.exp(0.5 * eps * sigma**2 * B) / K


# ============================================================
# Bellman equation method
# ============================================================

def bellman_expected_stopping_time(
    eps,
    B,
    sigma,
    grid_size=401,
    max_iter=100000,
    tol=1e-8,
    two_sided=True
):
    """
    Solve Bellman equation:
        H(l) = 1 + P_+(l) H(l + Delta_+(l))
                 + (1-P_+(l)) H(l + Delta_-(l))

    Boundary:
        H(l)=0 if |l|>=B (two-sided report stopping time tau_B).
    """
    if not two_sided:
        raise ValueError("This script implements only the two-sided stopping rule.")
    if B <= 0 or sigma <= 0 or eps <= 0 or grid_size < 3 or max_iter < 1 or tol <= 0:
        raise ValueError("Require positive B, sigma, eps, max_iter, tol and grid_size >= 3.")
    l_grid = np.linspace(-B, B, grid_size)
    interior = np.abs(l_grid) < B

    p_plus, delta_plus, delta_minus = llr_increments(l_grid, sigma, eps)

    next_plus = l_grid + delta_plus
    next_minus = l_grid + delta_minus

    H = np.zeros_like(l_grid)

    def interp_H(x, H_old):
        stopped = (x <= -B) | (x >= B)

        x_clip = np.clip(x, l_grid[0], l_grid[-1])
        val = np.interp(x_clip, l_grid, H_old)
        val[stopped] = 0.0
        return val

    for it in range(max_iter):
        H_old = H.copy()

        H_plus = interp_H(next_plus, H_old)
        H_minus = interp_H(next_minus, H_old)

        H_new = 1.0 + p_plus * H_plus + (1.0 - p_plus) * H_minus

        # Apply boundary condition
        H_new[~interior] = 0.0

        diff = np.max(np.abs(H_new - H_old))
        H = H_new

        if diff < tol:
            break
    else:
        raise RuntimeError(
            f"Bellman iteration did not converge for eps={eps:g}: "
            f"{max_iter} iterations, last change={diff:.3g}. "
            "Increase max_iter before using this result."
        )

    H0 = np.interp(0.0, l_grid, H)
    return H0, it + 1


# ============================================================
# Monte Carlo method
# ============================================================

def simulate_one_stopping_time(
    eps,
    B,
    sigma,
    theta_true,
    N_max,
    l_grid,
    delta_plus_grid,
    delta_minus_grid,
    two_sided=True
):
    """
    Simulate one trajectory and return stopping time.
    """
    if not two_sided:
        raise ValueError("This script implements only the two-sided stopping rule.")
    l = 0.0

    for n in range(1, N_max + 1):
        s = np.random.normal(loc=theta_true, scale=sigma)
        t = decision_threshold(l, sigma)

        a = 1 if s > t else -1

        flip_prob = 0.5 * np.exp(-eps * abs(s - t))
        x = -a if np.random.rand() < flip_prob else a

        l_clip = np.clip(l, l_grid[0], l_grid[-1])

        if x == 1:
            inc = np.interp(l_clip, l_grid, delta_plus_grid)
        else:
            inc = np.interp(l_clip, l_grid, delta_minus_grid)

        l += inc

        if abs(l) >= B:
            return n

    return np.nan


def mc_expected_stopping_time(
    eps,
    B,
    sigma,
    theta_true,
    n_runs,
    N_max,
    two_sided=True
):
    """
    Estimate E[tau_B | theta_true] using two-sided stopping.
    Trajectories are simulated in parallel using independent signals and flips.
    If any trajectory is unfinished at N_max, return NaN for the mean and
    report the hit rate: averaging completed trajectories would be biased.
    """
    if not two_sided:
        raise ValueError("This script implements only the two-sided stopping rule.")
    if n_runs < 1 or N_max < 1 or B <= 0 or sigma <= 0 or eps <= 0:
        raise ValueError("Require positive n_runs, N_max, B, sigma and eps.")
    if theta_true not in (-1, 1):
        raise ValueError("theta_true must be -1 or +1.")

    l_grid = np.linspace(-B, B, 1200)
    _, delta_plus_grid, delta_minus_grid = llr_increments(l_grid, sigma, eps)
    beliefs = np.zeros(n_runs)
    tau_values = np.full(n_runs, np.nan)
    active = np.ones(n_runs, dtype=bool)

    for n in range(1, N_max + 1):
        ids = np.flatnonzero(active)
        current = beliefs[ids]
        signals = np.random.normal(loc=theta_true, scale=sigma, size=ids.size)
        thresholds = decision_threshold(current, sigma)
        actions = np.where(signals > thresholds, 1, -1)
        flip_prob = 0.5 * np.exp(-eps * np.abs(signals - thresholds))
        reports = np.where(np.random.random(ids.size) < flip_prob, -actions, actions)
        increments = np.where(
            reports == 1,
            np.interp(current, l_grid, delta_plus_grid),
            np.interp(current, l_grid, delta_minus_grid),
        )
        beliefs[ids] += increments
        stopped_ids = ids[np.abs(beliefs[ids]) >= B]
        tau_values[stopped_ids] = n
        active[stopped_ids] = False
        if not np.any(active):
            break

    hit_rate = float(np.mean(~active))
    if np.any(active):
        warnings.warn(
            f"eps={eps:g}: {np.count_nonzero(active)}/{n_runs} trajectories "
            f"did not stop within {N_max} reports. Increase N_max and rerun; "
            "the unconditional MC mean is not reported.",
            RuntimeWarning,
        )
        return np.nan, hit_rate
    return float(np.mean(tau_values)), hit_rate


# ============================================================
# Run all methods
# ============================================================

def main():
    np.random.seed(123)
    bellman_vals = []
    bellman_iters = []

    mc_vals = []
    mc_hit_rates = []

    asym_vals = []

    for eps in eps_grid:
        print(f"epsilon = {eps:.3f}")

        H0, n_iter = bellman_expected_stopping_time(
            eps=eps,
            B=B,
            sigma=sigma,
            grid_size=bellman_grid_size,
            max_iter=bellman_max_iter,
            tol=bellman_tol,
            two_sided=two_sided_stop
        )

        mean_tau_mc, hit_rate = mc_expected_stopping_time(
            eps=eps,
            B=B,
            sigma=sigma,
            theta_true=theta_true,
            n_runs=n_runs,
            N_max=N_max,
            two_sided=two_sided_stop
        )

        tau_asym = asymptotic_stopping_time(eps, B, sigma)

        bellman_vals.append(H0)
        bellman_iters.append(n_iter)

        mc_vals.append(mean_tau_mc)
        mc_hit_rates.append(hit_rate)

        asym_vals.append(tau_asym)

    bellman_vals = np.array(bellman_vals)
    mc_vals = np.array(mc_vals)
    asym_vals = np.array(asym_vals)
    mc_hit_rates = np.array(mc_hit_rates)


    # ============================================================
    # Plot
    # ============================================================

    plt.figure(figsize=(9, 5.8))

    plt.plot(
        eps_grid,
        bellman_vals,
        color="black",
        linewidth=2.5,
        marker="o",
        markersize=4,
        label="Bellman equation"
    )

    plt.plot(
        eps_grid,
        mc_vals,
        color="blue",
        linewidth=2.2,
        marker="s",
        markersize=4,
        label="Monte Carlo simulation"
    )

    plt.plot(
        eps_grid,
        asym_vals,
        color="red",
        linewidth=2.2,
        linestyle="--",
        label=r"Asymptotic approximation (with $K$)"
    )

    plt.xlabel(r"Privacy budget $\varepsilon$", fontsize=16)

    plt.ylabel(r"Expected belief-threshold stopping time $\mathbb{E}[\tau_B]$", fontsize=13)

    plt.title(
        fr"Belief-Threshold Stopping Time vs. Privacy Budget, $B={B}$, $\sigma={sigma}$",
        fontsize=16
    )

    plt.grid(alpha=0.25)
    plt.legend(fontsize=12)
    plt.tight_layout()
    os.makedirs("figures", exist_ok=True)
    plt.savefig("figures/continuous_report_stopping_time_comparison.pdf", bbox_inches="tight")
    plt.show()

    # ============================================================
    # Print summary
    # ============================================================

    print("\nSummary")
    print("=======")
    print("eps\tBellman\t\tMC\t\tAsymptotic\tMC hit rate\tBellman iters")
    for eps, b, m, a, h, it in zip(
        eps_grid,
        bellman_vals,
        mc_vals,
        asym_vals,
        mc_hit_rates,
        bellman_iters
    ):
        print(f"{eps:.3f}\t{b:.2f}\t\t{m:.2f}\t\t{a:.2f}\t\t{h:.6f}\t\t{it}")


if __name__ == "__main__":
    main()
