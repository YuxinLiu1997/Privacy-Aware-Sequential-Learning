import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

OUTPUT_DIR = Path(__file__).resolve().parent / 'figures'
OUTPUT_DIR.mkdir(exist_ok=True)
import math

# ============================================================
# Settings
# ============================================================

p_values = [0.55, 0.7, 0.9]  # manuscript Figure 3
c = 0.01                 # stopping cost
use_participation = True # True: E[T] = E[tau]/eta(eps), False: E[tau]

# If use_participation = True, assume epsilon_i ~ Uniform[0, kappa]
kappa = 5.0

eps_grid = np.linspace(0.1, 4.95, 3000)


# ============================================================
# Participation model
# ============================================================

def eta_eps(eps):
    """
    eta(eps) = P(epsilon_i >= eps).
    If epsilon_i ~ Uniform[0, kappa], eta(eps)=1-eps/kappa.
    """
    if eps >= kappa:
        return np.nan
    return 1.0 - eps / kappa


# ============================================================
# Binary model functions
# ============================================================

def u_eps(eps):
    return 1.0 / (1.0 + np.exp(eps))


def rho_eps(eps, p):
    u = u_eps(eps)
    numerator = (1 - u) * (1 - p) + u * p
    denominator = u * (1 - p) + p * (1 - u)
    return numerator / denominator


def mu_eps(eps, p):
    u = u_eps(eps)
    return (2 * p - 1) * (1 - 2 * u)


def k_eps(eps, p):
    """
    k(eps,p) = floor(log_{rho(eps,p)}((1-p)/p)) + 1
    """
    r = rho_eps(eps, p)
    alpha = (1 - p) / p

    if r <= 0 or abs(r - 1) < 1e-10:
        return np.nan

    return math.floor(np.log(alpha) / np.log(r)) + 1


def prob_correct_final_decision(eps, p):
    """
    P_correct = (rho^k - 1) / (rho^(2k) - 1)
    """
    k = k_eps(eps, p)

    if np.isnan(k):
        return np.nan, np.nan

    r = rho_eps(eps, p)
    denominator = r ** (2 * k) - 1

    if abs(denominator) < 1e-12:
        return np.nan, k

    P_correct = (r ** k - 1) / denominator
    return P_correct, k


def expected_report_stopping_time(eps, p):
    """
    E[tau_{eps,p}] = k/mu * (1-rho^k)/(1+rho^k)
    """
    k = k_eps(eps, p)

    if np.isnan(k):
        return np.nan, np.nan

    r = rho_eps(eps, p)
    mu = mu_eps(eps, p)

    if abs(mu) < 1e-10:
        return np.nan, k

    E_tau = (k / mu) * ((1 - r ** k) / (1 + r ** k))
    return E_tau, k


def expected_stopping_time(eps, p):
    """
    If use_participation=True:
        E[T_{eps,p}] = E[tau_{eps,p}] / eta(eps)
    Otherwise:
        E[T_{eps,p}] = E[tau_{eps,p}]
    """
    E_tau, k = expected_report_stopping_time(eps, p)

    if np.isnan(E_tau):
        return np.nan, k

    if not use_participation:
        return E_tau, k

    eta = eta_eps(eps)

    if np.isnan(eta) or eta <= 0:
        return np.nan, k

    return E_tau / eta, k


def objective_value(eps, p, c):
    """
    J(eps) = P_correct(eps,p) - c * E[T_{eps,p}]
    """
    P_correct, k = prob_correct_final_decision(eps, p)
    E_T, _ = expected_stopping_time(eps, p)

    if np.isnan(P_correct) or np.isnan(E_T):
        return np.nan, np.nan, np.nan, np.nan

    J = P_correct - c * E_T
    return J, P_correct, E_T, k


# ============================================================
# Compute objective
# ============================================================

objective_curves = {
    p: np.array([objective_value(eps, p, c)[0] for eps in eps_grid])
    for p in p_values
}


# ============================================================
# Plot
# ============================================================

plt.figure(figsize=(10, 5.8))

for p, style in zip(p_values, ["-", "--", ":"]):
    plt.plot(eps_grid, objective_curves[p], linewidth=2,
             linestyle=style, label=f"p = {p}")

plt.xlabel(r"Privacy budget $\varepsilon$", fontsize=16)

if use_participation:
    plt.ylabel(
        r"Objective value",
        fontsize=14
    )
    plt.title(
        fr"Platform Objective vs. Privacy Budget, $c={c}$, $\varepsilon_i\sim U[0,{kappa}]$",
        fontsize=16
    )
else:
    plt.ylabel(
        r"Objective value ",
        fontsize=14
    )
    plt.title(
        fr"Platform Objective vs. Privacy Budget, $c={c}$",
        fontsize=16
    )

plt.grid(alpha=0.25)
plt.legend()
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "figure_03.pdf", bbox_inches="tight", metadata={"Author": ""})
plt.savefig(OUTPUT_DIR / "figure_03.png", bbox_inches="tight", dpi=200)
plt.close()
