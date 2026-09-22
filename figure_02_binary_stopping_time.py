import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

OUTPUT_DIR = Path(__file__).resolve().parent / 'figures'
OUTPUT_DIR.mkdir(exist_ok=True)
import math
from matplotlib.lines import Line2D
import os

# ============================================================
# Global plot style: larger fonts for paper figures
# ============================================================

plt.rcParams.update({
    "font.size": 16,
    "axes.titlesize": 22,
    "axes.labelsize": 20,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 14,
    "legend.title_fontsize": 15,
    "lines.linewidth": 2.6,
})

# ============================================================
# Participation setting
# ============================================================

kappa = 5.0  # epsilon_i ~ Uniform[0, kappa]

def eta_eps(eps):
    """
    Participation probability:
        eta(eps) = P(epsilon_i >= eps)

    If epsilon_i ~ Uniform[0, kappa], then:
        eta(eps) = 1 - eps / kappa
    """
    if eps >= kappa:
        return np.nan
    return 1.0 - eps / kappa


# ============================================================
# Binary-signal model functions
# ============================================================

def u_eps(eps):
    return 1.0 / (1.0 + np.exp(eps))


def rho_eps(eps, p):
    u = u_eps(eps)
    numerator = (1 - p) * (1 - u) + p * u
    denominator = p * (1 - u) + (1 - p) * u
    return numerator / denominator


def mu_eps(eps, p):
    u = u_eps(eps)
    return (2 * p - 1) * (1 - 2 * u)


def k_eps(eps, p):
    """
    k(epsilon,p) = floor(log_{rho(epsilon,p)} ((1-p)/p)) + 1
    """
    r = rho_eps(eps, p)
    alpha = (1 - p) / p

    if r <= 0 or abs(r - 1) < 1e-10:
        return np.nan

    return math.floor(np.log(alpha) / np.log(r)) + 1


def expected_report_stopping_time(eps, p):
    """
    Submitted-report expected stopping time:

        E[tau_{eps,p}]
        =
        k / mu * (1 - rho^k) / (1 + rho^k)
    """
    k = k_eps(eps, p)

    if np.isnan(k):
        return np.nan, np.nan

    r = rho_eps(eps, p)
    mu = mu_eps(eps, p)

    if abs(mu) < 1e-10:
        return np.nan, k

    E_tau = (k / mu) * ((1 - r**k) / (1 + r**k))
    return E_tau, k


def expected_calendar_stopping_time(eps, p):
    """
    Calendar expected stopping time:

        E[T_{eps,p}] = E[tau_{eps,p}] / eta(eps)
    """
    E_tau, k = expected_report_stopping_time(eps, p)

    eta = eta_eps(eps)

    if np.isnan(eta) or eta <= 0:
        return np.nan, k

    return E_tau / eta, k


def expected_nonprivate_stopping_time(p):
    """
    Non-private benchmark:
    epsilon = infinity, u = 0, k = 2, rho = (1-p)/p, mu = 2p - 1
    """
    k = 2
    rho = (1 - p) / p
    mu = 2 * p - 1

    return (k / mu) * ((1 - rho**k) / (1 + rho**k))


# ============================================================
# Plot function
# ============================================================

def plot_binary_stopping_time(use_participation, output_name):
    """
    use_participation = False: plot submitted-report stopping time E[tau]
    use_participation = True: plot calendar stopping time E[T] = E[tau]/eta
    """

    p_values = [0.55, 0.7, 0.9]

    line_styles = {
        0.55: ":",
        0.7: "--",
        0.9: "-"
    }

    k_colors = {
        2: "red",
        3: "limegreen",
        4: "blue",
        5: "purple",
    }

    eps_grid = np.linspace(0.05, 4.5, 3000)

    fig, ax = plt.subplots(figsize=(10.5, 6.2))

    # ========================================================
    # Plot private stopping time
    # ========================================================

    for p in p_values:
        E_values = []
        k_values = []

        for eps in eps_grid:
            if use_participation:
                E_val, k = expected_calendar_stopping_time(eps, p)
            else:
                E_val, k = expected_report_stopping_time(eps, p)

            E_values.append(E_val)
            k_values.append(k)

        E_values = np.array(E_values)
        k_values = np.array(k_values)

        valid_mask = ~np.isnan(k_values) & ~np.isnan(E_values)
        valid_k_values = sorted(set(k_values[valid_mask].astype(int)))

        for k_val in valid_k_values:
            if k_val not in k_colors:
                continue

            mask = valid_mask & (k_values == k_val)
            idx = np.where(mask)[0]

            if len(idx) == 0:
                continue

            # Split continuous segments so curves do not connect across threshold jumps
            segments = np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)

            for seg in segments:
                ax.plot(
                    eps_grid[seg],
                    E_values[seg],
                    color=k_colors[k_val],
                    linestyle=line_styles[p],
                    linewidth=2.8
                )

    # ========================================================
    # Non-private benchmark
    # ========================================================

    for p in p_values:
        E_np = expected_nonprivate_stopping_time(p)

        ax.axhline(
            y=E_np,
            color="black",
            linestyle=line_styles[p],
            linewidth=2.6
        )

    # ========================================================
    # Axis labels and title
    # ========================================================

    ax.set_xlabel(r"Privacy budget $\varepsilon$")

    if use_participation:
        ax.set_ylabel(r"Expected calendar stopping time")
        ax.set_title("Calendar Stopping Time vs. Privacy Budget")
    else:
        ax.set_ylabel(r"Expected submitted-report stopping time")
        ax.set_title("Submitted-Report Stopping Time vs. Privacy Budget")

    ax.grid(alpha=0.25)

    # ========================================================
    # Separate legends
    # ========================================================

    k_handles = [
        Line2D([0], [0], color=color, lw=2.8, linestyle="-", label=fr"$k={k_val}$")
        for k_val, color in k_colors.items()
    ]

    nonprivate_handle = Line2D(
        [0], [0],
        color="black",
        lw=2.8,
        linestyle="-",
        label="Non-private"
    )

    legend1 = ax.legend(
        handles=k_handles + [nonprivate_handle],
        title="Threshold",
        loc="upper right",
        bbox_to_anchor=(1.0, 1.0),
        frameon=True
    )

    ax.add_artist(legend1)

    p_handles = [
        Line2D([0], [0], color="black", lw=2.8, linestyle=line_styles[p], label=fr"$p={p}$")
        for p in p_values
    ]

    ax.legend(
        handles=p_handles,
        title="Signal Accuracy",
        loc="upper right",
        bbox_to_anchor=(1.0, 0.62),
        frameon=True
    )

    fig.tight_layout()

    OUTPUT_DIR.mkdir(exist_ok=True)

    fig.savefig(OUTPUT_DIR / f"{output_name}.pdf", bbox_inches="tight", metadata={"Author": ""})
    fig.savefig(OUTPUT_DIR / f"{output_name}.png", dpi=300, bbox_inches="tight")

    plt.close()


# ============================================================
# Generate two figures
# ============================================================

plot_binary_stopping_time(
    use_participation=False,
    output_name="figure_02a"
)

plot_binary_stopping_time(
    use_participation=True,
    output_name="figure_02b"
)
