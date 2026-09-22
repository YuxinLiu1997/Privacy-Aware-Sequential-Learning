import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.integrate import quad


# ============================================================
# Basic parameters
# ============================================================

sigma = 1.0
epsilon = 3.0       # fixed mDP upper bound
B = 10.0             # belief threshold
n_grid = 801
tol = 1e-8
max_iter = 20000


# ============================================================
# Threshold and flip probability families
# ============================================================

def threshold(l, sigma):
    """
    Gaussian decision threshold.
    Initial recommendation is +1 iff s > t(l).
    """
    return -0.5 * sigma**2 * l


def flip_exponential(d, beta):
    """
    Exponential family:
        psi_beta(d) = 0.5 exp(- beta d)
    Feasible if beta <= epsilon.
    """
    return 0.5 * np.exp(-beta * d)


def flip_polynomial(d, alpha, gamma):
    """
    Polynomial-decay family:
        psi_{alpha,gamma}(d) = 0.5 (1 + alpha d)^(-gamma)

    Sufficient feasibility condition:
        alpha * gamma <= epsilon.
    """
    return 0.5 * (1.0 + alpha * d) ** (-gamma)


def feasible_polynomial(alpha, gamma, epsilon):
    return alpha * gamma <= epsilon + 1e-12


# ============================================================
# Transition probabilities and Bellman equation
# ============================================================

def p_x_plus_given_theta(l, theta, sigma, flip_func):
    """
    Compute P(x=+1 | l, theta) under a threshold-based flip mechanism.
    """
    t = threshold(l, sigma)

    def density(s):
        return norm.pdf(s, loc=theta, scale=sigma)

    def integrand_left(s):
        # s < t: initial recommendation is -1, report +1 if flipped
        d = abs(s - t)
        u = flip_func(d)
        return u * density(s)

    def integrand_right(s):
        # s > t: initial recommendation is +1, report +1 if not flipped
        d = abs(s - t)
        u = flip_func(d)
        return (1.0 - u) * density(s)

    left_val, _ = quad(integrand_left, -np.inf, t, limit=100)
    right_val, _ = quad(integrand_right, t, np.inf, limit=100)

    return left_val + right_val


def transition_at_l(l, sigma, flip_func):
    """
    Return P_+(+1|l), P_+(-1|l), Delta(+1), Delta(-1).
    """
    p_plus_pos = p_x_plus_given_theta(l, theta=+1, sigma=sigma, flip_func=flip_func)
    p_minus_pos = p_x_plus_given_theta(l, theta=-1, sigma=sigma, flip_func=flip_func)

    eps_num = 1e-14
    p_plus_pos = np.clip(p_plus_pos, eps_num, 1.0 - eps_num)
    p_minus_pos = np.clip(p_minus_pos, eps_num, 1.0 - eps_num)

    p_plus_neg = 1.0 - p_plus_pos
    p_minus_neg = 1.0 - p_minus_pos

    delta_pos = np.log(p_plus_pos / p_minus_pos)
    delta_neg = np.log(p_plus_neg / p_minus_neg)

    return p_plus_pos, p_plus_neg, delta_pos, delta_neg


def interp_H(l_new, grid, H, B):
    """
    Boundary condition: H(l)=0 if |l| >= B.
    Otherwise linearly interpolate.
    """
    if abs(l_new) >= B:
        return 0.0
    return np.interp(l_new, grid, H)


def expected_stopping_time_bellman(sigma, B, flip_func,
                                   n_grid=801, tol=1e-8, max_iter=20000):
    """
    Solve:
        H(l) = 1 + P_+(+1|l)H(l+Delta_+)
                 + P_+(-1|l)H(l+Delta_-)

    under true state theta=+1.
    Return H(0).
    """
    grid = np.linspace(-B, B, n_grid)

    H = np.full(n_grid, B**2)
    H[0] = 0.0
    H[-1] = 0.0

    # Precompute transitions
    trans = []
    for l in grid:
        if abs(l) >= B:
            trans.append(None)
        else:
            trans.append(transition_at_l(l, sigma, flip_func))

    for it in range(max_iter):
        H_new = H.copy()

        for j, l in enumerate(grid):
            if abs(l) >= B:
                H_new[j] = 0.0
                continue

            p_pos, p_neg, d_pos, d_neg = trans[j]

            H_new[j] = (
                1.0
                + p_pos * interp_H(l + d_pos, grid, H, B)
                + p_neg * interp_H(l + d_neg, grid, H, B)
            )

        diff = np.max(np.abs(H_new - H))
        H = H_new

        if diff < tol:
            break

    return np.interp(0.0, grid, H), it + 1


# ============================================================
# Figure 1: Exponential family vs beta
# ============================================================

beta_values = np.linspace(0.05, epsilon, 20)
ET_exp = []

print("===== Exponential family: vary beta =====")
for beta in beta_values:
    flip_func = lambda d, beta=beta: flip_exponential(d, beta)
    ET, n_iter = expected_stopping_time_bellman(
        sigma=sigma,
        B=B,
        flip_func=flip_func,
        n_grid=n_grid,
        tol=tol,
        max_iter=max_iter
    )
    ET_exp.append(ET)
    print(f"beta={beta:.4f}, E[T_B]={ET:.4f}, iterations={n_iter}")

ET_exp = np.array(ET_exp)

plt.figure(figsize=(7, 4.5))
plt.plot(beta_values, ET_exp, marker="o")
plt.xlabel(r"$\beta$")
plt.ylabel(r"Expected stopping time $\mathbb{E}[T_B]$")
plt.title(r"Exponential family: $\psi_\beta(d)=\frac{1}{2}e^{-\beta d}$")
plt.tight_layout()
plt.show()


# ============================================================
# Figure 2: Polynomial family vs alpha, fixing gamma
# ============================================================

alpha_values = np.arange(0.05, 2.05, 0.05)
fixed_gammas = [0.5, 1.0, 2.0, 3.0, 4.0]

plt.figure(figsize=(7, 4.5))

print("\n===== Polynomial family: vary alpha, fixed gamma =====")
for gamma in fixed_gammas:
    xs, ys = [], []

    for alpha in alpha_values:
        if not feasible_polynomial(alpha, gamma, epsilon):
            continue

        flip_func = lambda d, alpha=alpha, gamma=gamma: flip_polynomial(d, alpha, gamma)
        ET, n_iter = expected_stopping_time_bellman(
            sigma=sigma,
            B=B,
            flip_func=flip_func,
            n_grid=n_grid,
            tol=tol,
            max_iter=max_iter
        )

        xs.append(alpha)
        ys.append(ET)
        print(f"gamma={gamma:.2f}, alpha={alpha:.2f}, alpha*gamma={alpha*gamma:.3f}, E[T_B]={ET:.4f}")

    if xs:
        plt.plot(xs, ys, marker="o", label=fr"$\gamma={gamma}$")

plt.xlabel(r"$\alpha$")
plt.ylabel(r"Expected stopping time $\mathbb{E}[T_B]$")
plt.title(r"Polynomial family vs. $\alpha$: $\psi_{\alpha,\gamma}(d)=\frac{1}{2}(1+\alpha d)^{-\gamma}$")
plt.legend()
plt.tight_layout()
plt.show()


# ============================================================
# Figure 3: Polynomial family vs gamma, fixing alpha
# ============================================================

gamma_values = np.arange(0.10, 5.10, 0.10)
fixed_alphas = [0.2, 0.4, 0.6, 0.8, 1.0]

plt.figure(figsize=(7, 4.5))

print("\n===== Polynomial family: vary gamma, fixed alpha =====")
for alpha in fixed_alphas:
    xs, ys = [], []

    for gamma in gamma_values:
        if not feasible_polynomial(alpha, gamma, epsilon):
            continue

        flip_func = lambda d, alpha=alpha, gamma=gamma: flip_polynomial(d, alpha, gamma)
        ET, n_iter = expected_stopping_time_bellman(
            sigma=sigma,
            B=B,
            flip_func=flip_func,
            n_grid=n_grid,
            tol=tol,
            max_iter=max_iter
        )

        xs.append(gamma)
        ys.append(ET)
        print(f"alpha={alpha:.2f}, gamma={gamma:.2f}, alpha*gamma={alpha*gamma:.3f}, E[T_B]={ET:.4f}")

    if xs:
        plt.plot(xs, ys, marker="o", label=fr"$\alpha={alpha}$")

plt.xlabel(r"$\gamma$")
plt.ylabel(r"Expected stopping time $\mathbb{E}[T_B]$")
plt.title(r"Polynomial family vs. $\gamma$: $\psi_{\alpha,\gamma}(d)=\frac{1}{2}(1+\alpha d)^{-\gamma}$")
plt.legend()
plt.tight_layout()
plt.show()