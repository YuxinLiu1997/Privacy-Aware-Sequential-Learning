"""Figure J.1: deterministic positive-report recurrence and its approximation."""
from pathlib import Path
import numpy as np
from scipy.special import ndtr, log_ndtr, expit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def probability_minus(belief, state, epsilon, sigma=1.0):
    threshold = -sigma*sigma*belief/2
    z = (threshold-state)/sigma
    a = expit(-epsilon)
    left = np.exp(epsilon*(state-threshold)+0.5*(epsilon*sigma)**2 + log_ndtr(z-epsilon*sigma))
    right = np.exp(epsilon*(threshold-state)+0.5*(epsilon*sigma)**2 + log_ndtr(-z-epsilon*sigma))
    return ndtr(z)-a*left+a*right


def trajectory(epsilon, n=200, sigma=1.0):
    values = np.empty(n)
    values[0] = 6
    for i in range(1, n):
        l = values[i-1]
        values[i] = l + np.log1p(-probability_minus(l, 1, epsilon, sigma)) - np.log1p(-probability_minus(l, -1, epsilon, sigma))
    return values


def approximate(n=200, epsilon=1.0, sigma=1.0):
    values = np.empty(n)
    values[0] = 6
    c = expit(-epsilon) * (np.exp(epsilon+epsilon**2*sigma**2/2)-np.exp(-epsilon+epsilon**2*sigma**2/2))
    for i in range(1, n):
        values[i] = values[i-1] + c*np.exp(-epsilon*sigma**2*values[i-1]/2)
    return values


def main():
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    x = np.arange(1, 201)
    for epsilon, color in [(0.1, 'blue'), (0.5, 'red'), (1, 'lime'), (100, 'black')]:
        axes[0].plot(x, trajectory(epsilon), color=color, ls='--' if epsilon == 100 else '-',
                     label='Non-private' if epsilon == 100 else rf'$\varepsilon={epsilon}$', lw=2)
    axes[1].plot(x, trajectory(1), color='blue', lw=2, label='Exact')
    axes[1].plot(x, approximate(), color='red', lw=2, label='Approximate')
    for i, ax in enumerate(axes):
        ax.set(xlabel='Step (n)', ylabel=r'$l_n$')
        ax.set_title(f'({chr(97+i)})', loc='left')
        ax.legend(loc='upper left', frameon=False)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    out = Path(__file__).resolve().parent / 'figures'
    out.mkdir(exist_ok=True)
    for ext in ['pdf', 'png']:
        fig.savefig(out / f'figure_J1.{ext}', dpi=200, metadata={'Author': ''} if ext == 'pdf' else None)
    plt.close(fig)


if __name__ == '__main__':
    main()
