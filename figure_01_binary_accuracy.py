"""Figure 1: probability of a correct final decision in the binary model."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def transition(k, p):
    if k == 2:
        return np.inf
    alpha = (1 - p) / p
    power = alpha ** ((k - 2) / (k - 1))
    w = (1 - power) / (1 - power + alpha ** (-1 / (k - 1)) - alpha)
    return np.log((1 - w) / w)


def accuracy(epsilon, k, p):
    u = 1 / (1 + np.exp(epsilon))
    q = u * (1 - p) + p * (1 - u)
    rho = (1 - q) / q
    return 1 / (1 + rho ** k)


def main():
    fig, ax = plt.subplots(figsize=(9, 6))
    grid = np.linspace(0, 5, 1000)
    colors = ['red', 'lime', 'blue', 'purple']
    styles = ['-', '--', ':']
    for p, style in zip([0.9, 0.7, 0.55], styles):
        for k, color in zip(range(2, 6), colors):
            lower, upper = transition(k + 1, p), transition(k, p)
            x = grid[(grid > lower) & (grid <= upper)]
            ax.plot(x, accuracy(x, k, p), color=color, ls=style, lw=2)
            ax.axvline(lower, color='gray', ls='--', lw=0.8)
        ax.axhline(p*p/(2*p*p-2*p+1), color='black', ls=style, lw=2)
    first = ax.legend([Line2D([], [], color=c, lw=2) for c in ['black'] + colors],
                      ['Non-private'] + [f'k = {k}' for k in range(2, 6)],
                      loc='center right', frameon=False)
    ax.add_artist(first)
    ax.legend([Line2D([], [], color='black', ls=s, lw=2) for s in styles],
              ['p = 0.9', 'p = 0.7', 'p = 0.55'], loc='lower right', frameon=False)
    ax.set(xlim=(0, 5), ylim=(0.5, 1), xlabel=r'$\varepsilon$',
           ylabel='Probability of Correct Cascade')
    fig.tight_layout()
    out = Path(__file__).resolve().parent / 'figures'
    out.mkdir(exist_ok=True)
    for ext in ['pdf', 'png']:
        fig.savefig(out / f'figure_01.{ext}', dpi=200, metadata={'Author': ''} if ext == 'pdf' else None)
    plt.close(fig)


if __name__ == '__main__':
    main()
