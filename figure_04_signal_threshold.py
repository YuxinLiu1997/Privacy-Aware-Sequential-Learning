"""Figure 4: reconstruct the signal/threshold schematic shown in the paper."""
from pathlib import Path
import numpy as np
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    fig, ax = plt.subplots(figsize=(8, 4.6))
    x = np.linspace(-5, 7, 1500)
    threshold, mean, sigma = -1.6, 1.0, 1.8
    density = norm.pdf(x, loc=mean, scale=sigma)
    # The vertical density scale is schematic, as in the paper illustration.
    density *= 0.85 / norm.pdf(mean, loc=mean, scale=sigma)
    blue, orange = '#006699', 'orangered'
    ax.plot(x, density, color=orange, lw=1)
    ax.plot([-6, threshold], [0, 0], color=blue, lw=1)
    ax.plot([threshold, 6], [1, 1], color=blue, lw=1)
    ax.vlines(threshold, 0, 1, colors=blue, linestyles=':', lw=0.8)
    ax.vlines(mean, 0, 0.85, colors=orange, linestyles=':', lw=0.8)
    ax.annotate('', xy=(7, 0), xytext=(-6, 0), arrowprops={'arrowstyle': '->', 'lw': 0.7})
    ax.annotate('', xy=(0, 1.35), xytext=(0, 0), arrowprops={'arrowstyle': '->', 'lw': 0.7})
    ax.text(-0.1, 1.02, '1', color=blue, ha='right')
    ax.text(1.3, 1.08, r'$P(a_n=+1\mid s_n,\theta=+1)$', color=blue)
    ax.text(3, 0.76, r'$f(s_n\mid\theta=+1)$', color=orange)
    ax.text(threshold-0.45, 0.07, r'$s_n$', ha='center')
    ax.text(threshold+0.48, 0.07, r"$s_n'$", ha='center')
    ax.vlines([threshold-0.45, threshold+0.48], 0, 0.025, colors='black', lw=0.6)
    ax.text(threshold, -0.065, 'threshold', ha='center')
    ax.text(mean, -0.065, '+1', ha='center')
    ax.text(7, -0.065, r'$s_n$', ha='center')
    ax.set(xlim=(-6.2, 7.4), ylim=(-0.16, 1.42))
    ax.axis('off')
    fig.tight_layout()
    out = Path(__file__).resolve().parent / 'figures'
    out.mkdir(exist_ok=True)
    for ext in ['pdf', 'png']:
        fig.savefig(out / f'figure_04.{ext}', dpi=200, metadata={'Author': ''} if ext == 'pdf' else None)
    plt.close(fig)


if __name__ == '__main__':
    main()
