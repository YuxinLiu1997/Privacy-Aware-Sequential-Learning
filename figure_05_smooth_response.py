"""Figure 5: signal distributions, smooth response and flipped-report densities."""
from pathlib import Path
import numpy as np
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def curves():
    x = np.linspace(-5, 5, 1001)
    negative, positive = norm.pdf(x, loc=-1), norm.pdf(x, loc=1)
    flip = 0.5 * np.exp(-np.abs(x + 3))
    return x, negative, positive, flip


def main():
    x, negative, positive, flip = curves()
    selected = x >= -3
    fig, axes = plt.subplots(1, 3, figsize=(12, 5))
    for i, ax in enumerate(axes):
        ax.set_title(f'({chr(97+i)})', loc='left')
        ax.axvline(-3, color='black', ls='--', lw=1)
        ax.set(xlabel=r'$s_n$', ylabel='Density', ylim=(0, 0.1 if i == 2 else 0.5))
        ax.grid(alpha=0.2)
        for spine in ax.spines.values():
            spine.set_visible(False)
        for values, color in [(negative, 'blue'), (positive, 'orange')]:
            if i == 2:
                ax.plot(x[selected], (values*flip)[selected], color=color)
                ax.fill_between(x[selected], 0, (values*flip)[selected], color=color, alpha=0.2)
            else:
                ax.plot(x, values, color=color)
                if i == 0:
                    ax.fill_between(x[selected], 0, values[selected], color=color, alpha=0.2)
        if i == 1:
            ax.plot(x, flip, color='black')
    fig.tight_layout()
    out = Path(__file__).resolve().parent / 'figures'
    out.mkdir(exist_ok=True)
    for ext in ['pdf', 'png']:
        fig.savefig(out / f'figure_05.{ext}', dpi=200, metadata={'Author': ''} if ext == 'pdf' else None)
    plt.close(fig)


if __name__ == '__main__':
    main()
