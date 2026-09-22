"""Figure 7: five privacy regimes, retaining the original seeded trajectories."""
from pathlib import Path
import argparse
import numpy as np
from scipy.special import ndtr, ndtri, log_ndtr
from numpy.polynomial.legendre import leggauss
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def random_inputs(seed, n, regime):
    """Reproduce the original 32-bit uniform stream and inversion-normal draws.

    Initialize NumPy's MT19937 state using the original seed recurrence. Each
    normal consumes two uniforms; retain the signal/budget/flip draw order.
    No external interpreter or saved random sample files are required.
    """
    seed = int(seed) & 0xffffffff
    state = []
    for i in range(50 + 625):
        seed = (69069 * seed + 1) & 0xffffffff
        if i >= 51:
            state.append(seed)
    generator = np.random.MT19937()
    generator.state = {'bit_generator': 'MT19937',
                       'state': {'key': np.array(state, dtype=np.uint32), 'pos': 624}}
    width = {'uniform': 4, 'fixed': 3, 'nonprivate': 2}[regime]
    u = generator.random_raw(n*width).reshape(n, width).astype(float) / 4294967296.0
    u = np.maximum(u, 0.5 / 4294967295.0)
    normal = ndtri((np.floor(134217728*u[:, 0])+u[:, 1])/134217728)
    return normal, u


def probability_minus(belief, sigma, epsilon, state):
    threshold = -sigma*sigma*belief/2
    z = (threshold-state)/sigma
    left = np.exp(epsilon*(state-threshold)+0.5*(epsilon*sigma)**2 + log_ndtr(z-epsilon*sigma))
    right = np.exp(epsilon*(threshold-state)+0.5*(epsilon*sigma)**2 + log_ndtr(-z-epsilon*sigma))
    return ndtr(z) - 0.5*left + 0.5*right


NODES, WEIGHTS = leggauss(64)
NODES, WEIGHTS = (NODES+1)/2, WEIGHTS/2


def simulate(n, seed, regime, epsilon=None, sigma=1.0):
    normals, uniforms = random_inputs(seed, n, regime)
    signals = 1 + sigma*normals
    values = np.empty(n)
    belief = 0.0
    states = np.array([1.0, -1.0])
    for i, signal in enumerate(signals):
        values[i] = belief  # Match l_n before the nth report, including l_1=0.
        threshold = -sigma*sigma*belief/2
        report = 1 if signal > threshold else -1
        if regime == 'nonprivate':
            probabilities = ndtr((threshold-states)/sigma)
        else:
            budget = uniforms[i, 2] if regime == 'uniform' else epsilon
            if uniforms[i, -1] < 0.5*np.exp(-budget*abs(signal-threshold)):
                report = -report
            if regime == 'uniform':
                probabilities = WEIGHTS @ probability_minus(belief, sigma, NODES[:, None], states)
            else:
                probabilities = probability_minus(belief, sigma, epsilon, states)
        if report == 1:
            increment = np.log1p(-probabilities[0]) - np.log1p(-probabilities[1])
        else:
            increment = np.log(probabilities[0]) - np.log(probabilities[1])
        belief += increment
    return values


GROUPS = [
    ('uniform', None, 6, 'purple', '-', 'Uniform[0,1]'),
    ('fixed', 0.1, 16, 'red', '--', 'Low'),
    ('fixed', 0.5, 36, 'blue', ':', 'Medium'),
    ('fixed', 1.0, 46, 'darkgreen', '-.', 'High'),
    ('nonprivate', None, 56, 'black', (0, (7, 3)), 'Non-private'),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true', help='Use 100 steps and two paths; not the paper figure.')
    args = parser.parse_args()
    n, runs = (100, 2) if args.quick else (10000, 5)
    fig, ax = plt.subplots(figsize=(9, 6))
    handles = []
    for regime, epsilon, base_seed, color, style, label in GROUPS:
        print(f'{label}: {runs} paths, {n} steps', flush=True)
        paths = np.array([simulate(n, base_seed+i, regime, epsilon) for i in range(1, runs+1)])
        row_means = paths.mean(axis=1)
        center = np.argmin(np.abs(row_means-row_means.mean()))
        for i, path in enumerate(paths):
            ax.plot(np.arange(1, n+1), path, color=color, ls=style,
                    lw=2.5 if i == center else 1, alpha=1 if i == center else 0.3)
        handles.append(Line2D([], [], color=color, ls=style, lw=2, label=label))
    # Preserve the original plot's 4% axis expansion.
    ax.set(xlim=(1-0.04*(n-1), n+0.04*(n-1)), ylim=(-9.42, 42.42),
           xlabel=r'$n$', ylabel=r'$l_n$')
    ax.grid(alpha=0.3)
    ax.legend(handles=handles, loc='lower right', frameon=False)
    fig.tight_layout()
    out = Path(__file__).resolve().parent / 'figures'
    out.mkdir(exist_ok=True)
    stem = 'figure_07_quick' if args.quick else 'figure_07'
    for ext in ['pdf', 'png']:
        fig.savefig(out / f'{stem}.{ext}', dpi=200, metadata={'Author': ''} if ext == 'pdf' else None)
    plt.close(fig)


if __name__ == '__main__':
    main()
