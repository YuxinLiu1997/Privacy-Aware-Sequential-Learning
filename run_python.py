"""Run the manuscript's Python figures from any working directory."""
import argparse
import importlib.util
import os
from pathlib import Path
import runpy

ROOT = Path(__file__).resolve().parent
TARGETS = {
    '2': 'figure_02_binary_stopping_time.py',
    '3': 'figure_03_binary_objective.py',
    '6a': 'figure_06a_continuous_report_time.py',
    '6b': 'figure_06b_continuous_calendar_time.py',
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('figures', nargs='+', choices=list(TARGETS) + ['all'])
    parser.add_argument('--quick', action='store_true', help='Reduced continuous settings; not manuscript reproduction')
    args = parser.parse_args()
    os.environ['MPLBACKEND'] = 'Agg'
    import matplotlib.pyplot as plt
    import numpy as np
    # Original scripts use show(); save all results without opening GUI windows.
    plt.show = lambda: None
    output = ROOT / 'outputs' / ('quick' if args.quick else 'full') / 'python'
    output.mkdir(parents=True, exist_ok=True)
    os.chdir(output)
    figures = list(TARGETS) if 'all' in args.figures else args.figures
    for key in figures:
        path = ROOT / 'python' / TARGETS[key]
        print(f'Figure {key}: {path.name}', flush=True)
        if key.startswith('6') and args.quick:
            spec = importlib.util.spec_from_file_location('figure_' + key, path)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            module.B = 2.0
            module.eps_grid = np.array([0.5, 1.0])
            module.n_runs = 100
            module.N_max = 5000
            module.bellman_grid_size = 81
            module.main()
        else:
            runpy.run_path(str(path), run_name='__main__')
        plt.close('all')
    print(f'Output directory: {output}')


if __name__ == '__main__':
    main()
