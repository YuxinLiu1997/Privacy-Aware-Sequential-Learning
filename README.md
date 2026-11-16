# Figure Reproduction Code

This repository contains the Python scripts and dependency list needed to reproduce the figures in the paper. No external data are required.

## Figure-to-Script Mapping

| Paper figure | Python file | Description | Output files (in `figures/`) |
| --- | --- | --- | --- |
| Figure 1 | `figure_01_binary_accuracy.py` | Binary signals: probability of a correct final decision versus privacy budget | `figure_01.pdf/png` |
| Figures 2(a) and 2(b) | `figure_02_binary_stopping_time.py` | Binary signals: stopping times measured in reports and in calendar time | `figure_02a.pdf/png`, `figure_02b.pdf/png` |
| Figure 3 | `figure_03_binary_objective.py` | Binary signals: the platform objective function at three signal accuracy levels | `figure_03.pdf/png` |
| Figure 4 | `figure_04_signal_threshold.py` | Illustration of the step decision function and Gaussian signal distribution near the decision threshold | `figure_04.pdf/png` |
| Figure 5 | `figure_05_smooth_response.py` | Three panels showing signal densities, smooth randomized response, and action-flip probability | `figure_05.pdf/png` |
| Figure 6(a) | `figure_06a_continuous_report_time.py` | Continuous signals: stopping time measured in reports, comparing dynamic programming, simulation, and an asymptotic approximation | `figure_06a.pdf/png` |
| Figure 6(b) | `figure_06b_continuous_calendar_time.py` | Continuous signals: calendar-time stopping time accounting for the participation rate | `figure_06b.pdf/png` |
| Figure 7 | `figure_07_privacy_regimes.py` | Log-likelihood ratio trajectories under five privacy regimes | `figure_07.pdf/png` |
| Appendix Figure J.1 | `figure_J1_likelihood_evolution.py` | Recurrence trajectories for successive positive reports and a comparison of exact and asymptotic recursions | `figure_J1.pdf/png` |

## Installation and Usage

Tested environment: Python 3.13. See `requirements.txt` for dependency versions. A separate virtual environment is recommended.

```sh
python -m pip install -r requirements.txt
python figure_01_binary_accuracy.py
python figure_02_binary_stopping_time.py
python figure_03_binary_objective.py
python figure_04_signal_threshold.py
python figure_05_smooth_response.py
python figure_06a_continuous_report_time.py
python figure_06b_continuous_calendar_time.py
python figure_07_privacy_regimes.py
python figure_J1_likelihood_evolution.py
```

Each script can be run independently. Output is always saved to the `figures/` directory alongside the scripts. The scripts automatically save both PDF and PNG files without opening plot windows. The two panels of Figures 2 and 6 are saved separately; Figure 5 and Appendix Figure J.1 each combine their panels in a single figure. The code archive and repository contain source code only, without pregenerated figures.

## Parameter Settings

- **Figure 1**: Signal accuracy `p = 0.55, 0.7, 0.9`; cascade threshold `k = 2, 3, 4, 5`; 1,000 equally spaced privacy budgets on `[0, 5]`.
- **Figure 2**: Signal accuracy `p = 0.55, 0.7, 0.9`; 3,000 equally spaced privacy budgets on `[0.05, 4.5]`.
  Figure 2(a) computes the stopping time measured in reports. Figure 2(b) uses privacy tolerances uniformly distributed on `[0, 5]`, giving a participation rate of `1 - epsilon/5` over this budget range.
- **Figure 3**: Signal accuracy `p = 0.55, 0.7, 0.9`; waiting cost `c = 0.01`; privacy tolerances uniformly distributed on `[0, 5]`; 3,000 equally spaced privacy budgets on `[0.1, 4.95]`.
- **Figure 4**: The illustrative Gaussian curve has mean `1` and standard deviation `1.8`, with decision threshold `-1.6`.
  The curve's peak is scaled to `0.85`, and the step function has height `1`. These parameters are used for the illustration's layout and are not simulation estimates.
- **Figure 5**: The Gaussian signals in the two states have means `-1` and `1`, with standard deviation `sigma = 1`; privacy budget `epsilon = 1`; decision threshold `-3`; signal grid `[-5, 5]` with step size `0.01`; flip probability `0.5 * exp(-abs(s + 3))`.
- **Figure 6**: Belief threshold `B = 10`, signal standard deviation `sigma = 1`, and true state `+1`; 25 equally spaced privacy budgets on `[0.1, 2]`.
  Each budget uses 300 simulated paths, a maximum of 50,000 reports, and random seed `123`. Dynamic programming uses a 401-point grid, at most 100,000 iterations, and convergence tolerance `1e-8`. Figure 6(b) uses privacy tolerances uniformly distributed on `[0, 5]`.
- **Figure 7**: All five simulation groups use signal standard deviation **`sigma = 1`**, true state `+1`, and initial log-likelihood ratio `0`, with 5 paths per group and 10,000 steps per path.
  Three groups use fixed privacy budgets of `0.1`, `0.5`, and `1`; the heterogeneous group uses budgets uniformly distributed on `[0, 1]`; the fifth group is the nonprivate benchmark. The randomized-response amplitude is `1/2`.
  Seeds are `7–11` for the heterogeneous group; `17–21`, `37–41`, and `47–51` for the three fixed-budget groups, respectively; and `57–61` for the nonprivate group. Within each group, the highlighted path is the one whose time average is closest to the mean time average across the group's paths.
  The `--quick` option uses 100 steps and 2 paths per group, leaves the other parameters unchanged, and saves separate output as `figure_07_quick`.
- **Appendix Figure J.1**: Signal standard deviation `sigma = 1`, initial log-likelihood ratio `6`, and recursion length `200`; privacy budgets `0.1`, `0.5`, and `1`. The nonprivate benchmark uses `epsilon = 100`; the randomized-response amplitude is `1/(1 + exp(epsilon))`.
  The right panel compares the exact and asymptotic recursions at `epsilon = 1`.
