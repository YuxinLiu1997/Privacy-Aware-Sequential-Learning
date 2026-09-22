# Privacy-Aware Sequential Learning

Code accompanying **Privacy-Aware Sequential Learning**, organized against the supplied
MOR submission `privacy_aware_sequential_learning_MOR_submission (6).pdf`.
The manuscript itself is not included. Figure numbers below refer to this submission,
not the older filenames in this repository's history.

## Manuscript figures

| Figure | Entry point | Content / status |
| --- | --- | --- |
| 1 | `r/figure_01_binary_accuracy.R` | Binary final-decision accuracy; p = 0.55, 0.7, 0.9. |
| 2 | `python/figure_02_binary_stopping_time.py` | Submitted-report and calendar stopping times; two separate panels. |
| 3 | `python/figure_03_binary_objective.py` | Binary platform objective; c = 0.01, tolerances U[0,5]. |
| 4 | No standalone source located | Threshold/Gaussian schematic. `reference_figures/normal_distribution.pdf` is a candidate existing asset, not a verified reproduction. |
| 5 | `r/figure_05_smooth_randomized_response.R` | Smooth randomized response and signal densities. |
| 6(a) | `python/figure_06a_continuous_report_time.py` | Continuous submitted-report time: Bellman, Monte Carlo and asymptotic comparison. |
| 6(b) | `python/figure_06b_continuous_calendar_time.py` | Continuous calendar time with participation. |
| 7 | `r/figure_07_heterogeneous_learning.R` | Five privacy regimes; final plot uses sigma = 1, N = 10,000, five paths per regime. |
| J.1 | `r/figure_J1_llr_evolution.R` | Legacy deterministic positive-report recurrence and asymptotic comparison; see caveat below. |
| Table 1 | No numerical script required | Summary of theoretical results. |

See [verification and limitations](docs/VALIDATION.md) before interpreting these files
as exact reproductions of the submission. This is an organized source release, not a
claim that every numerical result or theoretical assertion has been independently verified.

## Python

Tested environment: Python 3.13, NumPy 2.2.4, SciPy 1.15.2, Matplotlib 3.10.5.
From the repository root:

```sh
python -m venv .venv
# Activate .venv using your shell's normal command, then:
python -m pip install -r requirements.txt
python run_python.py 2 3
python run_python.py 6a 6b --quick
```

For the original continuous settings (B = 10, sigma = 1, 25 epsilon values,
300 Monte Carlo paths, report cap 50,000), use:

```sh
python run_python.py 6a 6b
```

Full continuous runs can be lengthy. `--quick` changes B to 2, epsilon to {0.5, 1},
Monte Carlo paths to 100 and the Bellman grid to 81; it tests execution only.
Continuous scripts and their default model parameters are unchanged. The runner saves
plots without requiring a graphical desktop. Output is under `outputs/full/python/`
or `outputs/quick/python/`; scripts that use `figures/` save inside that subfolder.
Figure 3 saves `binary_platform_objective.pdf` directly in the output directory.

## R

Use R 4.x with `Rscript` on PATH (or use its full installed path):

```sh
Rscript --vanilla r/install_packages.R
Rscript --vanilla run_r.R 1 5 J1
Rscript --vanilla run_r.R 7 --quick
Rscript --vanilla run_r.R 7
```

The runner opens PDF devices explicitly and writes to `outputs/full/r/` or
`outputs/quick/r/`. Figure 7's original script contains preliminary experiments,
so its PDF includes preliminary plots before the final five-regime comparison.
Quick mode reduces only the final Figure 7 simulation to N = 100 and two paths per
regime. Package installation is explicit in `r/install_packages.R`.

## Contents and provenance

- `python/`, `r/`: current manuscript-related scripts with consistent names.
- `supplementary/`: related exploratory Python/R programs and a Mathematica notebook;
  these are not certified generators of current manuscript figures.
- `reference_figures/`: selected existing local PDFs, preserved as reference assets;
  they may predate the submission or the current code.
- `archive/original_github/`: pre-update public scripts and README, preserved verbatim.
- `docs/source_manifest.json`: original relative locations, source SHA-256 hashes
  and copy/edit notes.
- `docs/VALIDATION.md`: checks performed and remaining reproducibility issues.

Original research directories were left intact. WhatsApp datasets, network-inference
code, unrelated simulations, duplicate nested folders, R session dumps and temporary
files are excluded. No private datasets are needed by the selected main scripts.

The local Figure 7 parameters differ from the old GitHub version (sigma = 3,
N = 2,000); the local sigma = 1 version matches the submission caption. The old
`Figure_4_Expectation.R` is a different experiment, not current Figure 4, and is
kept in `supplementary/r/`. Its local sigma = sqrt(2) also differs from GitHub's
sigma = 1. The J.1 script uses sigma = 1 and now initializes the previously missing
`ln_values_new`. Its flip amplitude remains 1/(1+exp(epsilon)), unlike the 1/2
amplitude in the current main smooth-RR scripts; this requires author review for
exact manuscript consistency. Mathematical formulas were not silently rewritten. Figure 3 now plots all three
manuscript signal accuracies using the existing objective function; the original
single-p script is preserved in `supplementary/python/binary_objective_single_p.py`.

## Reference

Privacy-Aware Sequential Learning. The pre-existing repository cites
https://arxiv.org/abs/2502.19525. Use the bibliographic details of the version you cite.
