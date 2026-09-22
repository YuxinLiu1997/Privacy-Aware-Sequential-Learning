# Verification and known limitations

## Checks completed

- Parsed all included Python scripts successfully (syntax only for supplementary scripts).
- Executed the Figure 2 and updated Figure 3 scripts at their default settings.
  Figure 3 now includes p = 0.55, 0.7, 0.9; the objective formulas were preserved.
- Executed both Figure 6 programs through the quick runner: B = 2, sigma = 1,
  epsilon = 0.5 and 1, 100 Monte Carlo paths, 5,000-report cap, 81-point Bellman grid.
  All paths hit the stopping boundary. Bellman report means were 8.34 and 4.25;
  Monte Carlo report means were 8.61 and 4.37. Calendar results correctly divide
  these values by participation probabilities 0.9 and 0.8.
- Compared both continuous scripts' closed-form Gaussian report probabilities with
  independent SciPy quadrature: 108 comparisons at beliefs -4, 0, 4; states -1, +1;
  epsilon 0.1, 0.5, 2; sigma = 1. Tolerance: atol=1e-9, rtol=1e-8. All passed.
- Executed R Figures 1, 5 and J.1 with their original settings, and Figure 7 with
  the final simulation reduced to N = 100 and two paths per regime. All completed
  and produced PDF files. Figure 7's preliminary experiments were left unchanged.
- Python environment: 3.13, NumPy 2.2.4, SciPy 1.15.2, Matplotlib 3.10.5.
  R environment: 4.4.2 with the required packages already installed. R 4.5.3 did
  not have the dependencies; no system packages were installed or changed.
- R reported locale warnings and a viridisLite build-version warning (4.4.3).
  These did not stop execution. The J.1 integration library prints many diagnostic
  messages about integration over infinite domains.
- Existing public scripts are preserved byte-for-byte in archive/original_github/.
  Source paths and source hashes are recorded in source_manifest.json.

## Limitations requiring attention

1. Current Figure 4: no standalone generating code was located in the scanned
   research code directories. normal_distribution.pdf is an existing candidate
   asset only. This release does not assert a complete reproduction of Figure 4.
2. J.1: the legacy code uses amplitude 1/(1+exp(epsilon)), whereas the current
   main smooth randomized-response scripts use 1/2. It iterates a deterministic
   positive-report recurrence, not the full stochastic learning process. We fixed
   its undefined ln_values_new variable but preserved its mathematical expressions.
   The author should check it against the current appendix before claiming exact
   consistency with the submission. Its local sigma = 1 also differs from the
   previous GitHub version's sqrt(2).
3. Full Figure 6 runs at B = 10 and full Figure 7 simulations at N = 10,000 were
   not rerun. Quick-mode results are execution checks, not replacement paper figures.
   The small-B asymptotic approximation is not expected to agree closely.
4. Existing files in reference_figures/ were copied from local outputs and may be
   older than the current scripts. They were not certified as identical to the
   submission. Newly produced verification plots are in the local outputs/ folder,
   which is excluded from Git and the code ZIP.
5. Supplementary and archived programs were not fully executed. In particular,
   some older stopping-time scripts use nanmean over completed paths; when paths
   are censored this is not an unconditional expected stopping time. They are kept
   as historical experiments rather than recommended reproduction entry points.
6. Figure 7's R script contains several preliminary experiments with older labels;
   the last PDF page is the manuscript-related five-regime comparison. The final
   sigma = 1, N = 10,000 settings match the supplied submission's regime/scale setup.

## Deliberate source changes

- Consistent folder names and current manuscript figure numbering.
- J.1: added ln_values_new <- simulate_ln_new(ln_start, n_steps).
- Figure 3: plot the three manuscript p values instead of only p = 0.7, with a legend;
  retain the original single-p script under supplementary/python/.
- Added headless Python and explicit-PDF R runners, dependency instructions, and
  source provenance. Quick settings are applied in memory by the runners.
- Prior public files remain in archive/original_github/ and in Git history.

The code ZIP includes SHA256SUMS for its payload files. Its outer checksum is
written to the adjacent .zip.sha256 file. .gitattributes preserves bytes across
platforms so these payload checksums also apply to the released source tree.
