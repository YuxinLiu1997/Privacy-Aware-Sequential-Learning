# 📚 Sequential Learning under Privacy Constraints — Figures

This repository contains R scripts used to generate the figures in our study on **sequential learning with privacy-preserving mechanisms**. Each script corresponds to a specific figure in the paper and can be executed independently.

---

## 📁 Repository Structure

| File | Description |
|---|---|
| `Figure_1_Probability_of_correct_cascade.R` | **Probability of correct cascade vs. privacy budget** $\varepsilon$ across sequence lengths \(k\in\{2,3,4,5\}\) and signal accuracy \(p\in\{0.90,0.70,0.55\}\).<br>**Key steps:** computes \(w_k(\alpha)\) and thresholds \(\varepsilon_k=\log\!\big(\frac{1-w_k}{w_k}\big)\) with \(\alpha=(1-p)/p\); evaluates the closed-form probability \(\Pr[\text{right cascade}]\) on each segment \((\varepsilon_{k+1},\varepsilon_k]\); overlays **non-private** baselines \(p^2/(2p^2-2p+1)\).<br>**Visual cues:** colored lines for \(k\), line types for \(p\); vertical dashed lines at \(\varepsilon_{k+1}\). |
| `Figure_3_Smooth_randomized_response.R` | **Smooth randomized response around a threshold** \(t_{\ln}\). Three side-by-side panels:<br>(a) Densities of signals.<br>(b) smooth randomized response.<br>(c) probability of changing action.<br>**Defaults:** \(\mu_{\pm}=\pm1\), \(\sigma=1\), \(\varepsilon=1\), \(t_{\ln}=-3\). |
| `Figure_4_Expectation.R` | **Expected stopping time** \(\mathbb{E}[\tau]\) as a function of \(\varepsilon\) (here \([0.5,1]\)).<br>**Model components:** \( \sigma=1 \); \( C(\varepsilon)=\tfrac{\varepsilon\sigma^2}{2}\!\left(e^{\varepsilon+\frac{\varepsilon^2\sigma^2}{2}}-e^{-\varepsilon+\frac{\varepsilon^2\sigma^2}{2}}\right) \); \( s=\tfrac{2}{\varepsilon\sigma^2} \). Approximates \(\zeta(s)\) by \( \sum_{n=1}^{N} n^{-s} \) with \(N=10{,}000\).<br>**Output:** line plot of \(\mathbb{E}[\tau](\varepsilon)=C_1\,C(\varepsilon)^{-s}\zeta(s)\) and the minimizing \(\varepsilon\) printed to console. |
| `Figure_5_Heterogeneous_likelihood_ratio.R` | **Finite-sample trajectories of the log-likelihood ratio** \(l_n\) under multiple privacy regimes, using exact Bayesian updates for action distributions.<br>**Simulators:**<br>• `simulate_seq_learning` (non-private baseline).<br>• `simulate_seq_learning_rand` (fixed \(\varepsilon\)).<br>• `simulate_seq_learning_epsU_sample_ab` (heterogeneous \(\varepsilon\sim U[\ell,h]\).<br>•**Comparisons included:** heterogenous \(U(0,1)\), fixed \(\varepsilon\in\{0.1,0.5,1\}\), and non-private. Multiple seeds overlayed; a helper highlights the “center” trajectory closest to the mean. **Default signal noise:** \(\sigma=3\). |
| `Figure_K1_Evolution_of_log_likelihood_ratio.R` | **Clean comparison of \(l_n\) growth** for non-private vs. private \(\varepsilon\) (e.g., 0.5, 1.0).

---

## 🔧 Requirements

- **R** ≥ 4.0  
- Packages: `ggplot2`, `viridis`, `gridExtra`, `cowplot`, `gsl`  
- Unless noted in the file, **signals are normal with \(\sigma=3\)** and fixed seeds are used for reproducibility.
- 
## 📄 Reference

If you use this repository, please cite our paper:

**Privacy-Aware Sequential Learning**  
arXiv:https://arxiv.org/abs/2502.19525
