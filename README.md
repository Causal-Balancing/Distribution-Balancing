# Distribution Balancing: Reproducibility Repository

## Paper information

**Title:** *Estimating Counterfactual Distribution Functions via Distribution Balancing with Applications*  
**Authors:** Zongwu Cai, Ying Fang, Ming Lin, Yaqian Wu  
**Paper link:** [Google Scholar entry](https://scholar.google.com/citations?view_op=view_citation&hl=zh-CN&user=oHPEDJYAAAAJ&sortby=pubdate&citation_for_view=oHPEDJYAAAAJ:a9-T7VOCCH8C)

This repository contains code to reproduce numerical results for the manuscript on Distribution Balancing methods.

## Repository structure

- `R/`: Core implementation of the proposed estimation and inference methods.
- `sim_utils/`: Data-generating processes, benchmark methods, and utility functions used in simulations.
- `Simulations/`: Scripts for simulation tables (Section 4 and appendix tables).
- `Applications/`: Scripts and datasets for empirical and semi-synthetic applications (Section 5 and appendix).
- `Main_Replicate_All.R`: Master script for preflight checks and batch replication.

## Requirements

- R >= 4.2 (tested with R 4.5 environment during development)
- Recommended packages:
  - `foreach`, `parallel`, `doParallel`, `doRNG`
  - `quantreg`, `drf`, `quantregForest`, `osqp`, `BMisc`, `ATE.ncb`
  - Some scripts may additionally require `ggplot2`, `Rmisc`, `plotrix`, `multcomp`

Install packages (example):

```r
install.packages(c(
  "foreach", "doParallel", "doRNG", "quantreg", "drf",
  "quantregForest", "osqp", "BMisc", 
  "ggplot2", "Rmisc", "plotrix", "multcomp"，"devtools"
))

devtools::install_github("raymondkww/ATE.ncb")
devtools::install_github("chadhazlett/KBAL")
```

## Quick start

From repository root:

```bash
Rscript Main_Replicate_All.R
```

By default, the master script performs:

1. Preflight checks for common path/dependency issues.
2. Script discovery in `Simulations/` and `Applications/`.
3. Optional batch execution with per-script logs.

To run only simulation scripts:

```r
# In R
source("Main_Replicate_All.R")
run_all(replicate_simulations = TRUE, replicate_applications = FALSE, dry_run = FALSE)
```

## Function reference 

Below is a concise reference for functions in `R/` and `sim_utils/`.

### `R/main_estimation.R`
- `fun_CCDF_balance_new`: Main estimator (new version) for distribution balancing QTE with cross-fitting or full-sample path.
- `compute_all_counterfactuals`: Computes fold-wise and aggregated counterfactual CDFs under balanced weights.
- `full_sample_moment_balance`: Runs full-sample moment balancing.
- `fun_CCDF_balance_old`: Computes full sample counterfactual CDFs under balanced weights..

### `R/ccdf_estimation.R`
- `estimate_all_ccdfs_cv`: Estimates conditional CDFs with cross-fitting (separate treatment/control models).
- `est_ccdf`: Dispatcher for CCDF estimators by method name.
- `est_ccdf_kernel`: Kernel-based CCDF estimator.
- `est_ccdf_llr`: Local-linear-regression CCDF estimator.
- `est_ccdf_quantile`: Quantile-regression-based CCDF estimator.
- `est_ccdf_quantregforest`: Quantile regression forest CCDF estimator.
- `est_ccdf_distforest`: Distributional random forest CCDF estimator.
- `estimate_all_ccdfs`: Full-sample CCDF estimation.

### `R/weights_estimation.R`
- `get_balanced_weights`: Solves for balancing weights for one treatment arm at one quantile level.
- `build_balance_design`: Builds balancing feature matrix (local CCDF features + moment features).
- `weight_w`: Wrapper selecting weighting objective (`entropy` or `variance`).
- `weight_entropy`: Entropy balancing solver via smooth unconstrained optimization on dual variables.
- `weight_var`: Variance-minimizing balancing weights via quadratic programming (`osqp`).
- `.clamp_prob`: Clips probabilities away from 0 and 1 for numerical stability.
- `.uniform_weights`: Creates uniform weights.
- `.drop_zero_var_cols`: Removes near-zero-variance balance features.

### `R/utils.R`
- `make_cvgroup`: Randomly assigns K-fold group labels.
- `make_cvgroup_balanced`: Stratified random folds balanced by treatment group.
- `make_ordered_group`: Ordered two-way split helper.
- `make_ordered_group_balanced`: Ordered split helper with treatment balance.
- `get_ate`: Computes weighted ATE.
- `get_cdf`: Computes weighted empirical CDF on a grid.
- `get_qte`: Computes weighted quantiles and QTE from CDF inversion.
- `quantile_lookup`: Inverts CDF to quantiles on a discrete outcome grid.
- `ker`: Univariate kernel function family.
- `ker_h`: Multivariate kernel matrix with mixed continuous/discrete support.
- `fun_K`: Higher-order polynomial kernel.
- `Kh`: Higher-order boundary kernel for mixed covariates.

### `R/inference_bootstrap.R`
- `fun_CCDF_balance_with_inference`: Main estimator + multiplier-bootstrap inference wrapper.
- `bootstrap_qte_multiplier`: Multiplier bootstrap for QTE standard errors and confidence intervals.
- `.multiplier_draw`: Draws multiplier bootstrap weights (`poisson` or `exponential`).
- `weight_entropy_multipiler`: Re-solves entropy balancing under multiplier perturbations.

### `R/test_bootstrap.R`
- `fun_CCDF_balance_with_test`: Main estimator + stochastic-dominance test wrapper.
- `sd_test_multiplier`: Multiplier-bootstrap test for distributional dominance.

### `sim_utils/sim_datagen.R`
- `data1`: Simulation DGP 1 (low-dimensional nonlinear setup).
- `data2`: Simulation DGP 2 (moderate-dimensional nonlinear setup).
- `data3`: Simulation DGP 3 (high-dimensional setup).
- `data_ks_1`: Simulation DGP 5.
- `data_ks_2`: Simulation DGP 6.
- `makedata.sim`: Kallus etal. (2024) simulation DGP.

### `sim_utils/sim_emp_datagen.R`
- `emp_data_gen01`: Semi-synthetic empirical-scenario generator 1.
- `emp_data_gen02`: Semi-synthetic empirical-scenario generator 2.

### `sim_utils/SLE.R`
- `SLE`: SLE benchmark estimator with polynomial propensity features.
- `powers_of_x`: Builds polynomial basis expansion up to degree `d`.
- `ipw_fit`: Standard IPW benchmark estimator for ATE/QTE.
- `kernel_balance`: Kernel balancing benchmark estimator.

### `sim_utils/energy_balance.R`
- `energy_balance`: Energy-distance balancing benchmark solved as quadratic optimization.

### `sim_utils/Other_SD_test.R`
- `DH_test`: DH-style stochastic-dominance test with multiplier bootstrap critical values.
- `SLE_test`: SLE-based stochastic-dominance test with bootstrap p-values.

### `sim_utils/output.R`
- `output`: Summarizes simulation results (bias, RMSE, etc.) against true QTE.

### `sim_utils/LDML.R`
- `make.cvgroup`: Random fold assignment (LDML module).
- `make.cvgroup.balanced`: Treatment-balanced fold assignment (LDML module).
- `make.cvgroup.balanced2`: Fold assignment with treatment/instrument balance.
- `cross.fit.propensities`: Cross-fitted propensity estimation utility.
- `solve.cumsum`: Inverse-CDF helper used by quantile routines.
- `density_`: Weighted density estimator helper.
- `check.data`: Input validation for treatment-effect routines.
- `check.data2`: Input validation for IV routines.
- `est.quantile.ipw`: Quantile treatment effect via IPW.
- `est.quantile.ldml`: LDML quantile treatment effect estimator.
- `est.quantile.dml`: DML quantile treatment effect estimator.
- `est.ivquantile.ipw`: IV quantile treatment effect via IPW.
- `est.ivquantile.ldml`: IV quantile treatment effect via LDML.
- `const`: Constant-model learner wrapper.
- `boost`: Boosting learner wrapper.
- `forest`: Random-forest learner wrapper.
- `neuralnet`: Neural-network learner wrapper.
- `lassor`: Lasso learner wrapper.
- `logistic`: Logistic learner wrapper.
- `reglm`: GLM learner wrapper.
- `forestcdf`: Forest-based conditional CDF learner.

## Citation

If you use this code, please cite:

Cai, Z., Fang, Y., Lin, M., and Wu, Y. (2026). *Estimating Counterfactual Distribution Functions via Distribution Balancing with Applications*.

## Correspondence

- Please contact the corresponding author listed in the manuscript for questions about methods, replication, or data details.

## Contact

For questions or bug reports, please open a GitHub issue.
