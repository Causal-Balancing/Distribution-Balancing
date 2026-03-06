# Distribution Balancing: Reproducibility Repository

This repository contains code to reproduce numerical results for the manuscript on Distribution Balancing methods for treatment effect estimation.

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
  "quantregForest", "osqp", "BMisc", "ATE.ncb",
  "ggplot2", "Rmisc", "plotrix", "multcomp"
))
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

## Current reproducibility notes

Some legacy scripts still include machine-specific settings from the original development environment, for example:

- hard-coded `setwd("Z:/...")` or `.libPaths("Z:/...")`
- references to non-existent files such as `sim_utils/datagen_new.R` or `sim_utils/datagen.R`

The master script flags these issues before execution. For public release, we recommend gradually standardizing script headers to use relative paths from repository root.

## Output

Most simulation/application scripts create timestamped output folders and CSV summaries. Output location is controlled inside each script.

## Citation

If you use this code, please cite the corresponding manuscript.

## Contact

For questions or bug reports, please open a GitHub issue.
